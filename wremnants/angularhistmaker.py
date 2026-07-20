import ROOT
ROOT.ROOT.DisableImplicitMT()
import narf
import hist
import subprocess
from wremnants.production import datasets 
from wremnants.utilities import binning
from wremnants.utilities import parsing 
from wremnants.utilities import common
from wremnants.production.datasets import dataset_tools
from wremnants.production import systematics, generator_level_definitions as gld
from wums import ioutils 
from wremnants.production.theory_corrections import make_qcd_uncertainty_helper_by_helicity as core_helper
from wremnants.production.theory_corrections import make_pdfs_uncertainties_helper_by_helicity
from wremnants.production import helicity_utils
import os 
import h5py 
import numpy as np 

helper_path = os.path.abspath("outputs/smoothed_helper_input_7D.hdf5")
# logger = logging.child_logger(__name__)


def reformat_helper_input(
    boson_key,
    in_path="outputs/smoothhistmaker.hdf5",
    dataset_key="Zmumu_2016PostVFP",
    hist_key="Zmumu_2016PostVFP_helicity_xsecs_scale",
):
    print("Reformatting smoothhistmaker output for the helper...")
    with h5py.File(in_path, "r") as f_in:
        payload = ioutils.pickle_load_h5py(f_in["results"])
        h_tensor = payload[dataset_key]["output"][hist_key].get()

    nominal_slice = h_tensor[{"muRfact": 1j, "muFfact": 1j}]
    if np.any(nominal_slice.values() == 0):
        print("WARNING: Zeroes detected in nominal slice! Normalization will fail.")
    
    # This is a temporary fix to allow the code to run
    h_tensor.view().value[h_tensor.view().value == 0] = 1e-9
        
    ax_mass = hist.axis.Regular(1, 60.0, 120.0, name="massVgen", flow=True)
    ax_charge = hist.axis.Regular(1, -2.0, 2.0, name="chargeVgen", flow=True)
    
    # Define the exact 7D axis layout expected by the C++ helper
    axes_7d = [
        ax_mass,
        h_tensor.axes["absYVgen"],
        h_tensor.axes["ptVgen"],
        ax_charge,
        h_tensor.axes["helicity"],
        h_tensor.axes["muRfact"],
        h_tensor.axes["muFfact"]
    ]
    
    #  Initialize empty 7D histograms for both nominal and LHE
    h_nominal_7d = hist.Hist(*axes_7d, storage=hist.storage.Weight())
    h_lhe_7d = hist.Hist(*axes_7d, storage=hist.storage.Weight())
    
    #  Transpose your 5D data to align with the new 7D axis order
    # h_tensor shape is (pt, absy, helicity, muR, muF) -> (40, 50, 9, 3, 3)
    # We need to map this to (absy, pt, helicity, muR, muF) -> (50, 40, 9, 3, 3)
    val_array = h_tensor.view(flow=False)
    val_array_transposed = np.transpose(val_array, (1, 0, 2, 3, 4))
    
    # 5. Populate the 7D LHE histogram
    h_lhe_7d.view(flow=False)[0, :, :, 0, :, :, :] = val_array_transposed
    
    # Create the 3D sliced nominal array and broadcast it to the 7D nominal histogram
    # This selects the nominal scale variation (muR=1, muF=1)
    nominal_slice = h_tensor[{"muRfact": 1j, "muFfact": 1j}].view(flow=False)
    nominal_slice_transposed = np.transpose(nominal_slice, (1, 0, 2)) # (absy, pt, helicity)
    
    # Broadcast nominal values across all scale variations in the nominal 7D hist
    for r in range(3):
        for f in range(3):
            h_nominal_7d.view(flow=False)[0, :, :, 0, :, r, f] = nominal_slice_transposed


    return {
        boson_key: h_nominal_7d ,       
        f"{boson_key}_lhe": h_lhe_7d,  
    }
       

def make_qcd_uncertainty_helper_by_helicity(
    is_z=False,
    filename=f"{common.data_dir}/angularCoefficients/w_z_helicity_xsecs.hdf5",  # look at this filename again before running code
    rebin_ptVgen=False,
    rebin_absYVgen=False,
    rebin_massVgen=False,
    return_tensor=False,
):
    
    return core_helper(
        is_z=is_z,
        filename=filename,
        rebin_ptVgen=rebin_ptVgen,
        rebin_absYVgen=rebin_absYVgen,
        rebin_massVgen=rebin_massVgen,
        return_tensor=return_tensor,
    )
    



def add_qcdScale_hist(results, df, axes, cols, dataset_name):
    """
    Books the multidimensional histogram applying the smoothed QCD weight tensor.
    """
    
    smoothed_hist = df.HistoBoost(
        f"{dataset_name}_smoothed_qcd_scale", #  The name goes first
        axes,                                 #  The hist.axis objects
        cols + ["pdfCT18ZWeights"]           #  Columns list + the new weight
    )
    
    # Append the RDataFrame action to the results list
    results.append(smoothed_hist)

print("Executing 7D reformatting step...")
formatted_data = reformat_helper_input(boson_key="Z")

new_filename = "outputs/smoothed_helper_input_7D.hdf5"
with h5py.File(new_filename, "w") as f_out:
    ioutils.pickle_dump_h5py("results", formatted_data, f_out)
print(f"Successfully saved 7D helper input to {new_filename}")

def build_graph(df, dataset): 
    results = [] 
    
    # Define 'info' dictionary so **info doesn't throw a NameError
    info = {"dataset_name": dataset.name}
    
    # Create new virtual column in dataset
    if dataset.is_data: 
        df = df.DefinePerSample("weight", "1.0")
        df = df.Define("nominal_weight_uncorr", "weight") 
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")  # extracts the sign of genWeight MC
        df = df.Define("nominal_weight_uncorr", "weight")  
        
    axis_pt = hist.axis.Regular(40, 0, 80, name="ptVgen", underflow=False, overflow=False)
    axis_eta = hist.axis.Regular(50, -2.5, 2.5, name="absYVgen", underflow=False, overflow=False)
    scale_axes = [axis_pt, axis_eta]
    scale_cols = ["ptVgen", "absYVgen"]
    
    df=df.Define("isEvenEvent", "event % 2 ==0")

    df = gld.define_prefsr_vars(df)

    is_z = bool("Z" in dataset.name) # Force cast to bool
    filename = "wremnants-data/data/TheoryCorrections/ByHelicity/PDFs/w_z_gen_dists_maxFiles_m1_ct18z_pdfByHelicity_skimmed.hdf5"


    # Initialize the pdf uncertainty helper
    is_z = "Z" in dataset.name
    helpers = make_pdfs_uncertainties_helper_by_helicity(
        # is_z=is_z,
        # filename=filename,
        # return_tensor=True,
        # rebin_ptVgen=False,
        proc="Z",
        pdfs=["ct18z",]
        )
    helper=helpers["pdfCT18Z"]
   
    # Evaluate helper to generate the required weight tensor column
    generator = "pdf"
    df = df.Define(
        f"{generator}Weight_tensor",
        helper,
        [
            "massVgen",
            "absYVgen",
            "ptVgen",
            "chargeVgen",
            "csSineCosThetaPhigen",
            "nominal_weight_uncorr",
        ],
    )
    
    df = df.Define("pdfCT18ZWeights", "pdfWeight_tensor")
   #Store the FULL vector/tensor of PDF weights, not just the 0th element
    df = df.Define("nominal_weight", "pdfWeight_tensor[0]")


    # Book the scale histograms (called AFTER defining the qcdWeight_tensor)
    add_qcdScale_hist(results, df, scale_axes, scale_cols, **info) 
    
    weightsum = (df.Sum("genWeight"), df.Count())
    

    
    axis_pt_weight = hist.axis.Variable([0, 1, 2, 5, 10, 15, 20, 35, 60, 100], name="ptVgen")

    
    axis_weight = hist.axis.Regular(100, 0.95, 1.15, name="weight_val")

    
    df = df.Define("diagnostic_weight_val", "nominal_weight")

    
    results.append(
        df.HistoBoost(
            "weights_hist", 
            [axis_pt_weight, axis_weight], 
            ["ptVgen", "diagnostic_weight_val"]
        )
    )
    cols = df.GetColumnNames()
    print("--- Available Columns in DataFrame ---")
    for col in cols:
        
        if "weight" in str(col):
            print(f"Found weight column: {col}")
    return results, weightsum


#    Execution Framework 
parser, initargs = parsing.common_parser()
args = parser.parse_args() 

dataset_list = dataset_tools.getDatasets(
    maxFiles=-1,
    filt="Zmumu_2016PostVFP",
    nanoVersion="v9",
    base_path='/scratch/submit/cms/wmass/NanoAOD',
    era="2016PostVFP",
) 

formatted_data = reformat_helper_input(boson_key="Z")

os.makedirs("outputs", exist_ok=True)
with h5py.File("outputs/smoothed_helper_input.hdf5", "w") as f_out:
    for key, value in formatted_data.items():
        ioutils.pickle_dump_h5py(key, value, f_out)

out = narf.build_and_run(dataset_list, build_graph)

os.makedirs("outputs", exist_ok=True)
with h5py.File("outputs/angularhistmaker.hdf5", "w") as f:
    ioutils.pickle_dump_h5py("results", out, f)
    
print("Histograms with smoothed QCD scales successfully created!") 