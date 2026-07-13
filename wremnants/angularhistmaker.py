import ROOT
import narf
import hist
from wremnants.production import datasets 
from wremnants.utilities import binning
from wremnants.utilities import parsing 
from wremnants.utilities import common
from wremnants.production.datasets import dataset_tools
from wremnants.production import systematics 
from wremnants.
from wums import ioutils 
import os 
import h5py  


def make_qcd_uncertainty_helper_by_helicity(
    is_z=False,
    filename=f"{common.data_dir}/angularCoefficients/w_z_helicity_xsecs.hdf5",
    rebin_ptVgen=binning.ptV_binning,
    rebin_absYVgen=False,
    rebin_massVgen=False,
    return_tensor=True,
):s
    # Place holder or actual implementation
    pass 

def build_graph(df, dataset): 
    results = [] 
    
    # Define 'info' dictionary so **info doesn't throw a NameError
    info = {"base_name": dataset.name}
    
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
    
    if "antilepton" not in df.GetColumnNames():
        df = df.Define("antilepton", "ROOT::Math::PtEtaPhiMVector(GenDressedLepton_pt[0], GenDressedLepton_eta[0], GenDressedLepton_phi[0], GenDressedLepton_mass[0])")
    
    if "lepton" not in df.GetColumnNames():
        df = df.Define("lepton", "ROOT::Math::PtEtaPhiMVector(GenDressedLepton_pt[1], GenDressedLepton_eta[1], GenDressedLepton_phi[1], GenDressedLepton_mass[1])")
    
    if "csSineCosThetaPhigen" not in df.GetColumnNames():
        df = df.Define("csSineCosThetaPhigen", "wrem::csSineCosThetaPhi(antilepton, lepton)")

    if "genV" not in df.GetColumnNames():
        df = df.Define("genV", "antilepton + lepton")

    if "ptVgen" not in df.GetColumnNames():
        df = df.Define("ptVgen", "genV.pt()")

    if "absYVgen" not in df.GetColumnNames():
        df = df.Define("absYVgen", "std::abs(genV.Rapidity())")

    if "massVgen" not in df.GetColumnNames():
        df = df.Define("massVgen", "genV.mass()")

    if "chargeVgen" not in df.GetColumnNames():
        is_z = "Z" in dataset.name
        df = df.Define("chargeVgen", "0" if is_z else "GenDressedLepton_charge[0]")

    # Initialize the QCD uncertainty helper
    is_z = "Z" in dataset.name
    helper = make_qcd_uncertainty_helper_by_helicity(is_z=is_z)

    # Evaluate helper to generate the required weight tensor column
    generator = "qcd"
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

    # Book the scale histograms (called AFTER defining the qcdWeight_tensor)
    add_qcdScale_hist(results, df, scale_axes, scale_cols, **info) 
    
    weightsum = (df.Sum("genWeight"), df.Count())
    return results, weightsum


# --- Execution Framework ---
parser, initargs = parsing.common_parser()
args = parser.parse_args() 

dataset_list = dataset_tools.getDatasets(
    maxFiles=-1,
    filt="Zmumu_2016PostVFP",
    nanoVersion="v9",
    base_path='/scratch/submit/cms/wmass/NanoAOD',
    era="2016PostVFP",
) 

out = narf.build_and_run(dataset_list, build_graph)

os.makedirs("outputs", exist_ok=True)
with h5py.File("outputs/smoothhistmaker.hdf5", "w") as f:
    ioutils.pickle_dump_h5py("results", out, f)
    
print("Histograms with smoothed QCD scales successfully created!") 