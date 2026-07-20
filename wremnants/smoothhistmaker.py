import os
import ROOT
import narf
import hist
from wremnants.production import datasets
from wremnants.utilities import binning 
from wremnants.utilities import parsing
from wremnants.utilities import common
from wremnants.production.datasets import dataset_tools 
from wremnants.production import systematics, generator_level_definitions as gld 
from wums import ioutils 
import h5py
import subprocess  


def build_graph(df, dataset):
    results = []
    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    axis_pt = hist.axis.Regular(40, 0, 80, name="ptVgen", underflow=False, overflow=False)
    axis_eta = hist.axis.Regular(50, -2.5, 2.5, name="absYVgen", underflow=False, overflow=False)
    axes = [axis_pt, axis_eta]
    cols = ["ptVgen", "absYVgen"]


    scale_tensor_axes = systematics.scale_tensor_axes
    storage = hist.storage.Weight()
    base_name = dataset.name

    print("Available columns:", [str(c) for c in df.GetColumnNames() if "lepton" in str(c).lower()])

    df = df.Define("isEvenEvent", "event % 2 == 0")

    df = gld.define_prefsr_vars(df)



    print(df.GetColumnType("csSineCosThetaPhigen"))


    # Handle scale weights and helicity tensors (reco-level, unchanged)
    if "scaleWeights_tensor" not in df.GetColumnNames():
        df = df.Define("scaleWeights_tensor", "wrem::makeScaleTensor(LHEScaleWeight, 1.0)")

    if "helicity_xsecs_scale_tensor" not in df.GetColumnNames():
        df = df.Define(
            "helicity_xsecs_scale_tensor",
            "wrem::makeHelicityMomentScaleTensor(csSineCosThetaPhigen, scaleWeights_tensor, weight)"
        )


    helicity_xsecs_scale = df.HistoBoost(
        f"{base_name}_helicity_xsecs_scale",
        axes,
        [*cols, "helicity_xsecs_scale_tensor"],
        tensor_axes=[binning.axis_helicity, *scale_tensor_axes],
        storage=storage,
    )



    results.append(helicity_xsecs_scale)
    weightsum = (df.Sum("weight"), df.Count())
    return results, weightsum

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

# Make sure the outputs directory exists
os.makedirs("outputs", exist_ok=True)

# Save the results to the HDF5 file your plotting script looks for
with h5py.File("outputs/smoothhistmaker.hdf5", "w") as f:
    ioutils.pickle_dump_h5py("results", out, f) 

#     print("HDF5 file successfully written! Launching plotting script...")

# # Now launch the plotting script automatically
# subprocess.run(["python", "wremnants/smoothplothist.py"], check=True)
# print("✅ Data processing complete! Output saved to outputs/smoothhistmaker.hdf5")

