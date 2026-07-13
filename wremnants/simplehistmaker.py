import ROOT
import hist 
import narf 
import numpy as np
import subprocess  
import pickle  
import os
from wremnants.utilities import parsing  
from wremnants.production.theory_corrections import make_theory_corr_weight_info
from wremnants.production import theory_corrections
from wremnants.production.histmaker_tools import (write_analysis_output,
)
from wremnants.production.datasets.dataset_tools import getDatasets
import matplotlib.pyplot as plt
import mplhep as mh

parser, initargs = parsing.common_parser()
args = parser.parse_args()


def build_graph(df, dataset):
    results = []
    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")
    weightsum=df.SumAndCount("weight") 
    df = df.Define("theory_weight_truncate", "10.")
    df = df.Filter("nMuon == 2", "Select events with exactly 2 muons")
    # create histograms, append to results
    df = df.Define("Muon_pt_vector_sum", 
                   "std::sqrt(Muon_pt[0]*Muon_pt[0] + Muon_pt[1]*Muon_pt[1] + 2*Muon_pt[0]*Muon_pt[1]*std::cos(Muon_phi[0] - Muon_phi[1]))"
)
    weight_info = make_theory_corr_weight_info("ct18z", alphas=True, renorm=True)
    values = weight_info["weights"]
    df = df.Define(
            f"pdfCT18ZASWeights_tensor",
            f"Eigen::TensorFixedSize<double, Eigen::Sizes<{len(values)}>> res; "
            + "; ".join(
                [
                    f"res({i}) = {entry}"
                    for i, entry in enumerate(values)
                ]
            )
            + "; return res;",
        )
    df= df.Define("pdfCT18ZWeights","pdfCT18ZASWeights_tensor[0]")
    hist_muon_pt_sum = df.Histo1D(
    ("muon_pt_vector_sum", "Magnitude of Muon Vector p_{T} Sum;p_{T} [GeV];Events", 100, 0, 100.),
    "Muon_pt_vector_sum" 
)
    
    axis_alphaS = hist.axis.StrCategory(["central", "alphaSDown", "alphaSUp"],name="alphaS")
    axis_muon_pt_sum = hist.axis.Regular(100, 0, 100, overflow=False, underflow=False, name="A")
    axis_weights = hist.axis.Regular(100, 0, 2, overflow=False, underflow=False, name="weight")
    hist_Z_pt_alphaS = df.HistoBoost("muon_pt_vector_sum_alphas_up",[axis_muon_pt_sum],["Muon_pt_vector_sum", "pdfCT18ZASWeights_tensor"],tensor_axes=[axis_alphaS])
# Weight histogram creation 
    hist_Z_pt_weight = df.HistoBoost("weights_hist",[axis_muon_pt_sum, axis_weights],["Muon_pt_vector_sum", "pdfCT18ZWeights"])

    results.append(hist_muon_pt_sum) 
    results.append(hist_Z_pt_alphaS) 
    results.append(hist_Z_pt_weight)
   
    
    return results, weightsum



datasets = getDatasets(
    maxFiles=100,
    filt="Zmumu_2016PostVFP",
    nanoVersion="v9",
    base_path='/scratch/submit/cms/wmass/NanoAOD',
    era="2016PostVFP",
)
print("process datasets", datasets)
results = narf.build_and_run(datasets, build_graph)
print("data processing finished! Now plotting ...")


dataset_name = datasets[0].name  

# plot muon pt vector sum with alphaS weights
print(results[dataset_name]["output"]["muon_pt_vector_sum_alphas_up"].get())

# plot muon pt vector sum and variations
alphaS_hist_2d = results[dataset_name]["output"]["muon_pt_vector_sum_alphas_up"].get()
muon_hist = alphaS_hist_2d[:, 0]
alphaS_up = alphaS_hist_2d[:, 1]
alphaS_down = alphaS_hist_2d[:, 2]    

# plot weight hist
weight_plot = results[dataset_name]["output"]["weights_hist"].get()

mh.style.use()  
muon_hist_clean = hist.Hist(alphaS_up.axes[0], storage=hist.storage.Weight())
muon_hist_clean.view().value = muon_hist.values()
muon_hist_clean.view().variance = muon_hist.variances() 

hist_data_dump = {
    "muon_hist_clean": muon_hist_clean,
    "alphaS_up": alphaS_up,
    "alphaS_down": alphaS_down,
    "weights_hist": results[dataset_name]["output"]["weights_hist"].get() # Saving your 2D weight histogram too!
}

# Dump them into a static file on your disk
with open("raw_analysis_histograms.pkl", "wb") as f:
    pickle.dump(hist_data_dump, f)

print("\n[Histmaker] Raw histograms saved successfully to raw_analysis_histograms.pkl")
print("[Histmaker] Launching plotter script...")

current_dir = os.path.dirname(os.path.abspath(__file__))


plotter_path = os.path.join(current_dir, "plothist.py")

# Automatically run your new plotting script right now
subprocess.run(["python", plotter_path], check=True)




