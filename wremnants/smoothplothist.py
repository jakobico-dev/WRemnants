import matplotlib.pyplot as plt
import mplhep as hep
from wums import ioutils 
import h5py 
import numpy as np

# Apply standard CMS styling layout
hep.style.use("CMS")

#  Load the raw high-dimensional tensor data
with h5py.File("outputs/smoothhistmaker.hdf5", "r") as f:
    #  Pass the open file object 'f' into the loader
    payload = ioutils.pickle_load_h5py(f["results"]) 

    h_tensor = payload["Zmumu_2016PostVFP"]["output"]["Zmumu_2016PostVFP_helicity_xsecs_scale"].get()  
    
    #dataset_output = payload["Zmumu_2016PostVFP"]["output"]
    
   # print("--- Available Histograms Found ---")
  #  for key in dataset_output.keys():
  #  print(f" -> {key}")
   # print("----------------------------------")
    
    # locate the new histogram booked by add_qcdScale_hist
    #qcd_keys = [k for k in dataset_output.keys() if "qcd" in k.lower() or "scale" in k.lower()]
    #if not qcd_keys:
     #   raise KeyError("Could not find the new scale histogram! Check the keys printed above.")
    
    #target_key = qcd_keys[0]
    #print(f"Extracting tensor for: {target_key}")
    #h_tensor = dataset_output[target_key].get()
    

# Slice the tensor directly down to a 1D pT distr.
h_1d = h_tensor[{"helicity": 0, "muRfact": 1j, "muFfact": 1j, "absYVgen": sum}] 


# the canvas and draw the histogram step line
#fig, ax = plt.subplots(figsize=(8, 8))
#hep.histplot(h_1d, 
             #ax=ax,
             #histtype="step", 
             #color="darkblue", 
             #label=r"Smoothed QCD Scale $Z \rightarrow \mu\mu$")

fig,ax = plt.subplots(figsize=(8,8))
hep.histplot(h_1d, ax=ax, histtype="step", color="darkblue", label=r"Smoothed $Z \rightarrow \mu^+\mu^-$" )

# add standard axis labels 
ax.set_xlabel(r"$p_{T}^{V} \text{ (Gen)} \text{ [GeV]}$")
ax.set_ylabel("Events / Bin")
ax.set_xlim(0, 80)
hep.cms.label(ax=ax, llabel="Preliminary Simulation", rlabel="13.6 TeV")


plt.savefig("smoothed_helicity_distribution.png", bbox_inches="tight", dpi=300)
print("Plot successfully saved as smoothed_helicity_distribution.png!") 
