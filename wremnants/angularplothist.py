import h5py
import matplotlib.pyplot as plt
import mplhep as hep
import hist
import numpy as np
from wums import ioutils  



#This plots the smoothed histogram used for alpha s extraction



 #Apply the standard CMS plotting styles
hep.style.use("CMS")

# Filepath matching line 278 of your angularhistmaker.py
filepath = "outputs/angularhistmaker.hdf5"

print(f"Opening {filepath}...")
with h5py.File(filepath, "r") as f:
    # 1. Safely unpack the pickled dictionary payload
    payload = ioutils.pickle_load_h5py(f["results"])

    # 2. Dynamically locate the dataset key (e.g., "Zmumu_2016PostVFP")
    dataset_key = list(payload.keys())[0]
    print(f"Detected Dataset Key: '{dataset_key}'")

    # 3. Extract the output dictionary and find the booked scale histogram
    output_dict = payload[dataset_key]["output"]

    # This matches f"{dataset_name}_smoothed_qcd_scale" from line 198
    hist_key = next((k for k in output_dict.keys() if "smoothed_qcd_scale" in k), None)

    if not hist_key:
        # Fallback to look for any key containing 'scale' if naming changed
        hist_key = next((k for k in output_dict.keys() if "scale" in k), list(output_dict.keys())[0])

    print(f"Loading Histogram Key: '{hist_key}'")

    # 4. RESOLVE THE PROXY HERE (while the file 'f' is still open)
    h_obj = output_dict[hist_key]
    h_tensor = h_obj.get() if hasattr(h_obj, "get") else h_obj

# The HDF5 file is now safely closed, but h_tensor is fully loaded in memory!

# Let's inspect the exact physical structure of your loaded histogram
print("\nSuccessfully loaded histogram! Detected axes:")
for ax in h_tensor.axes:
    print(f"  - Axis name: '{ax.name}', size: {ax.size}")

# 5. Dynamic Multi-Dimensional Slicing
# We want to keep ONLY the physical axes ('ptVgen' and 'absYVgen').
# Any variation/tensor axis (like 'pdfCT18ZWeights') will be collapsed to 0 (nominal).
slice_dict = {}
for ax in h_tensor.axes:
    if ax.name not in ["ptVgen", "absYVgen"]:
        slice_dict[ax.name] = 0

print(f"\nApplying nominal slice selection: {slice_dict}")
h_2d = h_tensor[slice_dict]

# 6. Project down to 1D distributions
h_1d_absy = h_2d[{"ptVgen": sum}]
h_1d_pt = h_2d[{"absYVgen": sum}]

# --- Plot 1: Rapidity (absYVgen) ---
fig, ax = plt.subplots(figsize=(8,8))
hep.histplot(
    h_1d_absy, 
    ax=ax, 
    histtype="step", 
    color="darkblue", 
    label=r"$Z \rightarrow \mu^+\mu^-$ (Nominal)"
)

ax.set_xlabel(r"$|Y^{V}| \text{ (Gen)}$")
ax.set_ylabel("Events / Bin")

# Dynamically set X-limits to match the exact binning edges defined in your axis
ax.set_xlim(h_1d_absy.axes[0].edges[0], h_1d_absy.axes[0].edges[-1])
hep.cms.label(ax=ax, llabel="Preliminary Simulation", rlabel="13.6 TeV")

plt.savefig("rapidity_plot.png", bbox_inches="tight", dpi=300)
print("Plot successfully saved as rapidity_plot.png!") 

# --- Plot 2: Transverse Momentum (ptVgen) ---
fig, ax = plt.subplots(figsize=(8,8))
hep.histplot(
    h_1d_pt, 
    ax=ax, 
    histtype="step", 
    color="darkblue", 
    label=r"$Z \rightarrow \mu^+\mu^-$ (Nominal)"
)

ax.set_xlabel(r"$p_T^{V} \text{ (Gen) [GeV]}$")
ax.set_ylabel("Events / Bin")

# Dynamically set X-limits to match the exact binning edges defined in your axis
ax.set_xlim(h_1d_pt.axes[0].edges[0], h_1d_pt.axes[0].edges[-1])
hep.cms.label(ax=ax, llabel="Simulation ", rlabel="13.6 TeV")

plt.savefig("pt_plot.png", bbox_inches="tight", dpi=300)
print("Plot successfully saved as pt_plot.png!")



# This plots the weight distribution



# # Apply standard CMS plotting style
# mh.style.use("CMS")

# filepath = "outputs/angularhistmaker.hdf5"

# # --- 1. Load the data from the HDF5 file ---
# print(f"Opening {filepath}...")
# with h5py.File(filepath, "r") as f:
#     payload = ioutils.pickle_load_h5py(f["results"])
#     dataset_key = list(payload.keys())[0]
#     print(f"Detected Dataset Key: '{dataset_key}'")
    
#     output_dict = payload[dataset_key]["output"]
    
#     hist_key = next((k for k in output_dict.keys() if "weights_hist" in k or "weight" in k), None)
#     if not hist_key:
#         raise KeyError("Could not find a weight histogram in the output file!")
        
#     print(f"Loading Histogram Key: '{hist_key}'")
#     h_obj = output_dict[hist_key]
#     weights_hist = h_obj.get() if hasattr(h_obj, "get") else h_obj

# print("Successfully loaded weights histogram!")

# # --- 2. Projection & Direct Bin Slicing ---
# # Project the entire 2D histogram down to a 1D weight distribution
# weight_hist_1d = weights_hist.project("weight_val") 

# # Slice using direct bin indices on the first axis (ptVgen)
# # Bin 0 corresponds to [0, 1] GeV
# weight_hist_0_1 = weights_hist[0, :]

# # Bin 1 corresponds to [1, 2] GeV
# weight_hist_1_2 = weights_hist[1, :]

# # --- 3. Plotting Canvas 1: Full Event Weights ---
# fig_w, ax_w = plt.subplots(figsize=(8, 6))  
# mh.histplot(
#     weight_hist_1d, 
#     ax=ax_w, 
#     color="darkviolet", 
#     histtype="step", 
#     linewidth=2, 
#     label="All Event Weights"
# )
# ax_w.set_xlabel("PDF Weight Value")
# ax_w.set_ylabel("Events")
# ax_w.set_title("Diagnostic: Total Event Weight Distribution", fontsize=12, pad=12)

# # Official CMS experiment headers (removes overlapping titles)
# mh.cms.label(ax=ax_w, data=False, label="Preliminary", rlabel="13.6 TeV", loc=1)

# ax_w.legend(loc='upper right', frameon=False) 
# plt.tight_layout()  # Automatically adjusts margins to prevent label clipping
# plt.savefig("Total_Weight_Distribution.png", bbox_inches='tight', dpi=300)
# # --- 4. Plotting Canvas 2: Sliced Comparisons (Density Normalized) ---
# fig_slice, ax_slice = plt.subplots(figsize=(8, 6))

# mh.histplot(
#     weight_hist_0_1,
#     ax=ax_slice,
#     color="teal",
#     histtype="step",
#     linewidth=2, 
#     density=True,
#     label=r"$0 < p_T < 1$ GeV"
# )

# mh.histplot(
#     weight_hist_1_2,
#     ax=ax_slice,
#     color="crimson",
#     histtype="step",
#     linewidth=2,
#     linestyle="--", 
#     density=True,
#     label=r"$1 < p_T < 2$ GeV"
# )

# ax_slice.set_xlabel("PDF Weight Value")
# ax_slice.set_ylabel("Probability Density")
# ax_slice.set_xlim(0.95, 1.15)  # Zooms into the high-resolution region

# # Official CMS experiment headers
# mh.cms.label(ax=ax_slice, data=False, label="Preliminary", rlabel="13.6 TeV")

# ax_slice.legend(loc='upper left', frameon=False, fontsize=12) 
# plt.tight_layout()  # Automatically adjusts margins to prevent label clipping
# plt.savefig("Weight_Kinematic_Comparison.png", bbox_inches='tight', dpi=300)

# # --- 5. Statistical Analysis Engine ---
# def print_stats(target_hist, label_name):
#     bin_centers = target_hist.axes[0].centers
#     bin_contents = target_hist.values()
#     total_events = np.sum(bin_contents)
    
#     if total_events == 0:
#         print(f"--- STATS FOR {label_name} --- Empty Histogram!")
#         return

#     mean_weight = np.sum(bin_contents * bin_centers) / total_events
#     variance_weight = np.sum(bin_contents * (bin_centers - mean_weight)**2) / total_events
#     std_dev_weight = np.sqrt(variance_weight) 

#     print("\n" + "="*45)
#     print(f"--- STATS FOR {label_name} ---")
#     print(f"Total Events:       {total_events:,.2f}")
#     print(f"Mean Weight:        {mean_weight:.6f}")
#     print(f"Standard Deviation: {std_dev_weight:.6f}")
#     print("="*45)

# # Run diagnostics
# print_stats(weight_hist_0_1, "pT: 0 to 1 GeV")
# print_stats(weight_hist_1_2, "pT: 1 to 2 GeV")