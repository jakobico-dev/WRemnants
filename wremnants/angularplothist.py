import h5py
import matplotlib.pyplot as plt
import mplhep as hep
import hist
import numpy as np
from wums import ioutils  # Matches line 12 of your angularhistmaker.py

# Apply the standard CMS plotting styles
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
hep.cms.label(ax=ax, llabel="Preliminary Simulation", rlabel="13.6 TeV")

plt.savefig("smoothed_hist", bbox_inches="tight", dpi=300)
print("Plot successfully saved as pt_plot.png!")