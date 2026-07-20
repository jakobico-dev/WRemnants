import mplhep as mh
import pickle
import matplotlib.pyplot as plt
import numpy as np
import hist


# # plot the histogram (aesthetically)

mh.style.use("CMS")

# Load the saved raw data from the histmaker
with open("raw_analysis_histograms.pkl", "rb") as f:
    loaded_data = pickle.load(f)


muon_hist_clean = loaded_data["muon_hist_clean"]
alphaS_up = loaded_data["alphaS_up"]
alphaS_down = loaded_data["alphaS_down"]
weights_hist = loaded_data["weights_hist"]

weight_hist_1d = weights_hist.project(1) 
fig_w, ax_w = plt.subplots(figsize=(8, 6))  

weight_hist_0_1 = weights_hist[hist.loc(0):hist.loc(1), :].project(1)
weight_hist_1_2 = weights_hist[hist.loc(1):hist.loc(2), :].project(1)

fig_slice, ax_slice = plt.subplots(figsize=(8, 6))

mh.histplot(
    weight_hist_1d, 
    ax=ax_w, 
    color="darkviolet", 
    histtype="step", 
    linewidth=2, 
    label="Event Weights"
)
mh.histplot(
    weight_hist_0_1,
    ax=ax_slice,
    color="teal",
    histtype="step",
    linewidth=2, 
    density=True,
    label="Event Weights: 0<pT<1 GeV"
)

# 4. Plot the 1 to 2 GeV slice on the same canvas 
mh.histplot(
    weight_hist_1_2,
    ax=ax_slice,
    color="crimson",
    histtype="step",
    linewidth=2,
    linestyle="--", 
    density=True,
    label="Event Weights : 1<pT<2 GeV"
)



# 3. Clean up the legend frame
ax_w.legend(loc='upper left', frameon=False, fontsize='small')

ax_slice.set_xlabel("PDF Weight Value")
ax_slice.set_ylabel("Events (Normalized)")

# 2. Provide 30% headroom so the peak completely clears the CMS labels
ax_slice.set_ylim(0, ax_slice.get_ylim()[1] * 1.3)

# 3. Apply the official CMS Simulation Preliminary and 13.6 TeV headers
mh.cms.label(ax=ax_slice, data=False, label="Preliminary", rlabel="13.6 TeV")

# 4. Draw the clean, frameless legend for the slice categories
ax_slice.legend(loc='upper right', frameon=False)

# statistical analysis 

target_hist = weight_hist_0_1

bin_centers = target_hist.axes[0].centers
bin_contents = target_hist.values()

total_events = np.sum(bin_contents)

mean_weight = np.sum(bin_contents * bin_centers) / total_events

variance_weight = np.sum(bin_contents * (bin_centers - mean_weight)**2) / total_events
std_dev_weight = np.sqrt(variance_weight) 

print("\n" + "="*40)
print(f"--- STATS FOR MUON pT: 0 to 1 GeV ---")
print(f"Total Events:       {total_events:,.2f}")
print(f"Mean Weight:        {mean_weight:.6f}")
print(f"Standard Deviation: {std_dev_weight:.6f}")
print("="*40 + "\n")

print("\n" + "="*40)
print(f"--- STATS FOR MUON pT: 1 to 2 GeV ---")
print(f"Total Events:       {total_events:,.2f}")
print(f"Mean Weight:        {mean_weight:.6f}")
print(f"Standard Deviation: {std_dev_weight:.6f}")
print("="*40 + "\n")


plt.savefig("Weight_Event_Distribution.png", bbox_inches='tight', dpi=300)
print("[Plotter] Saved weight distribution plot successfully.")


# This is all for the muon_pt_sum_comparison_WITH_RATIO.png


 
# plot muon pt vectyor sum with alphaS weights
#print(results[dataset_name]["output"]["muon_pt_vector_sum_alphas_up"].get()) 

# mh.style.use("CMS")

# results = loaded_data
# dataset_name = list(results.keys())[0]



# axis_pt = loaded_data["alphaS_up"].axes[0]
# axis_vars = hist.axis.Integer(0, 3, underflow=False, overflow=False, name="vars")
# alphaS_hist_2d_real = hist.Hist(axis_pt, axis_vars, storage=hist.storage.Weight())

# # 2. Assign your original 1D data directly into the 2D slices (Data is 100% untouched)
# alphaS_hist_2d_real.view(flow=True).value[:, 0] = loaded_data["muon_hist_clean"].view(flow=True).value
# alphaS_hist_2d_real.view(flow=True).variance[:, 0] = loaded_data["muon_hist_clean"].view(flow=True).variance

# alphaS_hist_2d_real.view(flow=True).value[:, 1] = loaded_data["alphaS_up"].view(flow=True).value
# alphaS_hist_2d_real.view(flow=True).variance[:, 1] = loaded_data["alphaS_up"].view(flow=True).variance

# alphaS_hist_2d_real.view(flow=True).value[:, 2] = loaded_data["alphaS_down"].view(flow=True).value
# alphaS_hist_2d_real.view(flow=True).variance[:, 2] = loaded_data["alphaS_down"].view(flow=True).variance

# # 3. Create a standard, native python nested dictionary structure
# class FrameworkWrapper:
#     def __init__(self, obj): self.obj = obj
#     def get(self): return self.obj

# dataset_name = "wremnants_analysis"
# results = {
#     dataset_name: {
#         "output": {
#             "muon_pt_vector_sum_alphas_up": FrameworkWrapper(alphaS_hist_2d_real)
#         }
#     }
# }

# # plot muon pt vector sum and variations
# alphaS_hist_2d = results[dataset_name]["output"]["muon_pt_vector_sum_alphas_up"].get()
# muon_hist = alphaS_hist_2d[:, 0]
# alphaS_up = alphaS_hist_2d[:, 1]
# alphaS_down = alphaS_hist_2d[:, 2]    

# # plot weight hist
# weights_hist = loaded_data["weights_hist"]


# muon_hist_clean = hist.Hist(alphaS_up.axes[0], storage=hist.storage.Weight())
# muon_hist_clean.view().value = muon_hist.values()
# muon_hist_clean.view().variance = muon_hist.variances() 

# # Initialize the canvas with your preferred top-panel order (Nominal first)
# fig, ax_main, ax_comp = mh.comp.hists(
#    muon_hist_clean,
#    alphaS_up,
#   comparison='ratio',
#     xlabel=r"Sum of Muon $p_{T}$ [GeV]",
#     ylabel='Events',
#     h1_label='Nominal Data',
#     h2_label=r"$\alpha_S$ Up"
# )

# # 96. Plot the Down variation on the main upper panel
# mh.histplot(alphaS_down, ax=ax_main, label=r"$\alpha_S$ Down", linestyle=":", color="blue")


# # Clear the automatic, mismatched ratio markers from the bottom panel
# ax_comp.clear()
# ax_comp.axhline(1.0, color='black', linestyle='--', alpha=0.5) # Baseline at 1.0

# # Extract values and variances for the clean manual ratio math
# val_nom = muon_hist_clean.values()
# var_nom = muon_hist_clean.variances()

# # 1. Manually calculate the Up Ratio (alphaS_up / Nominal)
# ratio_up = hist.Hist(alphaS_up.axes[0], storage=hist.storage.Weight())
# val_up = alphaS_up.values()
# var_up = alphaS_up.variances()

# ratio_val_up = np.divide(val_up, val_nom, out=np.zeros_like(val_up), where=val_nom != 0)
# ratio_var_up = np.divide(
#     var_up + (ratio_val_up**2) * var_nom,
#     val_nom**2,
#     out=np.zeros_like(var_up),
#     where=val_nom != 0
# )
# ratio_up.view().value = ratio_val_up
# ratio_up.view().variance = ratio_var_up

# # 2. Manually calculate the Down Ratio (alphaS_down / Nominal)
# ratio_down = hist.Hist(alphaS_down.axes[0], storage=hist.storage.Weight())
# val_down = alphaS_down.values()
# var_down = alphaS_down.variances()

# ratio_val_down = np.divide(val_down, val_nom, out=np.zeros_like(val_down), where=val_nom != 0)
# ratio_var_down = np.divide(
#     var_down + (ratio_val_down**2) * var_nom,
#     val_nom**2,
#     out=np.zeros_like(var_down),
#     where=val_nom != 0
# )
# ratio_down.view().value = ratio_val_down
# ratio_down.view().variance = ratio_var_down

# # 3. Plot both matching ratio curves onto the cleared bottom panel
# mh.histplot(ratio_up, ax=ax_comp, linestyle="--", color="orange")
# mh.histplot(ratio_down, ax=ax_comp, linestyle=":", color="blue")

# # 4. Re-apply bottom panel labels wiped out by ax_comp.clear()
# ax_comp.set_xlabel(r"Sum of Muon $p_{T}$ [GeV]")
# ax_comp.set_ylabel(r"$\text{Variation} / \text{Nominal}$")
# ax_comp.set_ylim(0.95, 1.05) # Cleans up the ratio scale focal area

# #  LABELS, LEGEND, AND EXCLUSIVE SAVE 
# # txt_obj = mh.add_text("Default", loc='over left', ax=ax_main)
# mh.cms.label(ax=ax_main, data=False, label="Preliminary", rlabel="13.6 TeV") 

# # Updates the main legend to display your custom labels in clean order
# ax_main.legend(loc='upper right') 

# # Crop out the final empty overflow bin at 100 GeV to remove the drop-to-zero line
# ax_main.set_xlim(0, 98) 

# # Save the final singular double-panel plot canvas safely
 
# plt.savefig("muon_pt_sum_comparison_WITH_RATIO.png", bbox_inches='tight')
# print("[Plotter] Saved weight distribution plot.")



