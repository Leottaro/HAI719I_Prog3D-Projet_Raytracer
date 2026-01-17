import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV file with proper decimal handling
file_path = "avg_mesh_benchmark.csv"
data = pd.read_csv(file_path, sep=';', index_col=0, decimal=',')

# Convert column names to integers (MaxLeafSize)
data.columns = [int(col) if col.isdigit() else col for col in data.columns]

# Ensure all data is numeric
data = data.apply(pd.to_numeric, errors='coerce')

print("Data loaded from CSV:")
print(data)

max_leaf_sizes = [col for col in data.columns if isinstance(col, int)]

# Plot 1: AvgElapsedMs vs MaxLeafSize
plt.figure(figsize=(10, 6), tight_layout=True)
adjusted_max_leaf_sizes = [x if x != 0 else 0.1 for x in max_leaf_sizes]
for mesh in data.index:
    # Filter out zero values to avoid vertical lines
    valid_indices = data.loc[mesh, max_leaf_sizes] > 0
    plt.plot(
        np.array(adjusted_max_leaf_sizes)[valid_indices],
        data.loc[mesh, max_leaf_sizes][valid_indices],
        label=mesh
    )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("MaxLeafSize")
plt.ylabel("AvgElapsedMs (log scale)")
xticks = [0.1, 1] + [x for x in max_leaf_sizes if x > 1][::2]
plt.gca().set_xticks(xticks)
plt.gca().set_xticklabels(["0", "1"] + [str(x) for x in max_leaf_sizes if x > 1][::2])
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.savefig("avg_elapsed_vs_maxleafsize.png", bbox_inches="tight", dpi=600)
plt.close()

# Plot 2: Relative Speedup vs MaxLeafSize
plt.figure(figsize=(10, 6), tight_layout=True)
adjusted_max_leaf_sizes = [x if x != 0 else 0.1 for x in max_leaf_sizes]
reference = data.loc["NONE", 0]
for mesh in data.index:
    if mesh == "NONE":
        continue
    speedup = reference / data.loc[mesh, max_leaf_sizes]
    plt.plot(adjusted_max_leaf_sizes, speedup, label=mesh)

plt.xscale("log")
plt.yscale("linear")
plt.xlabel("MaxLeafSize")
plt.ylabel("Speedup (relative to NONE)")
xticks = [0.1, 1] + [x for x in max_leaf_sizes if x > 1][::2]
plt.gca().set_xticks(xticks)
plt.gca().set_xticklabels(["0", "1"] + [str(x) for x in max_leaf_sizes if x > 1][::2])
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.savefig("relative_speedup_vs_maxleafsize.png", bbox_inches="tight", dpi=600)
plt.close()