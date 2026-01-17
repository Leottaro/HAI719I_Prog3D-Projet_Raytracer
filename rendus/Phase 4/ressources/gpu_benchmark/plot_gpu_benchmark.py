import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file with proper decimal handling
file_path = "avg_gpu_benchmark.csv"
data = pd.read_csv(file_path, sep=';', decimal=',')

# Extract columns
nb_objects = data['NbObjects']
cpu_elapsed = data['CPUElapsedMs']
gpu_elapsed = data['GPUElapsedMs']

# Plot NbObjects vs CPUElapsedMs and GPUElapsedMs
plt.figure(figsize=(10, 6), tight_layout=True)
plt.plot(nb_objects, cpu_elapsed, label="CPU Elapsed Time (ms)", marker='o')
plt.plot(nb_objects, gpu_elapsed, label="GPU Elapsed Time (ms)", marker='s')

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Number of Objects")
plt.ylabel("Elapsed Time (ms)")
plt.title("CPU vs GPU Elapsed Time")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.xlim(left=1)  # Set x-axis to start at 1
plt.ylim(bottom=1)  # Set y-axis to start at 1
plt.savefig("cpu_vs_gpu_elapsed_time.png", bbox_inches="tight", dpi=300)
plt.close()

print("Plot saved as cpu_vs_gpu_elapsed_time.png")