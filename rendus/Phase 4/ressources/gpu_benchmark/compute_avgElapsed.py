import csv

# Open the CSV file
raw_gpu_benchmark = csv.reader(open('raw_gpu_benchmark.csv', mode='r'), delimiter=';')

TotalElapsedMs = {}

# Process the rows
# Skip the headers
next(raw_gpu_benchmark)

for row in raw_gpu_benchmark:
    nb_objects = int(row[0])
    CPUElapsedMs = int(row[1])
    GPUElapsedMs = int(row[2])
    total_cpu_ms = 0
    total_gpu_ms = 0
    n = 0

    if nb_objects in TotalElapsedMs:
        (total_cpu_ms, total_gpu_ms, n) = TotalElapsedMs[nb_objects]

    TotalElapsedMs[nb_objects] = (total_cpu_ms + CPUElapsedMs, total_gpu_ms + GPUElapsedMs, n + 1)

avg_gpu_benchmark = csv.writer(open('avg_gpu_benchmark.csv', mode='w'), delimiter=';')
avg_gpu_benchmark.writerow(['NbObjects', 'CPUElapsedMs', 'GPUElapsedMs'])
for nb_objects, (total_cpu_ms, total_gpu_ms, n) in TotalElapsedMs.items():
    avgCPUElapsedMs = total_cpu_ms / n
    avgGPUElapsedMs = total_gpu_ms / n
    avg_gpu_benchmark.writerow([nb_objects, f"{avgCPUElapsedMs:.2f}".replace('.', ','), f"{avgGPUElapsedMs:.2f}".replace('.', ',')])