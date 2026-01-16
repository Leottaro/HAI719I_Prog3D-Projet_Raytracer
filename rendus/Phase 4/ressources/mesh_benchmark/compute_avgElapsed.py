import csv

# Open the CSV file
raw_mesh_benchmark = csv.reader(open('raw_mesh_benchmark.csv', mode='r'), delimiter=';')

TotalElapsedMs = {}

# Process the rows
# Skip the headers
next(raw_mesh_benchmark)

for row in raw_mesh_benchmark:
    mesh = row[0]
    maxLeafSize = int(row[1])
    elapsedMs = int(row[2])
    total_ms = 0
    n = 0

    if (mesh, maxLeafSize) in TotalElapsedMs:
        (total_ms, n) = TotalElapsedMs[(mesh, maxLeafSize)]

    TotalElapsedMs[(mesh, maxLeafSize)] = (total_ms + elapsedMs, n + 1)

avg_mesh_benchmark = csv.writer(open('avg_mesh_benchmark.csv', mode='w'), delimiter=',')
avg_mesh_benchmark.writerow(['Mesh', 'MaxLeafSize', 'AvgElapsed'])
for (mesh, maxLeafSize), (totalElapsedMs, n) in TotalElapsedMs.items():
    avgElapsedMs = totalElapsedMs / n
    avg_mesh_benchmark.writerow([mesh, maxLeafSize, round(avgElapsedMs, 2)])