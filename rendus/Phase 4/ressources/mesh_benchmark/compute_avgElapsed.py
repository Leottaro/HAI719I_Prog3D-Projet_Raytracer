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

    if mesh not in TotalElapsedMs:
        TotalElapsedMs[mesh] = {}

    if maxLeafSize not in TotalElapsedMs[mesh]:
        total_ms = 0
        for mesh2 in TotalElapsedMs.keys():
            TotalElapsedMs[mesh2][maxLeafSize] = (0, 0)

    (total_ms, n) = TotalElapsedMs[mesh][maxLeafSize]

    TotalElapsedMs[mesh][maxLeafSize] = (total_ms + elapsedMs, n + 1)

avg_mesh_benchmark = csv.writer(open('avg_mesh_benchmark.csv', mode='w'), delimiter=';')
avg_mesh_benchmark.writerow([''] + [maxLeafSize for maxLeafSize in TotalElapsedMs[next(iter(TotalElapsedMs))].keys()])
for mesh in TotalElapsedMs.keys():
    avg_mesh_benchmark.writerow([mesh] + [f'{round(TotalElapsedMs[mesh][maxLeafSize][0]/TotalElapsedMs[mesh][maxLeafSize][1], 2)}'.replace('.', ',') if TotalElapsedMs[mesh][maxLeafSize][1] != 0 else 0 for maxLeafSize in TotalElapsedMs[mesh]])