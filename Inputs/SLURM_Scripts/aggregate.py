import glob
import pandas as pd
import re
import sys
import os

job_id = sys.argv[1]


pattern = f"/mnt/home/monkhayd/Simulations/Ejection_modeling/Mass_Fraction/Outputs/Ejection_Results/{job_id}-*_ejection_results.txt"

files = glob.glob(pattern)
results = []
for file in files:
    with open(file, "r") as f:
        result = f.read().strip()
        results.append(result)
    os.remove(file)

output_file = f"/mnt/home/monkhayd/Simulations/Ejection_modeling/Mass_Fraction/Outputs/Ejection_Results/{job_id}_aggregated_results.txt"
with open(output_file, "w") as f:
    for result in results:
        f.write(f"{result}\n")
