import numpy as np
import pandas as pd
import os
import math
import glob

grid_dim=[1,1]

m_over_m_earth_values = np.logspace(2, 3, 3*grid_dim[0])
a_p_values = np.logspace(0,1.6,3*grid_dim[1])

rows = []

for m_p in m_over_m_earth_values:
    for a_p in a_p_values:
        rows.append({
            "m/m_earth": m_p,
            "a_p": a_p,
        })

param_df = pd.DataFrame(rows)

param_df.to_csv("Parameter_CSVs/m_a_parameter_values.csv", index=False)