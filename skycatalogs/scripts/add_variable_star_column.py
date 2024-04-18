import os
import glob
import numpy as np
import pandas as pd

files = sorted(glob.glob("skyCatalogs/pointsource_*.parquet"))
dest_dir = "skyCatalogs_variable_stars"
os.makedirs(dest_dir, exist_ok=True)

for i, item in enumerate(files):
    print(i, item, len(files))
    df = pd.read_parquet(item)

    # Get random indexes of 10% of rows to set the varible flag as true.
    nrows = len(df)
    index = np.arange(nrows)
    np.random.shuffle(index)
    index = index[:nrows // 10]
    index.sort()
    nvars = len(index)

    variable = np.array([False for _ in range(nrows)])
    variable[index] = True
    df['is_variable'] = variable

    # Draw a random period from a log-uniform distribution from 0.1 to 10 days.
    periods = np.zeros(nrows)
    periods[index] = 10.**np.random.uniform(-1, 1, size=nvars)
    df['period'] = periods

    # Draw a random amplitude in magnitudes from range 0.1 to 1.
    amplitudes = np.zeros(nrows)
    amplitudes[index] = np.random.uniform(0.1, 1, size=nvars)
    df['mag_amplitude'] = amplitudes

    # Random phase at MJD=60400.
    phases = np.zeros(nrows)
    phases[index] = np.random.uniform(0, 2.*np.pi, size=nvars)
    df['phase'] = phases

    # Write out the revised data frame to a new parquet file.
    dest = os.path.join(dest_dir, os.path.basename(item))
    df.to_parquet(dest)
