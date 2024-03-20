import os
import pandas as pd
from tqdm import tqdm
import glob
import pathlib

brick = pathlib.Path('brick')
brick.mkdir(exist_ok=True)
rsmi = glob.glob('raw/**/*.rsmi', recursive=True)

for path in tqdm(rsmi, desc="Processing files", unit="file"):
    df = pd.read_csv(path, sep='\t', dtype={'ParagraphNum': str})
    out = brick / os.path.basename(path).replace('.rsmi', '.parquet')
    df.to_parquet(out)
    
# Get the number of lines of first rsmi file
numlines = sum(1 for line in open(rsmi[-1]))
assert (numlines - 1) == len(df), f"Number of lines in {rsmi[-1]} does not match the number of rows in the DataFrame"

    

    