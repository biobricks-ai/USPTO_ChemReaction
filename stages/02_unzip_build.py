import os, shutil, mimetypes, subprocess, pandas as pd, tqdm, glob, pathlib, glob
from py7zr import unpack_7zarchive

#%% EXTRACT DATA TO TEMPORARY DIRECTORY ===================================================
# Register 7z file format if not already registered
if '7zip' not in [x[0] for x in shutil.get_unpack_formats()]:
    shutil.register_unpack_format('7zip', ['.7z'], unpack_7zarchive)

if os.path.exists('raw'): shutil.rmtree('raw')
os.makedirs('raw', exist_ok=True)
shutil.unpack_archive('download/chemical_reactions.zip', 'raw')

for zip in glob.glob('raw/**/*.zip', recursive=True):
    shutil.unpack_archive(zip, os.path.dirname(zip))

for z7 in glob.glob('raw/**/*.7z', recursive=True):
    shutil.unpack_archive(z7, os.path.dirname(z7), '7zip')

#%% TRANSFORM DATA INTO BRICK ==============================================================
brick = pathlib.Path('brick') 
brick.mkdir(exist_ok=True)
rsmi = glob.glob('raw/**/*.rsmi', recursive=True)
for path in tqdm.tqdm(rsmi, desc="Processing files", unit="file"):
    df = pd.read_csv(path, sep='\t', dtype={'ParagraphNum': str})
    df.to_parquet(brick / os.path.basename(path).replace('.rsmi', '.parquet'))

#%% TESTS =====================================================================================
# Get the number of lines of first rsmi file
numlines = sum(1 for line in open(rsmi[-1]))
df = pd.read_parquet(brick / os.path.basename(rsmi[-1]).replace('.rsmi', '.parquet'))
assert (numlines - 1) == len(df), f"Number of lines in {rsmi[-1]} does not match the number of rows in the DataFrame"

#%% CLEANUP ===================================================================================
shutil.rmtree('raw')
# %%
