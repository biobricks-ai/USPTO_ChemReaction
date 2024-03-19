import os
import pandas as pd
from tqdm import tqdm

def read_rsmi_file(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t', dtype={'ParagraphNum': str})
        return df
    except pd.errors.ParserError as e:
        print(f"Error reading file {file_path}: {e}")
        return None

def find_rsmi_files(start_path):
    for dirpath, dirnames, filenames in os.walk(start_path):
        for filename in filenames:
            if filename.endswith('.rsmi'):
                yield os.path.join(dirpath, filename)

def build_dataframe(rsmi_files):
    frames = []
    for file in tqdm(rsmi_files, desc="Reading files", unit="file"):
        df = read_rsmi_file(file)
        if df is not None:
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

def main():
    start_path = 'raw'  # Update this path as needed
    rsmi_files = list(find_rsmi_files(start_path))
    df = build_dataframe(rsmi_files)
    
    if not df.empty:
        df.to_parquet('output.parquet')
        print("Data saved to output.parquet")
    else:
        print("No data to save")

if __name__ == '__main__':
    main()



# Read the Parquet file
df = pd.read_parquet('output.parquet')

# Display the first few rows of the DataFrame
print(df.head())
# ==================================================
from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions

def visualize_reaction(reaction_smiles):
    # Parse the reaction SMILES
    reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles)
    
    # Generate an image of the reaction
    img = Draw.ReactionToImage(reaction)
    img.show()

# Example Reaction SMILES
reaction_smiles = 'CCO.O>[Pt]>OCC'  # Replace with your Reaction SMILES
visualize_reaction(reaction_smiles)

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Frequency of patents over the years
patents_per_year = df.groupby('Year')['PatentNumber'].nunique()
plt.figure(figsize=(10, 6))
sns.lineplot(x=patents_per_year.index, y=patents_per_year.values)
plt.title('Number of Unique Patents per Year')
plt.xlabel('Year')
plt.ylabel('Number of Patents')
plt.show()