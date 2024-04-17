import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# Read the Parquet file
df = pd.read_parquet('brick/1976_Sep2016_USPTOgrants_smiles.parquet')

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
# ====================================================================


# Frequency of patents over the years
patents_per_year = df.groupby('Year')['PatentNumber'].nunique()
plt.figure(figsize=(10, 6))
sns.lineplot(x=patents_per_year.index, y=patents_per_year.values)
plt.title('Number of Unique Patents per Year')
plt.xlabel('Year')
plt.ylabel('Number of Patents')
# save the plot
plt.savefig("unique_patents_per_year.png", format='png', dpi=300)

plt.show()
# =====================================================
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

def visualize_molecules(molecules):
    img = Draw.MolsToGridImage(molecules, subImgSize=(200, 200))
    return img

# Assuming interpret_reaction_smiles function as previously defined

# Example Reaction SMILES
reaction_smiles = "CCO.O>[Pt]>OCC"
reactants, reagents, products = interpret_reaction_smiles(reaction_smiles)

# Visualize each component
print("Reactants:")
img = visualize_molecules(reactants)
img.show()  # This will open the image in the default image viewer

if reagents:
    print("Reagents:")
    img = visualize_molecules(reagents)
    img.show()

print("Products:")
img = visualize_molecules(products)
img.show()
# ================================================
from rdkit import Chem
from rdkit.Chem import Descriptors, rdChemReactions



def compute_properties(mol):
    """ Compute molecular properties for a given molecule """
    # Sanitize the molecule to ensure proper calculation of properties
    Chem.SanitizeMol(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    return mw, logp, rotatable_bonds

def process_reaction_smiles(reaction_smiles):
    """ Process reaction SMILES and compute properties """
    reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles)
    reactants, products = reaction.GetReactants(), reaction.GetProducts()

    reactant_properties = [compute_properties(mol) for mol in reactants]
    product_properties = [compute_properties(mol) for mol in products]

    return reactant_properties, product_properties

# Example Reaction SMILES
reaction_smiles = "CCO.O>[Pt]>OCC"  # Replace with your Reaction SMILES
reactant_props, product_props = process_reaction_smiles(reaction_smiles)

# Output the properties
print("Reactant Properties (MW, LogP, Rotatable Bonds):", reactant_props)
print("Product Properties (MW, LogP, Rotatable Bonds):", product_props)


# =============================================================
from rdkit import Chem
from rdkit.Chem import rdChemReactions, rdMolDescriptors

def sanitize_molecules(mols):
    for mol in mols:
        if mol is not None:
            Chem.SanitizeMol(mol)  # Sanitize each molecule

def classify_reaction(reaction_smiles):
    reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles)
    reactants, products = reaction.GetReactants(), reaction.GetProducts()

    # Sanitize molecules
    sanitize_molecules(reactants)
    sanitize_molecules(products)

    # Counting atoms in reactants and products
    reactant_atoms = sum([rdMolDescriptors.CalcNumAtoms(mol) for mol in reactants])
    product_atoms = sum([rdMolDescriptors.CalcNumAtoms(mol) for mol in products])

    if reactant_atoms < product_atoms:
        return "Addition Reaction"
    elif reactant_atoms > product_atoms:
        return "Elimination Reaction"
    else:
        return "Substitution or Other Reaction Type"

# Example usage
reaction_smiles = "CCO.O>[Pt]>OCC"  # Example reaction
reaction_type = classify_reaction(reaction_smiles)
print("The reaction is classified as:", reaction_type)




