import pandas as pd
import os
from tqdm import tqdm
from rdkit import Chem
import re
import warnings

class NotCanonicalizableSmilesException(ValueError):
    pass

def canonicalize_smi(smi, remove_atom_mapping=False):
    r"""
    Canonicalize SMILES
    """
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        raise NotCanonicalizableSmilesException("Molecule not canonicalizable")
    if remove_atom_mapping:
        for atom in mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                atom.ClearProp("molAtomMapNumber")
    return Chem.MolToSmiles(mol)



def process_reaction(row):
    try:
        rxn = row['canonical_rxn']
        reactants, products = rxn.split(">>")
        reactants = canonicalize_smi(reactants)
        products = canonicalize_smi(products)
        catalyst = canonicalize_smi(row['catalyst']) if pd.notna(row['catalyst']) else ''
        solvent = canonicalize_smi(row['solvent']) if pd.notna(row['solvent']) else ''
        reagent = str(row['reagent']) if pd.notna(row['reagent']) else ''
        temp = str(row['temperature']) if pd.notna(row['temperature']) else ''

        # Concatenating the reaction components
        reaction_components = [comp for comp in [reactants, catalyst, solvent, reagent] if comp]
        reactant_side = '.'.join(reaction_components)
        if temp:
            reactant_side += f"|{temp}"

        return f"{reactant_side}>>{products}"  # Assuming the product is reactants processed or another column should be specified
    except NotCanonicalizableSmilesException:
        return ""

if __name__ == '__main__':
    df = pd.read_csv('../../dataset/pretrain_data/pistachio_reaction_condition_identify.csv')
    df['final_reaction'] = df.apply(process_reaction, axis=1)
    df = df[df['final_reaction'] != '']
    
    # Shuffle the dataframe
    df_shuffled = df.sample(frac=1, random_state=42)  # Random state for reproducibility
    final_reactions = df_shuffled['final_reaction'].tolist()

    # Splitting the data into training and evaluation sets (80/20 split)
    split_index = int(0.8 * len(final_reactions))
    train_reactions = final_reactions[:split_index]
    eval_reactions = final_reactions[split_index:]

    # Write training reactions to a text file
    with open('../../dataset/pretrain_data/final_reactions_train.txt', 'w') as f_train:
        for reaction in tqdm(train_reactions, desc="Writing training data"):
            f_train.write(reaction + '\n')

    # Write evaluation reactions to a text file
    with open('../../dataset/pretrain_data/final_reactions_eval.txt', 'w') as f_eval:
        for reaction in tqdm(eval_reactions, desc="Writing evaluation data"):
            f_eval.write(reaction + '\n')
