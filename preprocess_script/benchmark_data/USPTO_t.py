import pandas as pd
from rdkit import Chem
from tqdm import tqdm
import random
import json
from rdkit import RDLogger
from sklearn.model_selection import train_test_split

RDLogger.DisableLog('rdApp.*')

class NotCanonicalizableSmilesException(ValueError):
    pass

def canonicalize_smi(smi, remove_atom_mapping=False):
    """
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

def process_reaction(rxn, reagent):
    """
    Process and canonicalize reaction SMILES 
    """
    reactants, products = rxn.split(">>")
    try:
        precursors = [canonicalize_smi(r, False) for r in reactants.split(".")]
        if reagent:
            reagent_molecules = reagent.split(".")
            reagent_molecules = reagent_molecules[:-1]
            for r in reagent_molecules:
                if '|' in r:
                    # Split on the '|' and only process the first part for canonicalization
                    molecule, special_suffix = r.split('|', 1)
                    processed_molecule = canonicalize_smi(molecule, False) + '|' + special_suffix
                    precursors.append(processed_molecule)
                else:
                    precursors.append(canonicalize_smi(r, False))

        # Process products
        products = [canonicalize_smi(p, False) for p in products.split(".")]

    except NotCanonicalizableSmilesException:
        return ""

    # Reformat the reaction SMILES with all components
    return ".".join(precursors) + ">>" + ".".join(products)


if __name__ == "__main__":
    random.seed(123)

    uspto_df = pd.read_csv('../../dataset/source_dataset/USPTO_yield_t/USPTO_gram_t.csv', encoding='utf-8')
    uspto_df['catalyst'] = uspto_df['catalyst'].astype(str)
    uspto_df['reagent'] = uspto_df['reagent'].astype(str)
    uspto_df['solvent'] = uspto_df['solvent'].astype(str)
    uspto_df['yield'] = uspto_df['yield'].round(2)

    # List to store the combined conditions
    condition_merge = []

    # Iterate through each row
    for cat_, rea_, sol_, temp_ in zip(uspto_df['catalyst'].tolist(),
                                   uspto_df['reagent'].tolist(),
                                   uspto_df['solvent'].tolist(),
                                   uspto_df['temperature'].tolist()):
        # Split each string by '.' and filter out empty and 'nan' entries
        cat_s = [x for x in cat_.split('.') if x and x.lower() != 'nan']
        rea_s = [x for x in rea_.split('.') if x and x.lower() != 'nan']
        sol_s = [x for x in sol_.split('.') if x and x.lower() != 'nan']

        # Combine the lists and join with '.'
        condition_ls = cat_s + rea_s + sol_s
        conditions = '.'.join(condition_ls)
        complete_condition = conditions + '|' + str(temp_)
        condition_merge.append(complete_condition)


    uspto_df['merge_condition'] = condition_merge
    assert len(condition_merge) == len(uspto_df)

    final_rxn_smi_ls = []
    for rxn_smi, conditon in tqdm(zip(uspto_df['canonical_rxn'].tolist(),
                                      uspto_df['merge_condition'].tolist()),
                                  total=len(uspto_df)):
        fin_rxn_smiles = process_reaction(rxn_smi, conditon)
        final_rxn_smi_ls.append(fin_rxn_smiles)

    uspto_df['final_rxn_smi'] = final_rxn_smi_ls    
    uspto_df_f = uspto_df[['final_rxn_smi', 'yield']]


    train_df, temp_df = train_test_split(uspto_df_f, test_size=0.3, random_state=123)  # 70% for training, 30% for temp

    # Split the temporary dataset into validation and test
    val_df, test_df = train_test_split(temp_df, test_size=0.5, random_state=123)  # Split 30% of original into 15% val, 15% test

    # Save to CSV files
    train_df.to_csv('../../dataset/source_dataset/USPTO_yield_t/USPTO_gram_train_t.csv', index=False)
    val_df.to_csv('../../dataset/source_dataset/USPTO_yield_t/USPTO_gram_val_t.csv', index=False)
    test_df.to_csv('../../dataset/source_dataset/USPTO_yield_t/USPTO_gram_test_t.tsv', sep='\t', index=False)