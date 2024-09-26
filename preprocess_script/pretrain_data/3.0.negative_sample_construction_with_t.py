import pandas as pd
from rdkit import Chem
from tqdm import tqdm
import random
import json
from rdkit import RDLogger

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

    pistachio_df = pd.read_csv('../../dataset/pretrain_data/pistachio_reaction_condition_identify.csv', encoding='utf-8')
    with open('../../dataset/pretrain_data/pistachio_condition_dict.json', encoding='utf-8') as f:
        condition_dict = json.load(f)

    catalyst_ls = condition_dict['catalyst']
    reagent_ls = condition_dict['reagent']
    solvent_ls = condition_dict['solvent']

    negative_data = {
        'ID': [],
        'rxn_smi': [],
        'catalyst': [],
        'reagent': [],
        'solvent': [],
        'temperature': []  # Added temperature handling
    }

    for (id, canonical_rxn, cat, rea, sol, temp) in zip(pistachio_df['_ID'].tolist(),
                                                        pistachio_df['canonical_rxn'].tolist(), 
                                                        pistachio_df['catalyst'].tolist(), 
                                                        pistachio_df['reagent'].tolist(), 
                                                        pistachio_df['solvent'].tolist(),
                                                        pistachio_df['temperature'].tolist()):  # Assuming temperature column exists

        negative_data['ID'].append(id)
        negative_data['rxn_smi'].append(canonical_rxn)

        # Catalyst handling
        if not pd.isna(cat):
            cat_ = cat.split(';')
            if len(cat_) != 1:
                cat_ = '.'.join(cat_)
                negative_data['catalyst'].append(cat_)
            else:
                negative_data['catalyst'].append(cat_[0])
        else:
            negative_data['catalyst'].append('')

        # Reagent handling
        if not pd.isna(rea):
            rea_ = rea.split(';')
            if len(rea_) != 1:
                rea_ = '.'.join(rea_)
                negative_data['reagent'].append(rea_)
            else:
                negative_data['reagent'].append(rea_[0])
        else:
            negative_data['reagent'].append('')

        # Solvent handling
        if not pd.isna(sol):
            sol_ = sol.split(';')
            if len(sol_) != 1:
                sol_ = '.'.join(sol_)
                negative_data['solvent'].append(sol_)
            else:
                negative_data['solvent'].append(sol_[0])
        else:
            negative_data['solvent'].append('')

        # Temperature handling
        if not pd.isna(temp):
            negative_data['temperature'].append(temp)
        else:
            negative_data['temperature'].append('')

        # Generate negative samples (2 samples for each reaction)
        for i in range(1):
            negative_data['ID'].append(id)
            negative_data['rxn_smi'].append(canonical_rxn)

            # Catalyst random sampling
            if not pd.isna(cat):
                cat_ls = cat.split(';')
                cat_len = len(cat_ls)
                random_cat = random.sample(catalyst_ls, cat_len)
                if cat_len != 1:
                    random_cat = '.'.join(random_cat)
                    negative_data['catalyst'].append(random_cat)
                else:
                    negative_data['catalyst'].append(random_cat[0])
            else:
                random_cat = random.sample(catalyst_ls, 1)
                negative_data['catalyst'].append(random_cat[0])

            # Reagent random sampling
            if not pd.isna(rea):
                rea_ls = rea.split(';')
                rea_len = len(rea_ls)
                random_rea = random.sample(reagent_ls, rea_len)
                if rea_len != 1:
                    random_rea = '.'.join(random_rea)
                    negative_data['reagent'].append(random_rea)
                else:
                    negative_data['reagent'].append(random_rea[0])
            else:
                random_rea = random.sample(reagent_ls, 1)
                negative_data['reagent'].append(random_rea[0])

            # Solvent random sampling
            if not pd.isna(sol):
                sol_ls = sol.split(';')
                sol_len = len(sol_ls)
                random_sol = random.sample(solvent_ls, sol_len)
                if sol_len != 1:
                    random_sol = '.'.join(random_sol)
                    negative_data['solvent'].append(random_sol)
                else:
                    negative_data['solvent'].append(random_sol[0])
            else:
                random_sol = random.sample(solvent_ls, 1)
                negative_data['solvent'].append(random_sol[0])

            # Temperature random sampling
            if not pd.isna(temp):
                varied_temp = float(temp)
                negative_data['temperature'].append(varied_temp)
            else:
                negative_data['temperature'].append('')

    assert len(negative_data['ID']) == len(negative_data['rxn_smi'])
    assert len(negative_data['ID']) == len(negative_data['catalyst'])
    assert len(negative_data['ID']) == len(negative_data['reagent'])
    assert len(negative_data['ID']) == len(negative_data['solvent'])
    assert len(negative_data['ID']) == len(negative_data['temperature'])

    condition_aug_df = pd.DataFrame.from_dict(negative_data)

    # Merging conditions into a single string (excluding temperature)
    condition_merge = []
    for cat_, rea_, sol_,temp_ in zip(condition_aug_df['catalyst'].tolist(),
                                condition_aug_df['reagent'].tolist(),
                                condition_aug_df['solvent'].tolist(),
                                condition_aug_df['temperature'].tolist()):
        cat_s = cat_.split('.')
        rea_s = rea_.split('.')
        sol_s = sol_.split('.')
        cat_s_ = [x for x in cat_s if x != '']
        rea_s_ = [x for x in rea_s if x != '']
        sol_s_s = [x for x in sol_s if x != '']
        condition_ls = cat_s_ + rea_s_ + sol_s_s
        conditions = '.'.join(condition_ls)
        complete_condition = conditions + '|' + str(temp_)
        condition_merge.append(complete_condition)

    condition_aug_df['merge_condition'] = condition_merge
    assert len(condition_merge) == len(condition_aug_df)

    # Generate final reaction SMILES including variations
    final_rxn_smi_ls = []
    for rxn_smi, conditon in tqdm(zip(condition_aug_df['rxn_smi'].tolist(),
                                      condition_aug_df['merge_condition'].tolist()),
                                  total=len(condition_aug_df)):
        fin_rxn_smiles = process_reaction(rxn_smi, conditon)
        final_rxn_smi_ls.append(fin_rxn_smiles)

    condition_aug_df['final_rxn_smi'] = final_rxn_smi_ls
    condition_aug_df = condition_aug_df.loc[condition_aug_df['final_rxn_smi'] != '']
    print(len(condition_aug_df))
    condition_aug_df.to_csv('../../dataset/pretrain_data/pistachio_negative_sample_construction_t.csv', encoding='utf-8', index=False)
