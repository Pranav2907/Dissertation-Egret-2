import pandas as pd
from rdkit import Chem
from tqdm import tqdm
import random
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

def get_random_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol, doRandom=True)

def get_random_rxn(rxn):
    react, prod = rxn.split('>>')
    pipe_index = react.rfind('|')  # Find the last occurrence of '|' in case it's in the molecular structure
    if pipe_index != -1:
       special_value = react[pipe_index:]  # Get the '|25' and everything after (if there's more)
       react = react[:pipe_index]  # Get everything before '|25'
    else:
        special_value = ''
    react_list = react.split('.')
    react_list = [get_random_smiles(smi) for smi in react_list]
    random.shuffle(react_list)
    new_react = '.'.join(react_list)
    new_react_f = new_react + str(special_value)
    prod_list = prod.split('.')
    prod_list = [get_random_smiles(smi) for smi in prod_list]
    random.shuffle(prod_list)
    new_prod = '.'.join(prod_list)
    return f'{new_react_f}>>{new_prod}'

def construct_positive_rxn_smiles(rxn_smi):
    results = []
    for i in range(2):
        results.append(get_random_rxn(rxn_smi))
    posi_smi_1, posi_smi_2 = results
    return posi_smi_1, posi_smi_2

if __name__ == '__main__':
    import random
    random.seed(123)

    # Load the dataset and extract temperature information
    df = pd.read_csv('../../dataset/pretrain_data/pistachio_negative_sample_construction_t.csv', encoding='utf-8')
    
    # Assuming there is a 'temperature' column in your dataset
    rxn_id_ls = list(set(df['ID'].tolist()))
    rxn_id_ls.sort()
    val_id_ls = random.sample(rxn_id_ls, 1000)

    df_dict = {
        'ID': [],
        'org_smi': [],
        'positive_smi_1': [],
        'positive_smi_2': [],
    }

    for rxn_id, final_rxn_smi in tqdm(zip(df['ID'].tolist(), df['final_rxn_smi'].tolist()), total=len(df)):

        
        # Construct positive SMILES with the appended temperature
        posi_smi_1, posi_smi_2 = construct_positive_rxn_smiles(final_rxn_smi)
        
        df_dict['ID'].append(rxn_id)
        df_dict['org_smi'].append(final_rxn_smi)
        df_dict['positive_smi_1'].append(posi_smi_1)
        df_dict['positive_smi_2'].append(posi_smi_2)

    # Create the dataframe from the dictionary
    contrstive_df = pd.DataFrame.from_dict(df_dict)

    # Split into train and validation sets
    train_idx = []
    val_idx = []
    for i, reaction_id in enumerate(contrstive_df['ID'].tolist()):
        if reaction_id in val_id_ls:
            val_idx.append(i)
        else:
            train_idx.append(i)
    
    assert len(train_idx) + len(val_idx) == len(contrstive_df)
    
    pistachio_contrasive_learning_train_df = contrstive_df.loc[train_idx]
    pistachio_contrasive_learning_eval_df = contrstive_df.loc[val_idx]
    
    # Save the CSV files
    contrstive_df.to_csv('../../dataset/pretrain_data/aug_pis_pretraining_t.csv', encoding='utf-8', index=False)
    pistachio_contrasive_learning_train_df.to_csv('../../dataset/pretrain_data/aug_pis_train_t.csv', encoding='utf-8', index=False)
    pistachio_contrasive_learning_eval_df.to_csv('../../dataset/pretrain_data/aug_pis_eval_t.csv', encoding='utf-8', index=False)
