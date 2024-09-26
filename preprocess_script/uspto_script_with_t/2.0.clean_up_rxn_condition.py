import re
from collections import OrderedDict
from joblib import Parallel, delayed
import pandas as pd
import multiprocessing
import os
import argparse
import torch
from tqdm import tqdm
from utils import canonicalize_smiles, get_writer
from rxnmapper import RXNMapper

debug = False

def remap_and_reassign_condition_role(org_rxn, org_solvent, org_catalyst, org_reagent, temperature):
    if org_rxn.split('>') == 1:
        return None
    if '|' in org_rxn:
        rxn, frag = org_rxn.split(' ')
    else:
        rxn, frag = org_rxn, ''

    org_solvent, org_catalyst, org_reagent = [
        canonicalize_smiles(x) for x in [org_solvent, org_catalyst, org_reagent]
    ]
    try:
        results = rxn_mapper.get_attention_guided_atom_maps([rxn])[0]
    except Exception as e:
        print('\n'+rxn+'\n')
        print(e)
        return None

    remapped_rxn = results['mapped_rxn']
    confidence = results['confidence']

    new_precursors, new_products = remapped_rxn.split('>>')

    pt = re.compile(r':(\d+)]')
    new_react_list = []
    new_reag_list = []
    for precursor in new_precursors.split('.'):
        if re.findall(pt, precursor):
            new_react_list.append(precursor)
        else:
            new_reag_list.append(precursor)

    new_reactants = '.'.join(new_react_list)
    react_maps = sorted(re.findall(pt, new_reactants))
    prod_maps = sorted(re.findall(pt, new_products))
    if react_maps != prod_maps:
        return None

    new_reagent_list = []
    c_list = org_catalyst.split('.')
    s_list = org_solvent.split('.')
    r_list = org_reagent.split('.')
    for r in new_reag_list:
        if (r not in c_list + s_list) and (r not in r_list):
            new_reagent_list.append(r)
    new_reagent_list += [x for x in r_list if x != '']

    catalyst = org_catalyst
    solvent = org_solvent
    reagent = '.'.join(new_reagent_list)
    can_react = canonicalize_smiles(new_reactants, clear_map=True)
    can_prod = canonicalize_smiles(new_products, clear_map=True)
    can_rxn = '{}>>{}'.format(can_react, can_prod)
    
    results = OrderedDict()
    results['remapped_rxn'] = remapped_rxn
    results['frag'] = frag
    results['confidence'] = confidence
    results['can_rxn'] = can_rxn
    results['catalyst'] = catalyst
    results['solvent'] = solvent
    results['reagent'] = reagent
    results['temperature'] = temperature  # include temperature in the results

    return results

def run_tasks(task):
    idx, rxn, solvent, catalyst, reagent, source, temperature = task  # include temperature in the task tuple
    if pd.isna(solvent):
        solvent = ''
    if pd.isna(catalyst):
        catalyst = ''
    if pd.isna(reagent):
        reagent = ''
    results = remap_and_reassign_condition_role(
        rxn, solvent, catalyst, reagent, temperature
    )

    return idx, results, source

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--split_group', type=int, default=4)
    parser.add_argument('--group', type=int, default=1)
    args = parser.parse_args()

    assert args.group <= args.split_group - 1

    print('Debug:', debug)
    print('Split group:', args.split_group)
    print('Group number:', args.group)
    print('GPU index:', args.gpu)

    device = torch.device('cuda:{}'.format(args.gpu) if args.gpu >= 0 else 'cpu')
    rxn_mapper = RXNMapper()
    source_data_path = '../../dataset/source_dataset/'
    rxn_condition_fname = 'uspto_cleaned_dataset_t.csv'
    new_database_fpath = os.path.join(
        source_data_path, 'uspto_rxn_condition_remapped_and_reassign_condition_role_group_{}.csv'.format(args.group))

    database = pd.read_csv(os.path.join(source_data_path, rxn_condition_fname))
    print('All data number:', len(database))

    group_size = len(database) // args.split_group
    database = database.iloc[args.group * group_size:] if args.group == args.split_group - 1 else database.iloc[args.group * group_size:(args.group + 1) * group_size]
    print('Calculate index {} to {}'.format(database.index.min(), database.index.max()))

    header = [
        'remapped_rxn',
        'fragment',
        'confidence',
        'canonical_rxn',
        'catalyst',
        'solvent',
        'reagent',
        'temperature',  # Added temperature
        'source',
        'org_rxn'
    ]
    fout, writer = get_writer(new_database_fpath, header=header)
    
    for row in tqdm(database.itertuples(), total=len(database)):
        task = (row.Index, row.rxn_smiles, row.solvent, row.catalyst, row.reagent, row.source, row.temperature)  # Include temperature
        try:
            idx, results, source = run_tasks(task)
            if results:
                results['source'] = source
                results['org_rxn'] = row.rxn_smiles  # Correctly refer to the original reaction SMILES
                row = list(results.values())
                assert len(row) == len(header)
                writer.writerow(row)
                fout.flush()
        except Exception as e:
            print(e)
            continue  # Changed from pass to continue to better signify the intention

    fout.close()
    new_database = pd.read_csv(new_database_fpath)
    reset_header = [  # Adjusting the order if necessary or keeping as is
        'source',
        'org_rxn',
        'fragment',
        'remapped_rxn',
        'confidence',
        'canonical_rxn',
        'catalyst',
        'solvent',
        'reagent',
        'temperature'  # Ensure temperature is included in the final header
    ]
    new_database = new_database[reset_header]
    new_database.to_csv(new_database_fpath, index=False)
    print('Done!')
