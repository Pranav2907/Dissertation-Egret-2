import pandas as pd
import json

if __name__ == "__main__":
    df = pd.read_csv('../../../dataset/pretrain_data/pistachio_reaction_condition_identify.csv', encoding='utf-8')
    condition_dict = {
        'catalyst': [],
        'reagent': [],
        'solvent': [],
        'temperature': [],  # Added to handle temperature data
    }
    # Loop through each row's relevant columns
    for (cat, rea, sol, temp) in zip(df['catalyst'].tolist(), 
                                     df['reagent'].tolist(), 
                                     df['solvent'].tolist(), 
                                     df['temperature'].tolist()):  # Adjusted to include temperature
        if not pd.isna(cat):
            cat = cat.split(';')
            condition_dict['catalyst'].extend(cat)
        if not pd.isna(rea):
            rea = rea.split(';')
            condition_dict['reagent'].extend(rea)
        if not pd.isna(sol):
            sol = sol.split(';')
            condition_dict['solvent'].extend(sol)
        if not pd.isna(temp):  # Directly append the temperature if it's not NaN
            condition_dict['temperature'].append(temp)

    # Create lists from sets to remove duplicates and sort them
    catalyst_ls = list(set(condition_dict['catalyst']))
    catalyst_ls.sort()
    reagent_ls = list(set(condition_dict['reagent']))
    reagent_ls.sort()
    solvent_ls = list(set(condition_dict['solvent']))
    solvent_ls.sort()
    temperature_ls = list(set(condition_dict['temperature']))
    temperature_ls.sort()

    # Final dictionary to be stored in JSON
    condition_dict_final = {
        'catalyst': catalyst_ls,
        'reagent': reagent_ls,
        'solvent': solvent_ls,
        'temperature': temperature_ls,
    }

    # Write the dictionary to a JSON file
    with open('../../../dataset/pretrain_data/pistachio_condition_dict.json', 'w') as f:
        json.dump(condition_dict_final, f)
