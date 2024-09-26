import pandas as pd


def selcet_n_sample(dataframe):
    select_dataframe = dataframe.iloc[::2]
    return select_dataframe

if __name__ == "__main__":
    all_trn_df = pd.read_csv('../../dataset/pretrain_data/aug_pis_train_t.csv', encoding='utf-8')
    all_eval_df = pd.read_csv('../../dataset/pretrain_data/aug_pis_eval_t.csv', encoding='utf-8')
    select_trn_df = selcet_n_sample(all_trn_df)
    select_eval_df = selcet_n_sample(all_eval_df)
    select_trn_df.to_csv('../../dataset/pretrain_data/aug_pis_train_dataset_t.csv', encoding='utf-8', index=False)
    select_eval_df.to_csv('../../dataset/pretrain_data/aug_pis_eval_dataset_t.csv', encoding='utf-8', index=False)        
    print('Done!')


