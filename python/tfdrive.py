import numpy as _np
import pandas as _pd
import os as _os
import pyarrow as _pyarrow
import scipy.stats as _st
from collections import Counter as _Counter
from sklearn.linear_model import LogisticRegression as _LogisticRegression
from sklearn.ensemble import RandomForestClassifier as _RandomForestClassifier
from joblib import load as _load


def tfpred(gene_ids, species = 'hsa'):
    '''
    Predicts transcription factors driving differential gene expression.
    Parameters:
    ----------
    gene_ids : array_like, int
        1-D array of Entrez Gene IDs of differentially expressed genes.

    species : str
        Three-letter species ID. Currently, function works only with human genes. 
        Therefore, 'hsa' is the default value and should not be changed.

    Returns : DataFrame
        Returns a list of transcription factors ranked according to their probability 
        scores.
    '''
    # Check the input
    if type(gene_ids) is str:
        raise ValueError("'gene_ids' should be a list of IDs")

    try:
        gene_ids = list(gene_ids)
    except TypeError:
        print("'gene_ids' argument is not iterable.")
        raise

    if len(_np.array(gene_ids).shape) > 1:
        raise ValueError('The list of gene IDs should be flat')

    if len(gene_ids)<10:
        raise ValueError('The list of genes is too short')

    for i in range(len(gene_ids)):
        if type(gene_ids[i]) is not int:
            try:
                gene_ids[i] = int(gene_ids[i])
            except ValueError:
                print("An element of 'gene_ids' cannot be converted to an integer.")
                raise

    if species!='hsa':
        raise ValueError("Only 'hsa' option is currently supported.")
    
    # Load data
    file_list = ['tfdrive_data/pathways_hsa_full_hg38.parquet',
                 'tfdrive_data/go_terms_hsa_processes.parquet',
                 'tfdrive_data/tf_names_hsa.parquet',
                 'tfdrive_data/tf_families_hsa.parquet',
                 'tfdrive_data/tf_matrix_transfac.parquet',
                 'tfdrive_data/tf_matrix_chea.parquet',
                 'tfdrive_data/rf_mod.joblib',
                 'tfdrive_data/log_mod.joblib']
                 
    if not all([_os.path.isfile(f) for f in file_list]):
        raise OSError("Some required data files are missing.")

    # Table of gene ID to pathway associations
    pathways_db = _pd.read_parquet('tfdrive_data/pathways_hsa_full_hg38.parquet')
    pathways_db.drop_duplicates(['gene_id','pathway'],keep='first', inplace=True)

    gene_ids = list(set(gene_ids))
    if len([g for g in gene_ids if g in pathways_db['gene_id'].values]) < 10:
        raise ValueError("Too few or none of provided gene IDs were found in the database.")

    # Table of gene ID to GO terms (from 'processes' cathegory) associations
    go_terms_db = _pd.read_parquet('tfdrive_data/go_terms_hsa_processes.parquet')
    go_terms_db.drop_duplicates(['gene_id','GO_term'],keep='first', inplace=True)

    # Table of TF to TF ID associations
    tf_names_db = _pd.read_parquet('tfdrive_data/tf_names_hsa.parquet')

    # Load TF classification table
    tf_families_df = _pd.read_parquet('tfdrive_data/tf_families_hsa.parquet')

    # Tables of TF-gene interactions
    transfac_df = _pd.read_parquet('tfdrive_data/tf_matrix_transfac.parquet')
    chea_df = _pd.read_parquet('tfdrive_data/tf_matrix_chea.parquet')

    # Arbitrary distribution for calculating 'pathway importance'
    dist = _st.beta(2,2)

    # Prediction models
    rf_model = _load('tfdrive_data/rf_mod.joblib')
    log_model = _load('tfdrive_data/log_mod.joblib') 

    full_data_df = _pd.DataFrame(columns=['TF', 'TF_id', 'p_share', 
                                         'all_p', 'p_score', 'go_share', 
                                         'all_go', 'go_score'])

    for tf_id in tf_names_db.index:
        tf_name = tf_names_db.at[tf_id, 'gene']  

        # Get unique pathways and terms related to the TF and their counts
        tf_ptws = list(pathways_db[pathways_db['gene_id'].isin([tf_id])]['pathway'])
        tf_terms = list(go_terms_db[go_terms_db['gene_id'].isin([tf_id])]['GO_term'])
        tf_pw_count = _Counter(tf_ptws)
        tf_term_count = _Counter(tf_terms)
        all_p = len(tf_pw_count)
        all_go = len(tf_term_count)

        # Get unique pathways and terms related to differentially expressed genes
        # and their counts
        case_ptws = list(pathways_db[pathways_db['gene_id'].isin(gene_ids)]['pathway'])
        case_terms = list(go_terms_db[go_terms_db['gene_id'].isin(gene_ids)]['GO_term'])
        case_ptw_count = _Counter(case_ptws)
        case_term_count = _Counter(case_terms)

        # Calculate relative pathway/GO term importance
        max_pw_count = max(_Counter(pathways_db['pathway'].values).values())
        max_go_count = max(_Counter(go_terms_db['GO_term'].values).values())

        pw_df = _pd.DataFrame.from_dict(case_ptw_count, orient='index',columns=['count'])
        pw_df['importance'] = pw_df['count'].apply(lambda x: x/max_pw_count)
        pw_df['importance']  = dist.pdf(pw_df['importance'])

        go_df = _pd.DataFrame.from_dict(case_term_count, orient='index',columns=['count'])
        go_df['importance'] = go_df['count'].apply(lambda x: x/max_go_count)
        go_df['importance']  = dist.pdf(go_df['importance'])

        # Calculate relative pathway/term overlap and scores based on 'importance'
        common_p = 0
        common_go = 0
        p_score = 0.0
        go_score = 0.0
        coef = 100
        
        if all_p > 0:
            for i in tf_pw_count:
                if i in case_ptw_count:
                    common_p += 1
                    p_score += coef * pw_df.at[i, 'importance'] * pw_df.at[i, 'count'] / len(gene_ids)

        if all_go > 0:
            for t in tf_term_count:
                if t in case_term_count:
                    common_go += 1
                    go_score += coef * go_df.at[t, 'importance'] * go_df.at[t, 'count'] / len(gene_ids)
  
        p_share = common_p/all_p if all_p > 0 else 0.0
        go_share = common_go/all_go if all_go > 0 else 0.0
        
        # Add to DataFrame
        full_data_df.loc[len(full_data_df)] = [tf_name, tf_id, p_share, all_p,  
                                               p_score, go_share, all_go, go_score]

    #Drop TFs which have less than 3 pathways associated with them
    full_data_df = full_data_df[full_data_df['all_p']>2].copy(deep=True)
    full_data_df = full_data_df[full_data_df['all_go']>2].copy(deep=True)

    # Calculate relative interaction frequencies for each family for each group
    tf_group_scores_df = _pd.DataFrame(columns=['family', 'transfac', 'chea'])
    gene_ids_str = [str(g) for g in gene_ids]
    for f in tf_families_df.index:
        tf_ids = tf_families_df.at[f, 'ids']
        # Transfac
        t_value = 0.0
        t_tf_ids = [t for t in tf_ids if t in transfac_df.columns]
        t_g_ids = [h for h in gene_ids_str if h in transfac_df.index]    
        if (len(t_tf_ids) > 0) & (len(t_g_ids) > 0):
            t_value = (transfac_df.loc[t_g_ids, t_tf_ids].sum().sum() / len(t_g_ids)) / \
                      (transfac_df[t_tf_ids].sum().sum() / len(transfac_df))
        # ChEA
        c_value = 0.0
        c_tf_ids = [c for c in tf_ids if c in chea_df.columns]
        c_g_ids = [h for h in gene_ids_str if h in chea_df.index]    
        if (len(c_tf_ids) > 0) & (len(c_g_ids) > 0):
            c_value = (chea_df.loc[c_g_ids, c_tf_ids].sum().sum() / len(c_g_ids)) / \
                      (chea_df[c_tf_ids].sum().sum() / len(chea_df))
        tf_group_scores_df.loc[len(tf_group_scores_df)] = [f, t_value, c_value]
            
    # Create TF to TF family association table
    tf_family_pairs = _pd.DataFrame(columns=['family'])
    for i in tf_families_df.index:
        for tf in tf_families_df.at[i, 'ids']:
            tf_family_pairs.loc[tf] = [i]
            
    # Assign values to full_data_df
    full_data_df['transfac'] = 0.0
    full_data_df['chea'] = 0.0
    for i in full_data_df.index:
        tf_id = str(full_data_df.at[i, 'TF_id'])
        if tf_id in tf_family_pairs.index:
            family = tf_family_pairs.at[tf_id, 'family']
            full_data_df.at[i, 'transfac'] = tf_group_scores_df[tf_group_scores_df['family'] == family]['transfac'].values[0]
            full_data_df.at[i, 'chea'] = tf_group_scores_df[tf_group_scores_df['family'] == family]['chea'].values[0]
        
        # If TF is not in the TF families, calculate individual values
        else:
            t_value = 0.0
            t_g_ids = [h for h in gene_ids_str if h in transfac_df.index]    
            if (tf_id in transfac_df.columns) & (len(t_g_ids) > 0):
                t_value = (transfac_df.loc[t_g_ids, tf_id].sum().sum() / len(t_g_ids)) / \
                    (transfac_df[tf_id].sum().sum() / len(transfac_df))
                full_data_df.at[i, 'transfac'] = t_value 

            c_value = 0.0
            c_g_ids = [h for h in gene_ids_str if h in chea_df.index]    
            if (tf_id in chea_df.columns) & (len(c_g_ids) > 0):
                c_value = (chea_df.loc[c_g_ids, tf_id].sum().sum() / len(c_g_ids)) / \
                    (chea_df[tf_id].sum().sum() / len(chea_df))
                full_data_df.at[i, 'chea'] = c_value 

    # Fill zeroes if only one ('chea' or 'transfac') value is available
    for i in full_data_df.index:
        if full_data_df.at[i, 'chea'] == 0.0:
            full_data_df.at[i, 'chea'] = full_data_df.at[i, 'transfac']
        if full_data_df.at[i, 'transfac'] == 0.0:
            full_data_df.at[i, 'transfac'] = full_data_df.at[i, 'chea'] 

    # Predictions of Random Forest model
    full_data_df['RF_prob'] = rf_model.predict_proba(full_data_df[['TF_id',
                                                                    'p_share','all_p',
                                                                    'p_score', 'go_share',
                                                                    'all_go', 'go_score',
                                                                    'transfac','chea']].values)[:,[1]]

    # Predictions of Logistic Regression model
    full_data_df['LR_prob'] = log_model.predict_proba(full_data_df[['p_share','p_score', 
                                                'go_share', 'go_score', 'RF_prob']].values)[:,[1]]


    full_data_df.drop(columns=['RF_prob', 'p_share', 'all_p', 'p_score', 'go_share', 'all_go', 
                               'go_score', 'transfac', 'chea'], inplace=True)
    intersection = [g for g in tf_names_db.index if g in gene_ids]
    full_data_df = full_data_df[~full_data_df['TF_id'].isin(intersection)].copy(deep=True)
    full_data_df.sort_values('LR_prob', ascending=False, inplace=True)
    full_data_df.reset_index(drop=True, inplace=True)

    return full_data_df