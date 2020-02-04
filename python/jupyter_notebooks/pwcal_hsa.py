import numpy as np
import pandas as pd
import pyarrow
from multiprocessing import Pool
from collections import Counter
from itertools import repeat
import scipy.stats as st

def __get_ptws(ids, pathways_db, go_terms_db):
    '''Get list of patways/terms based in gene IDs'''
    pathways = list(pathways_db[pathways_db['gene_id'].isin(ids)]['pathway'])
    go_terms = list(go_terms_db[go_terms_db['gene_id'].isin(ids)]['GO_term'])
    return pathways, go_terms


def __pcal(n, tf_id, tf_name, cases_db, pathways_db, go_terms_db, dist, exclude_deg_tf):
    ''' Find an overlap between KEGG pathways/GO terms related to differentially 
    expressed genes in each case and each TF in the library and calculate corresponding
    scores based on this overlap and the relative pathway/term importance'''
    full_data_df = pd.DataFrame(columns=['case_id', 'case', 'TF', 'TF_ids', 'p_share', 
                                        'all_p', 'p_score', 'go_share', 'all_go', 
                                        'go_score'])
    indx = 0
    print(f'Processing TF {n+1}: {tf_name} {tf_id}')  

    # Get unique pathways and terms related to the TF and their counts 
    tf_ptws, tf_terms = __get_ptws([tf_id], pathways_db, go_terms_db)
    tf_pw_count = Counter(tf_ptws)
    tf_term_count = Counter(tf_terms)
    all_p = len(tf_pw_count)
    all_go = len(tf_term_count)

    for c in cases_db['case'].index: 
        # Get unique pathways and terms related to differentially expressed genes in each
        # case and their counts. If TF from the library is among the genes, its pathways
        # and terms can be excluded to avoid data leakage (for the testing purposes).
        if exclude_deg_tf:
            case_ids = [g for g in cases_db.at[c,'degs_entrez'] if g != tf_id]
        else:
            case_ids = cases_db.at[c,'degs_entrez'] 
        case_ptws, case_terms = __get_ptws(case_ids, pathways_db, go_terms_db)
        case_ptw_count = Counter(case_ptws)
        case_term_count = Counter(case_terms)

        # Calculate relative pathway importance for this case
        # The score is calculated based on beta distribution
        pw_df = pd.DataFrame.from_dict(case_ptw_count, orient='index',columns=['count'])
        pw_df['importance'] = 0.0
        for i in pw_df.index:
            pw_df.at[i, 'importance'] = pw_df.at[i, 'count'] /np.max(pw_df['count'].values)
        pw_df['importance']  = dist.pdf(pw_df['importance'])

        # Calculate relative GO term importance for this case
        go_df = pd.DataFrame.from_dict(case_term_count, orient='index',columns=['count'])
        go_df['importance'] = 0.0
        for i in go_df.index:
            go_df.at[i, 'importance'] = go_df.at[i, 'count'] /np.max(go_df['count'].values)
        go_df['importance']  = dist.pdf(go_df['importance'])

        # Calculate relative pathway/term overlap and scores based on 'importance'
        common_p = 0
        common_go = 0
        p_score = 0.0
        go_score = 0.0
        coef = 1
        if all_p > 0:
            for i in tf_pw_count:
                if i in case_ptw_count:
                    common_p += 1
                    p_score += coef*pw_df.at[i, 'importance'] *case_ptw_count[i] 
        if all_go > 0:
            for t in tf_term_count:
                if t in case_term_count:
                    common_go +=1
                    go_score += coef*go_df.at[t, 'importance'] *case_term_count[t]

        if all_p > 0:
            p_share = common_p/all_p
        else: p_share = 0.0

        if all_go > 0:
            go_share = common_go/all_go
        else: go_share = 0.0
        
        # Add to DataFrame
        full_data_df.loc[indx] = [c, cases_db.at[c,'case'], tf_name, tf_id, p_share, all_p,  
                                  p_score, go_share, all_go, go_score]
        indx +=1

    return full_data_df
        

def pwcal (cases_db, exclude_deg_tf = False, n_jobs = 2):
    '''
    Calculate pathway and GO term similarity between TFs and DEGs
    Parameters:
    ----------
    cases_db : DataFrame object
        Preprocessed table containing lists of differentially expressed genes

    exclude_deg_tf : bool, default False
        If True, KEGG pathways and GO terms related to TFs, which appear in the list of DEGs
        differentially expressed genes, are not included.

    n_jobs : int, default 2
        Number of parallel processes
    
    Returns : DataFrame
        Returns DataFrame with corresponding pathway/term scores
    '''

    # Table of TF to TF ID associations
    tf_names_db = pd.read_parquet('data/tf_names_hsa.parquet')
    
    # Table of gene ID to pathway associations
    pathways_db = pd.read_parquet('data/pathways_hsa_full_hg38.parquet')
    pathways_db.drop_duplicates(['gene_id','pathway'],keep='first', inplace=True)
    
    # Table of gene ID to GO terms (from 'processes' cathegory) associations
    go_terms_db = pd.read_parquet('data/go_terms_hsa_processes.parquet')
    go_terms_db.drop_duplicates(['gene_id','GO_term'],keep='first', inplace=True)
    
    # Arbitrary distribution for calculating 'pathway importance'
    dist = st.beta(3,5)

    # Create a process pool 
    pool = Pool(processes=n_jobs) 
    print(" {} processes started...".format(n_jobs))

    results = pool.starmap(__pcal, zip(
        range(len(tf_names_db.index)),
        list(tf_names_db.index),
        list(tf_names_db['gene'].values),
        repeat(cases_db),
        repeat(pathways_db),
        repeat(go_terms_db),
        repeat(dist),
        repeat(exclude_deg_tf)
        ))

    pool.close()
    pool.join()
    results = pd.concat(results)
    results.reset_index(drop=True, inplace=True)
    print('Finished')
    return results



