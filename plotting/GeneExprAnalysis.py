'''
Created on Jun 6, 2017

@author: husensofteng
'''
import matplotlib
import math
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import sys, os
from scipy import stats
from multiprocessing import Pool
import seaborn as sns
from utils import draw_text
from decimal import Decimal
import argparse
#sns.set_context("paper", font_scale=1)                                                  
import matplotlib.ticker as ticker
    
def get_sample_data(meata_data_input):#, search_by='icgc_donor_id', col_to_get='rna_seq_aliquot_id', value_to_search_for=''):
    df = pd.read_table(meata_data_input, sep='\t')
    df = df[df['wgs_white_black_gray']=="Whitelist"]
    df = df[df['wgs_exclusion_white_gray']=="Whitelist"]
    return  df#df.loc[df[search_by] == value_to_search_for][col_to_get].values

def read_elements(elements_input):
    elements = pd.read_table(elements_input, sep='\t', skiprows=6, header=0)
    elements = elements[(elements['#Samples(RegMuts)']>=1)]
    return elements

def read_genes_elements(genes_mutated_input):
    genes = pd.read_table(genes_mutated_input, sep='\t', header=None, names='Gene_symbol,GeneID,#RegMuts,#RegMutsSamples,#Muts,#MutsSamples,#Elements,Elements,RegSamples,Samples'.split(','))
    #genes = genes[(genes['#MutsSamples']>30)]# & (genes['#Elements']==1)]
    genes = genes[(genes['#RegMutsSamples']>=1)]# & (genes['#Elements']==1)]
    return genes
    
def read_gene_expr(gene_exp_input, mutated_genes = []):
    print(gene_exp_input+'_filtered')
    if os.path.exists(gene_exp_input+'_filtered'+str(len(mutated_genes))):
        return pd.read_table(gene_exp_input+'_filtered'+str(len(mutated_genes)), sep='\t'), gene_exp_input+'_filtered'+str(len(mutated_genes))
    gene_expr = pd.read_table(gene_exp_input, sep='\t')
    if len(mutated_genes)>0:
        gene_expr = gene_expr[gene_expr['feature'].isin(mutated_genes)]
    gene_expr.to_csv(gene_exp_input+'_filtered'+str(len(mutated_genes)), sep='\t')
    return gene_expr, gene_exp_input+'_filtered'+str(len(mutated_genes))

def get_exper_per_sample_per_gene(gene_counts, gene_info, sample_col_to_use, donor_id_col, aliquot_id, specimen_type, all_samples, meta_data):
    gene_counts_list = []
    samples_added = []
    for mutated_donor in gene_info[sample_col_to_use].split(','):
        tumor_gene_expr_this_donor = {}
        tumor_gene_expr_this_donor_info = {}
        for j, sample_info in meta_data.loc[meta_data[donor_id_col] == mutated_donor].iterrows():
            #get expr for the current gene in the current sample
            gene_exp = gene_counts.loc[gene_counts['feature'] == gene_info['GeneID']][sample_info[aliquot_id]].values[0]
            samples_added.append(sample_info[aliquot_id])
            
            if 'Normal'.lower() in sample_info[specimen_type].lower():
                gene_counts_list.append([gene_info['GeneID'], gene_info['Gene_symbol'], gene_exp, mutated_donor,
                                     sample_info[aliquot_id], sample_info['histology_abbreviation'],
                                     sample_info['is_tumour'], sample_info['is_tumour']])
                continue
            else:
                k = ''
                if 'Primary tumour' in sample_info[specimen_type]:
                    k  = 'Primary tumour'
                elif 'Recurrent tumour' in sample_info[specimen_type]:
                    k = 'Recurrent tumour'
                elif 'Metastatic tumour' in sample_info[specimen_type]:
                    k = 'Metastatic tumour'
                else:
                    print('Unrecognized specimen type', sample_info[specimen_type])
                try:
                    tumor_gene_expr_this_donor[k].append(gene_exp)
                    tumor_gene_expr_this_donor_info[k].append(sample_info)
                except KeyError:
                    tumor_gene_expr_this_donor[k] = [gene_exp]
                    tumor_gene_expr_this_donor_info[k] = [sample_info]
        if len(tumor_gene_expr_this_donor.keys())>0:
            to_use = ''
            if 'Primary tumour' in tumor_gene_expr_this_donor.keys():
                to_use = 'Primary tumour'
            elif 'Recurrent tumour' in tumor_gene_expr_this_donor.keys():
                to_use='Recurrent tumour'
            elif 'Metastatic tumour' in tumor_gene_expr_this_donor.keys():
                to_use='Metastatic tumour'
            else:
                print('Unrecognized specimen types: ', tumor_gene_expr_this_donor.keys())
            gene_expr_avg = np.mean(tumor_gene_expr_this_donor[to_use])
            samples = tumor_gene_expr_this_donor_info[to_use]
            aliquot_ids = []
            histology = []
            is_tumour = []
            is_mut = []
            for sample in samples:
                aliquot_ids.append(sample[aliquot_id])
                histology.append(sample['histology_abbreviation'])
                is_tumour.append(sample['is_tumour'])
                is_mut.append('yes')
            #print gene_info['GeneID'], gene_expr_avg, mutated_donor, gene_info['Gene_symbol'], aliquot_ids, histology, is_tumour
            
            gene_counts_list.append([gene_info['GeneID'], gene_info['Gene_symbol'], gene_expr_avg, mutated_donor, 
                                     ','.join(list(set(aliquot_ids))), ','.join(list(set(histology))), 
                                     ','.join(list(set(is_tumour))), ','.join(list(set(is_mut)))])
    
    #add the remaining samples as not mutated in this genes
    for not_mut_sample in all_samples:
        #check if the samples of this donor have already been added
        if not_mut_sample in samples_added:
            continue
        #print 'donor id: ', meta_data.loc[meta_data[aliquot_id]==not_mut_sample]['icgc_donor_id'].values[0]
        try:
            mutated_donor = meta_data.loc[meta_data[aliquot_id]==not_mut_sample][donor_id_col].values[0]
        except IndexError:
            #print 'sample not found in metadata', not_mut_sample
            continue
        #in case of RegMotifs, skip samples that have any mutation (CFRM or others)
        if mutated_donor in gene_info['Samples'].split(','):
            continue 
        tumor_gene_expr_this_donor = {}
        tumor_gene_expr_this_donor_info = {}
        for j, sample_info in meta_data.loc[meta_data[donor_id_col] == mutated_donor].iterrows():
            #get expr for the current gene in the current sample
            gene_exp = gene_counts.loc[gene_counts['feature'] == gene_info['GeneID']][sample_info[aliquot_id]].values[0]
            samples_added.append(sample_info[aliquot_id])
            
            if 'Normal'.lower() in sample_info[specimen_type].lower():
                gene_counts_list.append([gene_info['GeneID'], gene_info['Gene_symbol'], gene_exp, mutated_donor,
                                     sample_info[aliquot_id], sample_info['histology_abbreviation'],
                                     sample_info['is_tumour'], sample_info['is_tumour']])
                continue
            else:
                k = ''
                if 'Primary tumour' in sample_info[specimen_type]:
                    k  = 'Primary tumour'
                elif 'Recurrent tumour' in sample_info[specimen_type]:
                    k = 'Recurrent tumour'
                elif 'Metastatic tumour' in sample_info[specimen_type]:
                    k = 'Metastatic tumour'
                else:
                    print('Unrecognized specimen type', sample_info[specimen_type])
                try:
                    tumor_gene_expr_this_donor[k].append(gene_exp)
                    tumor_gene_expr_this_donor_info[k].append(sample_info)
                except KeyError:
                    tumor_gene_expr_this_donor[k] = [gene_exp]
                    tumor_gene_expr_this_donor_info[k] = [sample_info]
        if len(tumor_gene_expr_this_donor.keys())>0:
            to_use = ''
            if 'Primary tumour' in tumor_gene_expr_this_donor.keys():
                to_use = 'Primary tumour'
            elif 'Recurrent tumour' in tumor_gene_expr_this_donor.keys():
                to_use='Recurrent tumour'
            elif 'Metastatic tumour' in tumor_gene_expr_this_donor.keys():
                to_use='Metastatic tumour'
            else:
                print('Unrecognized specimen types: ', tumor_gene_expr_this_donor.keys())
            gene_expr_avg = np.mean(tumor_gene_expr_this_donor[to_use])
            samples = tumor_gene_expr_this_donor_info[to_use]
            aliquot_ids = []
            histology = []
            is_tumour = []
            is_mut = []
            for sample in samples:
                aliquot_ids.append(sample[aliquot_id])
                histology.append(sample['histology_abbreviation'])
                is_tumour.append(sample['is_tumour'])
                is_mut.append('no')
            #print gene_info['GeneID'], gene_expr_avg, mutated_donor, gene_info['Gene_symbol'], aliquot_ids, histology, is_tumour
            
            gene_counts_list.append([gene_info['GeneID'], gene_info['Gene_symbol'], gene_expr_avg, mutated_donor, 
                                     ','.join(list(set(aliquot_ids))), ','.join(list(set(histology))), 
                                     ','.join(list(set(is_tumour))), ','.join(list(set(is_mut)))])
    return gene_counts_list

def get_expr_per_sample(mutated_genes, meta_data, gene_counts, gene_counts_file, sample_col_to_use='Samples'):
    donor_id_col = 'icgc_donor_id'
    aliquot_id = 'rna_seq_aliquot_id'
    specimen_type = 'specimen_type'
    
    if os.path.exists(gene_counts_file+'_counts'+'_'+sample_col_to_use):
        return pd.read_table(gene_counts_file+'_counts'+'_'+sample_col_to_use, sep='\t'), gene_counts_file+'_counts'+'_'+sample_col_to_use 
    
    all_samples = gene_counts.columns[2:]
    gene_counts_list = []
    p = Pool(10)
    for i, gene_info in mutated_genes.iterrows():
        if gene_info['GeneID']=='None':
            continue
        print(gene_info['GeneID'], gene_info['Gene_symbol'])
        #generated_list_per_sample_per_gene = get_exper_per_sample_per_gene(gene_info, sample_col_to_use, donor_id_col, aliquot_id, specimen_type, all_samples)
        #gene_counts_list.extend(generated_list_per_sample_per_gene)
        p.apply_async(get_exper_per_sample_per_gene, args=(gene_counts, gene_info, sample_col_to_use, donor_id_col, aliquot_id, specimen_type, all_samples, meta_data), callback=(gene_counts_list.extend))
    p.close()
    p.join()
    
    gene_counts_info = pd.DataFrame(gene_counts_list, columns='GeneID,Gene_symbol,Expr,DonorID,AliquotID,Cancer_type,is_tumor,is_mutated'.split(','))
    gene_counts_info.to_csv(gene_counts_file+'_counts'+'_'+sample_col_to_use, sep='\t')
    return gene_counts_info, gene_counts_file+'_counts'+'_'+sample_col_to_use

def get_wilcoxon_rank_sum_pval(l1, l2):
    pval = stats.wilcoxon(l1,l2)
    return pval

def get_pval(element_score, avg, sd):
    if sd==0:
        sd = 1
    z_score = (element_score - avg)/sd
    p_val= stats.norm.sf(abs(z_score))*2.0#gives the same result as cdf
    #p_val = (1-stats.norm.cdf(x=element_score, loc=avg, scale=sd))*2.0
    return p_val

def get_fold_change(x, y):
    if y>0.0:
        return (x)/y*1.0
    return x

def compute_pval_by_permutation(stat_val, sample1_values, sample2_values, num_permuations=100000):
    
    combined_values = np.concatenate((sample1_values, sample2_values), axis=0)
    permutated_t_values = []
    for i in np.arange(num_permuations):
        permuted_values = np.random.permutation(combined_values)
        perm_t, perm_p = stats.ttest_ind(permuted_values[0:len(sample1_values)], permuted_values[len(sample1_values):], equal_var = False)
        permutated_t_values.append(perm_t)
    pval = get_pval(stat_val, np.mean(permutated_t_values), sd=np.std(permutated_t_values))
    return pval

def process_gene_counts_per_gene(gene_df, gene_id):
    results = []
    for cancer_type, cancer_type_df in gene_df.groupby('Cancer_type'):
        #for this gene in this cancer type collect the following info
        mutated_samples = cancer_type_df[(cancer_type_df['is_tumor']=='yes') & (cancer_type_df['is_mutated']=='yes')]['DonorID'].values
        mutated_aliquots = cancer_type_df[(cancer_type_df['is_tumor']=='yes') & (cancer_type_df['is_mutated']=='yes')]['AliquotID'].values
        mutated_values = cancer_type_df[(cancer_type_df['is_tumor']=='yes') & (cancer_type_df['is_mutated']=='yes')]['Expr'].values
        notmutated_values = cancer_type_df[(cancer_type_df['is_tumor']=='yes') & (cancer_type_df['is_mutated']=='no')]['Expr'].values
        normal_values = cancer_type_df[(cancer_type_df['is_tumor']=='no') & (cancer_type_df['is_mutated']=='no')]['Expr'].values
        normal_notmut_combined_values = cancer_type_df[(cancer_type_df['is_mutated']=='no')]['Expr'].values
        
        #tumor mutated versus matching normal
        mut_donors_with_matching_normal = []
        matching_normals_values = []
        matchin_tumor_values = []
        for donor, donor_df in cancer_type_df.groupby('DonorID'):
            if len(donor_df)>=2 and 'no' in donor_df['is_tumor'].values and 'yes' in donor_df['is_mutated'].values:
                matching_normals_values.extend(donor_df[(donor_df['is_tumor']=='no') & (donor_df['is_mutated']=='no')]['Expr'].values)
                matchin_tumor_values.extend(donor_df[(donor_df['is_tumor']=='yes') & (donor_df['is_mutated']=='yes')]['Expr'].values)
                mut_donors_with_matching_normal.append(donor)
        fc_notmut = None
        fc_normal = None
        fc_normal_notmut_combined = None
        if len(mutated_values)>0:
            fc_notmut = get_fold_change(np.mean(mutated_values), np.mean(notmutated_values))
            fc_normal = get_fold_change(np.mean(mutated_values), np.mean(normal_values))
            fc_normal_notmut_combined = get_fold_change(np.mean(mutated_values), np.mean(normal_notmut_combined_values))
        fc_matching = None
        fc_per_matching_sample = []
        if len(matchin_tumor_values)>0:
            fc_matching = get_fold_change(np.mean(matchin_tumor_values), np.mean(matching_normals_values))
            fc_per_matching_sample = []
            for i, v in enumerate(matchin_tumor_values):
                fc_per_matching_sample.append(get_fold_change(v, matching_normals_values[i]))
        
        #calculate stats
        permuted_pval_notmut = None
        permuted_pval_normal = None
        permuted_pval_normal_notmut_combined = None
        if len(mutated_values)>=10:
            t_val_notmut, p_val_notmut = stats.ttest_ind(mutated_values, notmutated_values, equal_var = False)
            permuted_pval_notmut = compute_pval_by_permutation(t_val_notmut, mutated_values, notmutated_values)
            
            t_val_normal, p_val_normal = stats.ttest_ind(mutated_values, normal_values, equal_var = False)
            permuted_pval_normal = compute_pval_by_permutation(t_val_normal, mutated_values, normal_values)
            
            t_val_normal_notmut_combined, p_val_normal_notmut_combined = stats.ttest_ind(mutated_values, normal_notmut_combined_values, equal_var = False)
            permuted_pval_normal_notmut_combined = compute_pval_by_permutation(t_val_normal_notmut_combined, mutated_values, normal_notmut_combined_values)
            
            #print gene_df['Gene_symbol'].values[0], cancer_type
        if len(mutated_samples)>0:
            results.append([gene_id, gene_df['Gene_symbol'].values[0], cancer_type, mutated_genes.loc[mutated_genes['GeneID'] == gene_id]['#RegMuts'].values[0], mutated_genes.loc[mutated_genes['GeneID'] == gene_id]['#Muts'].values[0],
                            len(mutated_values), len(matchin_tumor_values),len(normal_values), len(notmutated_values), len(normal_notmut_combined_values), 
                            fc_normal, fc_notmut, fc_normal_notmut_combined, fc_matching, ','.join(str(x) for x in fc_per_matching_sample),
                            permuted_pval_normal, permuted_pval_notmut, permuted_pval_normal_notmut_combined,
                            ','.join(mut_donors_with_matching_normal), ','.join(str(x) for x in matchin_tumor_values), ','.join(str(x) for x in matching_normals_values),
                            ','.join(mutated_samples), ','.join(mutated_aliquots),
                            ','.join(str(x) for x in mutated_values), ','.join(str(x) for x in normal_values), ','.join(str(x) for x in notmutated_values), ','.join(str(x) for x in normal_notmut_combined_values)
                            ])
    return results

def process_gene_counts(gene_counts_info, mutated_genes, gene_counts_info_file):
    '''
    Analysis:
    1. for each gene check its expression in  the mutated samples against non mutated samples in the same cancer type: 
        process: groupby gene then by cancer type
        report: fold-change differences and a z-score 
    2. for each gene check its expression in  the mutated samples against all normal samples in the same cancer type: 
        process: groupby gene then by cancer type
        report: fold-change differences and a z-score
    2. for each gene check its expression in  the mutated samples against all normal samples and not-mutated tumors in the same cancer type: 
        process: groupby gene then by cancer type
        report: fold-change differences and a z-score 
    3. for each gene check expression in mutated samples and matching normal samples (if any)
        process: group by gene, groupby donorID; if len >=2 and it had normal, then compare the normal value with the tumor value
        report: fold-change | z-score for each donor
        
    '''
    if os.path.exists(gene_counts_info_file+'_results.tsv'):
        return pd.read_table(gene_counts_info_file+'_results.tsv', sep='\t'), gene_counts_info_file+'_results.tsv'
    names = ['gene_id', 'Gene_symbol', 'cancer_type', '#RegMuts', '#Muts', 
        'num_mutated_values', 'num_matchin_tumor_values', 'num_normal_values', 'num_notmutated_values', 'num_normal_notmut_combined_values', 
        'fc_normal', 'fc_notmut', 'fc_normal_notmut_combined', 
        'fc_matching', 'fc_per_matching_sample',
        'permuted_pval_normal', 'permuted_pval_notmut', 'permuted_pval_normal_notmut_combined',
        'mut_donors_with_matching_normal', 
        'matchin_tumor_values', 'matching_normals_values',
        'mutated_samples', 'mutated_aliquots',
        'mutated_values', 'normal_values', 'notmutated_values', 'normal_notmut_combined_values']
    results = []
    p = Pool(10)
    #mutated versus nonmutated tumors per gene per cancer type
    for gene_id, gene_df in gene_counts_info.groupby('GeneID'):
        #if mutated_genes.loc[mutated_genes['GeneID'] == gene_id]['#RegMuts'].values[0]<10:
        #    continue
        p.apply_async(process_gene_counts_per_gene, args=(gene_df, gene_id), callback=results.extend)
    p.close()
    p.join()
    results_combined = [x for x in results if len(x)==len(names)]
    df = pd.DataFrame(results_combined, columns=names)
    df.to_csv(gene_counts_info_file+'_results.tsv', sep='\t')
    return df, gene_counts_info_file+'_results.tsv'

def process_results(gene_stats, gene_stats_file):
    min_expr_to_consider = 0.1
    min_normal_samples = 5
    min_notmut_samples = 5
    min_matching_samples = 0
    min_mutated_values = 3
    min_percentage_matching_samples = 0.75
    pval_df = [['GeneID', 'Gene_symbol', 'Cancer_type', 'num_mutated_values', 'num_matchin_tumor_values', 'Avg FC - Not Mutated (log10)', 'P-val (-log10)', 'DiffCheck']]
    fc_matching_df = [['GeneID', 'Gene_symbol', 'Cancer_type', 'num_mutated_values', 'num_matchin_tumor_values', 'Avg FC - Not Mutated (log10)', 'Avg FC - Matching Normal (log10)', 'DiffCheck']]
    fc_pop_df = [['GeneID', 'Gene_symbol', 'Cancer_type', 'num_mutated_values', 'num_matchin_tumor_values', 'Avg FC - Not Mutated (log10)', 'Avg FC - Normal (log10)', 'DiffCheck']]
    
    for gene,gene_df in gene_stats.groupby('gene_id'):
        gene_symbol = gene_df['Gene_symbol'].values[0]
        for cancer_type, cancer_type_df in gene_df.groupby('cancer_type'):
            #if cancer_type=="Lymph-BNHL":
            #    continue
            for i, row in cancer_type_df.iterrows():

                avg_fc_matching = 'nan'#no diff
                percentage_fc_matching = 'nan'
                if (np.mean([float(x) for x in str(row['matchin_tumor_values']).split(',')])>min_expr_to_consider or np.mean([float(x) for x in str(row['matching_normals_values']).split(',')])>min_expr_to_consider):
                    if row['num_matchin_tumor_values']>min_matching_samples:
                        avg_fc_matching = row['fc_matching']
                        n = 0
                        for i, v in enumerate(row['fc_per_matching_sample'].split(',')):
                            if v>1.5 or v<0.66:
                                n+=1
                        percentage_fc_matching = n/row['num_matchin_tumor_values']*1.0
                            
                permuted_pval_notmut = 'nan'
                fc_notmut = 'nan'
                if ((np.mean([float(x) for x in str(row['mutated_values']).split(',')])>min_expr_to_consider or np.mean([float(x) for x in str(row['notmutated_values']).split(',')])>min_expr_to_consider) and
                    row['num_notmutated_values']>=min_notmut_samples): 
                    if row['num_mutated_values']>=min_mutated_values:
                        permuted_pval_notmut = row['permuted_pval_notmut']
                    
                    fc_notmut =  row['fc_notmut']
                
                permuted_pval_normal = 'nan'
                fc_normals = 'nan'
                if ((np.mean([float(x) for x in str(row['mutated_values']).split(',')])>min_expr_to_consider or np.mean([float(x) for x in str(row['normal_values']).split(',')])>min_expr_to_consider) and 
                    row['num_normal_values']>=min_normal_samples):
                    if row['num_mutated_values']>=min_mutated_values:
                        permuted_pval_normal = row['permuted_pval_normal']
                    fc_normals = row['fc_normal']
                
                '''
                Check wether the mutated samples are different that the 
                    1) not_mutated tumors (pval<0.05) or nan, 
                    2) fc_matching>1.5 or <0.66 in at least min_percentage_matching_samples percentage or nan
                    3) permuted_pval_normal<0.05 or nan
                '''
                if permuted_pval_notmut!='nan':
                    #print 'not nan permuted_pval_notmut', permuted_pval_notmut
                    if (permuted_pval_notmut<0.05 and 
                        (avg_fc_matching>1.5 or avg_fc_matching<0.66 or avg_fc_matching=='nan') and (percentage_fc_matching>min_percentage_matching_samples or percentage_fc_matching=='nan') and
                        (permuted_pval_normal<0.05 or permuted_pval_normal=='nan')
                        ):
                        pval_df.append([gene, gene_symbol, cancer_type, cancer_type_df['num_mutated_values'].values[0], cancer_type_df['num_matchin_tumor_values'].values[0], fc_notmut, permuted_pval_notmut, 'Diff'])
                    else:
                        pval_df.append([gene, gene_symbol, cancer_type, cancer_type_df['num_mutated_values'].values[0], cancer_type_df['num_matchin_tumor_values'].values[0], fc_notmut, permuted_pval_notmut, 'NotDiff'])
                if fc_notmut=='nan':
                    fc_notmut=1
                if avg_fc_matching!='nan' and percentage_fc_matching!='nan' and permuted_pval_notmut=='nan': 
                    if ((avg_fc_matching>1.5 or avg_fc_matching<0.66) and (percentage_fc_matching>min_percentage_matching_samples)):
                        fc_matching_df.append([gene, gene_symbol, cancer_type, cancer_type_df['num_mutated_values'].values[0], cancer_type_df['num_matchin_tumor_values'].values[0], fc_notmut, avg_fc_matching, 'Diff'])
                    else:
                        fc_matching_df.append([gene, gene_symbol, cancer_type, cancer_type_df['num_mutated_values'].values[0], cancer_type_df['num_matchin_tumor_values'].values[0], fc_notmut, avg_fc_matching, 'NotDiff'])
                
                if fc_normals!='nan' and fc_notmut!='nan':
                    if (avg_fc_matching=='nan' and  percentage_fc_matching=='nan' and permuted_pval_notmut=='nan' and
                          (fc_normals>1.5 or fc_normals<0.66) and (fc_notmut>1.5 or fc_notmut<0.66)):
                        fc_pop_df.append([gene, gene_symbol, cancer_type, cancer_type_df['num_mutated_values'].values[0], cancer_type_df['num_matchin_tumor_values'].values[0], fc_notmut, fc_normals, 'Diff'])
                    else:
                        fc_pop_df.append([gene, gene_symbol, cancer_type, cancer_type_df['num_mutated_values'].values[0], cancer_type_df['num_matchin_tumor_values'].values[0], fc_notmut, fc_normals, 'NotDiff'])
    
    df_pval = pd.DataFrame(pval_df[1:], columns=pval_df[0]) 
    df_matching_normal = pd.DataFrame(fc_matching_df[1:], columns=fc_matching_df[0])
    df_fc_pop = pd.DataFrame(fc_pop_df[1:], columns=fc_pop_df[0])
    return [df_pval]#, df_matching_normal, df_fc_pop]

def plot_gene_expression(dfs,output_dir):
    sns.set_style(style='white')
    fig = plt.figure(figsize=(8,6), linewidth=1.0)#design a figure with the given size
    gs = gridspec.GridSpec(len(dfs), 1, wspace=0.0, hspace=0.0)#create 4 rows and three columns with the given ratio for each
    axes = []
    x_accept = 6
    y_accept = 6
    for i,df in enumerate(dfs):
        if i>0:
            ax = fig.add_subplot(gs[i,0], sharex=axes[0])
        else:
            ax = fig.add_subplot(gs[i,0])
        axes.append(ax)
        x_col = 'Avg FC - Not Mutated (log10)'
        y_col = 'P-val (-log10)'
        if 'Avg FC - Matching Normal (log10)' in df.columns:
            y_col = 'Avg FC - Matching Normal (log10)'
        elif 'Avg FC - Normal (log10)' in df.columns:
            y_col = 'Avg FC - Normal (log10)'
        
        df.to_csv(output_dir + '/{y_label}_expr_df.tsv'.format(y_label=y_col), sep='\t')
        log_value = 10
        if 'log10' in y_col:
            log_value = 10
        
        df[x_col] = np.where(df[x_col]==0, 1e-300, df[x_col])
        df[y_col] = np.where(df[y_col]==0, 1e-300, df[y_col])
        df[x_col] = df[x_col].apply(lambda x: math.log(x,10))
        df[y_col] = df[y_col].apply(lambda x: math.log(x,log_value)*-1) 
        
        min_x = int(df[x_col].min())
        if min_x < x_accept*-1:
            min_x = x_accept*-1
            df[x_col] = np.where(df[x_col]<x_accept*-1, x_accept*-1, df[x_col])
            
        max_x = df[x_col].max()
        if max_x>x_accept:
            max_x=x_accept
            df[x_col] = np.where(df[x_col]>x_accept, x_accept, df[x_col])
        
        min_y = int(df[y_col].min())
        if min_y < -1*y_accept:
            min_y = -1*y_accept
            df[y_col] = np.where(df[y_col] < -1*y_accept, -1*y_accept, df[y_col])
        
        max_y = df[y_col].max()
        if max_y>y_accept:
            max_y=y_accept
            df[y_col] = np.where(df[y_col]>y_accept, y_accept, df[y_col])
        df['col'] = np.where(df['DiffCheck']=='Diff', 'green', 'grey')
        ax.scatter(x=x_col, y=y_col, data=df, color=df['col'])
        
        #ax.plot([i for i in np.arange(min_x, df[x_col].max()+1)], [-1*np.log10(0.05) for i in np.arange(min_x, df[x_col].max()+1)], color='red', linewidth=0.5)
        
        #plot gene labels
        for i, r in df.iterrows():
            if (r['DiffCheck']=='Diff' and ('P-val' in y_col and (r[y_col]>2 or r[y_col]<-2))):#((r[x_col]<=min_x or r[x_col]>=max_x) and (r[y_col]<=min_y or r[y_col]>=max_y)) or 
                x_p = r[x_col]
                if x_p<min_x:
                    x_p = min_x
                elif x_p>max_x:
                    x_p = max_x
                y_pos = r[y_col] 
                if y_pos<min_y:
                    y_pos = min_y
                elif y_pos>max_y:
                    y_pos = max_y
                if 'P-val' in y_col:
                    draw_text(ax=ax, x=x_p+0.3, y=y_pos+0.5, text=r['Gene_symbol'], rotation=0)#, color, fontsize, horizontalalignment, rotation, verticalalignment)
        #set axis limits
        if min_y>0:
            min_y=0
        ax.set_xlim(min_x-0.5, max_x+2)
        if min_y<0:
            min_y = min_y-2
        ax.set_ylim(min_y, max_y+2)
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        
        '''ax.plot([math.log(15/10.0, 2) for i in np.arange(min_y,max_y)], [i for i in np.arange(min_y,max_y)], linestyle='--', color='red', linewidth=0.5)
        ax.plot([math.log(10/15.0, 2) for i in np.arange(min_y,max_y)], [i for i in np.arange(min_y,max_y)], linestyle='--', color='red', linewidth=0.5)
        
        if 'P-val' in y_col:
            print min_x
            ax.plot([i for i in np.arange(min_x,max_x)], [-1*np.log10(0.05) for i in np.arange(min_x,max_x)], linestyle='-', color='red', linewidth=0.5)
        else:
            ax.plot([i for i in np.arange(min_x,max_x)], [math.log(15/10.0,2) for i in np.arange(min_x,max_x)], linestyle='-', color='red', linewidth=0.5)
            ax.plot([i for i in np.arange(min_x,max_x)], [math.log(10/15.0, 2) for i in np.arange(min_x,max_x)], linestyle='-', color='red', linewidth=0.5)
        '''
    for ax in axes[:-1]:
        ax.get_xaxis().set_visible(False)
        #sns.despine(bottom=True, ax=ax)
    sns.despine(ax=ax)
    gs.tight_layout(fig, pad=2, h_pad=0.0, w_pad=0.0)
    
    plt.savefig(output_dir + 'sig.pdf')
    plt.clf()
    plt.close()
    print('done')

#plot RP11-731F5.1
def plot_gene_expr(dfs, output_dir):

    sns.set_style(style='white')
    fig = plt.figure(figsize=(8,6), linewidth=1.0)#design a figure with the given size
    gs = gridspec.GridSpec(len(dfs), 1, wspace=0.0, hspace=0.0)#create 4 rows and three columns with the given ratio for each
    axes = []
    for i,df in enumerate(dfs):
        ax = fig.add_subplot(gs[i,0])
        axes.append(ax)
        x_col = 'Avg FC - Not Mutated (log10)'
        y_col = 'P-val (-log10)'
        if 'Avg FC - Matching Normal (log10)' in df.columns:
            y_col = 'Avg FC - Matching Normal (log10)'
        elif 'Avg FC - Normal (log10)' in df.columns:
            y_col = 'Avg FC - Normal (log10)'
        
        df.to_csv('/home/huum/projs/regMotifs/analysis/RNA-seq/{y_label}_expr_df.tsv'.format(y_label=y_col), sep='\t')
        
        if 'log2' in y_col:
            log_value = 2
            df[y_col] = np.where(df[y_col]<=1e-2, 1e-2, df[y_col])
        else:
            df[y_col] = np.where(df[y_col]==0, 1e-100, df[y_col])
            log_value = 10
        df[x_col] = np.where(df[x_col]<=1e-2, 1e-2, df[x_col])
        df[x_col] = df[x_col].apply(lambda x: math.log(x,2))
        #except TypeError:
        #    print x_col, y_col, df[x_col]
        df[y_col] = df[y_col].apply(lambda x: math.log(x,log_value)*-1) 
        
        df['col'] = np.where(df['DiffCheck']=='Diff', 'green', 'grey')
        ax.scatter(x=x_col, y=y_col, data=df, color=df['col'])
        
        #plot gene labels
        for i, r in df.iterrows():
            if (('P-val' in y_col) and (r['DiffCheck']=='Diff' or r[x_col]<-1.5 or r[x_col]>1.5)):#((r[x_col]<=min_x or r[x_col]>=max_x) and (r[y_col]<=min_y or r[y_col]>=max_y)) or 
                if 'P-val' in y_col:
                    draw_text(ax=ax, x=r[x_col]+0.2, y=r[y_col]+0.8, text=r['Gene_symbol'], rotation=0)#, color, fontsize, horizontalalignment, rotation, verticalalignment)
        #set axis limits
        #ax.set_xlim(min_x-0.5, max_x+2)
        #ax.set_ylim(min_y, max_y+2)
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        
        '''ax.plot([math.log(15/10.0, 2) for i in np.arange(min_y,max_y)], [i for i in np.arange(min_y,max_y)], linestyle='--', color='red', linewidth=0.5)
        ax.plot([math.log(10/15.0, 2) for i in np.arange(min_y,max_y)], [i for i in np.arange(min_y,max_y)], linestyle='--', color='red', linewidth=0.5)
        
        if 'P-val' in y_col:
            print min_x
            ax.plot([i for i in np.arange(min_x,max_x)], [-1*np.log10(0.05) for i in np.arange(min_x,max_x)], linestyle='-', color='red', linewidth=0.5)
        else:
            ax.plot([i for i in np.arange(min_x,max_x)], [math.log(15/10.0,2) for i in np.arange(min_x,max_x)], linestyle='-', color='red', linewidth=0.5)
            ax.plot([i for i in np.arange(min_x,max_x)], [math.log(10/15.0, 2) for i in np.arange(min_x,max_x)], linestyle='-', color='red', linewidth=0.5)
        '''
    for ax in axes[:-1]:
        #ax.get_xaxis().set_visible(False)
        #sns.despine(bottom=True, ax=ax)
        ax.set_xlabel('')
    sns.despine(ax=ax)
    gs.tight_layout(fig, pad=2, h_pad=0.0, w_pad=0.0)
    
    plt.savefig(output_dir + 'sig.pdf')
    plt.clf()
    plt.close()
    print('done')
    
def get_sig_expr_events(gene_stats, gene_stats_file):
    pval_df = [['GeneID', 'Gene_symbol', 'Cancer_type', 'num_mutated_values', 'num_matchin_tumor_values', 'Avg FC - Not Mutated (log10)', 'P-val (-log10)', 'DiffCheck']]
    min_mutated_values = 10
    min_notmutated_values = 5
    min_expr_to_consider = 0.1
    
    for gene,gene_df in gene_stats.groupby('gene_id'):
        gene_symbol = gene_df['Gene_symbol'].values[0]
        for cancer_type, cancer_type_df in gene_df.groupby('cancer_type'):
            #if cancer_type=="Lymph-BNHL":
            #    continue
            for i, row in cancer_type_df.iterrows():
                permuted_pval_notmut = 'nan'
                fc_notmut = 'nan'
                if ((np.mean([float(x) for x in str(row['mutated_values']).split(',')])>min_expr_to_consider or np.mean([float(x) for x in str(row['notmutated_values']).split(',')])>min_expr_to_consider) and
                    row['num_notmutated_values']>=min_notmutated_values): 
                    if row['num_mutated_values']>=min_mutated_values:
                        permuted_pval_notmut = row['permuted_pval_notmut']
                    
                    fc_notmut =  row['fc_notmut']
                
                if permuted_pval_notmut!='nan':
                    if (permuted_pval_notmut<0.05):
                        pval_df.append([gene, gene_symbol, cancer_type, cancer_type_df['num_mutated_values'].values[0], cancer_type_df['num_matchin_tumor_values'].values[0], fc_notmut, permuted_pval_notmut, 'Diff'])
                    else:
                        pval_df.append([gene, gene_symbol, cancer_type, cancer_type_df['num_mutated_values'].values[0], cancer_type_df['num_matchin_tumor_values'].values[0], fc_notmut, permuted_pval_notmut, 'NotDiff'])
    df = pd.DataFrame(pval_df[1:], columns=pval_df[0])
    return df

def plot_scatter_geneexpr(df, output_dir):
    
    sns.set_style(style='white')
    fig = plt.figure(figsize=(6,4), linewidth=1.0)#design a figure with the given size
    gs = gridspec.GridSpec(1, 1, wspace=0.0, hspace=0.0)#create 4 rows and three columns with the given ratio for each
    ax = fig.add_subplot(gs[0,0])
    x_col = 'Avg FC - Not Mutated (log10)'
    y_col = 'P-val (-log10)'
    df.to_csv(output_dir + '/{y_label}_expr_df.tsv'.format(y_label=y_col), sep='\t')
    
    df[x_col] = np.where(df[x_col]<=1e-2, 1e-2, df[x_col])
    df[x_col] = df[x_col].apply(lambda x: math.log(x,10))
    df[y_col] = np.where(df[y_col]==0, 1e-100, df[y_col])
    df[y_col] = df[y_col].apply(lambda x: math.log(x,10)*-1) 
    print(df['num_mutated_values'])
    df['col'] = np.where(df['DiffCheck']=='Diff', 'green', 'grey')
    ax.scatter(x=x_col, y=y_col, data=df, color=df['col'])
        
    #plot gene labels
    for i, r in df.iterrows():
        if (r['DiffCheck']=='Diff'):#((r[x_col]<=min_x or r[x_col]>=max_x) and (r[y_col]<=min_y or r[y_col]>=max_y)) or 
            x_shift = 0.2
            y_shift = 0.3
            if 'RP11' in r['Gene_symbol']:
                x_shift = 1.0   
                y_shift = -0.3
            draw_text(ax=ax, x=r[x_col]+x_shift, y=r[y_col]+y_shift, text=r['Gene_symbol'], rotation=0)#, color, fontsize, horizontalalignment, rotation, verticalalignment)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_ylim(0,df[y_col].max()+1)
    sns.despine(ax=ax)
    gs.tight_layout(fig, pad=2, h_pad=0.0, w_pad=0.0)
    plt.savefig(output_dir + '/sig_min10samples.pdf')
    plt.close()
    print('done')
    
def box_plot_per_gene_cancertype(fig, gs, row_num, gene_counts_info_stats, genes_cancertypes):
    axes = []
    col_start = 0
    col_end = 2
    for i, gene in enumerate(genes_cancertypes):
        ax = fig.add_subplot(gs[row_num,col_start:col_end])
        col_start+=2
        col_end+=2
        axes.append(ax)
        df = gene_counts_info_stats[(gene_counts_info_stats['Gene_symbol']==gene.split(':')[0]) & (gene_counts_info_stats['cancer_type']== gene.split(':')[1])]   
        gene_values = []
        p_val = 'nan'
        for ind, r in df.iterrows():
            p_val = ('%.0E' % Decimal(r['permuted_pval_notmut'])).replace('E','e')
            gene_values.extend([[float(i), 'Mutatant', gene.split(':')[0], r['cancer_type']] for i in r['mutated_values'].split(',')])
            gene_values.extend([[float(i), 'WT', gene.split(':')[0], r['cancer_type']] for i in r['notmutated_values'].split(',')])
            #try:
            #    gene_values.extend([[float(i), 'Normal', gene, r['cancer_type']] for i in r['normal_values'].split(',')])
            #except AttributeError:
            #    continue
        gene_df = pd.DataFrame(gene_values, columns=['Expr (FPKM-UQ)', 'Mutation Status','Gene', 'Cancer Type'])
        print(gene_df)
        tick_spacing = 5
        if gene_df['Expr (FPKM-UQ)'].max() > 100:
            tick_spacing =50
        if gene_df['Expr (FPKM-UQ)'].max() > 300:
            tick_spacing =150
        ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        #ax.set_ylim(np.arange(0, gene_df['Expr (FPKM-UQ)'].max(), step))
        sns.boxplot(x='Mutation Status', y='Expr (FPKM-UQ)', data=gene_df, ax=ax, palette={'Mutatant':'red', 'WT':'blue'},
                    fliersize=1.0, linewidth=1, notch=False)
        #draw_text(ax, x, y, text, color, fontsize, horizontalalignment, rotation, verticalalignment)
        ax.set_title(gene.split(':')[0] + "\n(P = {pval})".format(pval=p_val))
        for i,artist in enumerate(ax.artists):
            # Set the linecolor on the artist to the facecolor, and set the facecolor to None
            col = artist.get_facecolor()
            artist.set_edgecolor(col)
            artist.set_facecolor('None')
        
            # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
            # Loop over them here, and use the same colour as above
            for j in range(i*6,i*6+6):
                line = ax.lines[j]
                line.set_color(col)
                line.set_mfc(col)
                line.set_mec(col)
            
    for ax in axes[1:]:
        ax.set_ylabel('')
        ax.set_xlabel('')
    
    for ax in axes:
        ax.set_xlabel('')
        #ax.spines['left'].set_color('grey')
        #ax.spines['left'].set_linewidth(1)
        sns.despine(ax=ax, bottom=True)
    
    return


def parse_args():
    '''Parse command line arguments'''
    
    parser = argparse.ArgumentParser(description='Plot Expression')
    parser.add_argument('--genes_mutated_input', default='', help='')
    parser.add_argument('--meta_data', default='', help='')
    parser.add_argument('--gene_expr_intput', default='', help='')
    parser.add_argument('--output_dir', default='', help='')
    

    
    return parser.parse_args(sys.argv[1:])

if __name__ == '__main__':
    #merged200bp_extended200bp_nofullexon_pancan , PancanElements
    args = parse_args()
    
    genes_mutated_input = args.genes_mutated_input
    meta_data = args.meta_data
    gene_expr_intput = args.gene_expr_intput
    
    meta_data = get_sample_data(meta_data)
    print('metadata')
    mutated_genes = read_genes_elements(genes_mutated_input)
    print('mutated_genes')
    gene_counts, gene_counts_file =  read_gene_expr(gene_expr_intput, mutated_genes['GeneID'])
    print('expr loaded')
    gene_counts_info, gene_counts_info_file = get_expr_per_sample(mutated_genes, meta_data, gene_counts, gene_counts_file, sample_col_to_use='Samples')
    print('gene counts')
    gene_counts_info_stats, gene_counts_info_stats_file = process_gene_counts(gene_counts_info, mutated_genes, gene_counts_info_file)
    print('stats done')
    #make a scatter plot for genes that are mutated in at least 10 samples with expr data (pval and avg FC (WT)
    df = get_sig_expr_events(gene_counts_info_stats, gene_counts_info_stats_file)
    plot_scatter_geneexpr(df, args.output_dir)
    
    #box_plot_per_gene_cancertype(gene_counts_info_stats, genes_cancertypes=['VHL:Kidney-RCC', 'BCL2:Lymph-BNHL', 'MYC:Lymph-BNHL', 'RP11-731F5.1:Lymph-BNHL'])#, 'TERT':['Skin-Melanoma', 'Bladder-TCC','CNS-Oligo','Thy-AdenoCA']})
    #box_plot_per_gene_cancertype(gene_counts_info_stats, genes_cancertypes=['VHL:Kidney-RCC'], out_ext='vhlexpr', fig_width=4, fig_hieght=2)
    
    
    #dfs = process_results(gene_counts_info_stats, gene_counts_info_stats_file)
    #plot_gene_expression(dfs, output_dir)
    #plot_gene_expr(dfs, output_dir)
    
