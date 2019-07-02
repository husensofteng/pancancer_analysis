'''
Created on Jun 6, 2017

@author: husensofteng
'''
import pandas as pd
import numpy as np
import sys, os
from scipy import stats
from multiprocessing import Pool

def get_sample_data(meata_data_input):#, search_by='icgc_donor_id', col_to_get='rna_seq_aliquot_id', value_to_search_for=''):
    df = pd.read_table(meata_data_input, sep='\t')
    df = df[df['wgs_white_black_gray']=="Whitelist"]
    df = df[df['wgs_exclusion_white_gray']=="Whitelist"]
    return  df#df.loc[df[search_by] == value_to_search_for][col_to_get].values

def read_elements(elements_input):
    elements = pd.read_table(elements_input, sep='\t', skiprows=6, header=0)
    elements = elements[(elements['#Samples(RegMuts)']>1)]
    return elements

def read_genes_elements(genes_mutated_input):
    genes = pd.read_table(genes_mutated_input, sep='\t', header=None, names='Gene_symbol,GeneID,#RegMuts,#RegMutsSamples,#Muts,#MutsSamples,#Elements,Elements,RegSamples,Samples'.split(','))
    genes = genes[(genes['#RegMutsSamples']>=2)]# & (genes['#Elements']==1)]
    return genes
    
def read_gene_expr(gene_exp_input, mutated_genes = []):
    print gene_exp_input+'_filtered'
    if os.path.exists(gene_exp_input+'_filtered'+str(len(mutated_genes))):
        return pd.read_table(gene_exp_input+'_filtered'+str(len(mutated_genes)), sep='\t'), gene_exp_input+'_filtered'+str(len(mutated_genes))
    gene_expr = pd.read_table(gene_exp_input, sep='\t')
    if len(mutated_genes)>0:
        gene_expr = gene_expr[gene_expr['feature'].isin(mutated_genes)]
    gene_expr.to_csv(gene_exp_input+'_filtered'+str(len(mutated_genes)), sep='\t')
    return gene_expr, gene_exp_input+'_filtered'+str(len(mutated_genes))

def get_exper_per_sample_per_gene(gene_info, sample_col_to_use, donor_id_col, aliquot_id, specimen_type, all_samples):
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
                    print 'Unrecognized specimen type', sample_info[specimen_type]
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
                print 'Unrecognized specimen types: ', tumor_gene_expr_this_donor.keys()
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
                    print 'Unrecognized specimen type', sample_info[specimen_type]
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
                print 'Unrecognized specimen types: ', tumor_gene_expr_this_donor.keys()
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
    p = Pool(15)
    for i, gene_info in mutated_genes.iterrows():
        if gene_info['GeneID']=='None':
            continue
        print gene_info['GeneID'], gene_info['Gene_symbol']
        #generated_list_per_sample_per_gene = get_exper_per_sample_per_gene(gene_info, sample_col_to_use, donor_id_col, aliquot_id, specimen_type, all_samples)
        #gene_counts_list.extend(generated_list_per_sample_per_gene)
        p.apply_async(get_exper_per_sample_per_gene, args=(gene_info, sample_col_to_use, donor_id_col, aliquot_id, specimen_type, all_samples), callback=(gene_counts_list.extend))
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

def compute_pval_by_permutation(stat_val, sample1_values, sample2_values, num_permuations=1000):
    
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
        if len(mutated_values)>=3:
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
    p = Pool(15)
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
            #print gene_id, cancer_type, matching_normals_values, matchin_tumor_values, mut_donors_with_matching_normal
            #print gene_df['Gene_symbol'].values[0], cancer_type, len(mutated_values), np.mean(mutated_values), np.std(mutated_values), len(notmutated_values), np.mean(notmutated_values), np.std(notmutated_values)
    
if __name__ == '__main__':
    #merged200bp_extended200bp_nofullexon_pancan , PancanElements
    genes_mutated_input = 'analysis/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv_Genes.tsv'
    meta_data = 'analysis/RNA-seq/extended_meatadata_syn7416381'
    gene_expr_intput = 'analysis/RNA-seq/tophat_star_fpkm_uq.v2_aliquot_gl.tsv'
    
    meta_data = get_sample_data(meta_data)
    print 'metadata'
    mutated_genes = read_genes_elements(genes_mutated_input)
    print 'mutated_genes'
    gene_counts, gene_counts_file =  read_gene_expr(gene_expr_intput, mutated_genes['GeneID'])
    print 'expr loaded'
    gene_counts_info, gene_counts_info_file = get_expr_per_sample(mutated_genes, meta_data, gene_counts, gene_counts_file, sample_col_to_use='RegSamples')
    print 'gene counts'
    process_gene_counts(gene_counts_info, mutated_genes, gene_counts_info_file)
    print 'stats done'
    