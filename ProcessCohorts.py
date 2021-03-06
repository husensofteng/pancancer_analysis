'''
Created on Feb 9, 2017

@author: husensofteng
'''
import os
import sys
import argparse
from multiprocessing import Pool
import pandas as pd
import numpy as np
import subprocess
from collections import Counter
from pybedtools import BedTool, set_tempdir, cleanup
from decimal import Decimal
from multiprocessing import Pool
from itertools import product
from scipy import stats
import math
import glob


import Utilities

def generate_cohorts(mutation_input_files, cohorts, 
                     mutations_cohorts_dir, stats_ext, num_cores):
    created_cohorts_dict = {}
    if num_cores>1:
        p = Pool(num_cores)
    
    for mutations_input_file in mutation_input_files:
        for cohort_value in cohorts:
            cohort_name = cohort_value.split('=')[0].split('::')[0]
            if cohort_name not in created_cohorts_dict.keys():
                created_cohorts_dict[cohort_name] = []
            cohort_file = mutations_cohorts_dir+'/'+cohort_value.split('=')[0].split('::')[0] + '_' + mutations_input_file.split('/')[-1]
            cohort_file_statpvalues = cohort_file + stats_ext
            created_cohorts_dict[cohort_name].append(cohort_file)
            if os.path.exists(cohort_file) or os.path.exists(cohort_file_statpvalues):
                continue
            print("Generating {}...".format(cohort_file))
            if num_cores>1:
                p.apply_async(generate_cohorts_per_mut_file, args = (
                    cohort_value, cohort_file, mutations_input_file))
            else:
                generate_cohorts_per_mut_file(cohort_value, cohort_file, 
                                              mutations_input_file)
            #read the mutations file line by line and check if the sample id is inthe lis of samples
    if num_cores>1:
        p.close()
        p.join()
        
    return created_cohorts_dict

def generate_cohorts_per_mut_file(cohort_value, cohort_file, mutations_input_file):
    if '=' not in cohort_value:
        if '::' not in cohort_value:
            awk_stm = """awk 'BEGIN{FS=OFS="\t"}{if($6=="%s") print $0 >> "%s"}' %s""" %(cohort_value, cohort_file, mutations_input_file)
            os.system(awk_stm)
        elif '::' in cohort_value:
            cohort_keywords_to_exlude = cohort_value.split('::')[1].split(",")
            cond = """ $6!~"%s" """ %(cohort_keywords_to_exlude[0])
            if len(cohort_keywords_to_exlude)>1:
                for keyword in cohort_keywords_to_exlude[1::]:
                    cond += """ && $6!~"%s" """ %(keyword)
            awk_stm = """awk 'BEGIN{FS=OFS="\t"}{if(%s) print $0 >> "%s"}' %s""" %(cond, cohort_file, mutations_input_file)
            os.system(awk_stm)
        elif '*' in cohort_value:
            awk_stm = """awk 'BEGIN{FS=OFS="\t"}{print $0 >> "%s"}' %s""" %(cohort_file, mutations_input_file)
            os.system(awk_stm) 
    elif '=' in cohort_value:
        cohort_sample_ids = cohort_value.split('=')[1].split(',')
        if len(cohort_sample_ids)>0:
            n=0
            with open(mutations_input_file, 'r') as mutations_infile:
                #write mutations that belong to the samples from cohort_value to a temp file
                with open(cohort_file + "_mutations_temp", 'w') as mutations_temp_outfile:
                    line = mutations_infile.readline()
                    while line!="":
                        if line.strip().split('\t')[8].strip() in cohort_sample_ids or line.strip().split('\t')[7].strip() in cohort_sample_ids:
                            mutations_temp_outfile.write(line)
                            n+=1
                        line = mutations_infile.readline()
        print("written lines to {}: {}".format(cohort_file + "_mutations_temp", n))
        #write the temp file to cohort_file (safe writting)
        awk_stm =r"""awk 'BEGIN{{FS=OFS="\t"}}{{ printf ("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32)}}' {cohort_temp} > {cohort_ofile} """.format(cohort_temp = cohort_file + "_mutations_temp", cohort_ofile =cohort_file)

        #awk_stm = """awk 'BEGIN{FS=OFS="\t"}{print $0 >> "%s"}' %s""" %(cohort_file, cohort_file + "_mutations_temp")
        os.system(awk_stm)
        print("written lines to {}: {}".format(cohort_file, n))
        if os.path.exists(cohort_file + "_mutations_temp"):
            os.remove(cohort_file + "_mutations_temp")            
    
    return cohort_file

def get_cohorts(cohort_names_input):
    cohorts = []
    lines = [cohort_names_input]
    if os.path.isfile(cohort_names_input) and os.path.exists(cohort_names_input):
        with open(cohort_names_input, 'r') as cohort_names_infile:
            lines = cohort_names_infile.readlines()
    for line in lines:
        if line.strip()!='' and not line.startswith('//') and not line.startswith('#'):
            if '=' not in line and '::' not in line and '*' not in line:
                cohorts.append(line.strip())
            elif '=' in line:
                if len(line.strip().split('='))==2:
                    cohort_name = line.strip().split('=')[0]
                    cohort_samples_file = line.strip().split('=')[1]
                    cohort_samples = []
                    if os.path.exists(cohort_samples_file):
                        with open(cohort_samples_file, 'r') as cohort_samples_readfile:
                            samples_lines = cohort_samples_readfile.readlines()
                            for sample_line in samples_lines:
                                if sample_line!='' and len(sample_line.split('\t'))==2:
                                    cohort_samples.append(sample_line.strip().split('\t')[0])
                        cohorts.append(cohort_name + '=' + ','.join(cohort_samples))
            elif '::' in line:
                if len(line.strip().split('::'))==2:
                    cohort_name = line.strip().split('::')[0]
                    cohort_keyword_to_exlude = line.strip().split('::')[1].split('_')[1]#::NOT_Lymph-
                cohorts.append(cohort_name + '::' + cohort_keyword_to_exlude)
            elif '*' in line:
                if len(line.strip().split('*'))==2:
                    cohort_name = line.strip().split('*')[0]
                    cohort_value_name = line.strip().split('*')[1]
                cohorts.append(cohort_name + '::' + cohort_value_name)
    return cohorts

def get_sig_merged_elements(unified_mutation_input_files, cohort_full_name, 
                            output_extension, distance_to_merge, 
                            merged_mut_sig_threshold, local_domain_window, 
                            chr_lengths_file, sig_elements_output_file, 
                            sim_sig_thresh, p_value_on_score):
    
    if os.path.exists(sig_elements_output_file):
        return sig_elements_output_file
    
    
    
    combined_simulated_muts_merged_output_file = cohort_full_name + output_extension + '_unified_combined' + '_merged{distance_to_merge}bp'.format(distance_to_merge=distance_to_merge) + '_combined'
#     merged_muts_output_ext = "_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)
# #    "Merge the simluated mutations into elements"
#     merged_simulated_element_files = []
#     for unified_muts_file_wihtmotifinfo in unified_mutation_input_files[1:]:
#         merged_muts_output_file = unified_muts_file_wihtmotifinfo+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)
#         Utilities.merge_muts(muts_input_file=unified_muts_file_wihtmotifinfo, 
#                               merged_muts_output_ext=merged_muts_output_ext, 
#                               filter_mut_motifs=False, filter_col_index=16, 
#                               filter_value=0.05, mut_score_index=9, 
#                               motifs_col_index =10, ref_alt_col_index=11, 
#                               mutpos_col_index=12, motifname_col_index=13, 
#                               motif_col_index=14, distance_to_merge=distance_to_merge)
#         merged_simulated_element_files.append(merged_muts_output_file)
     
    
    merged_simulated_element_files = []
    merged_muts_output_ext = "_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)
    filter_mut_motifs=False
    filter_col_index=16
    filter_value=0.05
    mut_score_index=9
    motifs_col_index =10
    ref_alt_col_index=11 
    mutpos_col_index=12
    motifname_col_index=13 
    motif_col_index=14
    distance_to_merge=distance_to_merge

    pm = Pool(15)
    merged_simulated_element_files = pm.starmap(Utilities.merge_muts, product(
        unified_mutation_input_files[1:], [merged_muts_output_ext], [filter_mut_motifs], 
        [filter_col_index], [filter_value], [mut_score_index], [motifs_col_index], [ref_alt_col_index], [mutpos_col_index], 
        [motifname_col_index], [motif_col_index], [distance_to_merge]))
    pm.close()
    pm.join()
    
        
    if not os.path.exists(combined_simulated_muts_merged_output_file):
        with open(combined_simulated_muts_merged_output_file, 'w') as combined_simulated_muts_merged_outfile:
            for merged_muts_output_file in merged_simulated_element_files:
                with open(merged_muts_output_file, 'r') as merged_muts_read_file:
                    combined_simulated_muts_merged_outfile.write(merged_muts_read_file.read())
    
    "Merge the observed mutations into elements"
    unified_observed_muts_file_wihtmotifinfo = unified_mutation_input_files[0]
    #merged_muts_output_file = unified_observed_muts_file_wihtmotifinfo+"_mergedmuts{distance_to_merge}bp".format(
    #                                        distance_to_merge=distance_to_merge)
    merged_muts_output_file= Utilities.merge_muts(muts_input_file=unified_observed_muts_file_wihtmotifinfo, 
                         merged_muts_output_ext=merged_muts_output_ext, 
                         filter_mut_motifs=False, filter_col_index=15, 
                         filter_value=merged_mut_sig_threshold, mut_score_index=9, 
                         motifs_col_index =10, ref_alt_col_index=11, mutpos_col_index=12, 
                         motifname_col_index=13, motif_col_index=14, 
                         distance_to_merge=distance_to_merge)
    
    '''Calcuate pval for each element by comparing its score to 
       the std and mean of scores in all elements across the genome
    '''
    merged_elements_statspvalues = merged_muts_output_file+"_statspvaluesSimSig" + str(sim_sig_thresh)
    merged_elements_statspvaluesonlysig = merged_muts_output_file+"_statspvaluesSimSig{simulated_mut_motif_sig_val}onlysig{merged_mut_sig_threshold}".format(simulated_mut_motif_sig_val=sim_sig_thresh, merged_mut_sig_threshold=merged_mut_sig_threshold)
#     Utilities.assess_stat_elements(observed_input_file=merged_muts_output_file, 
#                                    simulated_input_file=combined_simulated_muts_merged_output_file, 
#                                    merged_elements_statspvalues=merged_elements_statspvalues, 
#                                    merged_elements_statspvaluesonlysig=merged_elements_statspvaluesonlysig, 
#                                    merged_mut_sig_threshold=merged_mut_sig_threshold, 
#                                    score_index_observed_elements=3, score_index_sim_elements=3, p_value_on_score=p_value_on_score)
    
    
    #tmp to omit comparison to the background
    with open(merged_muts_output_file, 'r') as observed_infile, open(merged_elements_statspvalues, 'w') as merged_elements_statspvalues_outfile:
        l = observed_infile.readline()
        while l:
            sl = l.strip().split('\t')
            #print(dict_pvals[l_number])
            #print(l.strip())
            merged_elements_statspvalues_outfile.write(l.strip() + '\t' + '1.0' + '\t' + '1.0' + '\n')
            l = observed_infile.readline()
    
    #based on a local domain distribution of scores
    '''Calcuate pval for each element by comparing its score to 
       the std and mean of scores in elements within local_domain_window
    '''
    merged_elements_statspvalues_local = merged_elements_statspvalues+"_statspvalueslocalw{local_domain_window}".format(local_domain_window=local_domain_window)
    sig_elements_output_file= merged_elements_statspvalues_local+"onlysig"+str(merged_mut_sig_threshold)
    #sig_elements_output_file = merged_elements_statspvalues+"_statspvalueslocalw{local_domain_window}onlysig{merged_mut_sig_threshold}".format(local_domain_window=local_domain_window, merged_mut_sig_threshold=merged_mut_sig_threshold)
    Utilities.assess_stat_elements_local_domain(
        observed_input_file=merged_elements_statspvalues, 
        simulated_input_file=combined_simulated_muts_merged_output_file, 
        merged_elements_statspvalues=merged_elements_statspvalues_local, 
        merged_elements_statspvaluesonlysig=sig_elements_output_file, 
        chr_lengths_file=chr_lengths_file, local_domain_window=local_domain_window, 
        merged_mut_sig_threshold=merged_mut_sig_threshold, 
        score_index_observed_elements=3, score_index_sim_elements=3, p_value_on_score=p_value_on_score)
   
    #if os.path.exists(combined_simulated_muts_merged_output_file):
    #    os.remove(combined_simulated_muts_merged_output_file)
    
    return sig_elements_output_file




def run_cohort(cohort, created_cohorts, mutation_input_files, mutations_cohorts_dir, motif_name_index, 
               f_score_index, motif_breaking_score_index,chromatin_cat_index,
               background_window, background_window_size, 
               filter_on_qval, sig_category, sig_thresh, sim_sig_thresh,
               sim_output_extension,
               filter_cond, operation_on_unify, output_extension, 
               distance_to_merge, merged_mut_sig_threshold,
               local_domain_window, chr_lengths_file,
               sig_elements_output_file, sig_tfs_file, sig_tfpos_file, tmp_dir, n_cores_fscore, p_value_on_score, active_driver_script_dir,
               active_driver_min_mut, n_cores):    
    
    "get the cohort name to use for output file names"
    cohort_full_name = created_cohorts[cohort][0].split('_')[0]
    if '/' in created_cohorts[cohort][0]:
        cohort_full_name = '/'.join(created_cohorts[cohort][0].split('/')[0:-1]) + "/" + created_cohorts[cohort][0].split('/')[-1].split('_')[0]
    
    print('Processing: ', cohort_full_name)
    
    '''Calculate std, nummotifs and mean of the scores per TF-motif in the simulation sets
       The first file in each created_cohorts[cohort] is the observed set, so skip it
    '''
    #dict_simulated_mean_sd_per_TF_motif_output_file = cohort_full_name + "_meansdrand{}sets.dict".format(len(mutation_input_files)-1)
    #print(dict_simulated_mean_sd_per_TF_motif_output_file)
    
    '''As background consider simulated mutations in provided bacground window size around mutation
    '''

    
#     if background_window:
#         #print(background_window_size)
#         dict_type_mean_std_scores = Utilities.get_simulated_mean_sd_per_TF_motif_background_window(
#             cohort_full_name = cohort_full_name,
#             annotated_input_file = created_cohorts[cohort][0],
#             simulated_annotated_input_files=created_cohorts[cohort][1:],
#             mutations_cohorts_dir = mutations_cohorts_dir,
#             cohort_mean_sd_per_tf_overall_output_dict_file= dict_simulated_mean_sd_per_TF_motif_output_file, 
#             chr_lengths_file = chr_lengths_file,
#             background_window_size = background_window_size, 
#             motif_name_index = motif_name_index, f_score_index = f_score_index, 
#             motif_breaking_score_index = motif_breaking_score_index,
#             chromatin_cat_index = chromatin_cat_index, tmp_dir = tmp_dir, n_cores_fscore=n_cores_fscore)
#     else:
#         '''As background consider whole genome
#         '''
#         dict_type_mean_std_scores = Utilities.get_simulated_mean_sd_per_TF_motif(
#             simulated_annotated_input_files=created_cohorts[cohort][1:], 
#             cohort_mean_sd_per_tf_overall_output_dict_file= dict_simulated_mean_sd_per_TF_motif_output_file, 
#             motif_name_index = motif_name_index, f_score_index = f_score_index, 
#             motif_breaking_score_index = motif_breaking_score_index)
#     
#     '''For each mutation in the observed set created_cohorts[cohort][0]
#        calculate pval and qval by comparing its score to the std and mean
#        scores in the corresponding TF motif.
#        Filter out mutations that don't have a sig score or dont' pass other filters
#     '''
#     muts_sig_per_TF_file = Utilities.get_muts_sig_per_TF(
#         annoted_input_file=created_cohorts[cohort][0], 
#         dict_type_mean_std_scores=dict_type_mean_std_scores, 
#         annoted_output_file_extension="_rand{}setsTF".format(len(mutation_input_files)-1), 
#         annoted_output_file_extension_onlysig="_rand{}setsTFsigQval{}".format(
#             len(mutation_input_files)-1, sig_thresh),
#         background_window = background_window,
#         motif_name_index = motif_name_index, f_score_index = f_score_index, 
#         motif_breaking_score_index = motif_breaking_score_index, 
#         filter_on_qval=filter_on_qval, sig_cat=sig_category, 
#         sig_thresh=sig_thresh,
#         filter_on_signal = True, dnase_index = 24, fantom_index = 25, num_other_tfs_index = 27)
#     sig_muts_per_tf_mutation_input_files = [muts_sig_per_TF_file]
#     
#     '''repeat the same process to keep only sig muts from the simulated, but use sim_sig_thresh
#     if sim_sig_thresh >=1.0 no sig is performed and all muts are written to the output'''
#     for mutations_input_file in created_cohorts[cohort][1:]: 
#         muts_sig_per_TF_file = Utilities.get_muts_sig_per_TF(
#             annoted_input_file=mutations_input_file, 
#             dict_type_mean_std_scores=dict_type_mean_std_scores,
#             annoted_output_file_extension="_rand{}setsTF".format(len(mutation_input_files)-1), 
#             annoted_output_file_extension_onlysig=sim_output_extension,
#             background_window = background_window,
#             motif_name_index = motif_name_index, f_score_index = f_score_index, 
#             motif_breaking_score_index = motif_breaking_score_index,
#             filter_on_qval=filter_on_qval, sig_cat=sig_category, 
#             sig_thresh=sim_sig_thresh
#             )
#         sig_muts_per_tf_mutation_input_files.append(muts_sig_per_TF_file)
#     
    '''Based on the mutations that have a significant score (specify if qval should be used)
       Count number of mutations per TF motif and per TF motif position
       For each calcualate a pvalue and qvalue based on number of mutations
       in the correponding motif in the simulated sets  
    '''
    #sig_tfs_file, sig_tfpos_file = Utilities.get_tf_pval(
    #    cohort, sig_muts_per_tf_mutation_input_files, p_value_on_score, motif_name_index, 
    #            f_score_index, motif_breaking_score_index, 
    #            filter_cond, fsep='\t', sig_tfs_file=sig_tfs_file, 
    #            sig_tfpos_file=sig_tfpos_file,
    #            filter_on_signal = True, dnase_index = 24, fantom_index = 25, 
   #             num_other_tfs_index = 27)
   
    #
    
    '''replace an observed regulatory mutation file with an all annotated observed mutation file
    '''
    #muts_per_tf_mutation_input_files  = [created_cohorts[cohort][0]+"_rand{}setsTF".format(len(mutation_input_files)-1)]
    
    #muts_per_tf_mutation_input_files.extend(sig_muts_per_tf_mutation_input_files[1:])
    #print(muts_per_tf_mutation_input_files)
    #print(len(muts_per_tf_mutation_input_files))
    
    '''Combine nearby mutations accross the cohort into one element'''
        
    
        
     
    '''Unify the mutations that have significant scores accross the cohorts
       Make one record for mutations that overlap multiple motifs
    '''
    
    unified_muts_files_obj=[]
    p=Pool(15)
    for mutations_input_file in created_cohorts[cohort]:
        unified_muts_file=p.apply_async(Utilities.unify_muts, args=(mutations_input_file,output_extension, True, filter_cond, operation_on_unify))
        unified_muts_files_obj.append(unified_muts_file)
    p.close()
    p.join()
    unified_muts_files = [p.get() for p in unified_muts_files_obj]
    
    
    unified_mutation_input_files_obj = []
    p=Pool(15)
    for unified_muts_file in unified_muts_files:
        unified_muts_file_wihtmotifinfo=p.apply_async(Utilities.get_max_motif_in_grouped_muts, args=(unified_muts_file,))
        unified_mutation_input_files_obj.append(unified_muts_file_wihtmotifinfo)

    
    p.close()
    p.join()
    unified_mutation_input_files = [p.get() for p in unified_mutation_input_files_obj]
    
    #for unified_muts_file in unified_muts_files:
    #        if os.path.exists(unified_muts_file):
    #            os.remove(unified_muts_file)
    
#     unified_mutation_input_files = []
#     for mutations_input_file in created_cohorts[cohort]:
#         unified_muts_file = mutations_input_file + output_extension + "_groupedbymut" 
#         unified_muts_file_wihtmotifinfo = unified_muts_file+"withmotifinfo"
#         if not os.path.exists(unified_muts_file_wihtmotifinfo):
#             print("Unifying: ", mutations_input_file)
#             Utilities.unify_muts(mutations_input_file, unified_muts_file_ext, 
#                                  filter_mut_motifs=True, filter_cond=filter_cond, 
#                                  operation_on_unify=operation_on_unify)
#             Utilities.get_max_motif_in_grouped_muts(
#                 annotated_mutations_grouped_file=unified_muts_file, 
#                 annotated_mutations_grouped_output_file=unified_muts_file_wihtmotifinfo)
#             os.remove(unified_muts_file)
#         unified_mutation_input_files.append(unified_muts_file_wihtmotifinfo)
#    print('Unified mutations input files: ', unified_mutation_input_files)
    '''   Evaluate the significance of each element based on: 
       - the element score (sum of the score of its mutations)
       - number of mutations in the element 
    '''
    get_sig_merged_elements(unified_mutation_input_files, cohort_full_name, 
                            output_extension,
                            #sim_output_extension+output_extension, 
                            distance_to_merge, merged_mut_sig_threshold, 
                            local_domain_window, chr_lengths_file, 
                                sig_elements_output_file, sim_sig_thresh, p_value_on_score=p_value_on_score)
    
    
    
    
    #ActiveDriverWGS
    
    active_driver_results = sig_elements_output_file + '_ActiveDriver_results'
    active_driver_results_sig = active_driver_results + '_sig'
    active_driver_output_file = sig_elements_output_file + '_ActiveDriver'
    active_driver_output_file_sig = active_driver_output_file + '_sig'
    
    if not os.path.exists(active_driver_output_file_sig):
            print(['Rscript', active_driver_script_dir, sig_elements_output_file,  created_cohorts[cohort][0], active_driver_min_mut, active_driver_output_file, active_driver_output_file_sig, active_driver_results,  n_cores])
                 
            try:
                subprocess.call(['Rscript', active_driver_script_dir, sig_elements_output_file,  created_cohorts[cohort][0], str(active_driver_min_mut), active_driver_output_file, active_driver_output_file_sig,  active_driver_results, str(n_cores)])
            except KeyError:
                open(active_driver_output_file_sig, 'a').close()
         
    
    sig_muts_file=created_cohorts[cohort][0]+"_sig"
    
    sig_muts_file = Utilities.get_sig_muts(
        elements_input_file=active_driver_output_file_sig,
        mutations_input_file=created_cohorts[cohort][0],
        sig_muts_file=sig_muts_file,
        motif_breaking_score_index=motif_breaking_score_index,
        tf_binding_index=30, 
        dnase_index=24)
    
        
         
    sig_muts_per_tf_mutation_input_files=[]
    sig_muts_per_tf_mutation_input_files=[sig_muts_file]
    for  mutations_input_file in created_cohorts[cohort][1:]: 
       sig_muts_per_tf_mutation_input_files.append(mutations_input_file) 
             
             
             
    sig_tfs_file, sig_tfpos_file = Utilities.get_tf_pval(
        cohort, sig_muts_per_tf_mutation_input_files, p_value_on_score, motif_name_index, 
                f_score_index, motif_breaking_score_index, 
               filter_cond, fsep='\t', sig_tfs_file=sig_tfs_file, 
                sig_tfpos_file=sig_tfpos_file,
                filter_on_signal = True, dnase_index = 24, fantom_index = 25, 
                num_other_tfs_index = 27)
    

    
    return sig_elements_output_file, active_driver_output_file_sig, sig_tfs_file, sig_tfpos_file
    
def process_cohorts(cohort_names_input, mutations_cohorts_dir, 
                    observed_input_file, simulated_input_dir,
                    chr_lengths_file, num_cores, 
                    background_window, background_window_size, 
                    filter_on_qval, sig_category, sig_thresh, sim_sig_thresh_pval,
                    distance_to_merge, 
                    merged_mut_sig_threshold, local_domain_window, tmp_dir, n_cores_fscore, p_value_on_score, active_driver_script_dir,
                    active_driver_min_mut, n_cores):
    
    simulated_input_files = [simulated_input_dir+'/'+x for x in os.listdir(simulated_input_dir) if '_annotated.bed9' in x]
    mutation_input_files = [observed_input_file]
    mutation_input_files.extend(simulated_input_files)
    print(mutation_input_files)
    motif_name_index = 17
    f_score_index = 9
    motif_breaking_score_index = 10
    chromatin_cat_index = 22
    sim_output_extension = "_rand{}setsTFsigPval{}".format(len(mutation_input_files)-1, sim_sig_thresh_pval)#len(mutation_input_files)-1 
    operation_on_unify = 'mean'#'max'
    
    output_extension = "_{operation_on_unify}TFExprMotifBreaking03Filters".format(operation_on_unify=operation_on_unify)
    
    "filter on gene expression level of the corresponding TF ($32) and on motif-break score (>=0.3)"
    filter_cond = 'if(($32>0 || $32~"nan") && $11>=0.3)'# && ($31>0 || $31~"nan")
    
    #output_extension = "_{operation_on_unify}TFExprTFBindMotifBreaking03Filters".format(operation_on_unify=operation_on_unify)
    #filter_cond = 'if(($31>0 || $31~"nan")  && ($32>0 || $32~"nan") && $11>=0.3)'
    
    generated_sig_merged_element_files = []
    generated_sig_merged_elements_active_files = []
    sig_tfs_files = []
    sig_tfpos_files = []
    if num_cores>1:
        p = Pool(num_cores)
    
    stats_ext = "_rand{}setsTF".format(len(mutation_input_files)-1)
    
    '''
    # to make it possible re-run later steps without requiring earlier ones when premilinary files are removed to save space.
    cohorts_initial = get_cohorts(cohort_names_input)
    cohorts_initial = [x.split("=")[0] for x in cohorts_initial]
    cohorts = cohorts_initial[0::]
    print("Cohorts:", cohorts)
    for i, cohort in enumerate(cohorts_initial):
        sig_elements_output_file = mutations_cohorts_dir + '/' + cohort + "_{}_rand{}setsTFsigQval{}".format(observed_input_file.split('/')[-1], len(mutation_input_files)-1, sig_thresh_pval) + output_extension + "_groupedbymut"+"withmotifinfo"+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)+"_statspvaluesSimSig"+str(sim_sig_thresh_pval)+"_statspvalueslocalw{local_domain_window}onlysig{merged_mut_sig_threshold}".format(local_domain_window=local_domain_window, merged_mut_sig_threshold=merged_mut_sig_threshold)
        sig_tfs_file = mutations_cohorts_dir + '/' + cohort + "_{}_rand{}setsTFsigQval{}".format(observed_input_file.split('/')[-1], len(mutation_input_files)-1, sig_thresh_pval) + '_sigTFs_{}'.format(sig_thresh_pval) 
        sig_tfpos_file = mutations_cohorts_dir + '/' + cohort + "_{}_rand{}setsTFsigQval{}".format(observed_input_file.split('/')[-1], len(mutation_input_files)-1, sig_thresh_pval) + '_sigTFpos_{}'.format(sig_thresh_pval)
        print(sig_elements_output_file)
        if os.path.exists(sig_elements_output_file):
            generated_sig_merged_element_files.append(sig_elements_output_file)
            sig_tfs_files.append(sig_tfs_file)
            sig_tfpos_files.append(sig_tfpos_file)
            print("Exists: ", sig_elements_output_file)
            del cohorts[cohorts.index(cohort)]
    '''
    cohorts = get_cohorts(cohort_names_input)
    print("Cohorts:", cohorts)
    print('process cohort')
    created_cohorts = generate_cohorts(mutation_input_files, cohorts, 
                                       mutations_cohorts_dir, stats_ext,
                                       num_cores)
    for cohort in created_cohorts.keys():
        sig_elements_output_file = mutations_cohorts_dir + '/' + cohort + "_{}".format(observed_input_file.split('/')[-1]) + output_extension + "_groupedbymut"+"withmotifinfo"+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)+"_statspvaluesSimSig"+str(sim_sig_thresh_pval)+"_statspvalueslocalw{local_domain_window}onlysig{merged_mut_sig_threshold}".format(local_domain_window=local_domain_window, merged_mut_sig_threshold=merged_mut_sig_threshold)
        #sig_elements_output_file = created_cohorts[cohort][0] + "_rand{}setsTFsigQval{}".format(len(mutation_input_files)-1, sig_thresh) + output_extension + "_groupedbymut"+"withmotifinfo"+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)+"_statspvaluesonlysig{merged_mut_sig_threshold}".format(
        #    merged_mut_sig_threshold=merged_mut_sig_threshold)
        sig_elements_active_output_file=sig_elements_output_file+'_ActiveDriver_sig'
        sig_tfs_file = created_cohorts[cohort][0] + "_rand{}setsTFsigQval{}".format(len(mutation_input_files)-1, sig_thresh) + '_sigTFs_mean_sd_onlySim_{}'.format(sig_thresh) 
        sig_tfpos_file = created_cohorts[cohort][0] + "_rand{}setsTFsigQval{}".format(len(mutation_input_files)-1, sig_thresh) + '_sigTFpos_min_mean_sd_onlySim_{}'.format(sig_thresh)
        if os.path.exists(sig_elements_output_file) and os.path.exists(sig_tfs_file) and os.path.exists(sig_tfpos_file) and os.path.exists(sig_elements_active_output_file):
            generated_sig_merged_element_files.append(sig_elements_output_file)
            generated_sig_merged_elements_active_files.append(sig_elements_active_output_file)
            sig_tfs_files.append(sig_tfs_file)
            sig_tfpos_files.append(sig_tfpos_file)
            continue
        if len(created_cohorts[cohort])<2:
            print('WARNING: no observed/simulated mutations file for ', cohort, created_cohorts[cohort])
            continue
        
        if num_cores>1:
            p.apply_async(run_cohort, args=(cohort, created_cohorts, 
                    mutation_input_files, mutations_cohorts_dir, motif_name_index, 
                    f_score_index, motif_breaking_score_index, chromatin_cat_index,
                    background_window, background_window_size,
                    filter_on_qval, sig_category, sig_thresh, sim_sig_thresh_pval,
                    sim_output_extension,
                    filter_cond, operation_on_unify, output_extension, 
                    distance_to_merge, merged_mut_sig_threshold,
                    local_domain_window, chr_lengths_file, sig_elements_output_file, 
                    sig_tfs_file, sig_tfpos_file, tmp_dir, n_cores_fscore, p_value_on_score,
                    active_driver_script_dir, active_driver_min_mut, n_cores))#, callback=generated_sig_merged_element_files.append)
        else:
            run_cohort(cohort, created_cohorts, 
                       mutation_input_files, mutations_cohorts_dir, motif_name_index, 
                       f_score_index, motif_breaking_score_index, chromatin_cat_index,
                       background_window, background_window_size, 
                       filter_on_qval, sig_category, sig_thresh, sim_sig_thresh_pval,
                       sim_output_extension,
                       filter_cond, operation_on_unify, output_extension, 
                       distance_to_merge, merged_mut_sig_threshold,
                       local_domain_window, chr_lengths_file, sig_elements_output_file, 
                       sig_tfs_file, sig_tfpos_file, tmp_dir, n_cores_fscore, p_value_on_score,
                       active_driver_script_dir, active_driver_min_mut, n_cores)
        
        generated_sig_merged_element_files.append(sig_elements_output_file)
        generated_sig_merged_elements_active_files.append(sig_elements_active_output_file)

        sig_tfs_files.append(sig_tfs_file)
        sig_tfpos_files.append(sig_tfpos_file)
        print("Generated for cohort ({}): {}".format(cohort, sig_elements_output_file))
        
    if num_cores>1:
        p.close()
        p.join()
    
    return generated_sig_merged_element_files, generated_sig_merged_elements_active_files, sig_tfs_files, sig_tfpos_files
    
#awk 'BEGIN{FS=OFS="\t"}{if($18<0.05) sig="element_sig"; else sig="element_notsig";  print $4,"Observed", sig}' mutations_cohorts_output/All-tumors_observed_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmutsonlysig0.1_mergedmuts20bp_statspvalues >> ../analysis/All-tumors_6setsobs_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_scores.col3

#awk 'BEGIN{FS=OFS="\t"}{sig="mut_sig"; if($16<0.1) sig="mut_sig"; else sig="mut_notsig";  print $10,"Rand3", sig}' All-tumors_randomised98f100kw50kbnonparametric_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts >> ../analysis/All-tumors_6setsobs_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_scores.col3
def parse_args():
    '''Parse command line arguments'''
    print('parse')
    parser = argparse.ArgumentParser(description='Process PCAWG Cohorts')
    parser.add_argument('-c', '--cohort_names_input', default='', help='')
    parser.add_argument('-o', '--mutations_cohorts_outdir', default='', help='')
    parser.add_argument('-m', '--observed_input_file', default='', help='')
    parser.add_argument('-s', '--simulated_input_dir', default='', help='')
    parser.add_argument('-l', '--chr_lengths_file', default='', help='')

    parser.add_argument('--background_window', action='store_const', const=True, help='Check mutation functional score significance by comparing to background window around mutation in simulated mutations, if the flag is missing it would use the whole genome as background')
    parser.add_argument('--background_window_size', type=int, default=50000, help='Background window around mutation for capturing simulated mutation to compare mutation functional score')
    parser.add_argument('--sig_thresh', type=float, default=0.05, help='Sig level threshold on mutation score level')
    parser.add_argument('--sim_sig_thresh', type=float, default=1.0, help='Sig level threshold for simulated mutations on score level')
    parser.add_argument('--merged_mut_sig_threshold', type=float, default=0.05, help='P-value threshold for simulated mutations on score level')
    parser.add_argument('--p_value_on_score', action='store_const', const=True, help='P-value of elements computed based on the score distribution, if the flag is missing it would use the probability density function to obtain p-value')
    parser.add_argument('--distance_to_merge', type=int, default=200, help='Window size (number of base-pairs) to merge nearby mutations within')
    parser.add_argument('--local_domain_window', type=int, default=25000, help='Window width for capturing simulated elements to compare mutation frequency ')
    parser.add_argument('--filter_on_qval', action='store_const', const=True, help='Filter on FDR (adjusted p-values), if the flag is missing it would filter on p-value')
    parser.add_argument('--sig_category', default = 'overallTFs', choices=['overallTFs', 'perTF', 'perChromatinCat', 'perTF_perChromatinCat'], help='')
    parser.add_argument('--num_cores', type=int, default=10, help='number of cores (cpus) to use in parallel')
    parser.add_argument('--tmp_dir', default='$SNIC_TMP', help='')
    parser.add_argument('--n_cores_fscore', type=int, default=10, help='number of cores (cpus) to use in parallel to obtain fscore')
    parser.add_argument('--active_driver_script_dir', default='/proj/snic2020-16-50/nobackup/pancananalysis/pancancer_analysis_KS/pancancer_analysis/ActiveDriverWGS_to_run.R', help='')
    parser.add_argument('--active_driver_min_mut', default=1, help='n')
    parser.add_argument('--num_cores_activedriver',  default=15, help='Number of cores to run ActiveDriverWGS')

    
    return parser.parse_args(sys.argv[1:])
    
if __name__ == '__main__':
    
    

    args = parse_args()
    
    temp_dir = args.tmp_dir
    set_tempdir(temp_dir)
    
    if not os.path.exists(args.mutations_cohorts_outdir):
        os.makedirs(args.mutations_cohorts_outdir)
    
    generated_sig_merged_element_files, generated_sig_merged_elements_active_files, sig_tfs_files, sig_tfpos_files = process_cohorts(
        args.cohort_names_input, args.mutations_cohorts_outdir, args.observed_input_file, 
        args.simulated_input_dir, args.chr_lengths_file, args.num_cores, 
        args.background_window, args.background_window_size, 
        args.filter_on_qval, args.sig_category, args.sig_thresh, args.sim_sig_thresh,  
        args.distance_to_merge, args.merged_mut_sig_threshold,
        args.local_domain_window, args.tmp_dir, args.n_cores_fscore, args.p_value_on_score,
        args.active_driver_script_dir, args.active_driver_min_mut, args.num_cores_activedriver)
    #print("Generated Sig. Element Sets: \n", '\n'.join(generated_sig_merged_element_files))
    
    