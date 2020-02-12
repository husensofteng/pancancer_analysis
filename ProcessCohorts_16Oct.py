'''
Created on Feb 9, 2017

@author: husensofteng
'''
import os
import sys
from multiprocessing import Pool

from Utilities import unify_muts, get_max_motif_in_grouped_muts, merge_muts, assess_stat_muts, assess_stat_elements, get_mean_and_sd_from_file, get_muts_sig_per_TF, get_simulated_mean_sd_per_TF_motif, assess_stat_elements_local_domain, get_tf_pval

def generate_cohorts(mutations_input_files, cohorts, mutations_cohorts_dir, run_in_parallel=False, num_process=8, stats_ext=""):
    created_cohorts_dict = {}
    if run_in_parallel:
        p = Pool(num_process)
        
    for mutations_input_file in mutations_input_files:
        for cohort_value in cohorts:
            cohort_name = cohort_value.split('=')[0].split('::')[0]
            if cohort_name not in created_cohorts_dict.keys():
                created_cohorts_dict[cohort_name] = []
            cohort_file = mutations_cohorts_dir+'/'+cohort_value.split('=')[0].split('::')[0] + '_' + mutations_input_file.split('/')[-1]
            cohort_file_statpvalues = cohort_file + stats_ext
            #print cohort_file
            created_cohorts_dict[cohort_name].append(cohort_file)
            #print cohort_file_statpvalues
            if os.path.exists(cohort_file) or os.path.exists(cohort_file_statpvalues):
                continue
            print "Generating {}...".format(cohort_file)
            if run_in_parallel:
                p.apply_async(generate_cohorts_per_mut_file, args = (cohort_value, cohort_file, mutations_input_file))
            else:
                generate_cohorts_per_mut_file(cohort_value, cohort_file, mutations_input_file)
            #read the mutations file line by line and check if the sample id is inthe lis of samples
    if run_in_parallel:
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
                with open(cohort_value.split('=')[0] + "_mutations_temp", 'w') as mutations_temp_outfile:
                    line = mutations_infile.readline()
                    while line!="":
                        if line.strip().split('\t')[8].strip() in cohort_sample_ids or line.strip().split('\t')[7].strip() in cohort_sample_ids:
                            mutations_temp_outfile.write(line)
                            n+=1
                        line = mutations_infile.readline()
        print "written lines to {}: {}".format(cohort_file, n)
        awk_stm = """awk 'BEGIN{FS=OFS="\t"}{print $0 >> "%s"}' %s""" %(cohort_file, cohort_value.split('=')[0] + "_mutations_temp")
        os.system(awk_stm)
        if os.path.exists(cohort_value.split('=')[0] + "_mutations_temp"):
            os.remove(cohort_value.split('=')[0] + "_mutations_temp")            
    
    return cohort_file

def get_cohorts(cohort_names_input):
    cohorts = []
    lines = [cohort_names_input]
    if os.path.isfile(cohort_names_input) and os.path.exists(cohort_names_input):
        with open(cohort_names_input, 'r') as cohort_names_infile:
            lines = cohort_names_infile.readlines()
    for line in lines:
        if line!='' and not line.startswith('//') and not line.startswith('#'):
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

"""
def run_cohort_old(cohort, created_cohorts,output_extension, distance_to_merge, mut_sig_threshold):
    
    print 'Processing: ', cohort
    
    combined_simulated_muts_output_file = '/'.join(created_cohorts[cohort][0].split('/')[0:-1]) + '/' + cohort + output_extension + '_unified_combined'
    combined_simulated_muts_merged_output_file = '/'.join(created_cohorts[cohort][0].split('/')[0:-1]) + '/' + cohort + output_extension + '_unifiedsig{mut_sig_threshold}_merged{distance_to_merge}bp'.format(mut_sig_threshold=mut_sig_threshold, distance_to_merge=distance_to_merge) + '_combined'
    stats_dict = {}
    if os.path.exists(combined_simulated_muts_output_file) and os.path.exists(combined_simulated_muts_merged_output_file):
        stats_dict = get_mean_and_sd_from_file(combined_simulated_muts_output_file, scores_index=9, index_mut_type = 6)
    else:
        if not os.path.exists(combined_simulated_muts_output_file):
            with open(combined_simulated_muts_output_file, 'w') as combined_simulated_muts_outfile:
                #combine all unified-simulated sets
                #if not os.path.exists(combined_simulated_muts_output_file):
                for unified_muts_file_wihtmotifinfo in created_cohorts[cohort][1:]:
                    with open(unified_muts_file_wihtmotifinfo, 'r') as unified_muts_readfile:
                        combined_simulated_muts_outfile.write(unified_muts_readfile.read())
                #get sig-muts for each unified-simulated set using the combined set and then merge the sig-mut sets for each simulated set             
        stats_dict = get_mean_and_sd_from_file(combined_simulated_muts_output_file, scores_index=9, index_mut_type = 6)
        
        if not os.path.exists(combined_simulated_muts_merged_output_file):
            with open(combined_simulated_muts_merged_output_file, 'w') as combined_simulated_muts_merged_outfile:
                for unified_muts_file_wihtmotifinfo in created_cohorts[cohort][1:]:
                    calculated_pvalues_unified_muts_file_wihtmotifinfo = unified_muts_file_wihtmotifinfo+"_statmuts{}".format(mut_sig_threshold)
                    sig_calculated_pvalues_unified_muts_file_wihtmotifinfo = calculated_pvalues_unified_muts_file_wihtmotifinfo+"onlysig"
                    calculated_pvalues_unified_muts_file_wihtmotifinfo, sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, n_sig_muts = assess_stat_muts(muts_input_file=unified_muts_file_wihtmotifinfo, simulated_input_file=combined_simulated_muts_output_file, observed_output_file=calculated_pvalues_unified_muts_file_wihtmotifinfo, observed_onlysig_output_file=sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, score_index_observed_elements=9, score_index_sim_elements=9, mut_sig_threshold=mut_sig_threshold, stats_dict=stats_dict)
                    
                    merged_muts_output_file = sig_calculated_pvalues_unified_muts_file_wihtmotifinfo+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)
                    
                    merged_muts_output_file = merge_muts(muts_input_file=sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, merged_muts_output_file=merged_muts_output_file, filter_mut_motifs=False, filter_col_index=16, filter_value=0.05, mut_score_index=9, motifs_col_index =10, ref_alt_col_index=11, mutpos_col_index=12, motifname_col_index=13, motif_col_index=14, distance_to_merge=distance_to_merge)
                    with open(merged_muts_output_file, 'r') as merged_muts_read_file:
                        combined_simulated_muts_merged_outfile.write(merged_muts_read_file.read())
    
    unified_observed_muts_file = created_cohorts[cohort][0]
    calculated_pvalues_unified_muts_file_wihtmotifinfo = unified_observed_muts_file +"_statmuts{}".format(mut_sig_threshold)
    sig_calculated_pvalues_unified_muts_file_wihtmotifinfo = calculated_pvalues_unified_muts_file_wihtmotifinfo+"onlysig"
    calculated_pvalues_unified_muts_file_wihtmotifinfo, sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, n_sig_muts = assess_stat_muts(muts_input_file=unified_observed_muts_file, simulated_input_file=combined_simulated_muts_output_file, observed_output_file=calculated_pvalues_unified_muts_file_wihtmotifinfo, observed_onlysig_output_file=sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, score_index_observed_elements=9, score_index_sim_elements=9, mut_sig_threshold=mut_sig_threshold, stats_dict=stats_dict)
    
    merged_muts_output_file = calculated_pvalues_unified_muts_file_wihtmotifinfo+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)
    merged_muts_output_file = merge_muts(muts_input_file=calculated_pvalues_unified_muts_file_wihtmotifinfo, merged_muts_output_file=merged_muts_output_file, filter_mut_motifs=False, filter_col_index=15, filter_value=mut_sig_threshold, mut_score_index=9, motifs_col_index =10, ref_alt_col_index=11, mutpos_col_index=12, motifname_col_index=13, motif_col_index=14, distance_to_merge=distance_to_merge)
    merged_elements_statspvalues = merged_muts_output_file+"_statspvalues"
    merged_elements_statspvaluesonlysig = merged_muts_output_file+"_statspvaluesonlysig"
    merged_elements_statspvalues, merged_elements_statspvaluesonlysig, n_sig = assess_stat_elements(observed_input_file=merged_muts_output_file, simulated_input_file=combined_simulated_muts_merged_output_file, merged_elements_statspvalues=merged_elements_statspvalues, merged_elements_statspvaluesonlysig=merged_elements_statspvaluesonlysig, score_index_observed_elements=3, score_index_sim_elements=3)
    #print merged_elements_statspvaluesonlysig, '\t', n_sig, '\t', n_sig_muts
    #merged_elements_statspvaluesonlysig = created_cohorts[cohort][0]+"_statmutsonlysig"+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)+"_statspvaluesonlysig"
    return merged_elements_statspvaluesonlysig
    

"""
def get_sig_merged_elements(unified_mutation_input_files, cohort_full_name, output_extension, 
                            distance_to_merge, merged_mut_sig_threshold, local_domain_window, 
                            chr_lengths_file, sig_elements_output_file, simulated_mut_motif_sig_val):
    
    combined_simulated_muts_merged_output_file = cohort_full_name + output_extension + '_unified_combined' + '_merged{distance_to_merge}bp'.format(distance_to_merge=distance_to_merge) + '_combined'
    merged_simulated_element_files = []
    for unified_muts_file_wihtmotifinfo in unified_mutation_input_files[1:]:
        merged_muts_output_file = unified_muts_file_wihtmotifinfo+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)
        merge_muts(muts_input_file=unified_muts_file_wihtmotifinfo, merged_muts_output_file=merged_muts_output_file, filter_mut_motifs=False, filter_col_index=16, filter_value=0.05, mut_score_index=9, motifs_col_index =10, ref_alt_col_index=11, mutpos_col_index=12, motifname_col_index=13, motif_col_index=14, distance_to_merge=distance_to_merge)
        merged_simulated_element_files.append(merged_muts_output_file)
    if not os.path.exists(combined_simulated_muts_merged_output_file):
        with open(combined_simulated_muts_merged_output_file, 'w') as combined_simulated_muts_merged_outfile:
            for merged_muts_output_file in merged_simulated_element_files:
                with open(merged_muts_output_file, 'r') as merged_muts_read_file:
                    combined_simulated_muts_merged_outfile.write(merged_muts_read_file.read())
    
    unified_observed_muts_file_wihtmotifinfo = unified_mutation_input_files[0]
    merged_muts_output_file = unified_observed_muts_file_wihtmotifinfo+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)
    merge_muts(muts_input_file=unified_observed_muts_file_wihtmotifinfo, merged_muts_output_file=merged_muts_output_file, filter_mut_motifs=False, filter_col_index=15, filter_value=merged_mut_sig_threshold, mut_score_index=9, motifs_col_index =10, ref_alt_col_index=11, mutpos_col_index=12, motifname_col_index=13, motif_col_index=14, distance_to_merge=distance_to_merge)
    
    #based on genome wide distribution
    
    merged_elements_statspvalues = merged_muts_output_file+"_statspvaluesSimSig" + str(simulated_mut_motif_sig_val)
    merged_elements_statspvaluesonlysig = merged_muts_output_file+"_statspvaluesSimSig{simulated_mut_motif_sig_val}onlysig{merged_mut_sig_threshold}".format(simulated_mut_motif_sig_val=simulated_mut_motif_sig_val, merged_mut_sig_threshold=merged_mut_sig_threshold)
    assess_stat_elements(observed_input_file=merged_muts_output_file, simulated_input_file=combined_simulated_muts_merged_output_file, merged_elements_statspvalues=merged_elements_statspvalues, merged_elements_statspvaluesonlysig=merged_elements_statspvaluesonlysig, merged_mut_sig_threshold=merged_mut_sig_threshold, score_index_observed_elements=3, score_index_sim_elements=3)
    
    #based on a local domain distribution of scores
    merged_elements_statspvalues_local = merged_elements_statspvalues+"_statspvalueslocalw{local_domain_window}".format(local_domain_window=local_domain_window)
    #sig_elements_output_file = merged_elements_statspvalues+"_statspvalueslocalw{local_domain_window}onlysig{merged_mut_sig_threshold}".format(local_domain_window=local_domain_window, merged_mut_sig_threshold=merged_mut_sig_threshold)
    assess_stat_elements_local_domain(observed_input_file=merged_elements_statspvalues, simulated_input_files=merged_simulated_element_files, 
                         merged_elements_statspvalues=merged_elements_statspvalues_local, merged_elements_statspvaluesonlysig=sig_elements_output_file, 
                         chr_lengths_file=chr_lengths_file, local_domain_window=local_domain_window, 
                         merged_mut_sig_threshold=merged_mut_sig_threshold, score_index_observed_elements=3, score_index_sim_elements=3)
    if os.path.exists(combined_simulated_muts_merged_output_file):
        os.remove(combined_simulated_muts_merged_output_file)
    return sig_elements_output_file


def run_cohort(cohort, created_cohorts, mutation_input_files, motif_name_index, f_score_index, motif_breaking_score_index,
               sig_level_per_TF_thresh, sim_output_extension, sim_sig_level_per_TF_thresh,
               filter_cond, operation_on_unify, output_extension, distance_to_merge, merged_mut_sig_threshold,
               local_domain_window, chr_lengths_file,
               sig_elements_output_file, sig_tfs_file, sig_tfpos_file):    
    
    cohort_full_name = created_cohorts[cohort][0].split('_')[0]
    
    if '/' in created_cohorts[cohort][0]:
        cohort_full_name = '/'.join(created_cohorts[cohort][0].split('/')[0:-1]) + "/" + created_cohorts[cohort][0].split('/')[1].split('_')[0]
    print 'Processing: ', cohort_full_name
    
    dict_simulated_mean_sd_per_TF_motif_output_file = cohort_full_name + "_meansdrand{}sets.dict".format(len(mutation_input_files)-1)
    dict_simulated_mean_sd_per_TF_motif = get_simulated_mean_sd_per_TF_motif(simulated_annotated_input_files=created_cohorts[cohort][1:], 
                                                                             cohort_mean_sd_per_tf_overall_output_dict_file= dict_simulated_mean_sd_per_TF_motif_output_file, 
                                                                             motif_name_index = motif_name_index, f_score_index = f_score_index, motif_breaking_score_index = motif_breaking_score_index)
    #get sig muts with fdr<0.2
    muts_sig_per_TF_file = get_muts_sig_per_TF(annoted_input_file=created_cohorts[cohort][0], dict_simulated_mean_sd_per_TF_motif=dict_simulated_mean_sd_per_TF_motif, 
                                               annoted_output_file_extension="_rand{}setsTF".format(len(mutation_input_files)-1), annoted_output_file_extension_onlysig="_rand{}setsTFsigQval{}".format(len(mutation_input_files)-1, sig_level_per_TF_thresh),
                                               motif_name_index = motif_name_index, f_score_index = f_score_index, motif_breaking_score_index = motif_breaking_score_index,
                                               filter_on_qval=False, sig_thresh=sig_level_per_TF_thresh, 
                                               filter_on_signal = True, dnase_index = 24, fantom_index = 25, num_other_tfs = 27)
    sig_muts_per_tf_mutation_input_files = [muts_sig_per_TF_file]
    #get sig muts from the simulated with pval<0.5
    
    for mutations_input_file in created_cohorts[cohort][1:]: 
        muts_sig_per_TF_file = get_muts_sig_per_TF(annoted_input_file=mutations_input_file, dict_simulated_mean_sd_per_TF_motif=dict_simulated_mean_sd_per_TF_motif,
                                                   annoted_output_file_extension="_rand{}setsTF".format(len(mutation_input_files)-1), annoted_output_file_extension_onlysig=sim_output_extension,
                                                   motif_name_index = motif_name_index, f_score_index = f_score_index, motif_breaking_score_index = motif_breaking_score_index,
                                                   filter_on_qval=False, sig_thresh=sim_sig_level_per_TF_thresh)
        sig_muts_per_tf_mutation_input_files.append(muts_sig_per_TF_file)
     
    #mutation_input_files[0]
    sig_tfs_file, sig_tfpos_file = get_tf_pval(cohort, sig_muts_per_tf_mutation_input_files, motif_name_index, f_score_index, motif_breaking_score_index,
               sig_level_per_TF_thresh, filter_cond, fsep='\t', sig_tfs_file=sig_tfs_file, sig_tfpos_file=sig_tfpos_file, filter_on_qval=False,
               filter_on_signal = True, dnase_index = 24, fantom_index = 25, num_other_tfs = 27)
    unified_mutation_input_files = []
    
    for mutations_input_file in sig_muts_per_tf_mutation_input_files:
        unified_muts_file = mutations_input_file + output_extension + "_groupedbymut" 
        unified_muts_file_wihtmotifinfo = unified_muts_file+"withmotifinfo"
        if not os.path.exists(unified_muts_file_wihtmotifinfo):
            print "Unifying: ", mutations_input_file
            unify_muts(mutations_input_file, unified_muts_file, filter_mut_motifs=True, filter_cond=filter_cond, operation_on_unify=operation_on_unify)
            get_max_motif_in_grouped_muts(annotated_mutations_grouped_file=unified_muts_file, annotated_mutations_grouped_output_file=unified_muts_file_wihtmotifinfo)
            os.remove(unified_muts_file)
        unified_mutation_input_files.append(unified_muts_file_wihtmotifinfo)
    #print 'Unified mutations input files: ', unified_mutation_input_files
    get_sig_merged_elements(unified_mutation_input_files, cohort_full_name, sim_output_extension+output_extension, distance_to_merge, merged_mut_sig_threshold, local_domain_window, chr_lengths_file, sig_elements_output_file, sim_sig_level_per_TF_thresh)
    
    return sig_elements_output_file, sig_tfs_file, sig_tfpos_file
    
def process_cohorts(cohort_names_input):
    
    mutations_cohorts_dir = 'mutations_cohorts_output_exclVEP'#_Alltumors_noexon'
    
    observed_input_file = 'mutations_files/obsann22May2017_exclVEP.bed9'
    simulated_input_dir = 'mutations_simulations_files_103sets_exclVEP'
    simulated_input_files = [simulated_input_dir+'/'+x for x in os.listdir(simulated_input_dir) if '_annotated.bed9' in x]
    mutation_input_files = [observed_input_file]
    mutation_input_files.extend(simulated_input_files)
    
    sig_level_per_TF_thresh = 0.05
    sim_sig_level_per_TF_thresh = 1.0#0.05
    motif_name_index = 17
    f_score_index = 9
    motif_breaking_score_index = 10
    sim_output_extension = "_rand{}setsTFsigPval{}".format(len(mutation_input_files)-1, sim_sig_level_per_TF_thresh)#len(mutation_input_files)-1 
    operation_on_unify = 'mean'#'max'
    
    distance_to_merge = 200#1000
    merged_mut_sig_threshold = 0.05
    
    local_domain_window = 25000#50000
    chr_lengths_file = 'datafiles/chr_lengths_hg19.txt'
    
    output_extension = "_{operation_on_unify}TFExprMotifBreaking03Filters".format(operation_on_unify=operation_on_unify)
    filter_cond = 'if(($32>0 || $32~"nan") && $11>=0.3)'# && ($31>0 || $31~"nan")
    
    #output_extension = "_{operation_on_unify}TFExprTFBindMotifBreaking03Filters".format(operation_on_unify=operation_on_unify)
    #filter_cond = 'if(($31>0 || $31~"nan")  && ($32>0 || $32~"nan") && $11>=0.3)'
    
    generated_sig_merged_element_files = []
    sig_tfs_files = []
    sig_tfpos_files = []
    run_parallel = False
    p = Pool(12)
    
    stats_ext = "_rand{}setsTF".format(len(mutation_input_files)-1)
    cohorts = get_cohorts(cohort_names_input)
    created_cohorts = generate_cohorts(mutation_input_files, cohorts, mutations_cohorts_dir, stats_ext=stats_ext)
    
    for cohort in created_cohorts.keys():
        sig_elements_output_file = created_cohorts[cohort][0] + "_rand{}setsTFsigQval{}".format(len(mutation_input_files)-1, sig_level_per_TF_thresh) + output_extension + "_groupedbymut"+"withmotifinfo"+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)+"_statspvaluesSimSig"+str(sim_sig_level_per_TF_thresh)+"_statspvalueslocalw{local_domain_window}onlysig{merged_mut_sig_threshold}".format(local_domain_window=local_domain_window, merged_mut_sig_threshold=merged_mut_sig_threshold)
        sig_tfs_file = created_cohorts[cohort][0] + "_rand{}setsTFsigQval{}".format(len(mutation_input_files)-1, sig_level_per_TF_thresh) + '_sigTFs_{}'.format(sig_level_per_TF_thresh) 
        sig_tfpos_file = created_cohorts[cohort][0] + "_rand{}setsTFsigQval{}".format(len(mutation_input_files)-1, sig_level_per_TF_thresh) + '_sigTFpos_{}'.format(sig_level_per_TF_thresh)
        if os.path.exists(sig_elements_output_file) and os.path.exists(sig_tfs_file) and os.path.exists(sig_tfpos_file):
            generated_sig_merged_element_files.append(sig_elements_output_file)
            sig_tfs_files.append(sig_tfs_file)
            sig_tfpos_files.append(sig_tfpos_file)
            continue
        if len(created_cohorts[cohort])<2:
            print 'WARNING: no observed/simulated mutations file for ', cohort, created_cohorts[cohort]
            continue
        
        if run_parallel:
            p.apply_async(run_cohort, args=(cohort, created_cohorts, mutation_input_files, motif_name_index, f_score_index, motif_breaking_score_index,
               sig_level_per_TF_thresh, sim_output_extension, sim_sig_level_per_TF_thresh,
               filter_cond, operation_on_unify, output_extension, distance_to_merge, merged_mut_sig_threshold,
               local_domain_window, chr_lengths_file, sig_elements_output_file, sig_tfs_file, sig_tfpos_file))#, callback=generated_sig_merged_element_files.append)
        else:
            run_cohort(cohort, created_cohorts, mutation_input_files, motif_name_index, f_score_index, motif_breaking_score_index,
               sig_level_per_TF_thresh, sim_output_extension, sim_sig_level_per_TF_thresh,
               filter_cond, operation_on_unify, output_extension, distance_to_merge, merged_mut_sig_threshold,
               local_domain_window, chr_lengths_file, sig_elements_output_file, sig_tfs_file, sig_tfpos_file)
        generated_sig_merged_element_files.append(sig_elements_output_file)
        sig_tfs_files.append(sig_tfs_file)
        sig_tfpos_files.append(sig_tfpos_file)
        print "Processed:" + sig_elements_output_file
        
    if run_parallel:
        p.close()
        p.join()
    
    return generated_sig_merged_element_files, sig_tfs_files, sig_tfpos_files
    
#awk 'BEGIN{FS=OFS="\t"}{if($18<0.05) sig="element_sig"; else sig="element_notsig";  print $4,"Observed", sig}' mutations_cohorts_output/All-tumors_observed_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmutsonlysig0.1_mergedmuts20bp_statspvalues >> ../analysis/All-tumors_6setsobs_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_scores.col3

#awk 'BEGIN{FS=OFS="\t"}{sig="mut_sig"; if($16<0.1) sig="mut_sig"; else sig="mut_notsig";  print $10,"Rand3", sig}' All-tumors_randomised98f100kw50kbnonparametric_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts >> ../analysis/All-tumors_6setsobs_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_scores.col3

if __name__ == '__main__':
    cohort_names_input = sys.argv[1]
    generated_sig_merged_element_files, sig_tfs_files, sig_tfpos_files = process_cohorts(cohort_names_input)
    #print "Generated Sig. Element Sets: \n", '\n'.join(generated_sig_merged_element_files)
    
