'''
Created on Mar 2, 2017

@author: husensofteng
'''
import matplotlib
from boto.file.key import Key
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns
#sns.set(style="ticks")
sns.set_style("whitegrid")
sns.set_context("paper")#talk

def plot_score_hist(df, x, y, hue, out_file):
    plt.figure(figsize=(8, 10))
    plt_results = sns.boxplot(y=y, x=x, hue=hue, data=df, palette="PRGn")
    sns.despine(offset=5, trim=True)
    fig = plt_results.get_figure()
    #fig.savefig(out_file+".svg")
    fig.savefig(out_file+".pdf")
    #plt.show()
    return
        
def plot_box_plots(data, x=None, y=None, hue=None, title='', out_file=''):
    width = 8
    height = 10
    
    plt.figure(figsize=(width, height))
    plt_results = sns.boxplot(x=x, y=y, hue=hue, data=data)#, palette="PRGn")
    #sns.despine(offset=5, trim=True)
    plt_results.axes.set_title(title)
    plt_results.set_xticklabels(plt_results.get_xticklabels(), rotation=90)
    fig = plt_results.get_figure()
    #fig.savefig(out_file+".svg")
    fig.savefig(out_file+".pdf")
    #plt.show()
    return

def plot_count_plots(data, x=None, y=None, hue=None, title="Candidate Regulatory Elements", out_file=''):
    
    plt.figure(figsize=(12, 6))
    plt_results = sns.barplot(x=x, y=y, hue=hue, data=data, ci=None, estimator=sum)#, palette="PRGn")
    plt_results.axes.set_title(title)
    plt_results.axes.set_ylabel(y)
    #sns.despine(offset=5, trim=True)
    plt_results.set_xticklabels(plt_results.get_xticklabels(), rotation=90)
    fig = plt_results.get_figure()
    #fig.savefig(out_file+".svg")
    fig.savefig(out_file+".pdf")
    #plt.show()
    return

def convert_value(x):
    if float(x)>0.05:
        return 'notsig'
    else:
        return 'sig'
    
def process_merged_muts(sig_elements_file, mut_info_cols=['chr', 'start', 'end', 'score', 'num_muts','num_samples', 'cancer_types', 'muts', 'motif_positions', 'motif_names', 'pvalues', 'scores', 'cancer_types_collapse', 'donors', 'donors_collapse', 'element_motif_info', 'pval', 'qval']):
    sig_elements = []
    
    with open(sig_elements_file, 'r') as sig_elements_infile:
        l = sig_elements_infile.readline()
        while l:
            ls = []
            for x in l.strip().split('\t'):
                try:
                    ls.append(float(x))
                except ValueError:
                    ls.append(x)
                #sig_elements.append([try: float(x) except: x for x in l.strip().split('\t')])
            sig_elements.append(ls)
            l = sig_elements_infile.readline()
    
    df = pd.DataFrame(sig_elements, columns=mut_info_cols)
    df['sig'] = df['qval'].apply(convert_value)
    #print df['score'].sort().head(10)
    df['score']=df['score'].apply(np.log2)
    plot_score_hist(df, x='sig', y='score', hue='sig', out_file='merged_muts.png')
    print 'plot is done'
    
def write_content_to_file(output_files_scores_combined, input_files, names_for_files, 
                          pval_index, score_index, sig_label, notsig_label, stat_value_threshold):
    
    with open(output_files_scores_combined, 'a') as output_files_scores_combined_outfile:
        for i,mut_file in enumerate(input_files):
            with open(mut_file, 'r') as mut_file_r:
                l = mut_file_r.readline()
                while l:
                    sl = l.strip().split('\t')
                    if pval_index is None:
                        output_files_scores_combined_outfile.write("{:.4f}".format(sl[score_index]) + '\t' + names_for_files[i] + '\t' + sig_label + '\n')
                    else:
                        if float(sl[pval_index])<stat_value_threshold:
                            output_files_scores_combined_outfile.write("{:.4f}".format(sl[score_index]) + '\t' + names_for_files[i] + '\t' + sig_label + '\n')
                        else:
                            output_files_scores_combined_outfile.write("{:.4f}".format(sl[score_index]) + '\t' + names_for_files[i] + '\t' + notsig_label + '\n')
                    l = mut_file_r.readline()
    return

def plot_score_boxplots(output_files_scores_combined="All-tumors_combined_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_statmuts0.1.col3"):
    
    simulation_mut_files = ['mutations_cohorts_output/All-tumors_simulation_broad_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1',
                            'mutations_cohorts_output/All-tumors_simulation_dkfz_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1',
                            'mutations_cohorts_output/All-tumors_simulation_Sangerneutral_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1',
                            'mutations_cohorts_output/All-tumors_randomised100f100kw50kbnonparametric_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1',
                            'mutations_cohorts_output/All-tumors_randomised99f100kw50kbnonparametric_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1',
                            'mutations_cohorts_output/All-tumors_randomised98f100kw50kbnonparametric_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1']
    names_for_simulation_files = ['Broad', 'DKFZ', 'Sanger', 'Random1', 'Random2', 'Random3']
    
    simulation_element_files = ['mutations_cohorts_output/All-tumors_simulation_broad_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1onlysig_mergedmuts20bp',
                            'mutations_cohorts_output/All-tumors_simulation_dkfz_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1onlysig_mergedmuts20bp',
                            'mutations_cohorts_output/All-tumors_simulation_Sangerneutral_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1onlysig_mergedmuts20bp',
                            'mutations_cohorts_output/All-tumors_randomised100f100kw50kbnonparametric_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1onlysig_mergedmuts20bp',
                            'mutations_cohorts_output/All-tumors_randomised99f100kw50kbnonparametric_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1onlysig_mergedmuts20bp',
                            'mutations_cohorts_output/All-tumors_randomised98f100kw50kbnonparametric_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1onlysig_mergedmuts20bp']
    
    observed_mut_file = ["mutations_cohorts_output/All-tumors_observed_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1"]
    observed_element_file = ["mutations_cohorts_output/All-tumors_observed_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1onlysig_mergedmuts20bp_statspvalues"]
    score_index_mut_files = 9
    score_index_element_files = 3
    pval_index_mut_files = 15
    qval_index_element_files = 17
    
    if not os.path.exists(output_files_scores_combined):
        #sim-mut
        print 'writing sim-muts'
        write_content_to_file(output_files_scores_combined=output_files_scores_combined, 
                input_files=simulation_mut_files, names_for_files=names_for_simulation_files, 
                pval_index=pval_index_mut_files, score_index=score_index_mut_files, sig_label='sig_mut', notsig_label='notsig_mut', stat_value_threshold=0.1)
        #sim-elements
        print 'writing sim-elements'
        write_content_to_file(output_files_scores_combined=output_files_scores_combined, 
                input_files=simulation_element_files, names_for_files=names_for_simulation_files, 
                pval_index=None, score_index=score_index_element_files, sig_label='sim_elements', notsig_label='sim_elements', stat_value_threshold=0.05)
        
        #obs-mut
        write_content_to_file(output_files_scores_combined=output_files_scores_combined, 
                input_files=observed_mut_file, names_for_files=['Observed'], 
                pval_index=pval_index_mut_files, score_index=score_index_mut_files, sig_label='sig_mut', notsig_label='notsig_mut', stat_value_threshold=0.1)
        #obs-elements
        write_content_to_file(output_files_scores_combined=output_files_scores_combined, 
                input_files=observed_element_file, names_for_files=['Observed'], 
                pval_index=qval_index_element_files, score_index=score_index_element_files, sig_label='sig_elements', notsig_label='notsig_elements', stat_value_threshold=0.05)
    
    #merged_mut_info_cols=['chr', 'start', 'end', 'score', 'num_muts','num_samples', 'cancer_types', 'muts', 'motif_positions', 'motif_names', 'pvalues', 'scores', 'cancer_types_collapse', 'donors', 'donors_collapse', 'element_motif_info', 'pval', 'qval'] 
    #process_merged_muts(sig_elements_file, merged_mut_info_cols)
    
    mut_info_cols=['Functional Score (log2)', 'Mutation Source', 'Significance']
    df = pd.read_csv(output_files_scores_combined, sep='\t', names=mut_info_cols, header=None)
    print df.head(5)
    df['Functional Score (log2)']=df['Functional Score (log2)'].apply(np.log2)
    plot_score_hist(df, x='Mutation Source', y='Functional Score (log2)', hue='Significance', out_file=output_files_scores_combined)
    print 'plot is done'
    
def plot_mutation_rate_per_cancer_type(output_files_mutation_rate, observed_mutation_file):
    
    index_cancer_type = 5
    index_donor_IDs = 8
    dict_cancer_types = {}
    
    with open(observed_mutation_file, 'r') as observed_element_ifile:
        l = observed_element_ifile.readline().strip().split('\t')
        while l and len(l)>=index_cancer_type:
            cancer_type = l[index_cancer_type]
            donor_ID = l[index_donor_IDs]
            try:
                dict_cancer_types[cancer_type][donor_ID] +=1
            except KeyError:
                try:
                    dict_cancer_types[cancer_type][donor_ID] = 1
                except KeyError:
                    dict_cancer_types[cancer_type] = {}
                    
            l = observed_element_ifile.readline().strip().split('\t')
    
    cancer_rows = []
    sample_rows = []
    number_muts = []
    for cancer_type in sorted(dict_cancer_types.keys()):
        for sample in sorted(dict_cancer_types[cancer_type].keys()):
            cancer_rows.append(cancer_type)
            sample_rows.append(sample)
            number_muts.append(dict_cancer_types[cancer_type][sample])
            
    df = pd.DataFrame()
    df['Cancer Types'] = cancer_rows
    df['Donor IDs'] = sample_rows
    df['Number of Mutations'] = number_muts
    plot_box_plots(df, x='Cancer Types', y='Number of Mutations', title="Rate of Observed Mutations Across Cancer Samples", out_file=output_files_mutation_rate+"persample")
    print df.head(3)
    plot_count_plots(df, x='Cancer Types', y='Number of Mutations', hue=None, title="Number of Observed Mutations Across Cancer Types", out_file=output_files_mutation_rate+"percancertype")
    
    '''
    col_names = ['chr','start','end','ref','alt','Cancer Types', 'Mutation Types', 'Sample IDs', 'Donor IDs']
    df = pd.read_csv(initial_observed_mutation_file, sep='\t', names=col_names, header=None)
    df['num'] = 1
    groupby_object = df[['Cancer Types', 'Donor IDs', 'num']].groupby(['Cancer Types', 'Donor IDs']).size()
    group_names = []
    group_values = []
    for i,group_index in enumerate(groupby_object.index.labels[0]):
        group_names.append(groupby_object.index.levels[0][group_index])
        group_values.append(groupby_object[i])
    new_df = pd.DataFrame()
    new_df['Cancer Types'] = group_names
    new_df['Number of Mutations per Sample'] = group_values
    plot_box_plots(new_df, x='Cancer Types', y='Number of Mutations per Sample', out_file=output_files_mutation_rate)
    #plot_box_plots(data=df, x='Cancer Types', out_file=output_files_mutation_rate)
    '''

def plot_element_mutation_rate_per_cancer_type(output_files_mutation_rate, observed_element_file):
    
    index_cancer_type = 11
    index_donor_IDs = 13
    dict_cancer_types = {}
    
    with open(observed_element_file, 'r') as observed_element_ifile:
        l = observed_element_ifile.readline().strip().split('\t')
        while l and len(l)>=index_cancer_type:
            
            cancer_types = l[index_cancer_type].split(',')
            donor_IDs = l[index_donor_IDs].split(',')
            for i, cancer_type in enumerate(cancer_types):
                try:
                    dict_cancer_types[cancer_type][donor_IDs[i]] +=1
                except KeyError:
                    try:
                        dict_cancer_types[cancer_type][donor_IDs[i]] = 1
                    except KeyError:
                        dict_cancer_types[cancer_type] = {}
                    
            l = observed_element_ifile.readline().strip().split('\t')
    
    cancer_rows = []
    sample_rows = []
    number_muts = []
    num_muts_per_cancer = {}
    num_samples_per_cancer = {}
    for cancer_type in sorted(dict_cancer_types.keys()):
        num_muts_per_cancer[cancer_type]=0
        num_samples_per_cancer[cancer_type]=[]
        for sample in sorted(dict_cancer_types[cancer_type].keys()):
            cancer_rows.append(cancer_type)
            sample_rows.append(sample)
            number_muts.append(dict_cancer_types[cancer_type][sample])
            num_samples_per_cancer[cancer_type].append(sample)
            num_muts_per_cancer[cancer_type]+=dict_cancer_types[cancer_type][sample]
    
    df = pd.DataFrame()
    df['Cancer Types'] = pd.Series(cancer_rows).values
    df['Donor IDs'] = pd.Series(sample_rows).values
    df['Number Candidate Regulatory Mutations per Sample'] = pd.Series(number_muts).values
    #print df
    plot_box_plots(df, x='Cancer Types', y='Number Candidate Regulatory Mutations per Sample', title="Distribution of Candidate Regulatory Mutations Across Cancer Sampels", 
                   out_file=output_files_mutation_rate+"number_muts_per_sample")
    
    #number samples with reg mutations per cancer type
    cancers_samples = []
    num_samples = []
    for c in sorted(num_samples_per_cancer.keys()):
        cancers_samples.append(c)
        num_samples.append(len(set(num_samples_per_cancer[c])))
    num_samples_per_cancer_df = pd.DataFrame()
    num_samples_per_cancer_df['Cancer Types'] = pd.Series(cancers_samples).values
    num_samples_per_cancer_df['Number of Samples'] = pd.Series(num_samples).values 
    print num_samples_per_cancer_df
    plot_count_plots(data=num_samples_per_cancer_df, x='Cancer Types', y='Number of Samples', title="Number of Samples with Candidate Regulatory Mutations Across Cancer Types", 
                     out_file=output_files_mutation_rate+"number_samples_per_cancer")
    
    #number mutations per cancer type
    cancers = []
    num_muts = []
    for c in sorted(num_muts_per_cancer.keys()):
        cancers.append(c)
        num_muts.append(num_muts_per_cancer[c])
    num_muts_per_cancer_df = pd.DataFrame()
    num_muts_per_cancer_df['Cancer Types'] = pd.Series(cancers).values
    num_muts_per_cancer_df['Number of Candidate Regulatory Mutations'] = pd.Series(num_muts).values
    plot_count_plots(num_muts_per_cancer_df, x='Cancer Types', y='Number of Candidate Regulatory Mutations',title="Distribution of Candidate Regulatory Mutations Across Cancer Types", 
                     out_file=output_files_mutation_rate+"number_muts_per_cancer")
    
def generate_dict_data(input_file):
    
    motif_mut_inifo_index = 14
    fsep = '\t'
    vsep = '#'
    msep = '*'
    MotifInfo, matching_motifs_sep, MaxMotif_sep, Pval_sep, Qval_sep = 'MotifInfo', 'MatchingMotifs', 'MaxMotif', 'Pval', 'Qval'
    
    cancer_type_index_mut = 5
    chromatin_region_index_motif = 14
    dnase1_index_motif = 16
    motif_name_index_motif = 9
    TF_binding_info_index_motif = 22
    
    num_muts_per_pos_per_TF = 0
    dict_data = {'Chromatin Regions':{}, 'Motifs':{}, 'Max Motifs': {}, 'Motifs with Binding Evidence':{},
                 'DNase1 Sites':{}, 'Chromatin Regions per Motif': {}}
    counts_dict = {'Motif Counts':{}, 'Chromatin Status Counts':{}}
    with open(input_file, 'r') as ifile:
        l =  ifile.readline().strip().split(fsep)
        while l and len(l)>=motif_mut_inifo_index:
            muts_in_this_element = l[motif_mut_inifo_index].split(',')
            for i,mut in enumerate(muts_in_this_element):
                mut_info = mut.split(matching_motifs_sep)[0].split(vsep)
                
                motifs_info = mut.split(matching_motifs_sep)[1].split(MaxMotif_sep)[0].split(MotifInfo)
                max_motif_info = mut.split(matching_motifs_sep)[1].split(MaxMotif_sep)[1].split(vsep)
                #Chromatin Regions
                
                try:
                    dict_data['Chromatin Regions'][mut_info[cancer_type_index_mut]][max_motif_info[chromatin_region_index_motif]] += 1
                except KeyError:
                    try:
                        dict_data['Chromatin Regions'][mut_info[cancer_type_index_mut]][max_motif_info[chromatin_region_index_motif]] = 1
                    except KeyError:
                        try:
                            dict_data['Chromatin Regions'][mut_info[cancer_type_index_mut]] = {max_motif_info[chromatin_region_index_motif]:1}
                        except KeyError:
                            dict_data['Chromatin Regions'] = {mut_info[cancer_type_index_mut]}
                            dict_data['Chromatin Regions'][mut_info[cancer_type_index_mut]] = {max_motif_info[chromatin_region_index_motif]:1}
                
                #DNase1 Sites
                if float(max_motif_info[dnase1_index_motif])>0:
                    try:
                        dict_data['DNase1 Sites'][mut_info[cancer_type_index_mut]]['DNase1 Sites'] += 1
                    except KeyError:
                        try:
                            dict_data['DNase1 Sites'][mut_info[cancer_type_index_mut]]['DNase1 Sites'] = 1
                        except KeyError:
                            try:
                                dict_data['DNase1 Sites'][mut_info[cancer_type_index_mut]] = {'DNase1 Sites':1}
                            except KeyError:
                                dict_data['DNase1 Sites'] = {mut_info[cancer_type_index_mut]}
                                dict_data['DNase1 Sites'][mut_info[cancer_type_index_mut]] = {'DNase1 Sites':1}
                  
                #Max Motifs
                try:
                    dict_data['Max Motifs'][mut_info[cancer_type_index_mut]][max_motif_info[motif_name_index_motif]] += 1
                except KeyError:
                    try:
                        dict_data['Max Motifs'][mut_info[cancer_type_index_mut]][max_motif_info[motif_name_index_motif]] = 1
                    except KeyError:
                        try:
                            dict_data['Max Motifs'][mut_info[cancer_type_index_mut]] = {max_motif_info[motif_name_index_motif]:1}
                        except KeyError:
                            dict_data['Max Motifs'] = {mut_info[cancer_type_index_mut]}
                            dict_data['Max Motifs'][mut_info[cancer_type_index_mut]] = {max_motif_info[motif_name_index_motif]:1}
                #Motifs
                motifs_counted = []
                chromatin_status_counted = []
                for motif in motifs_info:
                    motif = motif.split(vsep)
                    #in case a mutation was overlapping two motifs of the same TF just one instance should be counted 
                    if 'CTCF' in motif[motif_name_index_motif]:
                        #mut_motif_pos =motif[motif_name_index_motif-4]
#                         if '-' in mut_motif_pos:
#                             mut_motif_positions = range(int(mut_motif_pos.split('-')[0]), int(mut_motif_pos.split('-')[1])+1)
#                             for kk in mut_motif_positions:
#                                 if kk == 9:
#                                     num_muts_per_pos_per_TF+=1
#                         elif int(mut_motif_pos)==9:
                        pass#num_muts_per_pos_per_TF+=1
                        
                    if motif[motif_name_index_motif] in motifs_counted:
                        continue
                    motifs_counted.append(motif[motif_name_index_motif])
                    if 'SP1' in motif[motif_name_index_motif]:
                        num_muts_per_pos_per_TF+=1
                    #count chromatin status
                    if motif[chromatin_region_index_motif] not in chromatin_status_counted:
                        try:
                            counts_dict['Chromatin Status Counts'][motif[chromatin_region_index_motif]] +=1
                        except KeyError:
                            counts_dict['Chromatin Status Counts'][motif[chromatin_region_index_motif]] =1
                    chromatin_status_counted.append(motif[chromatin_region_index_motif])
                    
                    try:
                        dict_data['Chromatin Regions per Motif'][motif[chromatin_region_index_motif]][motif[motif_name_index_motif].split('_')[0]] += 1
                    except KeyError:
                        try:
                            dict_data['Chromatin Regions per Motif'][motif[chromatin_region_index_motif]][motif[motif_name_index_motif].split('_')[0]] = 1
                        except KeyError:
                            try:
                                dict_data['Chromatin Regions per Motif'][motif[chromatin_region_index_motif]] = {motif[motif_name_index_motif].split('_')[0]: 1}
                            except KeyError:
                                dict_data['Chromatin Regions per Motif'] = [motif[chromatin_region_index_motif]]
                                dict_data['Chromatin Regions per Motif'][motif[chromatin_region_index_motif]] = {motif[motif_name_index_motif].split('_')[0]: 1}
                            
                    #count num muts per TF motif (overall count)
                    try:
                        counts_dict['Motif Counts'][motif[motif_name_index_motif].split('_')[0]] +=1
                    except KeyError:
                        counts_dict['Motif Counts'][motif[motif_name_index_motif].split('_')[0]] =1
                    #cancer specific count
                    try:
                        dict_data['Motifs'][mut_info[cancer_type_index_mut]][motif[motif_name_index_motif]] += 1
                    except KeyError:
                        try:
                            dict_data['Motifs'][mut_info[cancer_type_index_mut]][motif[motif_name_index_motif]] = 1
                        except KeyError:
                            try:
                                dict_data['Motifs'][mut_info[cancer_type_index_mut]] = {motif[motif_name_index_motif]:1}
                            except KeyError:
                                dict_data['Motifs'] = {mut_info[cancer_type_index_mut]}
                                dict_data['Motifs'][mut_info[cancer_type_index_mut]] = {motif[motif_name_index_motif]:1}
                    
                    #Motifs with Binding Evidence
                    if motif[TF_binding_info_index_motif]!='nan':
                        if float(motif[TF_binding_info_index_motif])>0:
                            try:
                                dict_data['Motifs with Binding Evidence'][mut_info[cancer_type_index_mut]][motif[motif_name_index_motif]] += 1
                            except KeyError:
                                try:
                                    dict_data['Motifs with Binding Evidence'][mut_info[cancer_type_index_mut]][motif[motif_name_index_motif]] = 1
                                except KeyError:
                                    try:
                                        dict_data['Motifs with Binding Evidence'][mut_info[cancer_type_index_mut]] = {motif[motif_name_index_motif]:1}
                                    except KeyError:
                                        dict_data['Motifs with Binding Evidence'] = {mut_info[cancer_type_index_mut]}
                                        dict_data['Motifs with Binding Evidence'][mut_info[cancer_type_index_mut]] = {motif[motif_name_index_motif]:1}
                         
            l =  ifile.readline().strip().split(fsep)   
    print "num_muts_per_pos_per_TF: ", num_muts_per_pos_per_TF
    return dict_data, counts_dict

def plot_enrichment(dict_data, output_file):
    
    for k in dict_data.keys():
        main_items = []
        secondary_items = []
        number_items = []
        for main_item in sorted(dict_data[k].keys()):
            for secondary_item in sorted(dict_data[k][main_item].keys()):
                if dict_data[k][main_item][secondary_item]>=0:
                    main_items.append(main_item)
                    if "Motifs" in k:
                        secondary_items.append(secondary_item.split('_')[0])
                    else:
                        secondary_items.append(secondary_item)
                    number_items.append(int(dict_data[k][main_item][secondary_item]))
        
        print k, len(main_items), len(secondary_items), len(number_items) 
        
        df = pd.DataFrame()
        df['Cancer Types'] = pd.Series(main_items).values
        df[k] = pd.Series(secondary_items).values
        df['Values'] = pd.Series(number_items).values
        df = df.pivot('Cancer Types', k, 'Values').fillna(0)
        df = df[df.sum(axis=1) > 40]
        df_pivot = pd.DataFrame()
        for c in df.columns:
            if df[c].sum()>40:
                df_pivot[c] = df[c] 
        #df_pivot = df_pivot[df_pivot.sum(axis=0) > 20]
        if df_pivot.shape[0]>0 and df_pivot.shape[1]>0:
            print df_pivot
            
            plt.figure(figsize=(25, 10))
            plt_results = sns.heatmap(df_pivot, annot=True, fmt='g', linewidths=.5)
            plt_results.axes.set_title("Enrichment of mutations at {}".format(k))
            plt_results.axes.set_ylabel("Cancer Types")
            plt_results.axes.set_xlabel(k)
            fig = plt_results.get_figure()
            fig.savefig(output_file+ k +".pdf")
        
def plot_count_from_dict(dict_data, output_file):
    for k in dict_data.keys():
        main_items = []
        number_items = []
        for main_item in sorted(dict_data[k].keys()):
            if dict_data[k][main_item]>10:
                main_items.append(main_item)
                number_items.append(dict_data[k][main_item])
        
        print k, len(main_items), len(number_items) 
        
        df = pd.DataFrame()
        df[k] = pd.Series(main_items).values
        df['Values'] = pd.Series(number_items).values
        plot_count_plots(data=df, x=k, y='Values', title="Frequency of " + k, out_file= output_file+ k+"_freq")

#not used:
def stat_plot_scores(mut_file, ovcombined, ov, ovc, ovp):
    
    print mut_file
    print len(ov), np.mean(ov), np.std(ov)
    
    title = (mut_file.split('/')[-1]).split('.')[0]
    plt.hist(ov)
    plt.xlabel('Scores')
    plt.ylabel('#Mutations')
    plt.title('Distribution of scores \n ({})'.format(title))
    plt.savefig('{}.png'.format(mut_file))
    
    for c in ovc:
        print c, len(ovc[c]), np.mean(ovc[c]), np.std(ovc[c])
    for p in ovp:
        print p, len(ovp[p]), np.mean(ovp[p]), np.std(ovp[p])
    
    dfcombined = pd.DataFrame(ovcombined)
    plt.clf()
    plt.close()
    
    plt.figure()
    ax = dfcombined.boxplot(by=['chromatin_status'], rot=90)
    plt.title(title)
    fig = ax.get_figure()
    fig.savefig(mut_file+"_chromatin_status_boxplot"+".png")
    plt.clf()
    plt.close()
    
    plt.figure()
    ax = dfcombined.boxplot(by=['cancer_types'], rot=90)
    plt.title(title)
    fig = ax.get_figure()
    #fig = ax[0][0].get_figure()
    fig.savefig(mut_file+"_cancer_types_boxplot"+".png")
    plt.clf()
    plt.close()
    return

if __name__ == '__main__':
    
    #plot_score_boxplots(output_files_scores_combined="analysis/All-tumors_combined_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_statmuts0.1.col3")
    
    #plot_mutation_rate_per_cancer_type(output_files_mutation_rate="analysis/observed_muts_count_", observed_mutation_file = "mutations_files/observed.bed9")
    '''sig_elements_input_file = "mutations_cohorts_output/All-tumors_observed_annotated.bed9sigperTFOvarallQval0.2_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1onlysig_mergedmuts20bp_statspvaluesonlysig" #"mutations_cohorts_output/All-tumors_observed_annotated.bed9_simulated6setsmaxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_statmuts0.1onlysig_mergedmuts20bp_statspvaluesonlysig"
    plot_element_mutation_rate_per_cancer_type(output_files_mutation_rate="analysis/sginificant_elements_", observed_element_file = sig_elements_input_file)
    
    dict_data = generate_dict_data(input_file=sig_elements_input_file)
    plot_enrichment(dict_data, output_file="analysis/significant_elements_enrichment_of_")
    '''
    
    #sig_elements_input_file = "mutations_cohorts_output/Cervix-SCC_observed_annotated.bed9_TFsigQval0.2rand13sets_maxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts20bp_statspvaluesonlysig0.05"
    #sig_elements_input_file = "mutations_cohorts_output/All-tumors_observed_annotated.bed9_TFsigQval0.2rand13sets_maxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts20bp_statspvaluesonlysig0.05"#"mutations_cohorts_output/Breast-AdenoCa_observed_annotated.bed9_TFsigQval0.2rand13sets_maxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts20bp_statspvaluesonlysig0.05"
    #sig_elements_input_file = "mutations_cohorts_output/All-tumors-without-Lymphatic-system-tumors-without-Skin-Melanoma_observed_annotated.bed9_TFsigQval0.2rand16sets_maxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts20bp_statspvalues_statspvalueslocalw25000onlysig0.05"
    sig_elements_input_file = "mutations_cohorts_output/All-tumors_observed_annotated.bed9_TFsigQval0.2rand16sets_maxTFExprTFBindMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts1000bp_statspvalues_statspvalueslocalw25000onlysig0.05_onlyrec1_only"
    plot_element_mutation_rate_per_cancer_type(output_files_mutation_rate="analysis/{}_sginificant_elements_".format(sig_elements_input_file.split('/')[1].split('_')[0]), observed_element_file = sig_elements_input_file)
    
    dict_data, count_data = generate_dict_data(input_file=sig_elements_input_file)
    plot_enrichment(dict_data, output_file="analysis/{}_significant_elements_enrichment_of_".format(sig_elements_input_file.split('/')[1].split('_')[0]))
    print count_data
    plot_count_from_dict(count_data, output_file="analysis/{}_significant_elements_enrichment_of_".format(sig_elements_input_file.split('/')[1].split('_')[0]))
    
    