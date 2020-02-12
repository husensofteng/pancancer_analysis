'''
Created on Mar 2, 2017

@author: husensofteng
'''
import matplotlib
import pybedtools
from pybedtools.bedtool import BedTool
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
from decimal import Decimal
import os, sys
import seaborn as sns
import operator
sns.set(style="ticks")
#plt.style.use('ggplot')
sns.set_style("whitegrid")
sns.set_context("paper")#talk

from ProcessCohorts import process_cohorts
from ProcessSigElements import getSigElements
from Utilities import find_overlap_genesets_genelist

def plot_box_plots(data, x=None, y=None, hue=None, title='', out_file=''):
    width = 8
    height = 12
    
    plt.figure(figsize=(width, height))
    plt_results = sns.boxplot(x=x, y=y, hue=hue, data=data)#, palette="PRGn")
    #sns.despine(offset=5, trim=True)
    plt_results.axes.set_title(title)
    plt_results.set_xticklabels(plt_results.get_xticklabels(), rotation=90)
    fig = plt_results.get_figure()
    #fig.savefig(out_file+".svg")
    fig.savefig(out_file+".pdf")
    #plt.show()
    plt.close()
    return

def plot_count_plots(data, x=None, y=None, hue=None, title="Candidate Regulatory Elements", out_file='', estimator=sum):
    
    plt.figure(figsize=(10, 12))
    plt_results = sns.barplot(x=x, y=y, hue=hue, data=data, ci=None, estimator=sum)#, palette="PRGn")
    plt_results.axes.set_title(title)
    plt_results.axes.set_ylabel(y)
    #sns.despine(offset=5, trim=True)
    plt_results.set_xticklabels(plt_results.get_xticklabels(), rotation=90)
    fig = plt_results.get_figure()
    #fig.savefig(out_file+".svg")
    fig.savefig(out_file+".pdf")
    #plt.show()
    plt.close()
    return

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
    plot_box_plots(df, x='Cancer Types', y='Number of Mutations', title="Rate of Observed Mutations Across Cancer Samples", out_file=output_files_mutation_rate)
    plot_count_plots(df, x='Cancer Types', y='Number of Mutations', hue=None, title="Number of Observed Mutations Across Cancer Types", out_file=output_files_mutation_rate+"_barchart")


def plot_element_mutation_rate_per_cancer_type(box_plot_dict, output_files_mutation_rate):
    
    for k in sorted(box_plot_dict.keys()):
        main_items = []
        secondary_items = []
        main_items_for_count_plots = []
        number_secondary_items_for_count_plots = []
        type = box_plot_dict[k]['Type']
        plot_counts = False
        min_to_include = 0
        x_label = 'Cancer Types'
        if 'CountPlot' in box_plot_dict[k].keys():
            plot_counts = True
        if 'X-label' in box_plot_dict[k].keys():
            x_label = box_plot_dict[k]['X-label']
        if 'Minimum' in box_plot_dict[k].keys():
            min_to_include = box_plot_dict[k]['Minimum']
        for r in sorted(box_plot_dict[k].keys()):
            if r == 'Type' or r=='Minimum' or r=='X-label' or r=='CountPlot':
                continue
            values = box_plot_dict[k][r]
            if type == 'Samples':
                values = pd.Series(box_plot_dict[k][r]).value_counts()
            if len(values)>min_to_include:
                main_items_for_count_plots.append(r)
                number_secondary_items_for_count_plots.append(len(box_plot_dict[k][r]))
                for v in values:
                    main_items.append(r)
                    secondary_items.append(v)
        
        df = pd.DataFrame()
        df[x_label] = pd.Series(main_items).values
        df[k] = pd.Series(secondary_items).values
        plot_box_plots(df, x=x_label, y=k,title="Distribution of {k}".format(k=k), 
                         out_file=output_files_mutation_rate+k)
        if plot_counts:
        #count plots
            df = pd.DataFrame()
            df[x_label] = pd.Series(main_items_for_count_plots).values
            df[k] = pd.Series(number_secondary_items_for_count_plots).values
            df.sort_values(by = k, ascending=False, inplace=True)
            plot_count_plots(df, x=x_label, y=k,title="Frequency of {k}".format(k=k.replace('Scores in ', '').replace('Scores of', '').replace('Scores ', '')), 
                             out_file=output_files_mutation_rate+k+'_bar')
            

def plot_enrichment(dict_data, output_file, count_uniq_items):
    threshold_to_include_element = 40
    print dict_data.keys()
    for k in dict_data.keys():
        print 'Plotting ', k
        main_items = []
        secondary_items = []
        number_items = []
        x_label = 'Cancer Type'
        y_label = 'Number of Mutations'
        if 'X-label' in dict_data[k].keys():
            x_label = dict_data[k]['X-label']
        if 'Y-label' in dict_data[k].keys():
            y_label = dict_data[k]['Y-label']
        
        if 'Minimum' in dict_data[k].keys():
            threshold_to_include_element = dict_data[k]['Minimum']
        plot_counts = False
        if 'CountPlot' in dict_data[k].keys():
            plot_counts = True
        for main_item in sorted(dict_data[k].keys()):
            #if "Lymph" in main_item:
            #    continue
            if main_item=='X-label' or main_item=='Y-label' or main_item == 'Minimum' or main_item=='CountPlot':
                continue
            for secondary_item in sorted(dict_data[k][main_item].keys()):
                if len(set(dict_data[k][main_item][secondary_item]))>=0:
                    main_items.append(main_item)
                    secondary_items.append(secondary_item)
                    try:
                        number_items.append(int(dict_data[k][main_item][secondary_item]))
                    except TypeError:
                        if count_uniq_items:
                            number_items.append(len(set(dict_data[k][main_item][secondary_item])))
                        else:
                            number_items.append(len(dict_data[k][main_item][secondary_item]))
        print k, len(main_items), len(secondary_items), len(number_items)
        df = pd.DataFrame()
        df[x_label] = pd.Series(main_items).values
        df[k] = pd.Series(secondary_items).values
        df[y_label] = pd.Series(number_items).values
        df_pivot = df.pivot(x_label, k, y_label)#.fillna(0)
        df_pivot = df_pivot[df_pivot.sum(axis=1) > threshold_to_include_element]
        df_pivot_filtered = pd.DataFrame()
        for c in df_pivot.columns:
            if df_pivot[c].sum()>threshold_to_include_element:
                df_pivot_filtered[c] = df_pivot[c] 
        #df_pivot = df_pivot[df_pivot.sum(axis=0) > 20]
        if df_pivot_filtered.shape[0]>0 and df_pivot_filtered.shape[1]>0:
            plt.figure(figsize=(25, 10))
            plt_results = sns.heatmap(df_pivot_filtered, annot=True, fmt='g', linewidths=.5)
            plt_results.axes.set_title("Enrichment of mutations at {}".format(k))
            plt_results.axes.set_ylabel(x_label)
            plt_results.axes.set_xlabel(k)
            fig = plt_results.get_figure()
            fig.savefig(output_file+ k +".pdf")
        else:
            print 'Skipped: ', k, df_pivot_filtered.shape
        '''if plot_counts:
            print df.head(2)
            for main_item in sorted(set(df[x_label])):
                #df_main = df[(df[x_label] == main_item) & (df[y_label] > 10)].sort_values(by=y_label)
                df_main = df[(df[x_label] == main_item)].sort_values(by=y_label, ascending=False)
                if len(df_main[k])>1:
                    if len(df_main)>20:
                        df_main = df_main[0:20]
                    print len(df_main)
                    plot_count_plots(df_main, x=k, y=y_label,title="Frequency of {k} in {m}".format(k=k.replace(' per Chromatin Status', ' in '), m=main_item), 
                                 out_file=output_file+'_'+main_item.replace('/','-')+'_bar')
        '''
        
        
def plot_count_from_dict(dict_data, output_file, top_n=0):
    
    for k in dict_data.keys():
        type = 'Values'
        main_items = []
        number_items = []
        if 'Type' in dict_data[k].keys():
            type = dict_data[k]['Type']
        for main_item in sorted(dict_data[k].keys()):
            if main_item=='Type':
                continue
            values = dict_data[k][main_item]
            if type == 'Samples':
                values = pd.Series(dict_data[k][main_item]).value_counts()
            
            main_items.append(main_item)
            number_items.append(sum(values))#sum(dict_data[k][main_item]))
        
        print k, len(main_items), len(number_items)
        
        df = pd.DataFrame()
        df[k] = pd.Series(main_items).values
        df['No. Mutations'] = pd.Series(number_items).values
        df.sort_values(by='No. Mutations', ascending=False, inplace=True)
        if top_n>0:
            if len(df)>top_n:
                df = df[0:top_n]
        if len(df)>0:
            plot_count_plots(data=df, x=k, y='No. Mutations', title="Frequency of " + k, out_file= output_file+ k.replace('/','-')+"_freq")

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

def genearate_bar_char_dict(elements_input, cols_names_to_use_for_barchars):
    
    bar_chars_data = {}
    
    for col in cols_names_to_use_for_barchars:
        for cancer_types in elements_input[cols_names_to_use_for_barchars[col]]:
            for cancer_type in cancer_types.split(','):
                if cancer_type==".:1":#some regions might have exonic-regMutations only. Therefore no muts are reported
                    continue
                try:
                    bar_chars_data[col][cancer_type.split(':')[0]].append(int(cancer_type.split(':')[1]))
                except KeyError:
                    try:
                        bar_chars_data[col][cancer_type.split(':')[0]] = [int(cancer_type.split(':')[1])] 
                    except KeyError:
                        bar_chars_data[col] = {cancer_type.split(':')[0]:[int(cancer_type.split(':')[1])]}    
    return bar_chars_data
    
def generate_Muts_heatmap_dict(elements_input, cols_names_to_use, count_one_mut_per_sample_per_element, 
                               min_numnber_muts_to_consider_element=40, min_numnber_muts_to_consider_element_for_number_of_tfs=500):
    
    heatmap_dict_data = {'Other Annotations':{'Minimum':min_numnber_muts_to_consider_element}, 'TF per Chromatin State': {'Minimum':min_numnber_muts_to_consider_element_for_number_of_tfs}}
    box_plot_dict_data = {'Mutations in Sig Elements per Cancer Type': {'Type':'Samples'}}
    for col in cols_names_to_use:
        for element in elements_input[cols_names_to_use[col]]:
            heatmap_dict_data_samples = {'Other Annotations':{}}
            
            items = element.split(',')#muts in this element
            for item in items:
                annotations = item.split('#')[-1].split('|')
                cancer_type = item.split('#')[2]
                sample_id = item.split('#')[-2]
                if cancer_type=='.':
                    continue
                
                try:
                    box_plot_dict_data['Mutations in Sig Elements per Cancer Type'][cancer_type].append(sample_id)
                except KeyError:
                    box_plot_dict_data['Mutations in Sig Elements per Cancer Type'][cancer_type] = [sample_id]
                    
                if annotations[0]=='NaN':
                    continue
                tfbindings = []
                chromhmm_state = []
                for anno in annotations:
                    anno_assay = anno.split(':')[0]
                    anno_values = set(anno.split(':')[1].split(';'))
                    if anno_assay=="ContactingDomain" or anno_assay=='LoopDomain' or anno_assay=='FANTOM' or anno_assay=='DNase-seq':
                        anno_values = [anno_assay] 
                        anno_assay = 'Other Annotations'
                    if anno_assay=="TFBinding":
                        tfbindings = anno_values
                        min_numnber_muts_to_consider_element = min_numnber_muts_to_consider_element_for_number_of_tfs
                    elif anno_assay=="ChromHMM":
                        chromhmm_state=anno_values    
                    #all
                    for anno_val in anno_values:
                        sample_exists = False
                        try:
                            if sample_id in heatmap_dict_data_samples[anno_assay][cancer_type][anno_val]:
                                sample_exists = True
                            else:
                                heatmap_dict_data_samples[anno_assay][cancer_type][anno_val].append(sample_id)
                        except KeyError:
                            try: 
                                heatmap_dict_data_samples[anno_assay][cancer_type][anno_val] = [sample_id]
                            except KeyError:
                                try:
                                    heatmap_dict_data_samples[anno_assay][cancer_type] = {anno_val: [sample_id]}
                                except KeyError:
                                    heatmap_dict_data_samples[anno_assay] = {cancer_type: {anno_val: [sample_id]}}
                        
                        try:
                            if count_one_mut_per_sample_per_element:
                                if not sample_exists:
                                    heatmap_dict_data[anno_assay][cancer_type][anno_val].append(sample_id)
                            else:
                                heatmap_dict_data[anno_assay][cancer_type][anno_val].append(sample_id)
                        except KeyError:
                            try:
                                heatmap_dict_data[anno_assay][cancer_type][anno_val] = [sample_id]
                            except KeyError:
                                try:
                                    heatmap_dict_data[anno_assay][cancer_type] = {anno_val: [sample_id]}
                                except KeyError:
                                    heatmap_dict_data[anno_assay] = {cancer_type: {anno_val: [sample_id]}, 'Minimum':min_numnber_muts_to_consider_element}
                        
                if len(tfbindings)>0 and len(chromhmm_state)>0:
                    for state in chromhmm_state:
                        for tf in tfbindings:
                            try:
                                heatmap_dict_data['TF per Chromatin State'][state][tf].append(sample_id)
                            except KeyError:
                                try:
                                    heatmap_dict_data['TF per Chromatin State'][state][tf] = [sample_id]
                                except KeyError:
                                    try:
                                        heatmap_dict_data['TF per Chromatin State'][state] = {tf : [sample_id]}
                                    except KeyError:
                                        heatmap_dict_data['TF per Chromatin State'] = {state: {tf : [sample_id]}}
    
    return heatmap_dict_data, box_plot_dict_data


def generate_RegMuts_dicts(elements_input, cols_names_to_use, count_one_mut_per_sample_per_element):
    
    box_plot_dict_data = {'Reg. Mutations per Cancer Type': {'Type':'Samples'}, 'Reg. Mutation Scores per Cancer Type': {'Type':'Values'}}
    for col in cols_names_to_use:
        for element in elements_input[cols_names_to_use[col]]:
            samples_in_element = []
            
            items = element.split(',')#muts in this element
            for item in items:
                cancer_type = item.split('#')[5]
                sample_id = item.split('#')[8]
                mut_score = float(item.split('#')[9])
                sample_exists = False
                if sample_id in samples_in_element:
                    sample_exists = True
                else:
                    samples_in_element.append(sample_id)
                try:
                    if count_one_mut_per_sample_per_element: 
                        if not sample_exists:
                            box_plot_dict_data['Reg. Mutations per Cancer Type'][cancer_type].append(sample_id)
                            box_plot_dict_data['Reg. Mutation Scores per Cancer Type'][cancer_type].append(mut_score)
                    else:
                        box_plot_dict_data['Reg. Mutations per Cancer Type'][cancer_type].append(sample_id)
                        box_plot_dict_data['Reg. Mutation Scores per Cancer Type'][cancer_type].append(mut_score)
                except KeyError:
                    box_plot_dict_data['Reg. Mutations per Cancer Type'][cancer_type] = [sample_id]
                    box_plot_dict_data['Reg. Mutation Scores per Cancer Type'][cancer_type] = [mut_score]
    return box_plot_dict_data

def generate_RegMaxMotif_dicts(elements_input, cols_names_to_use, count_one_mut_per_sample_per_element, min_numnber_muts_to_consider_element=100):
    
    heatmap_dict_data = {'Reg. Mutations in Chromatin States per Cancer Type':{'X-label': 'Chromatin States', 'Y-label': 'Number of Mutations'}}
    for col in cols_names_to_use:
        for element in elements_input[cols_names_to_use[col]]:
            samples_in_element = []
            
            items = element.split(',')#muts in this element
            for item in items:
                cancer_type = item.split('#')[3]
                sample_id = item.split('#')[4]
                motif_name = item.split('#')[10].split('_')[0]
                chromatin_state = item.split('#')[15]
                sample_exists = False
                if sample_id+motif_name in samples_in_element:#count one sample per motif per element
                    sample_exists = True
                else:
                    samples_in_element.append(sample_id+motif_name)

                try:
                    if count_one_mut_per_sample_per_element:
                        if not sample_exists:
                            heatmap_dict_data['Reg. Mutations in Chromatin States per Cancer Type'][cancer_type][chromatin_state].append(sample_id)
                    else:
                        heatmap_dict_data['Reg. Mutations in Chromatin States per Cancer Type'][cancer_type][chromatin_state].append(sample_id)
                except KeyError:
                    try:
                        heatmap_dict_data['Reg. Mutations in Chromatin States per Cancer Type'][cancer_type][chromatin_state] = [sample_id]
                    except KeyError:
                        try:
                            heatmap_dict_data['Reg. Mutations in Chromatin States per Cancer Type'][cancer_type] = {chromatin_state: [sample_id]}
                        except KeyError:
                            heatmap_dict_data['Reg. Mutations in Chromatin States per Cancer Type'] = {cancer_type : {chromatin_state : [sample_id]}}
                
    return heatmap_dict_data

def generate_RegMotifs_dicts(elements_input, cols_names_to_use, count_one_mut_per_sample_per_element, min_numnber_muts_to_consider_element=100):
    
    heatmap_dict_data = {'Mutated-Motifs per Cancer Type': {'Minimum':min_numnber_muts_to_consider_element}, 
                         'Mutated-Motifs per Chromatin Status': {'X-label': 'Chromatin States', 'Y-label': 'Number of Mutations', 'Minimum':min_numnber_muts_to_consider_element, 'CountPlot': 'Yes'}}
    
    motif_names_to_plot_their_scores_per_cancer_type = ['CTCF', 'SP1', 'ELF1', 'Klf4', 'CEBPB', 'ZNF263']
    box_plot_dict_data = {'Scores of Mutated-Motifs in Sig. Elements per Cancer Type':{'Type':'Values'},
                          'Scores of Mutated-Motifs in Sig. Elements per TF':{'Type':'Values', 'Minimum': min_numnber_muts_to_consider_element, 'X-label': 'TF-Motifs'}}#,
    
    count_plots_dict_data = {'All Motifs': {'Type': 'Samples'}, 'All Motif Positions': {'Type': 'Samples'}}
    for motif_name in motif_names_to_plot_their_scores_per_cancer_type:
        box_plot_dict_data['Scores in Mutated-Motifs of ' + motif_name] = {'Type': 'Values', 'CountPlot':'Yes'}
    
    for motif_name in motif_names_to_plot_their_scores_per_cancer_type:
        box_plot_dict_data['Scores in Mutated-Motif Positions of ' + motif_name] = {'Type': 'Values', 'CountPlot':'Yes'}
    
    for col in cols_names_to_use:
        for element in elements_input[cols_names_to_use[col]]:
            samples_in_element = []
            
            items = element.split(',')#muts in this element
            for item in items:
                cancer_type = item.split('#')[3]
                sample_id = item.split('#')[4]
                mut_score = float(item.split('#')[0])
                motif_name = item.split('#')[10].split('_')[0]
                mut_pos = item.split('#')[6]
                mut_sig = item.split('#')[5]
                mut_breaking_score = float(item.split('#')[2])
                chromatin_state = item.split('#')[15]
                sample_exists = False
                if sample_id+motif_name in samples_in_element:#count one sample per motif per element
                    sample_exists = True
                else:
                    samples_in_element.append(sample_id+motif_name)
                
                try:
                    if count_one_mut_per_sample_per_element:
                        if not sample_exists:
                            heatmap_dict_data['Mutated-Motifs per Cancer Type'][cancer_type][motif_name].append(sample_id)
                    else:
                        heatmap_dict_data['Mutated-Motifs per Cancer Type'][cancer_type][motif_name].append(sample_id)
                except KeyError:
                    try:
                        heatmap_dict_data['Mutated-Motifs per Cancer Type'][cancer_type][motif_name] = [sample_id]
                    except KeyError:
                        try:
                            heatmap_dict_data['Mutated-Motifs per Cancer Type'][cancer_type] = {motif_name: [sample_id]}
                        except KeyError:
                            heatmap_dict_data['Mutated-Motifs per Cancer Type'] = {cancer_type : {motif_name : [sample_id]}}
                
                try:
                    if count_one_mut_per_sample_per_element:
                        if not sample_exists:
                            heatmap_dict_data['Mutated-Motifs per Chromatin Status'][chromatin_state][motif_name].append(sample_id)
                    else:
                        heatmap_dict_data['Mutated-Motifs per Chromatin Status'][chromatin_state][motif_name].append(sample_id)
                except KeyError:
                    try:
                        heatmap_dict_data['Mutated-Motifs per Chromatin Status'][chromatin_state][motif_name] = [sample_id]
                    except KeyError:
                        try:
                            heatmap_dict_data['Mutated-Motifs per Chromatin Status'][chromatin_state] = {motif_name: [sample_id]}
                        except KeyError:
                            heatmap_dict_data['Mutated-Motifs per Chromatin Status'] = {chromatin_state : {motif_name : [sample_id]}}
                
                try:
                    if count_one_mut_per_sample_per_element: 
                        if not sample_exists:
                            box_plot_dict_data['Scores of Mutated-Motifs in Sig. Elements per Cancer Type'][cancer_type].append(mut_score)
                    else:
                        box_plot_dict_data['Scores of Mutated-Motifs in Sig. Elements per Cancer Type'][cancer_type].append(mut_score)
                except KeyError:
                    box_plot_dict_data['Scores of Mutated-Motifs in Sig. Elements per Cancer Type'][cancer_type] = [mut_score]
                   
                try:
                    if count_one_mut_per_sample_per_element: 
                        if not sample_exists:
                            box_plot_dict_data['Scores of Mutated-Motifs in Sig. Elements per TF'][motif_name].append(mut_score)
                    else:
                        box_plot_dict_data['Scores of Mutated-Motifs in Sig. Elements per TF'][motif_name].append(mut_score)
                except KeyError:
                        box_plot_dict_data['Scores of Mutated-Motifs in Sig. Elements per TF'][motif_name] = [mut_score]
                
                if motif_name in motif_names_to_plot_their_scores_per_cancer_type:
                    try:
                        if count_one_mut_per_sample_per_element: 
                            if not sample_exists:
                                box_plot_dict_data['Scores in Mutated-Motif Positions of ' + motif_name][motif_name+'#'+mut_pos].append(mut_score)
                        else:
                            box_plot_dict_data['Scores in Mutated-Motif Positions of ' + motif_name][motif_name+'#'+mut_pos].append(mut_score)
                    except KeyError:
                        box_plot_dict_data['Scores in Mutated-Motif Positions of ' + motif_name][motif_name+'#'+mut_pos] = [mut_score]
                
                try:
                    if count_one_mut_per_sample_per_element:
                        if not sample_exists:
                            count_plots_dict_data['All Motif Positions'][motif_name+'#'+mut_pos].append(sample_id)
                    else:
                        count_plots_dict_data['All Motif Positions'][motif_name+'#'+mut_pos].append(sample_id)
                except KeyError:
                    count_plots_dict_data['All Motif Positions'][motif_name+'#'+mut_pos] = [sample_id]
                
                #MOTIF Positions per TF motif
                if motif_name in motif_names_to_plot_their_scores_per_cancer_type:
                    try:
                        if count_one_mut_per_sample_per_element: 
                            if not sample_exists:
                                box_plot_dict_data['Scores in Mutated-Motifs of ' + motif_name][cancer_type].append(mut_score)
                        else:
                            box_plot_dict_data['Scores in Mutated-Motifs of ' + motif_name][cancer_type].append(mut_score)
                    except KeyError:
                        box_plot_dict_data['Scores in Mutated-Motifs of ' + motif_name][cancer_type] = [mut_score]
                
                try:
                    if count_one_mut_per_sample_per_element:
                        if not sample_exists:
                            count_plots_dict_data['All Motifs'][motif_name].append(sample_id)
                    else:
                        count_plots_dict_data['All Motifs'][motif_name].append(sample_id)
                except KeyError:
                    count_plots_dict_data['All Motifs'][motif_name] = [sample_id]
                
                    
                try:
                    if count_one_mut_per_sample_per_element:
                        if not sample_exists:
                            count_plots_dict_data[chromatin_state][motif_name].append(sample_id)
                    else:
                        count_plots_dict_data[chromatin_state][motif_name].append(sample_id)
                except KeyError:
                    try:
                        count_plots_dict_data[chromatin_state][motif_name] = [sample_id]
                    except KeyError:
                        count_plots_dict_data[chromatin_state] = {motif_name: [sample_id], 'Type':'Samples'}
                        
                
                    
    return box_plot_dict_data, heatmap_dict_data, count_plots_dict_data

def elements_fdr(elements_input, output_file_ext):
    elements_input['FDR'] = np.where(elements_input['FDR']==0.0, 1e-300, elements_input['FDR'])#to replace zero with a number that can be converted to log10
    elements_input['FDR'] = elements_input['FDR'].apply(lambda x: np.log10(x)*-1)
    x = '#Muts'
    y = 'FDR'#'#RegMuts'
    key_name = 'Nearby-Genes(Name::ID::O|U|Ddistance::COSMIC|KCP|PCD)'#'Position'
    plt_results = sns.lmplot(x, y, data=elements_input, fit_reg=False, legend=False, legend_out=False, scatter=True, markers='o', scatter_kws={"s": 1})#, col='Cohorts'
    plt_results.set_axis_labels(x, y).set(ylim=[0,elements_input[y].max()], xlim=[elements_input[x].min(), elements_input[x].max()+10])
    sns.despine(right=True, top=True)
    plt.figure(figsize=(6, 6))
    regions_written = []
    num_labels_to_write = 20
    min_value_to_write_text = 10
    region_name_key = 'Position'
    filter_col_for_text = '#RegMuts'
    for i, row in elements_input.sort_values(by=y, ascending=False).iterrows():#[0:num_labels_to_write]
        key_value = '|'.join([v.split('::')[0] for v in row[key_name].split(',')])
        if row[filter_col_for_text]>=min_value_to_write_text:
            plt_results.ax.text(row[x] ,row[y], key_value, size='xx-small', color='darkslategrey', rotation=5, verticalalignment='bottom')
            regions_written.append(row[region_name_key])
            print key_value
    '''for i, row in elements_input.sort_values(by=x, ascending=False)[0:num_labels_to_write].iterrows():#len(regions_written)
        key_value = '|'.join([v.split('::')[0] for v in row[key_name].split(',')])
        if row[region_name_key] not in regions_written:
            plt_results.ax.text(row[x] ,row[y], key_value, size='xx-small', color='darkslategrey', rotation=25, verticalalignment='bottom')
            regions_written.append(row[region_name_key])
    '''
    plt_results.savefig(output_file_ext +".pdf")
    plt.close()

    
    return
def read_sig_elements_file(elements_input_file, output_file_ext="Fig"):
    elements_input = pd.read_table(elements_input_file, sep='\t', skiprows=0, header=0)
    elements_input = elements_input[(elements_input['#Samples(RegMuts)']>1)]
    print len(elements_input)
    print elements_input.head(2)
    elements_fdr(elements_input, output_file_ext+"_elementsFDRMuts")
    cols_names_to_use_for_barchars = {'Reg. Mutations in Sig Elements per Cancer Type': 'Cancer-Types:#RegMuts',
                                      'Unique Reg. Mutations in Sig Elements per Cancer Type': 'Cancer-Types:#Samples(RegMuts)',
                                      'Mutations in Sig. Elements per Cancer Type': 'Cancer-Types:#Muts',
                                      'Unique Mutations in Sig. Elements per Cancer Type': 'Cancer-Types:#Samples'
                                      }
    #plotting muts
    bar_chars_data = genearate_bar_char_dict(elements_input, cols_names_to_use_for_barchars)
    plot_count_from_dict(bar_chars_data, output_file=output_file_ext)
    
    cols_names_to_use_for_Muts = {'Muts Annotations':'Muts'}
    #Number of mutations
    heatmap_dict_data, box_plot_dict_data = generate_Muts_heatmap_dict(elements_input, cols_names_to_use=cols_names_to_use_for_Muts, count_one_mut_per_sample_per_element=False, min_numnber_muts_to_consider_element_for_number_of_tfs=1000)
    plot_enrichment(heatmap_dict_data, output_file=output_file_ext+"AllMuts_", count_uniq_items=False)
    
    #Number of unique mutations in each element
    heatmap_dict_data, box_plot_dict_data = generate_Muts_heatmap_dict(elements_input, cols_names_to_use=cols_names_to_use_for_Muts, count_one_mut_per_sample_per_element=True, min_numnber_muts_to_consider_element_for_number_of_tfs=1000)
    plot_enrichment(heatmap_dict_data, output_file=output_file_ext+"UniqMutsPerElement_", count_uniq_items=False)
    
    #Number of unique samples with mutation per item
    heatmap_dict_data, box_plot_dict_data = generate_Muts_heatmap_dict(elements_input, cols_names_to_use=cols_names_to_use_for_Muts, count_one_mut_per_sample_per_element=True)
    plot_enrichment(heatmap_dict_data, output_file=output_file_ext+"NumberSamples_", count_uniq_items=True)
    
    plot_element_mutation_rate_per_cancer_type(box_plot_dict_data, output_file_ext)
    
    cols_names_to_use_for_regMuts = {'Reg. Mutations':'RegMuts'}
    box_plot_dict_data = generate_RegMuts_dicts(elements_input, cols_names_to_use=cols_names_to_use_for_regMuts, count_one_mut_per_sample_per_element=False)
    plot_element_mutation_rate_per_cancer_type(box_plot_dict_data, output_file_ext+'AllRegMutsCounted_')
    
    box_plot_dict_data = generate_RegMuts_dicts(elements_input, cols_names_to_use=cols_names_to_use_for_regMuts, count_one_mut_per_sample_per_element=True)
    plot_element_mutation_rate_per_cancer_type(box_plot_dict_data, output_file_ext+'UniqueRegMutsPerElement_')
    
    cols_names_to_use_for_mutated_motifs = {'Reg. Mutations':'Mutated-Moitfs'}
    box_plot_dict_data, heatmap_dict_data, count_plots_dict_data = generate_RegMotifs_dicts(elements_input, cols_names_to_use=cols_names_to_use_for_mutated_motifs, count_one_mut_per_sample_per_element=True, min_numnber_muts_to_consider_element=200)
    plot_element_mutation_rate_per_cancer_type(box_plot_dict_data, output_file_ext)
    plot_enrichment(heatmap_dict_data, output_file=output_file_ext, count_uniq_items=False)
    plot_count_from_dict(count_plots_dict_data, output_file=output_file_ext, top_n=20)
    
    cols_names_to_use_for_mutated_motifs = {'Reg. Muts in Chromatin States per Cancer Type':'Max-RegMotif'}
    heatmap_dict_data = generate_RegMaxMotif_dicts(elements_input, cols_names_to_use=cols_names_to_use_for_mutated_motifs, count_one_mut_per_sample_per_element=True)
    plot_enrichment(heatmap_dict_data, output_file=output_file_ext, count_uniq_items=False)
    
    
def read_sig_tfs(sig_tfs_file, output_file_ext="Fig", key_name='TFs', threshold_to_include_element=100, fdr_log10neg_threshold=2):
    sig_tfs = pd.read_table(sig_tfs_file, sep='\t', header=0)#, names=[u'TFs', u'P-Val', u'FDR', u'#Mutated Motifs', u'MeanSim', u'#Mutated Motifs in Simulated Sets', u'Cohorts'])
    sig_tfs = sig_tfs[(sig_tfs['#Mutated Motifs']>10)]
    sig_tfs['FDR'] = np.where(sig_tfs['FDR']==0.0, 1e-300, sig_tfs['FDR'])#to replace zero with a number that can be converted to log10
    sig_tfs['FDR'] = sig_tfs['FDR'].apply(lambda x: np.log10(x)*-1)
    
    if 'Position' in key_name:
        sig_tfs[key_name] = sig_tfs[key_name].apply(lambda x: x.split('_')[0]+"#"+x.split('#')[-1])
    else:
        sig_tfs[key_name] = sig_tfs[key_name].apply(lambda x: x.split('_')[0])
    #sig_tfs = sig_tfs[(sig_tfs['#Mutated Motifs']>50)]
    print sig_tfs.columns
    print len(sig_tfs)
    print sig_tfs.head(2)
    df_pivot = sig_tfs.pivot('Cohorts', key_name, 'FDR')#.fillna(0)
    df_pivot_annot = sig_tfs.pivot('Cohorts', key_name, '#Mutated Motifs')#.fillna(0)
    df_pivot_annot_filtered = pd.DataFrame()
    df_pivot_filtered = pd.DataFrame()
    for c in df_pivot.columns:
        if df_pivot[c].max()>fdr_log10neg_threshold and df_pivot_annot[c].max()>threshold_to_include_element:
            df_pivot_filtered[c] = df_pivot[c]
            df_pivot_annot_filtered[c] = df_pivot_annot[c]
    print 'df_pivot_annot_filtered.shape: ', df_pivot_annot_filtered.shape
    #df_pivot = df_pivot[df_pivot.sum(axis=0) > 20]
    if df_pivot_filtered.shape[0]>0 and df_pivot_filtered.shape[1]>0:
        plt.figure(figsize=(25, 10))
        plt_results = sns.heatmap(df_pivot_filtered, annot=df_pivot_annot_filtered, fmt='g', linewidths=.5, annot_kws={"size": 12, "color":'white'}, vmin=0, vmax=10)
        plt_results.axes.set_title("Enrichment of mutations at {}".format(key_name))
        plt_results.axes.set_ylabel('Cohorts')
        plt_results.axes.set_xlabel(key_name)
        fig = plt_results.get_figure()
        fig.savefig(output_file_ext+ key_name +".pdf")
    plt.close()
    #plot_enrichment(dict_data, output_file=output_file_ext, count_uniq_items=False)

def plot_sig_tfs(sig_tfs_file, output_file_ext="Fig", key_name='TFs', threshold_to_include_element=100, fdr_log10neg_threshold=2, cohorts=['All-tumors-without-Lymphatic-system-Skin-Melanoma'],
                 num_labels_to_write=10):
    
    sig_tfs = pd.read_table(sig_tfs_file, sep='\t', header=0)
    #sig_tfs = sig_tfs[sig_tfs['#Mutated Motifs']>500]
    
    sig_tfs['Significance'] = np.where(sig_tfs['FDR']<0.05, 'Significant', 'Not_Significant')
    sig_tfs['FDR'] = np.where(sig_tfs['FDR']==0.0, 1e-300, sig_tfs['FDR'])#to replace zero with a number that can be converted to log10
    sig_tfs['FDR'] = sig_tfs['FDR'].apply(lambda x: math.log(x,10)*-1)
    if 'Position' in key_name:
        sig_tfs[key_name] = sig_tfs[key_name].apply(lambda x: x.split('_')[0]+"#"+x.split('#')[-1])
    else:
        sig_tfs[key_name] = sig_tfs[key_name].apply(lambda x: x.split('_')[0])
    
    for cohort in cohorts:
        plt.figure(figsize=(6, 6))
        sig_tfs_chort = sig_tfs[(sig_tfs['Cohorts']==cohort)]# | (sig_tfs['Cohorts']=='CNS-Medullo')]
        plt_results = sns.lmplot(x='#Mutated Motifs', y='FDR', data=sig_tfs_chort, fit_reg=False, legend=False, legend_out=False, scatter=True, markers=["*", "*"], hue='Significance', 
                                 palette=dict(Significant="r", Not_Significant="b"))#, col='Cohorts'
        plt_results.set_axis_labels("#Mutated Motifs", "FDR (-log10)").set(ylim=[0,sig_tfs_chort['FDR'].max()], xlim=[sig_tfs_chort['#Mutated Motifs'].min(), sig_tfs_chort['#Mutated Motifs'].max()+500])
        sns.despine(right=True, top=True)
        plt.figure(figsize=(6, 6))
        tfs_written = []
        for i, row in sig_tfs_chort.sort_values(by='FDR', ascending=False)[0:num_labels_to_write].iterrows():
            if row['Significance']=='Significant':
                plt_results.ax.text(row['#Mutated Motifs'] ,row['FDR'], row[key_name], size='xx-small', color='darkslategrey', rotation=25, verticalalignment='bottom')
                tfs_written.append(row[key_name])
        
        for i, row in sig_tfs_chort.sort_values(by='#Mutated Motifs', ascending=False)[0:num_labels_to_write].iterrows():
            if row[key_name] not in tfs_written:
                plt_results.ax.text(row['#Mutated Motifs'] ,row['FDR'], row[key_name], size='xx-small', color='darkslategrey', rotation=25, verticalalignment='bottom')
                tfs_written.append(row[key_name])
        #plt_results.axes.set_title("Significance of {}".format(key_name))
        #plt_results.axes.set_ylabel('#Mutated Motifs')
        #plt_results.axes.set_xlabel(key_name)
        
        #fig = plt_results.get_figure()
        plt_results.savefig(output_file_ext+ key_name +"regplot.pdf")
        plt.close()

def plot_motifs(input_file, x_label='TF-Motifs', y_label='Cancer Types', output_file_ext="motfsheatmap", threshold_to_include_element = 1000):
    
    sig_motifs = pd.read_table(input_file, sep='\t', header=None, usecols=[5,10,17,22,24,30,31], names=['Cancer Types', 'MBR', 'TF-Motifs', 'Chromatin States', 'DNase1', 'TFBinding', 'TFExpr'])
    sig_motifs = sig_motifs[((sig_motifs['MBR']>=0.3) & ((sig_motifs['TFExpr'].apply(math.isnan)) | (sig_motifs['TFExpr']>0)))]
    
    if 'Position' in x_label:
        sig_motifs[x_label] = sig_motifs[x_label].apply(lambda x: x.split('_')[0]+"#"+x.split('#')[-1])
    else:
        sig_motifs[x_label] = sig_motifs[x_label].apply(lambda x: x.split('_')[0])
    
    dict_data = {}
    for i, row in sig_motifs.iterrows():
        try:
            dict_data[row[y_label]][row[x_label]]+=1
        except KeyError:
            try:
                dict_data[row[y_label]][row[x_label]] = 1
            except KeyError:
                dict_data[row[y_label]] = {row[x_label]: 1}
    y_labels = []
    x_labels = []
    values = []
    values_label = '#Mutated Motifs'
    for y in sorted(dict_data.keys()):
        for x in sorted(dict_data[y].keys()):
            y_labels.append(y)
            x_labels.append(x)
            values.append(dict_data[y][x])
    df = pd.DataFrame({y_label: y_labels , x_label: x_labels, values_label: values})
    print df.columns, df.shape
    print df.head(10)
    df_pivot = df.pivot(y_label, x_label, values_label)#.fillna(0)
    df_pivot_filtered = pd.DataFrame()
    for c in df_pivot.columns:
        if df_pivot[c].sum()>threshold_to_include_element:
            df_pivot_filtered[c] = df_pivot[c]
            
    print df_pivot_filtered.head(4)
    print 'df_pivot_filtered.shape: ', df_pivot_filtered.shape
    #df_pivot = df_pivot[df_pivot.sum(axis=0) > 20]
    if df_pivot_filtered.shape[0]>0 and df_pivot_filtered.shape[1]>0:
        plt.figure(figsize=(25, 10))
        plt_results = sns.heatmap(df_pivot_filtered, annot=True, fmt='g', linewidths=.5, annot_kws={"size": 12, "color":'grey'})
        plt_results.axes.set_title("Enrichment of mutations at {}".format(x_label))
        plt_results.axes.set_ylabel(y_label)
        plt_results.axes.set_xlabel(x_label)
        plt.xticks(rotation=45)
        fig = plt_results.get_figure()
        fig.savefig(output_file_ext+ y_label +".pdf")
    plt.close()
    
    return

def get_info_from_sigelements(elements_input_file):
    elements_input = pd.read_table(elements_input_file, sep='\t', skiprows=6, header=0)
    elements_input_filtered = elements_input[(elements_input['#Samples(RegMuts)']>1)]
    sorted_elements_by_fdr = elements_input_filtered.sort_values(by='FDR', ascending=True)
    c = 0
    for i, r in sorted_elements_by_fdr.iterrows():
        if 'WDR74' in r['Nearby-Genes(Name::ID::O|U|Ddistance::COSMIC|KCP|PCD)']:
            print c,r
        c+=1

if __name__ == '__main__':
    
    sns.set_style("white", {"axes.linewidth": ".5"})
    #generated_sig_merged_element_files = process_cohorts()
    #sig_elements_file = getSigElements(generated_sig_merged_element_files)
    elements_input_file = '/home/huum/projs/regMotifs/analysis/Karolina_New_results/mutHIM_80_out.tsv'
#'/home/huum/projs/regMotifs/analysis/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv'
    sig_tfs_file = '/home/huum/projs/regMotifs/analysis/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_sigTFs_0.05.tsv'
    sig_tfpos_file = '/home/huum/projs/regMotifs/analysis/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_sigTFpos_0.05.tsv'
    input_file = '/home/huum/projs/regMotifs/mutations_cohorts_output/All-tumors-without-Lymphatic-system-Skin-Melanoma_observed_annotated_agreement_22May2017.bed9_rand103setsTFsigQval0.05'
    
    output_file_ext="/home/huum/projs/regMotifs/analysis/Fig_"
    read_sig_elements_file(elements_input_file, output_file_ext)
    
    output_file_ext = '/home/huum/projs/regMotifs/analysis/sigTFs_'
    read_sig_tfs(sig_tfs_file, output_file_ext, key_name='TFs', threshold_to_include_element=500, fdr_log10neg_threshold=3)
    plot_sig_tfs(sig_tfs_file, output_file_ext=output_file_ext+'_lmplot', key_name='TFs', threshold_to_include_element=100, fdr_log10neg_threshold=2, 
                 cohorts=['All-tumors-without-Lymphatic-system-Skin-Melanoma'], num_labels_to_write=10)
    
    output_file_ext = '/home/huum/projs/regMotifs/analysis/sigTFpos_'
    read_sig_tfs(sig_tfpos_file, output_file_ext, key_name='TF Positions', threshold_to_include_element=100, fdr_log10neg_threshold=3)
    
    plot_sig_tfs(sig_tfpos_file, output_file_ext=output_file_ext+"_lmplot", key_name='TF Positions', threshold_to_include_element=100, fdr_log10neg_threshold=2, 
                 cohorts=['All-tumors-without-Lymphatic-system-Skin-Melanoma'], num_labels_to_write=10)
    #plot_mutation_rate_per_cancer_type(output_files_mutation_rate="analysis/mutationratepersamplepercancertype", observed_mutation_file='mutations_files/observed.bed9')
    
    output_file_ext = '/home/huum/projs/regMotifs/analysis/FuncMotifsAllTumorsNoLymphSKin_'
    plot_motifs(input_file, output_file_ext=output_file_ext+"_motifsheatmap")
    plot_motifs(input_file, x_label='TF-Motifs', y_label='Chromatin States', output_file_ext=output_file_ext+"_motifsheatmap", threshold_to_include_element = 1000)
    
    '''
    #Heat-Maps for:
        - #muts|regmotifs|samples|regsamples per cancer type per chromatin state
        - #muts|regmotifs|samples|regsamples per cancer type per TF-motif
        - #muts|regmotifs|samples|regsamples per chromatin state per TF-motif
         
    #Box-Plots
        - Score distribution in sig elements (observed vs. simulated)
        - Score distribution in motifs of TFs with highest #mutated-motifs (observed vs. combined-simulated)
        - Number of regMuts|muts per sample per cancer types
        
    #Bar-Charts
        - #regMuts|muts|samples|regsamples per cancer type
        - #mutated-motifs per TF per chromatin state (multiple figures)
    
    #Stacked Bar-Charts:
        - #muts|regmotifs|samples|regsamples per cancer type per element (annotate them with associated COSMIC genes)
    
    Genes and Enrichments??
    '''
    #elements_input_file = 'analysis/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist100kb_within500kb_filteredrec1regmut.tsv'
    #get_info_from_sigelements(elements_input_file)
    
