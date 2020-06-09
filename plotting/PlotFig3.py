'''
Created on Jun 6, 2017

@author: husensofteng
'''
import matplotlib
matplotlib.use('Agg')
from numpy.lib.function_base import average
import math
from matplotlib.pyplot import tight_layout
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
import matplotlib.patches as patches
from matplotlib import transforms
import matplotlib.patheffects
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, ticker
import pandas as pd
import sys

from utils import *


def generate_RegMaxMotif_dicts(elements_input, cols_names_to_use, count_one_mut_per_sample_per_element, min_numnber_muts_to_consider_element=100):
    
    heatmap_dict_data = {'Reg. Mutations in Chromatin States per Cancer Type':{'X-label': 'Chromatin States', 'Y-label': 'Mutation Frequency'}}
    for col in cols_names_to_use:
        for element in elements_input[cols_names_to_use[col]]:
            samples_in_element = []
            
            items = element.split(',')#muts in this element
            for item in items:
                cancer_type = item.split('#')[3]
                sample_id = item.split('#')[4]
                motif_name = item.split('#')[10].split('_')[0]
                chromatin_state = replace_state(item.split('#')[15])
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

def plot_enrichment(ax, fig, dict_data, count_uniq_items, cancer_type_samples_dict):
    threshold_to_include_element = 0
    print dict_data.keys()
    for k in dict_data.keys():
        print 'Plotting ', k
        main_items = []
        secondary_items = []
        number_items = []
        x_label = 'Cancer Type'
        y_label = 'Mutations Frequency'
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
                            #print set(dict_data[k][main_item][secondary_item])
                            print main_item,secondary_item
                            print len(set(dict_data[k][main_item][secondary_item]))
                            print len(cancer_type_samples_dict[main_item])
                            print len(set(dict_data[k][main_item][secondary_item]))/ ((len(cancer_type_samples_dict[main_item])*1.0))
                            number_items.append( len(set(dict_data[k][main_item][secondary_item]))/ ((len(cancer_type_samples_dict[main_item])*1.0)))
                        else:
                            number_items.append(len(dict_data[k][main_item][secondary_item]) / (len(cancer_type_samples_dict[main_item])*1.0))
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
        states_order = 'Tx,Tss,Enh,Repr,Quies'.split(',')
        if df_pivot_filtered.shape[0]>0 and df_pivot_filtered.shape[1]>0:
            cbar_ax = fig.add_axes([0.1, 0.18, .12, .015])
            plt_results = sns.heatmap(df_pivot_filtered, annot=False,  linewidths=.5, ax=ax, square=True, 
                                      xticklabels=states_order, cbar_ax=cbar_ax, cbar=True, cbar_kws={"orientation": "horizontal"}, center=0.5, cmap='Greys')
            plt_results.axes.set_ylabel('')
            plt_results.axes.set_xlabel('')
            plt_results.set_yticklabels(plt_results.get_yticklabels(), rotation=0)
            plt_results.set_xticklabels(plt_results.get_xticklabels(), rotation=90)
            
        else:
            print 'Skipped: ', k, df_pivot_filtered.shape
            

def process_elements(ax, fig, elements_input_file, cancer_type_samples_dict):
    elements_input = pd.read_table(elements_input_file, sep='\t', skiprows=6, header=0)
    elements_input = elements_input[(elements_input['#Samples(RegMuts)']>1)]
    
    cols_names_to_use_for_mutated_motifs = {'Reg. Muts in Chromatin States per Cancer Type':'Max-RegMotif'}
    heatmap_dict_data = generate_RegMaxMotif_dicts(elements_input, cols_names_to_use=cols_names_to_use_for_mutated_motifs, count_one_mut_per_sample_per_element=False)
    plot_enrichment(ax, fig, heatmap_dict_data, count_uniq_items=True, cancer_type_samples_dict=cancer_type_samples_dict)
    

def get_samples_per_cancer_type(caner_samples_input):
    cancer_type_samples_dict = {}
    with open(caner_samples_input, 'r') as f:
        l = f.readline().strip().split('\t')
        while len(l)>1:
            
            try:
                if l[1] not in cancer_type_samples_dict[l[0]]:
                    cancer_type_samples_dict[l[0]].append(l[1])
            except:
                cancer_type_samples_dict[l[0]] = [l[1]]
            l = f.readline().strip().split('\t')
    return cancer_type_samples_dict

def draw_region_context(ax1):
    '''plot a stem line with the mutations in the region (x: mut position,y:frquency at that position)
     plot motif logo in the motif position
     draw rectangles for ChIP-seq peaks and annotate them with the name; match the box color with the matching-motif's box color
     plot a line for dnase1 using bedgraph format where x is the position of each base and y is the value in that position
     draw lines for genes and box for its exons
    '''
    return

def draw_plot():
    
    fig = plt.figure(figsize=(12,8), linewidth=1.0)#design a figure with the given size
    gs = gridspec.GridSpec(3, 2, height_ratios=[4,4,4,4], width_ratios=[4,4,4,4], wspace=0.0, hspace=0.0)#create 4 rows and three columns with the given ratio for each
    #heatmap
    ax0 = fig.add_subplot(gs[0:, 0])
    #context plots
    ax1 = fig.add_subplot(gs[0, 1])
    ax1 = fig.add_subplot(gs[1, 1])
    ax3 = fig.add_subplot(gs[2, 1])
    draw_region_context(ax1)
    #scatter plots
    
    
    gs.tight_layout(fig, pad=2, h_pad=0.0, w_pad=0.0)
    
    #sns.despine(top=True, right= True)
    '''
    elements_input_file = '../analysis/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv'
    caner_samples = '../analysis/cancer_types_samples.txt'
    cancer_type_samples_dict=get_samples_per_cancer_type(caner_samples_input=caner_samples)
    process_elements(ax0, fig, elements_input_file, cancer_type_samples_dict)
    '''
    plt.savefig("../analysis/Fig3.pdf")#, bbox_inches='tight')
    plt.savefig("../analysis/Fig3.svg")#, bbox_inches='tight')
    plt.close()
    
draw_plot()