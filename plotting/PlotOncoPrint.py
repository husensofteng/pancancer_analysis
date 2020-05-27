'''
Created on Jun 4, 2017

@author: husensofteng
'''
import matplotlib
from numpy.lib.function_base import average
from collections import Counter
import math
matplotlib.use('Agg')
#matplotlib.use('TkAgg')
from matplotlib.pyplot import tight_layout
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
import matplotlib.patches as patches
from matplotlib import transforms
import matplotlib.patheffects
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, ticker
import pandas as pd
from operator import itemgetter
import sys

from utils import *

def generate_dicts(input_file, groups_cols_file, min_num_to_include_row =50,
                   rows_index = 1, cols_index_reg = 6, cols_index = 7, 
                   measures_index = 10, rows_frequency_index= 5):
    
    #get the group name for each col
    groups_index = 0; col_id_in_groups_file = 1
    cols_groups_dict = {}
    with open(groups_cols_file, 'r') as f:
        l = f.readline().strip().split('\t')
        while l and len(l)>1:
            cols_groups_dict[l[col_id_in_groups_file]] = l[groups_index]
            l = f.readline().strip().split('\t')
    
    rows = []
    cols = []
    rows_cols_types_dict = {}
    rows_measures = {}
    groups_cols_dict = {}
    rows_cols_dict = {}
    
    with open(input_file, 'r') as f:
        l = f.readline().strip().split('\t')
        while l and len(l)>=10:
            if measures_index is not None:
                if float(l[measures_index]) > 0.01:
                    l = f.readline().strip().split('\t')
                    continue
            if int(l[rows_frequency_index]) < min_num_to_include_row:
                l = f.readline().strip().split('\t')
                continue
            all_samples = l[cols_index].split(',')
            all_samples.extend(l[cols_index_reg].split(','))
            
            for sample in set(all_samples):
                try:
                    if sample not in groups_cols_dict[cols_groups_dict[sample]]:
                        groups_cols_dict[cols_groups_dict[sample]].append(sample)
                except KeyError:
                    groups_cols_dict[cols_groups_dict[sample]] = [sample]
                
                try:
                    if sample not in rows_cols_dict[l[rows_index]]:
                        rows_cols_dict[l[rows_index]].append(sample)
                except KeyError:
                    rows_cols_dict[l[rows_index]] = [sample]
                
                color = '.'
                if sample in l[cols_index_reg].split(','):
                    color = '|'
                try:
                    rows_cols_types_dict[l[rows_index]][sample] = color
                except KeyError:
                    rows_cols_types_dict[l[rows_index]] = {sample: color}
                if measures_index is not None:
                    rows_measures[l[rows_index]] = l[measures_index]
            
            l = f.readline().strip().split('\t')
    groups_ordered = []
    #print sorted(groups_cols_dict.iteritems(), key=itemgetter(1))
    for group, group_cols in sorted(groups_cols_dict.viewitems(), key=lambda x: len(x[1]), reverse=True):
        groups_ordered.append(group)
        cols.extend(group_cols)
        cols.extend(['Gap' for i in range(0,5)])
        
    for row, row_cols in sorted(rows_cols_dict.viewitems(), key=lambda x: len(x[1])):
        #if len(row_cols)>=min_num_to_include_row and float(rows_measures[row])<0.01:
        rows.append(row)
        
    return rows, cols, rows_cols_types_dict, rows_measures, cols_groups_dict, groups_ordered, groups_cols_dict

def plot_oncoprint(ax, input_file, groups_cols_file, x_shift=400, min_num_to_include_row =265, rows_index = 1, cols_index_reg = 6, cols_index = 7, 
                   measures_index = 10, rows_frequency_index= 5, fontsize=10, plot_color_legend=False, extra_x_shift = 15, extra_x_shift_others = 40, min_num_samples_to_write_group=10,
                   max_num_pathways_to_draw=None):
    
    rows, cols, rows_cols_types_dict, rows_measures, cols_groups_dict, groups_ordered, groups_cols_dict = generate_dicts(
                                                                                       input_file, groups_cols_file, min_num_to_include_row =min_num_to_include_row,
                                                                                       rows_index = rows_index, cols_index_reg = cols_index_reg, cols_index = cols_index, 
                                                                                       measures_index = measures_index, rows_frequency_index= rows_frequency_index)
    #chaned the color of Lung-SCC: from #FDF5E6 to #473604  and Lung-AdenoCA: from #FFFFFF to #440447
    groups_colors_dict = {'Biliary-AdenoCA':'#00CD66','Bladder-TCC':'#EEAD0E','Bone-Osteosarc':'#FFD700','Bone-Leiomyo':'#FFEC8B','Bone-Epith':'#ADAC44','Breast-AdenoCa':'#CD6090','Cervix-SCC':'#79CDCD','CNS-Medullo':'#D8BFD8','CNS-PiloAstro':'#B0B0B0','CNS-GBM':'#3D3D3D','CNS-Oligo':'#787878','ColoRect-AdenoCA':'#191970','Eso-AdenoCa':'#1E90FF','Head-SCC':'#8B2323','Kidney-RCC':'#FF4500','Kidney-ChRCC':'#B32F0B','Liver-HCC':'#006400','Lung-SCC':'#473604','Lung-AdenoCA':'#440447','Lymph-BNHL':'#698B22','Lymph-CLL':'#F4A35D','Myeloid-MPN':'#FFC100','Myeloid-AML':'#CD6600','Ovary-AdenoCA':'#008B8B','Panc-AdenoCA':'#7A378B','Panc-Endocrine':'#E066FF','Prost-AdenoCA':'#87CEFA','Skin-Melanoma':'#000000','Stomach-AdenoCA':'#BFEFFF','Thy-AdenoCA':'#9370DB','Uterus-AdenoCA':'#FF8C69','Bone-Cart':'#DDCDCD','Breast-LobularCa':'#DDCDCD','Breast-DCIS':'#DDCDCD','Lymph-NOS':'#DDCDCD','Myeloid-MDS':'#DDCDCD','Cervix-AdenoCA':'#DDCDCD'}
    print (rows_cols_types_dict)
    print (len(cols))
    gap = 0
    if max_num_pathways_to_draw is not None:
        if len(rows) > max_num_pathways_to_draw:
            rows = rows[len(rows)-max_num_pathways_to_draw::]
    
    for y, row in enumerate(rows):
        #measure = rows_measures[row]
        #row_text = row.replace(' signaling pathway','') + ' ({:.1f}%)'.format((len(rows_cols_types_dict[row])/(len(cols)*1.0))*100)
        row_text = row.replace(' signaling pathway','').replace('misregulation in cancer','misregulation') + ' ({:.1f}%)'.format((len(rows_cols_types_dict[row])/(2515.0*1.0))*100)
        draw_text(ax, x=0, y=y+gap, text=row_text, horizontalalignment='left', fontsize=fontsize)
        for col in rows_cols_types_dict[row]:
            x = cols.index(col)
            color = groups_colors_dict[cols_groups_dict[col]]
            marker = rows_cols_types_dict[row][col]
            size = 10
            if marker=='.':
                size =6
            draw_marker(ax, x=x+x_shift, y=y+gap, color=color, marker=marker, markersize=size)
        gap+=1.0
    y = len(rows)+gap
    #draw_text(ax, x=0, y=y, text="Samples")
    if plot_color_legend:
        groups_written = []
        for x,col in enumerate(cols):
            #print cols_groups_dict[col]
            color = 'white'
            try:
                color =groups_colors_dict[cols_groups_dict[col]]
            except KeyError:
                #print 'sample: ' + col + 'not found in cancer type groups'
                color = 'white'
            draw_marker(ax, x=x+x_shift, y=y, color=color, marker='|', markersize=14)
            
            try:
                group = cols_groups_dict[col]
                color = groups_colors_dict[group]
                rotation = 90
                ext = 0
                if len(groups_cols_dict[group])<min_num_samples_to_write_group:
                    group = 'Others'
                    color = 'grey'
                    rotation = 0
                    ext = extra_x_shift_others
                if group not in groups_written:
                    draw_text(ax, x=x+x_shift+extra_x_shift+ext, y=y+1, text=group, color=color, fontsize=8, horizontalalignment='right', rotation=rotation, verticalalignment='bottom')
                    groups_written.append(group)
            except KeyError:
                pass
    
    legend_x = len(cols) - (len(cols)/4.0)
    legend_y = -3
    draw_text(ax, x=legend_x, y=legend_y, text="CFRMs", color='black', fontsize=10, horizontalalignment='left', rotation=0, verticalalignment='center')
    draw_marker(ax, x=legend_x-20, y=legend_y, marker='|', color='grey', markersize=12)
    draw_text(ax, x=legend_x+200, y=legend_y, text="Other Mutations", color='black', fontsize=10, horizontalalignment='left', rotation=0, verticalalignment='center')
    draw_marker(ax, x=legend_x+200-20, y=legend_y, marker='.', color='grey', markersize=12)
    #gap-=5.0
    #draw_text(ax, x=len(cols)+80, y=y, text="%, FDR", horizontalalignment='center')
    ''''x = 0
    y=-1
    if plot_color_legend:
        for i, group in enumerate(groups_ordered):
            if i % 6 == 0:
                y-=1
                x=0
            draw_text(ax, x=x, y=y, text=group, color=groups_colors_dict[group], fontsize=10, horizontalalignment='left')
            x+=300
    '''
    ax.set_xlim(0,len(cols)+x_shift)
    ax.set_ylim(legend_y-1,20)
    #ax.set_ylim(0,len(rows)+gap+10)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.minorticks_off()
    ax.set_frame_on(False)
    
    return

def draw_pathwyas(input_file, groups_cols_file, output_dir):
    
    #input_file = "/home/huum/projs/regMotifs/analysis_exclVEP/merged200bp_extended200bp_nofullexon_pancan_skip_exon/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv_GenesInclCDS.tsv_pathways_calculated_pval.tsv" 
#"/home/huum/projs/regMotifs/analysis/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv_GenesInclCDS.tsv_pathways_calculated_pval.tsv"
#"../analysis/data/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb_pancan.tsv_Genes.tsv_pathways_calculated_pval_sig.tsv"
    #sort pathways by FDR
    df = pd.read_table(input_file,
                       sep = '\t', header=None, names=
                       ['Pathways','Total Num. genes', 'Num. enriched genes', 'Num. mutated samples with CFRMs', 'Num. mutated samples', 'Enriched Genes', 'Samples', 'Genes', 'P-Value','FDR'])
    df['Enrichment size'] = df['Num. enriched genes']/df['Total Num. genes']
    df.sort_values(by='FDR', inplace=True)
    print ('/'.join(input_file.split('/')[:-1])+'/sigpathways_sorted.tsv')
    df.to_csv('/'.join(input_file.split('/')[:-1])+'/sigpathways_sorted.tsv', sep='\t')
    #groups_cols_file = "/home/huum/projs/regMotifs/cancer_types_samples.txt"
    sns.set_style('white', {'text.color': '.15'})
    rcParams['svg.fonttype'] = 'none'
    fig = plt.figure(figsize=(12.5,3.5))#design a figure with the given size
    gs = gridspec.GridSpec(1, 1, wspace=0.0, hspace=0.0)#create 2 rows and three columns with the given ratio for each
    ax0 = fig.add_subplot(gs[0, 0])
    gs.tight_layout(fig, pad=0, h_pad=0.0, w_pad=0.0)
    plot_oncoprint(ax0, input_file, groups_cols_file, plot_color_legend=True, max_num_pathways_to_draw=10)
    #plt.savefig("../analysis/Fig4PancanT.pdf", bbox_inches='tight')
    fig_dir = output_dir + "/Fig4PancanT_skip_exon.svg"
    plt.savefig(fig_dir, bbox_inches='tight')
    plt.close()     
       
def draw_genes(input_file,groups_cols_file, output_dir):
    
    #input_file = "/home/huum/projs/regMotifs/analysis_exclVEP/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv_GenesInclCDS.tsv"
#"/home/huum/projs/regMotifs/analysis/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv_GenesInclCDS.tsv"
#"analysis/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv_Genes.tsv"
    #groups_cols_file = "/home/huum/projs/regMotifs/cancer_types_samples.txt"
    sns.set_style('white', {'text.color': '.15'})
    rcParams['svg.fonttype'] = 'none'
    fig = plt.figure(figsize=(12,12))#design a figure with the given size
    gs = gridspec.GridSpec(1, 1, wspace=0.0, hspace=0.0)#create 2 rows and three columns with the given ratio for each
    ax0 = fig.add_subplot(gs[0, 0])
    
    fig_dir = output_dir + "FigGenesRegMutsmin10.pdf"
    
    plot_oncoprint(ax0, input_file, groups_cols_file, x_shift=300, min_num_to_include_row =10, rows_index = 0, cols_index_reg = 8, cols_index = 9, 
                   measures_index = None, rows_frequency_index= 3, fontsize=10, plot_color_legend=True, extra_x_shift = 20, extra_x_shift_others = 50, min_num_samples_to_write_group=10)
    plt.savefig(fig_dir, bbox_inches='tight')
    plt.clf()
    
    fig_dir = output_dir + "FigGenesMutsmin25.pdf"

    fig = plt.figure(figsize=(12,6))#design a figure with the given size
    gs = gridspec.GridSpec(1, 1, wspace=0.0, hspace=0.0)#create 2 rows and three columns with the given ratio for each
    ax0 = fig.add_subplot(gs[0, 0])
    
    plot_oncoprint(ax0, input_file, groups_cols_file, x_shift=300, min_num_to_include_row =25, rows_index = 0, cols_index_reg = 8, cols_index = 9, 
                   measures_index = None, rows_frequency_index= 5, fontsize=10, plot_color_legend=True, extra_x_shift = 20, extra_x_shift_others = 50, min_num_samples_to_write_group=10)
    
    plt.savefig(fig_dir, bbox_inches='tight')
    
    plt.close()


def parse_args():
    '''Parse command line arguments'''
    
    parser = argparse.ArgumentParser(description='Plot Onco Print')
    parser.add_argument('--input_file_path', default='', help='')
    parser.add_argument('--input_file_genes', default='', help='')

    parser.add_argument('--groups_cols_file', default='', help='')
    parser.add_argument('--output_dir', default='', help='')
    

    
    return parser.parse_args(sys.argv[1:])


if __name__ == '__main__':
    
    args = parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    
    print("Onco Print")


    draw_pathwyas(args.input_file_path, args.groups_cols_file, args.output_dir)
    draw_genes(args.input_file_genes, args.groups_cols_file, args.output_dir)
    
