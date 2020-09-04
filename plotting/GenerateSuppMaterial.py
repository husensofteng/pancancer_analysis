'''
Created on Jun 8, 2017

@author: husensofteng
'''
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import pybedtools
from pybedtools.bedtool import BedTool
from matplotlib.pyplot import tight_layout
import matplotlib.pyplot as plt
from pylab import gca
import pandas as pd
import math
import numpy as np
from decimal import Decimal
import os, sys
import seaborn as sns
import operator
import argparse
sns.set(style="ticks")
#plt.style.use('ggplot')
sns.set_style("white")
sns.set_context("paper")#talk
from utils import *

    
def get_mut_df(input='', x_col_index=5, y_col_index=8, x_col_name = 'Cancer types', y_col_name='Mutation Frequency'):
    if os.path.isfile(str(input)):
        names = [x_col_name, y_col_name]
        if x_col_index>y_col_index:
            names = [y_col_name, x_col_name]
            
        df = pd.read_table(input, sep='\t', header=None, usecols=[x_col_index, y_col_index], names=names)
        
        if x_col_name=="Chromatin States":
            df[x_col_name] = df[x_col_name].apply(get_unique_state).apply(replace_state)
    return df

def plot_boxplot(df, x_col_name = 'Cancer types', y_col_name='Mutation Frequency', title="",
                 groups_colors_dict=None, order=None, rotation=90, fig_width=8, fig_height=6, log=False
    ):
    plt.clf()
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(111)
        
    df['counts'] = [1 for i in range(0,len(df))]
    dfg = df.groupby(by=[x_col_name,y_col_name])
    counts = []
    if log:
        counts = dfg.size().apply(math.log10).rename('counts').tolist()
    else:
        counts = dfg.size().rename('counts').tolist()
    categories = []
    for c in dfg[x_col_name]:
        categories.append(c[0][0])
    sns.boxplot(y = counts, x=categories, palette=groups_colors_dict, ax=ax, order=order)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation)
    ax.set_ylabel(y_col_name)
    ax.set_title(label=title, loc='left')
    sns.despine()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    
    return fig
    
def plot_heatmap(df, x_col_name, y_col_name, fig_width=8, fig_height=6, title="", rotation=90, threshold_to_include_element=500):
    
    plt.clf()
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(111)
    df_pivot_filtered = pd.DataFrame()
    for c in df.columns:
            if df[c].sum()>threshold_to_include_element:
                df_pivot_filtered[c] = df[c] 
    cbar_ax = fig.add_axes([.75, 0.85, .2, .03])
    sns.heatmap(df_pivot_filtered, ax=ax, square=True, cbar_ax=cbar_ax, cbar=True, cbar_kws={"orientation": "horizontal"}, cmap = "YlGnBu" )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation)
    ax.set_ylabel(y_col_name)
    ax.set_title(label=title, loc='left')
    sns.despine()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    return fig
    

def get_df_from_elements(input_file='', col_to_use='RegMuts', sep='#', x_col_index=5, y_col_index=8, x_col_name = 'Cancer types', y_col_name='Mutation Frequency', 
                 col_to_check='#Samples(RegMuts)', threshold=1):
    elements_input = pd.read_table(input_file, sep='\t', skiprows=6, header=0, usecols=[col_to_check,col_to_use])
    elements_input = elements_input[(elements_input['#Samples(RegMuts)']>threshold)]
    box_plot_list = []
    heatmap_dict = {}
    for i, l in elements_input[col_to_use].iteritems():
        for m in l.split(','):
            x_col_value = m.split(sep)[x_col_index].split('_')[-1]
            if x_col_name=="Chromatin States":
                x_col_value = replace_state(x_col_value.split('_')[-1])
            y_col_value = m.split('#')[y_col_index].split('_')[0]
            box_plot_list.append([x_col_value,y_col_value])
            try:
                heatmap_dict[y_col_value][x_col_value] += 1
            except KeyError:
                try:
                    heatmap_dict[y_col_value][x_col_value] = 1
                except KeyError:
                    heatmap_dict[y_col_value] = {x_col_value: 1}
    box_plot_df = pd.DataFrame(box_plot_list, columns=[x_col_name, y_col_name])
    heatmap_df = pd.DataFrame(heatmap_dict)#.pivot(index=x_col_name, columns=y_col_name, values='counts')
    
    return box_plot_df, heatmap_df

def plot_barcharts(input_file, x_col_index = 17, col0_to_check = 10, col1_to_check = 24, col2_to_check=25, fig_width=6, fig_height=6,
                   tf_motifs_to_include = ['FOXP1', 'ZNF263', 'IRF1'], rotation=0, y_col_name='Mutation Frequency', title=''):
    plt.clf()
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(111)
    df = pd.read_table(input_file, sep='\t', header=None, usecols=[col0_to_check, x_col_index, col1_to_check, col2_to_check], names=['Motif diff score', 'TF motifs', 'DHS', 'TFBS'])
    df['TF motifs'] = df['TF motifs'].apply(lambda x: x.split('_')[0])
    df = df[df['TF motifs'].isin(tf_motifs_to_include)]
    df['TFBS'] = np.where(df['TFBS'].apply(math.isnan), 0.0, df['TFBS'])
    df['TFBS'] = df['TFBS'].apply(float)
    df['Potential Effect'] = np.where((df['Motif diff score']>=0.3) & (df['DHS']>0.0) | (df['TFBS']>0.0), 'Yes', 'No')
    df['count']=1
    sns.barplot(x='TF motifs', y='count', hue='Potential Effect', data=df, estimator=sum, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation)
    ax.set_ylabel(y_col_name)
    ax.set_title(label=title, loc='left')
    sns.despine()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    
    return fig





def draw_rec_sigregs(elements_input_file, title):
    sns.set_style('white', {'text.color': '.15', 'axes.linewidth': 0.5})
    fig = plt.figure(figsize=(6,4), linewidth=1.0)#design a figure with the given size
    ax = fig.add_subplot(111)
    elements_input = pd.read_table(elements_input_file, sep='\t', skiprows=6, header=0)
    df = elements_input[(elements_input['#Samples(RegMuts)']>1)]
    df['FDR'] = np.where(df['FDR']==0.0, 1e-300, df['FDR'])#to replace zero with a number that can be converted to log10
    df['FDR'] = df['FDR'].apply(lambda x: np.log10(x)*-1)
    y = 'Score'
    x = '#Samples'
    s = '#Samples(RegMuts)'
    df[x]=df[x].apply(np.log10)
    df[y]=df[y].apply(np.log10)
    colors = []
    colors_region_types_dict = {'CDS': 'green', 'UTR':'orange', 'intergenic': 'black', 'intronic': 'blue', 'proximal_promoter': 'red'}
    for r in df['Feature_type'].values:
        try:
            colors.append(colors_region_types_dict[r])
        except KeyError:
            colors.append('grey')
    ax.scatter(df[x], df[y], color=colors, s=df[s])
    #df['colors'] = colors
    #for i, r in df.iterrows():
    #    ax.scatter(r[x], r[y], color=r['colors'], s=r[s])
    
    for i,f in enumerate(sorted(colors_region_types_dict.keys())):
        ax.scatter(df[x].min()-0.8, df[y].max()-(i/3.0)-0.5, color=colors_region_types_dict[f])
        draw_text(ax, x=df[x].min()-0.6, y=df[y].max()-(i/3.0)-0.5, text=f.replace('proximal_','').replace('intronic', 'intron'), horizontalalignment='left', fontsize=8)
    labels_to_plot_df = df[(df['#Samples']>5) | (df['#Samples(RegMuts)']>=10)]# (df['Score']>5)]
    print(len(labels_to_plot_df))
    for i, r in labels_to_plot_df.iterrows():
        color = 'black'
        label = r['Nearby-Genes(Name::ID::O|U|Ddistance::COSMIC|KCP|PCD)'].split('::')[0]
        if r['#Samples(RegMuts)']<5:
            color = 'grey'
        rotation = 90
        y_shift = 0.2
        x_shift = -0.05
        if label == 'None':
            label = ''
        draw_text(ax, x=r[x]+x_shift, y=r[y]+y_shift, text=label, horizontalalignment='left', verticalalignment='bottom', color=color, fontsize=8, rotation=rotation)
    ax.set_xlabel('Number of mutated samples (log10)')
    ax.set_ylabel('Regulatory score (log10)')
    ax.set_xlim(0,df[x].max()+0.2)
    ax.set_ylim(0,df[y].max()+0.2)
    #ax.legend(loc='upper center')
    ax.set_title(label=title, loc='left')
    sns.despine(right=True, top=True, bottom=False, left=False)
    
    #ax.legend(loc='upper center')
    return fig



def parse_args():
    '''Parse command line arguments'''
    
    parser = argparse.ArgumentParser(description='Plot Supp Figure')
    parser.add_argument('--mut_input_file', default='', help='')
    parser.add_argument('--mut_anno_input_file', default='', help='')
    parser.add_argument('--motif_mut_input_file', default='', help='')
    parser.add_argument('--elements_input_file', default='', help='')
    parser.add_argument('--elements_input_file_Lymph', default='', help='')
    parser.add_argument('--output_dir', default='', help='')
    

    
    return parser.parse_args(sys.argv[1:])


if __name__ == '__main__':
    
    args = parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    
    print("Generate SuppFig")

    groups_colors_dict = {'Biliary-AdenoCA':'#00CD66','Bladder-TCC':'#EEAD0E','Bone-Osteosarc':'#FFD700','Bone-Leiomyo':'#FFEC8B','Bone-Epith':'#ADAC44','Breast-AdenoCa':'#CD6090','Cervix-SCC':'#79CDCD','CNS-Medullo':'#D8BFD8','CNS-PiloAstro':'#B0B0B0','CNS-GBM':'#3D3D3D','CNS-Oligo':'#787878','ColoRect-AdenoCA':'#191970','Eso-AdenoCa':'#1E90FF','Head-SCC':'#8B2323','Kidney-RCC':'#FF4500','Kidney-ChRCC':'#B32F0B','Liver-HCC':'#006400','Lung-SCC':'#FDF5E6','Lung-AdenoCA':'#FFFFFF','Lymph-BNHL':'#698B22','Lymph-CLL':'#F4A35D','Myeloid-MPN':'#FFC100','Myeloid-AML':'#CD6600','Ovary-AdenoCA':'#008B8B','Panc-AdenoCA':'#7A378B','Panc-Endocrine':'#E066FF','Prost-AdenoCA':'#87CEFA','Skin-Melanoma':'#000000','Stomach-AdenoCA':'#BFEFFF','Thy-AdenoCA':'#9370DB','Uterus-AdenoCA':'#FF8C69','Bone-Cart':'#DDCDCD','Breast-LobularCa':'#DDCDCD','Breast-DCIS':'#DDCDCD','Lymph-NOS':'#DDCDCD','Myeloid-MDS':'#DDCDCD','Cervix-AdenoCA':'#DDCDCD'}
    sfig_num = 1
    fig_dir = args.output_dir + "/SFigures.pdf" 
    with PdfPages(fig_dir) as pdf:
        #Box plot for all mutations per sample per cancer type
        df = get_mut_df(input=args.mut_input_file, x_col_index=5, y_col_index=8, x_col_name = 'Cancer types', y_col_name='Mutation Frequency (log10)')
        #df = get_mut_df(input='/home/huum/projs/regMotifs/mutations_files/obsagr22May2017_exclVEP.bed9', x_col_index=5, y_col_index=8, x_col_name = 'Cancer types', y_col_name='Mutation Frequency')

        #df = get_mut_df(input='/Users/karolinasg/Documents/pcawg/analysis/obsagr22May2017_exclVEP.bed9', x_col_index=5, y_col_index=8, x_col_name = 'Cancer types', y_col_name='Mutation Frequency')
        fig = plot_boxplot(df, groups_colors_dict=groups_colors_dict,  title="", log=True)
        #fig = plot_boxplot(df, groups_colors_dict=groups_colors_dict, title="SFig. {n} Mutation rate across cancer types".format(n=sfig_num), log=TRUE)

        pdf.savefig(fig)
        
        sfig_num+=1
        print(sfig_num)
        #Box plot for all mutations in motifs per sample per cancer type
        df = get_mut_df(input=args.motif_mut_input_file, x_col_index=5, y_col_index=8, x_col_name = 'Cancer types', y_col_name='Mutation Frequency (log10)')
        #df = get_mut_df(input='/home/huum/projs/regMotifs/analysis/motifmuts_all.bed12', x_col_index=5, y_col_index=8, x_col_name = 'Cancer types', y_col_name='Mutation Frequency')

        #df = get_mut_df(input='/Users/karolinasg/Documents/pcawg/analysis/motifmuts_all.bed12', x_col_index=5, y_col_index=8, x_col_name = 'Cancer types', y_col_name='Mutation Frequency')
        fig = plot_boxplot(df, groups_colors_dict=groups_colors_dict, title="", log=True)
        #fig = plot_boxplot(df, groups_colors_dict=groups_colors_dict, title="SFig. {n} Mutation rate in TF motifs across cancer types".format(n=sfig_num), log=TRUE)

        pdf.savefig(fig)
        
        #Skip the chromatin states figure
        #sfig_num+=1
        #print sfig_num
        #Box plot for all mutations in motifs per sample per chromatin state
        #order = 'Tx,Tss,Enh,Repr,NA,Quies'.split(',')
        #groups_colors_dict_chrom=['red', 'green', 'orange', 'brown', 'grey','lightgray']
        #df = get_mut_df(input='/home/huum/projs/regMotifs/analysis/motifmuts_all.bed12', x_col_index=9, y_col_index=8, x_col_name='Chromatin States', y_col_name='Mutation Frequency (log10)')

        #fig = plot_boxplot(df, order=order, groups_colors_dict=groups_colors_dict_chrom, rotation=0,  x_col_name='Chromatin States', y_col_name='Mutation Frequency (log10)', 
         #                  title="SFig. {n} Mutations rate in TF motifs across chromatin states".format(n=sfig_num), log=True)
        #pdf.savefig(fig)
        
        sfig_num+=3
        print(sfig_num)
        #Bar charts per TF-motif (with DHS|TFBS and without)
        motif_muts_file = args.mut_anno_input_file
        #motif_muts_file = '/Users/karolinasg/Documents/pcawg/analysis/obsann22May2017_exclVEP.bed9' 

        fig = plot_barcharts(motif_muts_file, x_col_index = 17, col1_to_check = 24, col2_to_check=25, fig_width=6, fig_height=6, title="")
        #, 
                             #title="SFig. {n} Number of active and inactive mutated-motifs".format(n=sfig_num)
                             #)
        pdf.savefig(fig)
        
        
        
        sfig_num+=1
        print(sfig_num)
        #Box plot for CFRMs in SF-MREs per sample per cancer type
        #elements_input_file = '/Users/karolinasg/Documents/pcawg/NEW_RESULTS_removig_VEP_23_october/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv'
        elements_input_file =args.elements_input_file
        #elements_input_file = '/home/huum/projs/regMotifs/analysis_exclVEP/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv'
        box_plot_df, heatmap_df = get_df_from_elements(elements_input_file, col_to_use='RegMuts', sep='#', x_col_index=5, y_col_index=8, x_col_name = 'Cancer types', y_col_name='Mutation Frequency (log10)', col_to_check='#Samples(RegMuts)', threshold=1)
        fig = plot_boxplot(box_plot_df, groups_colors_dict=groups_colors_dict, rotation=90, log=True, title="")
        
        #fig = plot_boxplot(box_plot_df, groups_colors_dict=groups_colors_dict, rotation=90, title="SFig. {n} regMuts across cancer types".format(n=sfig_num), log=TRUE)
        pdf.savefig(fig)
        
        
        
        
        sfig_num+=1
        print(sfig_num)
        #heat map for RegMotifs in SF-MREs per TF motif per chromatin state
        box_plot_df, heatmap_df = get_df_from_elements(elements_input_file, col_to_use='Mutated-Moitfs', sep='#', 
                                            x_col_index=15, y_col_index=10, x_col_name = 'Chromatin States', y_col_name='TF Motifs', 
                 col_to_check='#Samples(RegMuts)', threshold=1)
        fig = plot_heatmap(heatmap_df, x_col_name='Chromatin States', y_col_name='', title='')#, title="SFig. {n} Enrichment of mutated motifs in mutEs per chromatin state".format(n=sfig_num)
                          # )
        pdf.savefig(fig)
        
        sfig_num+=1
        print(sfig_num)
        #elements_input_file = '/Users/karolinasg/Documents/pcawg/NEW_RESULTS_removig_VEP_23_october/merged200bp_extended200bp_nofullexon_Lymph/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv'
        elements_input_file_Lymph =args.elements_input_file_Lymph
        #elements_input_file = '/home/huum/projs/regMotifs/analysis_exclVEP/merged200bp_extended200bp_nofullexon_Lymph/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv'
        fig = draw_rec_sigregs(elements_input_file_Lymph, title='') #, title='SFig. {n} mutEs identified from Lymphoma cohorts'.format(n=sfig_num)
                              # )
        pdf.savefig(fig)
        
        sfig_num+=1
        plt.close()
        
        
