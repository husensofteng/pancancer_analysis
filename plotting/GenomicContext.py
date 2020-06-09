'''
Created on Jun 16, 2017

@author: husensofteng
'''
import matplotlib
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import sys, os
from scipy import stats
from multiprocessing import Pool
import seaborn as sns
from utils import draw_text,stem_line, draw_marker
from decimal import Decimal
from matplotlib.patches import BoxStyle
sns.set_context("paper", font_scale=1)                                                  
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from pybedtools import BedTool
import os, sys
import argparse
from GeneExprAnalysis import get_sample_data, read_genes_elements, read_gene_expr, get_expr_per_sample, process_gene_counts, box_plot_per_gene_cancertype

def get_gene_features(infile, gene_name='',#start=10182000, end=10194000,
               features_to_show = ['gene', 'exon', 'proximal_promoter','UTR', 'start_codon', 'stop_codon', 'CDS'], 
               status_to_show=['KNOWN'], bio_types_to_show=['protein_coding', 'lincRNA'],
               anno_file_header = ['chr', 'start', 'end', 'feature', 'source', 'strand', 'ID', 'ID::name', 'biotype', 'status']):
    
    anno_df = pd.read_table(infile, sep='\t', header=None, names=anno_file_header)
    
    anno_df = anno_df[anno_df['feature'].isin(features_to_show) & anno_df['status'].isin(status_to_show) & (anno_df['biotype'].isin(bio_types_to_show))]
    anno_df['Gene_symbol'] = anno_df['ID::name'].apply(lambda x: x.split('::')[1])
    
    anno_df = anno_df[anno_df['Gene_symbol']==gene_name] 
    chr = anno_df['chr'].values[0]
    start= anno_df['start'].min()-2000
    end= anno_df['end'].max()+2000
    
    return anno_df, chr, start, end#.sort_values(by='start')

def draw_genes(ax, anno_df, start=10182000, end=10194000, 
               regionstart=10182000, regionend=10194000,
               draw_name=True, draw_type=True, 
               features_color_code = {'CDS': 'green', 'UTR':'orange', 'intergenic': 'black', 'intronic': 'blue', 'proximal_promoter': 'red',
                                      'gene':'lightgrey', 'start_codon':'yellow', 'stop_codon':'yellow', 'exon':'grey'},
               features_heights = {'gene':1, 'UTR':0.5, 'start_codon':0.5, 'stop_codon':0.5, 'CDS':1, 'exon':0.5, 'proximal_promoter':1}):
    count = -1
    for gene_id_name, gene_df in anno_df.groupby('ID::name'):
        count+=features_heights['CDS']
        
        for i, r in gene_df.iterrows():
            y_shift = 0
            ax.add_patch(patches.Rectangle((r['start'], count+y_shift), r['end']-r['start'], features_heights[r['feature']], 
                                           edgecolor = None, linewidth=1.0, fill=True, color=features_color_code[r['feature']]))#facecolor="green"
        count+=1
        strand = '> >'
        if gene_df['strand'].values[0]=='-':
            strand = '< <'
        draw_text(ax, x=gene_df['start'].min(), y = count, text=gene_id_name.split('::')[1], fontsize=10)
        ax.plot([i for i in np.arange(regionstart, regionend)], [count for i in np.arange(regionstart, regionend)], color='red', linewidth=1)
        ax.plot([regionstart for i in np.arange(0, 2)], [count+i for i in np.arange(0, 2)], color='red', linewidth=1.0)
        ax.plot([regionend for i in np.arange(0, 2)], [count+i for i in np.arange(0, 2)], color='red', linewidth=1.0)
    ax.set_xlim(start,end)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
    #ax.set_ylim(0, count+3)
    #ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.set_yticks([])
    ax.set_yticklabels([])
    if show_y_label:
        ax.set_ylabel("Overlapping\ngene", fontsize=10)
    else:
        ax.get_yaxis().set_visible(False)
    sns.despine(bottom=True,left=True,ax=ax)
    return

def get_mutations_to_plot(elements_infile, gene_to_draw='', regions_to_draw=[]):
    elements = pd.read_table(elements_infile, skiprows=6, header=0)
    elements['Genes'] = elements['Nearby-Genes(Name::ID::O|U|Ddistance::COSMIC|KCP|PCD)'].apply(lambda x: x.split(':')[0])
    if gene_to_draw!="":
        elements = elements[elements['Genes']==gene_to_draw]
    elif len(regions_to_draw)>0:
        elements = elements[elements['Position'].isin(regions_to_draw)]
    
    muts_to_plot = {}
    region_positions = []
    regions = {}
    motifs = {}
    for i, r in elements.iterrows():
        regions[r['Position'].split(':')[1]]=r['Feature_type']
        region_positions.append(int(r['Position'].split(':')[1].split('-')[0]))
        region_positions.append(int(r['Position'].split(':')[1].split('-')[1]))
        for mut_motif in r['Mutated-Moitfs'].split(','):
            mut_motif_info = mut_motif.split('#')
            motifs[mut_motif_info[8]+'-'+mut_motif_info[9]] = mut_motif_info[10]
        for mut in r['Muts'].split(','):
            mut_info = mut.split('#')
            muts_to_plot[mut_info[0].split(':')[1].split('-')[0] + "::"+mut_info[5]] = mut_info[2]
    return regions, muts_to_plot, motifs, sorted(region_positions)
    
def plot_muts(ax, regions, muts_to_plot, motifs, x_shift, start, end):
    
    regions_plotted = []
    features_color_code = {'CDS': 'green', 'UTR':'orange', 'intergenic': 'black', 'intronic': 'blue', 'proximal_promoter': 'red',
                                      'gene':'lightgrey', 'start_codon':'yellow', 'stop_codon':'yellow', 'exon':'grey'}
               
    #feature_type_colors = {'exon': 'green', 'intergenic': 'orange', 'intronic':'red'}
    for region in regions.keys():
        try:
            color = features_color_code[regions[region]]
        except KeyError:
            color = 'green'
        ax.add_patch(patches.Rectangle(
                                       (int(region.split('-')[0]), 0), int(region.split('-')[1])-int(region.split('-')[0]), 0.5, 
                                           linewidth=1.0, fill=False, edgecolor = color))
        regions_plotted.append(int(region.split('-')[0]))
        regions_plotted.append(int(region.split('-')[1]))
    
    for motif in motifs:
        ax.add_patch(patches.Rectangle(
                                       (int(motif.split('-')[0]), 0), int(motif.split('-')[1])-int(motif.split('-')[0]), 0.5, 
                                           linewidth=1.0, fill=True, 
                                           color = 'brown'))
    groups_colors_dict = {'Biliary-AdenoCA':'#00CD66','Bladder-TCC':'#EEAD0E','Bone-Osteosarc':'#FFD700','Bone-Leiomyo':'#FFEC8B','Bone-Epith':'#ADAC44','Breast-AdenoCa':'#CD6090','Cervix-SCC':'#79CDCD','CNS-Medullo':'#D8BFD8','CNS-PiloAstro':'#B0B0B0','CNS-GBM':'#3D3D3D','CNS-Oligo':'#787878','ColoRect-AdenoCA':'#191970','Eso-AdenoCa':'#1E90FF','Head-SCC':'#8B2323','Kidney-RCC':'#FF4500','Kidney-ChRCC':'#B32F0B','Liver-HCC':'#006400','Lung-SCC':'#FDF5E6','Lung-AdenoCA':'#FFFFFF','Lymph-BNHL':'#698B22','Lymph-CLL':'#F4A35D','Myeloid-MPN':'#FFC100','Myeloid-AML':'#CD6600','Ovary-AdenoCA':'#008B8B','Panc-AdenoCA':'#7A378B','Panc-Endocrine':'#E066FF','Prost-AdenoCA':'#87CEFA','Skin-Melanoma':'#000000','Stomach-AdenoCA':'#BFEFFF','Thy-AdenoCA':'#9370DB','Uterus-AdenoCA':'#FF8C69','Bone-Cart':'#DDCDCD','Breast-LobularCa':'#DDCDCD','Breast-DCIS':'#DDCDCD','Lymph-NOS':'#DDCDCD','Myeloid-MDS':'#DDCDCD','Cervix-AdenoCA':'#DDCDCD'}
    muts_plotted = {}
    cancer_types_showed = []
    for i, mut in enumerate(muts_to_plot.keys()):
        x_pos = int(mut.split('::')[0])
        try:
            muts_plotted[x_pos]+=1
        except KeyError:
            muts_plotted[x_pos]=1
        draw_marker(ax=ax,x=x_pos, y=muts_plotted[x_pos], marker='.', markersize=6,
                    color=groups_colors_dict[muts_to_plot[mut]])
        cancer_types_showed.append(muts_to_plot[mut])
    x_move = 10
    y_move = 6
    for i, cancer_type in enumerate(list(set(cancer_types_showed))):
        if 'Lymph' in cancer_type or 'Kidney' in cancer_type:
            draw_marker(ax=ax,x=start+x_move, y=y_move, marker='.', markersize=6, color=groups_colors_dict[cancer_type])
            draw_text(ax=ax, x=start+x_move+5, y=y_move, text=cancer_type, horizontalalignment='left')
        #x_move+=40
            y_move-=1
    
    '''for mut in muts_plotted.keys():
        stem_line(ax, x=[mut], y=[muts_plotted[mut]], marker='o', markerfacecolor='red', 
                  markeredgecolor=None, stemline_color='red', markersize=1.0, markeredgewidth=0.2)
    '''
    #get the max mutation count 
    ymax = sorted(muts_plotted.viewitems(), key=lambda x: x[1], reverse=True)[0][1]
    ax.set_ylim(-0.8, 6.2)#ymax+0.2)
    ax.plot([i for i in np.arange(start,end)], [-0.3 for i in np.arange(start,end)], color='grey', linewidth=1.0)
    draw_text(ax, x=start+((end-start)/2), y=-0.8, text=str(end-start)+' bp', fontsize=8)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    if show_y_label:
        ax.set_ylabel("Number of mutations", fontsize=10)
        sns.despine(ax=ax, bottom=True)
    else:
        ax.get_yaxis().set_visible(False)
        sns.despine(ax=ax, bottom=True, left=True)
    
    ax.set_xlim(start-x_shift,end+x_shift)#sorted(regions_plotted)[0], sorted(regions_plotted)[-1]+5)#start,end+1)#anno_df['start'].min(), anno_df['end'].max())
    ax.get_xaxis().set_visible(False)
    
    
    return

def get_boxes_to_plot(infile, chr, start, end, cell_names=[], factors=[], names=['chr', 'start', 'end', 'cell', 'factor', 'score']):
    
    df = pd.read_table(infile, header=None, names=names)
    if len(factors)>0:
        df = df[(df['chr'] == chr) & (df['start']>=start) & (df['end']<=end) (df['cell'].isin(cell_names)) & (df['factor'].isin(factors))]
    else:
        df = df[(df['chr'] == chr) & (df['cell'].isin(cell_names))  & 
                                      ( 
                                       ((df['start']>=start) & (df['end']<=end)) | 
                                       ((df['start']<=start) & (df['end']>start)) |
                                       ((df['start']<=end) & (df['end']>=end)) 
                                       )]
    df.sort_values(by=['cell', 'factor'])
    df.to_csv(infile+".csvtemp", sep='\t', index=False, header=None)
    cells_boxes = {}
    if len(df)==0:
        print('no peaks found overlapping this region:', chr, start, end)
        return cells_boxes
    merged = BedTool(infile+".csvtemp").sort().merge(c=[2,3,4,5,6], o=['collapse', 'collapse', 'collapse', 'collapse','collapse'])#.groupby(g=[4,5], c=[1,2,3,6], o=['distinct','min','max','distinct'], full=True)
    for r in merged:
        for i,cell in enumerate(r[5].split(',')):
            try:
                cells_boxes[cell].append((int(r[3].split(',')[i]), int(r[4].split(',')[i]), r[6].split(',')[i], i))
            except KeyError:
                cells_boxes[cell] = [(int(r[3].split(',')[i]), int(r[4].split(',')[i]), r[6].split(',')[i], i)]
    return cells_boxes
    
def plot_peaks(ax, start, end, cells_boxes):
    #ax.plot([i for i in range(0,100)], [1 for i in range(0,100)])
    max_num_peaks = 0
    for cell in cells_boxes.keys():
        for i, peak in enumerate(cells_boxes[cell]):
            if peak[3]>max_num_peaks:
                max_num_peaks = peak[3]
            ax.add_patch(patches.Rectangle(
                                       (peak[0], peak[3]*2), peak[1]-peak[0], 0.5, 
                                           linewidth=1.0, fill=True, color = 'grey'))
    #draw_text(ax, x=start+20, y=max_num_peaks*2-5, text="TF Peaks", fontsize=10)
    ax.set_ylim(0,(max_num_peaks*2)+1)
    if show_y_label:
        ax.set_ylabel("TF Peaks", fontsize=10)
        ax.set_yticks([])
        ax.set_yticklabels([])
    else:
        ax.get_yaxis().set_visible(False)

    ax.set_xlim(start,end)
    print(start, end, end-start)
    #ax.plot(x = [i for i in np.arange(start, end)], y=[0 for i in np.arange(start, end)])
    #ax.set_xticklabels([start/1000000.0, end/1000000.0])
    ax.xaxis.set_major_locator(ticker.FixedLocator(locs=[start, start+((end-start)/2.0), end], nbins=None))
    formatter = matplotlib.ticker.ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    #ax.set_xticklabels([float(str(x.get_text()))/1000000.0 for x in ax.get_xticklabels()])
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    #ax.get_xaxis().set_visible(False)
    sns.despine(ax=ax, left=True)#bottom=True, 
    
    return

def plot_genomic_context(ax1, ax2, ax3, x_shift, gene_name, cell_names, elements_infile, gene_infile, chip_seq_infile):
    regions, muts_to_plot, motifs, region_positions = get_mutations_to_plot(elements_infile, gene_to_draw=gene_name, regions_to_draw=[])
    plot_muts(ax1, regions, muts_to_plot, motifs, x_shift=x_shift, start=min(region_positions), end=max(region_positions))
    
    gene_df, chr, start, end = get_gene_features(infile=gene_infile, gene_name=gene_name, 
                                features_to_show=['exon'])
    draw_genes(ax2, gene_df, start=start, end=end, regionstart=min(region_positions), regionend=max(region_positions))
    
    cells_boxes = get_boxes_to_plot(chip_seq_infile, chr=chr, start=start, end=end, cell_names=cell_names, factors=[])
    plot_peaks(ax3, start=start, end=end, cells_boxes=cells_boxes)
    
def plot_gene_expr(fig, gs, row_num, genes_cancertypes,genes_mutated_input, meta_data,gene_expr_intput):

    #genes_mutated_input = '../analysis/PancanElements/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv_GenesInclExons.tsv'
    #genes_mutated_input = '../analysis/data/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb_pancan.tsv_GenesInclCDS.tsv'
    #meta_data = '../analysis/RNA-seq/extended_meatadata_syn7416381'
    #gene_expr_intput = '../analysis/RNA-seq/tophat_star_fpkm_uq.v2_aliquot_gl.tsv'
    meta_data = get_sample_data(meta_data)
    mutated_genes = read_genes_elements(genes_mutated_input)
    print('mutated_genes')
    gene_counts, gene_counts_file =  read_gene_expr(gene_expr_intput, mutated_genes['GeneID'])
    print('mutated genes extracted')
    gene_counts_info, gene_counts_info_file = get_expr_per_sample(mutated_genes, meta_data, gene_counts, gene_counts_file, sample_col_to_use='Samples')
    print('expr per sample extracted into:', gene_counts_info_file)
    gene_counts_info_stats, gene_counts_info_stats_file = process_gene_counts(gene_counts_info, mutated_genes, gene_counts_info_file)
    print('stats done')
    #make a scatter plot for genes that are mutated in at least 10 samples with expr data (pval and avg FC (WT)
    #df = get_sig_expr_events(gene_counts_info_stats, gene_counts_info_stats_file)
    #plot_scatter_geneexpr(df)
    box_plot_per_gene_cancertype( fig, gs, row_num, gene_counts_info_stats, genes_cancertypes)#, 'TERT':['Skin-Melanoma', 'Bladder-TCC','CNS-Oligo','Thy-AdenoCA']})
        
def parse_args():
    '''Parse command line arguments'''
    
    parser = argparse.ArgumentParser(description='Plot Fig 4')
    parser.add_argument('-e', '--elements_input_file', default='', help='')
    parser.add_argument('--gene_input_file', default='', help='')
    parser.add_argument('--chip_seq_input_file', default='', help='')
    parser.add_argument('--genes_elem_input_file', default='',  help='')
    parser.add_argument('--meta_data', default='',  help='')
    parser.add_argument('--gene_expr_intput', default='',  help='')
    parser.add_argument('--output_dir', default='', help='')
    

    
    return parser.parse_args(sys.argv[1:])


if __name__ == '__main__':
    
    args = parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    plt.clf()
    fig = plt.figure(figsize=(8, 6), linewidth=1.0)#design a figure with the given size
    genes_cancertypes = ['VHL:Kidney-RCC', 'BCL2:Lymph-BNHL', 'MYC:Lymph-BNHL', 'RP11-731F5.1:Lymph-BNHL']
    num_cols = len(genes_cancertypes)
    gs = gridspec.GridSpec(4, 8, height_ratios=[4,2,2,4], wspace=0.0, hspace=0.0)#create 4 rows and three columns with the given ratio for each
    sns.set_style('white', {'axes.linewidth': 1})
    #first genomic track
    ax1 = fig.add_subplot(gs[0,0:3])
    ax2 = fig.add_subplot(gs[1,0:3])
    ax3 = fig.add_subplot(gs[2,0:3])
    #second genomic track
    ax4 = fig.add_subplot(gs[0,3:8])
    ax5 = fig.add_subplot(gs[1,3:8])
    ax6 = fig.add_subplot(gs[2,3:8])
    
    gs.tight_layout(fig, pad=2, h_pad=2.0, w_pad=2.0)
    x_shift=100
    global show_y_label
    show_y_label = True
    
    #elements_infile = '/Users/husensofteng/Documents/workspace/ActiveMotifs/analysis/data/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb_pancan.tsv'
    #gene_infile = '/Users/husensofteng/Documents/workspace/ActiveMotifs/analysis/data/tracks/gencode.v19.annotation.gff3_extractedinfo'
    #chip_seq_infile='/Users/husensofteng/Documents/workspace/ActiveMotifs/analysis/data/tracks/all_chip-seq_data_CellInfo_combined.bed6'
    
    
    elements_infile = arsg.elements_input_file
    gene_infile = args.gene_input_file
    chip_seq_infile=args.chip_seq_input_file
    genes_mutated_input= args.genes_elem_input-file
    meta_data = args.meta_data
    gene_expr_intput=args.gene_expr_intput
    plot_genomic_context(ax1, ax2, ax3, x_shift=x_shift, gene_name='VHL', cell_names=['HEK293'], elements_infile=elements_infile, gene_infile=gene_infile, chip_seq_infile=chip_seq_infile)
    show_y_label = False
    plot_genomic_context(ax4, ax5, ax6, x_shift=x_shift, gene_name='MYC', cell_names=['GM12878'], elements_infile=elements_infile, gene_infile=gene_infile, chip_seq_infile=chip_seq_infile)
    
    plot_gene_expr(fig, gs, row_num=3, genes_cancertypes=genes_cancertypes, genes_mutated_input, meta_data,gene_expr_intput)
    fig4 = args.output_dir+'/Fig4'
    plt.savefig(fig4+'.pdf')
    plt.savefig(fig4+'.svg')
    