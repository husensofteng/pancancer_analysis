'''
Created on 21 Jul 2017

@author: husensofteng
'''
import matplotlib
from numpy.lib.function_base import average
import math
matplotlib.use('TkAgg')
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
from operator import itemgetter
import sys
import os
import argparse

from utils import *

def draw_sigregs(ax):
    x1 = 1
    
    draw_text(ax, x=x1+2.5, y=4.25, text="Identification of Regulatory Elements", fontsize=12)
    ax.plot([i for i in np.arange(0.5,x1+6.5)], [0 for i in np.arange(0.5,x1+6.5)], color='grey', linewidth=1)
    
    stem_line(ax, x=[x1,x1+1,x1+4], y=[2,2,2], marker='.', markerfacecolor='grey', markeredgecolor='grey', stemline_color='grey')
    stem_line(ax, x=[x1+2, x1+5], y=[3,2], marker='.', markerfacecolor='grey', markeredgecolor='red', stemline_color='grey')
    #ax.plot([i for i in np.arange(x1-1,x1+10)], [0 for i in np.arange(x1-1,x1+10)], color='grey')
    #draw_text(ax, x=x1+4, y=-0.75, text="Merge within 200bp")
    
    x1 = 10
    stem_line(ax, x=[x1,x1+1,x1+4], y=[2,2,2], marker='.', markerfacecolor='grey', markeredgecolor='grey', stemline_color='grey')
    stem_line(ax, x=[x1+2, x1+5], y=[3,2], marker='o', markerfacecolor='green', markeredgecolor='red', stemline_color='green')
    
    draw_marker(ax, x=x1+5, y=3.5, marker='*', color='black', markersize=7)
    draw_marker(ax, x=x1+5.5, y=3.5, marker='*', color='black', markersize=7)
    draw_marker(ax, x=x1+6, y=3.5, marker='*', color='black', markersize=7)
    
    #ax.plot([i for i in np.arange(x1-1,x1+9)], [0 for i in np.arange(x1-1,x1+9)], color='grey', linewidth=1)
    ax.plot([i for i in np.arange(0.5,x1+6.5)], [0 for i in np.arange(0.5,x1+6.5)], color='grey', linewidth=1)
    ax.plot([i for i in np.arange(x1,x1+6)], [-0.75 for i in np.arange(x1,x1+6)], color='orange', linewidth=2)
    draw_text(ax, x=x1+2.5, y=-1.75, text="Regulatory Element")
    
    ax.plot([i for i in np.arange(0,x1+8)], [-2.5 for i in np.arange(0,x1+8)], color='grey', linewidth=1)#bottom line
    ax.plot([0 for i in np.arange(-2.5,5)], [i for i in np.arange(-2.5,5)], color='grey', linewidth=1)#left
    ax.plot([17 for i in np.arange(-2.5,5)], [i for i in np.arange(-2.5,5)], color='grey', linewidth=1)#right
    ax.plot([i for i in np.arange(7,x1+8)], [4.5 for i in np.arange(7,x1+8)], color='grey', linewidth=1)#top line
    
    
    ax.set_xlim(0,20)
    #ax.set_ylim(-0.75,4)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.set_frame_on(False)

def draw_rec_sigregs(ax, elements_input_file):
    elements_input = pd.read_table(elements_input_file, sep='\t', skiprows=6, header=0)
    df = elements_input[(elements_input['#Samples(RegMuts)']>1)]
    #print(df)
    df['FDR'] = np.where(df['FDR']==0.0, 1e-300, df['FDR'])#to replace zero with a number that can be converted to log10
    df['FDR'] = df['FDR'].apply(lambda x: np.log10(x)*-1)
    y = 'Score'
    x = '#Samples'
    s = '#Samples(RegMuts)'
    df[x]=df[x].apply(np.log2)
    df[y]=df[y].apply(np.log2)
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
        y_shift = 0.1
        x_shift = -0.05
        if label == 'None':
            label = ''
        draw_text(ax, x=r[x]+x_shift, y=r[y]+y_shift, text=label, horizontalalignment='left', verticalalignment='bottom', color=color, fontsize=8, rotation=rotation)
    ax.set_xlabel('Number of mutated samples (log2)')
    ax.set_ylabel('Regulatory score (log2)')
    ax.set_xlim(0,df[x].max()+0.2)
    ax.set_ylim(0,df[y].max()+0.2)
    #ax.legend(loc='upper center')
    return ax

def draw_plot(elements_input_file, output_dir):
    sns.set_style('white', {'text.color': '.15', 'axes.linewidth': 0.5})
    #matplotlib.rc('axes',edgecolor='grey')
    #mpl.rcParams['font.family'] = fontfamily
    #rcParams['font.sans-serif'] = ['Verdana']
    #rcParams['svg.fonttype'] = 'none'
    
    fig = plt.figure(figsize=(8,6), linewidth=1.0)#design a figure with the given size
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,6], wspace=2.0, hspace=0.0)#create 2 rows and three columns with the given ratio for each
    ax0 = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, :])#take the entire first row for the first sub plot
    draw_sigregs(ax0)
    #elements_input_file = '../analysis/data/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb_ATELM.tsv'
    draw_rec_sigregs(ax1, elements_input_file)
    #gs.tight_layout(fig, pad=2.0, h_pad=2.0, w_pad=4.0)
    sns.despine(right=True, top=True, bottom=False, left=False)
    fig3=output_dir+'/Fig3max1'
    plt.savefig(fig3+".svg", bbox_inches='tight')
    plt.savefig(fig3+".pdf", bbox_inches='tight')
    plt.close()

def parse_args():
    '''Parse command line arguments'''
    
    parser = argparse.ArgumentParser(description='Plot Fig3')
    parser.add_argument('-e', '--elements_input_file', default='', help='')
    parser.add_argument('--output_dir', default='', help='')
    

    
    return parser.parse_args(sys.argv[1:])


if __name__ == '__main__':
    
    args = parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
       
    draw_plot(args.elements_input_file, args.output_dir)