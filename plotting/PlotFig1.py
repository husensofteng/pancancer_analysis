'''
Created on May 29, 2017

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

from utils import *

def draw_steam_lines(ax):
    #no motif no signal
    draw_text(ax, x=4, y=4.75, text="Identification of CFRMs", fontsize=12)
    y_text_checks = 4
    stem_line(ax, x=[1], y=[2], marker='.', markerfacecolor='grey', markeredgecolor='grey', stemline_color='grey')
    draw_marker(ax, x=1,y=y_text_checks, marker='$\\times$', color='black')
    draw_text(ax, x=1,y=-1.75, text="No peak\nNo motif")
    #ax.plot([i for i in range(0,3)], [0 for i in range(0,3)], color='grey')
    #effective motif but no signal
    stem_line(ax, x=[4.5], y=[2], marker='.', markerfacecolor='grey', markeredgecolor='red', stemline_color='grey')
    ax.add_patch(patches.FancyBboxPatch((4, 0.1), 1, 0.1, edgecolor = 'orange', boxstyle='round', fill=False, linewidth=2.0))#facecolor="green"
    draw_marker(ax, x=4.5,y=y_text_checks, marker='$\\times$', color='black')
    draw_text(ax, x=4.5,y=-1.25, text="No peak")
    #ax.plot([i for i in np.arange(3.5,6.5)], [0 for i in np.arange(3.5,6.5)], color='grey')
    #signal but no motif
    stem_line(ax, x=[10], y=[2], marker='.', markerfacecolor='grey', markeredgecolor='grey', stemline_color='grey')
    plot_sin(ax, shift_x=8)
    draw_marker(ax, x=10,y=y_text_checks, marker='$\\times$', color='black')
    draw_text(ax, x=10,y=-1.25, text="No motif")
    #ax.plot([i for i in range(8,13)], [0 for i in range(8,13)], color='grey')
    #signal and motif but not effective
    stem_line(ax, x=[15], y=[2], marker='.', markerfacecolor='grey', markeredgecolor='grey', stemline_color='grey')
    ax.add_patch(patches.FancyBboxPatch((14.5, 0.1), 1, 0.1, edgecolor = 'orange', boxstyle='round', fill=False, linewidth=2.0))#facecolor="green"
    plot_sin(ax, shift_x=13)
    draw_marker(ax, x=15,y=y_text_checks, marker='$\\times$', color='black')
    draw_text(ax, x=15,y=-1.75, text="No sig. effect\non motif")
    #ax.plot([i for i in range(13,18)], [0 for i in range(13,18)], color='grey')
    #signal and motif but not significant
    stem_line(ax, x=[20], y=[2], marker='.', markerfacecolor='grey', markeredgecolor='red', stemline_color='grey')
    ax.add_patch(patches.FancyBboxPatch((19.5, 0.1), 1, 0.1, edgecolor = 'orange', boxstyle='round', fill=False, linewidth=2.0))#facecolor="green"
    plot_sin(ax, shift_x=18)
    draw_marker(ax, x=20,y=y_text_checks, marker='$\\times$', color='black')
    draw_text(ax, x=20,y=-1.75, text="No sig. MFS")
    #ax.plot([i for i in range(18,23)], [0 for i in range(18,23)], color='grey')
    #effective motif, signal and significant
    stem_line(ax, x=[25.5], y=[2], marker='o', markerfacecolor='green', markeredgecolor='red', stemline_color='green')
    ax.add_patch(patches.FancyBboxPatch((25, 0.1), 1, 0.1, edgecolor = 'orange', boxstyle='round', fill=False, linewidth=2.0))#facecolor="green"
    plot_sin(ax, shift_x=23.5)
    draw_marker(ax, x=25.5,y=y_text_checks, marker='$\\checkmark$')
    draw_marker(ax, x=26.4,y=3, marker='*', color='black', markersize=7)
    
    draw_text(ax, x=26.4,y=-2.25, text="Sig. MFS\nSig. effect on motif   \nDNase1 peak           ")
    draw_marker(ax, x=23.15,y=-1, marker='*', color='black', markersize=7)
    draw_marker(ax, x=23.15,y=-1.55, marker='o', color='green', markersize=5, markeredgecolor='red', markeredgewidth=1)
    plot_sin(ax, shift_x=22.9, y_shift=-2.45, R=0.5, A=-0.5, B=-0.5, color='blue')
    #ax.plot([i for i in np.arange(23.5,28.5)], [0 for i in np.arange(23.5,28.5)], color='grey')
    #effective motif and same tf peak
    stem_line(ax, x=[32.5], y=[2], marker='.', markerfacecolor='grey', markeredgecolor='red', stemline_color='grey')
    ax.add_patch(patches.FancyBboxPatch((32, 0.1), 1, 0.1, edgecolor = 'orange', boxstyle='round', fill=False, linewidth=2.0))#facecolor="green"
    plot_sin(ax, shift_x=30.5, color='grey', linestyle='-')
    draw_marker(ax, x=32.5, y=y_text_checks, marker='$\\times$', color='black')
    draw_text(ax, x=32.5,y=-1.75, text="No motif-matcing\nTF-peak")
    #ax.plot([i for i in np.arange(30.5,35.5)], [0 for i in np.arange(30.5,35.5)], color='grey')
    
    #effective motif and same tf peak
    stem_line(ax, x=[38], y=[2], marker='o', markerfacecolor='green', markeredgecolor='red', stemline_color='green')
    ax.add_patch(patches.FancyBboxPatch((37.5, 0.1), 1, 0.1, edgecolor = 'orange', boxstyle='round', fill=False, linewidth=2.0))#facecolor="green"
    plot_sin(ax, shift_x=36, color='orange', linestyle='-')
    #plot_sin(ax3, shift_x=36, R=1, A=3,B=3,color='orange', linestyle='-')
    
    draw_marker(ax, x=38,y=y_text_checks, marker='$\\checkmark$')
    draw_text(ax, x=38.65,y=-2.25, text="Sig. effect on motif\nMotif-matching      \nTF-peak   ")
    draw_marker(ax, x=35.65,y=-1, marker='o', color='green', markersize=5, markeredgecolor='red', markeredgewidth=1)
    plot_sin(ax, shift_x=35.4, y_shift=-1.95, R=0.5, A=-0.5, B=-0.5, linestyle='-', color='orange')
    
    ##ax.plot([i for i in range(36,41)], [0 for i in range(36,41)], color='grey')
    ax.plot([i for i in np.arange(0,42)], [0 for i in np.arange(0,42)], color='grey', linewidth=1.0)
    
    ax.plot([i for i in np.arange(-1,43)], [-2.75 for i in np.arange(-1,43)], color='grey', linewidth=1.0)
    ax.plot([i for i in np.arange(9.5,43)], [5 for i in np.arange(9.5,43)], color='grey', linewidth=1.0)
    
    ax.plot([-1 for i in np.arange(-2.75,6)], [i for i in np.arange(-2.75,6)], color='grey', linewidth=1.0)
    ax.plot([42 for i in np.arange(-2.75,6)], [i for i in np.arange(-2.75,6)], color='grey', linewidth=1.0)
    
    ax.set_xlim(-1,42)
    ax.set_ylim(-4,5)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.set_frame_on(False)
    #ax.xaxis.set_ticks([])
        
#mutate rates    

#part B

def rate_per_chromatin_state(muts_input_file, ax, names='chr,start,Annotation'.split(','), usecols = [0,1,9],
              d = 300, chr_in=0, window_size=1000000):
    
    df = pd.read_table(muts_input_file, sep='\t', skiprows=0, header=None, usecols=usecols, names=names, nrows=1000)
    #plt.figure(figsize=(12, 12))
    df['chr'] = df['chr'].apply(getchrnum)
    if chr_in>0:
        df = df[df['chr']==chr_in]
    df['x'] = df['chr'] + df['start'].apply(float).apply(get_xaxis, args=(window_size, d,))
    df['State'] = df['Annotation'].apply(get_state)
    df['NumberMutsPerWindowPerState'] = df.groupby(by=['State','x'])['start'].transform('count')
    avg_per_state_per_mb = {}
    
    for i, r in df.iterrows():
        try:
            avg_per_state_per_mb[r['State']][r['x']] =  r['NumberMutsPerWindowPerState']
        except KeyError:
            avg_per_state_per_mb[r['State']] = {r['x']:  r['NumberMutsPerWindowPerState']}
    #plt.show()
    states = []
    counts = []
    for k in sorted(avg_per_state_per_mb.keys()):
        for i in range(len(avg_per_state_per_mb[k]), len(df['x'].unique())):
                avg_per_state_per_mb[k][str(i)+"nomut"]=0
        for w in sorted(avg_per_state_per_mb[k]):
            states.append(k)
            #if avg_per_state_per_mb[k][w]>0.0:
            #    counts.append(math.log10(avg_per_state_per_mb[k][w]/2515.0))
            #else:
            counts.append(avg_per_state_per_mb[k][w]/2515.0)
    
    #sns.set_style("white", {'axes.linewidth': 0.5})
    #sns.set_context("talk")
    #plt.figure(figsize=(4, 2))
    states_order = 'Tx,Tss,Enh,Repr,NA,Quies'.split(',')
    
    #hue_colors = {'TssA': , 'TssAFlnk': ,'TssBiv': , 'Tx': ,'TxFlnk': , 'TxWk': ,'BivFlnk': , 'Enh': ,'EnhBiv': ,'EnhG': , 'Het': ,'ZNF/Rpts': ,'ReprPCWk': ,'ReprPC': ,'Quies': ,'NoState':}
    plt_results = sns.barplot(states, counts, estimator=average, ci=95, ax=ax, 
                              order=states_order, hue_order=states_order, palette=['red', 'green', 'orange', 'brown', 'grey','lightgray'],
                              errwidth=0.5, capsize=0.15, orient='v')#'BuGn_r')
    plt_results.set_xticklabels(plt_results.get_xticklabels(), rotation=0)
    #plt.setp(ax.spines.values(), linewidth=1)
    ax.spines['left'].set_linewidth(1)#set line width fo the y axis
    ax.spines['left'].set_color('grey')
    #plt_results.set_yticklabels(range(0,10))
    #draw_text(plt_results, x = 0, y=3, text='c')
    ax.set_ylabel("Average mutation\nrate per megabse")
    tick_spacing=1
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    
    #ax.plot([-0.5 for i in np.arange(0,5)], [i for i in range(0,5)], linewidth=1.0, color='grey')
    #draw_text(ax, x=0, y=3, text='c', color='black', fontsize=10)
    
def plot_muts(muts_input_file, names='chr,start,DNase1,TFBinding'.split(','), usecols = [0,1,10,11],
              motifs=True, d = 6000.0, chr_in=0, window_size=50000, ymax_value = 1.0):
    
    df = pd.read_table(muts_input_file, sep='\t', skiprows=0, header=None, usecols=usecols, names=names, nrows=1000)
    #plt.figure(figsize=(12, 12))
    df['chr'] = df['chr'].apply(getchrnum)
    if chr_in>0:
        df = df[df['chr']==chr_in]
    df['x'] = df['chr'] + df['start'].apply(float).apply(get_xaxis, args=(window_size, d,))
    if motifs:
        df['Active'] = np.where((df['DNase1']>1e-300) | (df['TFBinding']>1e-300), True, False)#col='ChromatinState', col_wrap=5, 
    else:
        df['Active'] = df['Annotation'].str.contains('DNase1|TFBinding', na=False)
        #df['Active'] = np.where(('TFBinding' in df['Annotation']) | ('DNase1' in df['Annotation']), 'Yes', 'No')
    #df['x'] = df.start.apply(get_xaxis, args=(int(df['chr'].replace('X','23').replace('Y','24').replace('M','25').replace('chr','')),))
    #df['y'] = df.start.apply(get_yaxis)
    #df['ChromatinState'] = df['ChromatinState'].apply(get_unique_state)
    
    active_windows = []
    active_window_counts = []
    unactive_windows = []
    unactive_window_counts = []
    chromosomes = []
    max_window_numbers_per_chr = []
    hyper_mutated_windows = {}
    for active_label, df_activity in df.groupby('Active'):
        for chr_label, df_chr in df_activity.groupby('chr'):
            chromosomes.append(chr_label)
            max_window_numbers_per_chr.append(df_chr['x'].max())
            windows = []
            window_counts = []
            for window_label, df_window in df_chr.groupby('x'):
                num_muts_in_window = len(df_window)/2515.0
                if not motifs:
                    if num_muts_in_window>ymax_value:#2515
                        hyper_mutated_windows[window_label]=num_muts_in_window
                        num_muts_in_window = ymax_value#2500
                    
                if num_muts_in_window==0:
                    windows.append('nan')
                    window_counts.append('nan')
                windows.append(window_label)
                window_counts.append(num_muts_in_window)
                #window_counts.append(math.log(len(df_window), 2))
                
                if num_muts_in_window==0:
                    windows.append('nan')
                    window_counts.append('nan')
                    
            if active_label:
                active_windows.extend(windows)
                active_window_counts.extend(window_counts)
                active_windows.append('nan')
                active_window_counts.append('nan')
            else:
                unactive_windows.extend(windows)
                unactive_window_counts.extend(window_counts)
                unactive_windows.append('nan')
                unactive_window_counts.append('nan')
    return active_windows, active_window_counts, unactive_windows, unactive_window_counts, chromosomes, max_window_numbers_per_chr, hyper_mutated_windows

def plot_lines(windows, chromosomes, ax, max_window_numbers_per_chr, hyper_mutated_windows, ymax_value = 1.0):
    ax.plot(windows[0], windows[1], windows[2], 
             windows[3], windows[4], windows[5],
             windows[6], windows[7], windows[8],
             windows[9], windows[10], windows[11],)
    ax.spines['left'].set_color('grey')
    ax.spines['left'].set_linewidth(1)
    #ax.spines.values(), linewidth=1, color='grey')
    #ax.set_xticklabels(labels=chromosomes)
    chromosomes =  np.array(list(set(chromosomes)))
    chromosomes_pos = ((np.array(max_window_numbers_per_chr[0:len(chromosomes)])-chromosomes)/2.0)+chromosomes
    chromosomes = list(chromosomes)
    ax.set_xlabel("chromosomes (hg19)")
    ax.set_ylabel("Average mutation rate per 50Kb")
    
    #plot hyper mutated elements
    names_for_hyper_elements = {2:'IGK', 14:'IGH', 22:'IGL'}
    print hyper_mutated_windows
    for hyper_element in hyper_mutated_windows.keys():
        draw_text(ax, x=hyper_element, y=1+0.05, text='.', color='black', fontsize=8, horizontalalignment='left', rotation=0)
        draw_text(ax, x=hyper_element, y=1+0.10, text='.', color='black', fontsize=8, horizontalalignment='left', rotation=0)
        draw_text(ax, x=hyper_element, y=1+0.15, text='.', color='black', fontsize=8, horizontalalignment='left', rotation=0)
        draw_text(ax, x=hyper_element, y=ymax_value+0.3, text=names_for_hyper_elements[int(hyper_element)] + ' (%0.2f)' % hyper_mutated_windows[hyper_element], color='black', fontsize=8, horizontalalignment='left', rotation=45)
    
    for i,chr in enumerate(chromosomes):
        if chr==23:
            chromosomes[i] = 'X'
        elif chr==24:
            chromosomes[i] = 'Y'
    
    ax.set_xticks(list(chromosomes_pos))#,list(chromosomes))
    ax.set_xticklabels(chromosomes)
    ax.minorticks_off()
    ax.set_ylim(0.0,ymax_value)
    #ax.tick_params(axis='x', direction='out', length=1, width=1, colors='black', )
    #tick_spacing=1
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
def draw_piechart_muts(ax, radius=1.0, center=(0,0)):
    
    labels = ['No overlap', 'TFBS/DHS', 'TFBS/DHS - Matching cell line']
    sizes = [25277342-21698761, 21698761-937338, 937338]
    colors = ['lightgrey',  '#C0C0C0', 'grey']
    explode = (0.0, 0.0, 0.1)
    pie_wedge_collection = ax.pie(sizes, explode = explode, labels=labels, colors=colors, startangle=90, labeldistance=1.05, autopct='%1.1f%%',
                                  radius=radius, center=center)
    for pie_wedge in pie_wedge_collection[0]:
        pie_wedge.set_edgecolor('white')
    ax.axis('equal')
    return

def draw_piechart_motifs(ax, radius=1.0, center=(0,0)):
    
    labels = ['No matching TFBS/DHS', 'DHS - Matching cell line', 'Matching TFBS - Matching cell line']
    sizes = [3935511-222800-72317, 222800, 72317]
    colors = ['lightgrey',  '#C0C0C0', 'grey']
    explode = (0.0, 0.1, 0.1)
    pie_wedge_collection = ax.pie(sizes, explode = explode, labels=labels, colors=colors, startangle=90, labeldistance=1.05, autopct='%1.1f%%',
                                  radius=radius, center=center)
    for pie_wedge in pie_wedge_collection[0]:
        pie_wedge.set_edgecolor('white')
    ax.axis('equal')
    return

def draw_plot():
    
    sns.set_style('white', {'text.color': '.15'})
    #matplotlib.rc('axes',edgecolor='grey')
    #mpl.rcParams['font.family'] = fontfamily
    #rcParams['font.sans-serif'] = ['Verdana']
    #rcParams['svg.fonttype'] = 'none'
    
    fig = plt.figure(figsize=(12,8))#design a figure with the given size
    gs = gridspec.GridSpec(3, 3, height_ratios=[3 , 6, 3.5], wspace=0.0, hspace=0.0)#create 2 rows and three columns with the given ratio for each
    ax0 = fig.add_subplot(gs[1, :])
    ax1 = fig.add_subplot(gs[2, :])#take the entire first row for the first sub plot
    ax2 = fig.add_subplot(gs[0, -1:]) #take the last col of the second row for sub fig3
    #ax2 = fig.add_subplot(gs[0, :-1]) #the col1 and 2 of the second row for sub fig2
    
    ax3 = fig.add_subplot(gs[0, 0]) #the col1 and 2 of the second row for sub fig2
    ax4 = fig.add_subplot(gs[0, 1])
    #ax4 = fig.add_subplot(gs[0, 1]) #the col1 and 2 of the second row for sub fig2
    #ax7 = fig.add_subplot(gs[1, 0])
    #ax8 = fig.add_subplot(gs[1, 1])
    draw_piechart_muts(ax3)
    draw_piechart_motifs(ax4)
    #ax2.get_xaxis().set_visible(False)
    #ax2.get_yaxis().set_visible(False)
    #f, (ax1, ax2) = plt.subplots(2, figsize=(12,6))
    
    #Plot all muts for all chromosomes
    muts_input_file = '../analysis/data/observed_agreement_22May2017_annotated.bed10'
    motif_muts_input_file = '../analysis/data/motifmuts_all.bed12'
    active_windows, active_window_counts, unactive_windows, unactive_window_counts, chromosomes, max_window_numbers_per_chr, hyper_mutated_windows = plot_muts(muts_input_file, names='chr,start,Annotation'.split(','),
              usecols = [0,1,9], motifs=False, chr_in=0, ymax_value = 1.0)
    print len(active_windows)
    
    windows = [unactive_windows, unactive_window_counts, '#666666', active_windows, active_window_counts, '#0000ff']
    active_windows, active_window_counts, unactive_windows, unactive_window_counts, chromosomes, max_window_numbers_per_chr, hyper_mutated_windows_motifs = plot_muts(motif_muts_input_file, ymax_value = 1.0)
    print len(active_windows)
    windows.extend([unactive_windows, unactive_window_counts, '#4dff4d', active_windows, active_window_counts, '#ff3333'])
    plot_lines(windows, chromosomes, ax0, max_window_numbers_per_chr, hyper_mutated_windows)
    
    #Plot CFRMs
    draw_steam_lines(ax1)
    #Plot SF-MREs
    #draw_sigregs(ax3)
    
    #Plot Barchart for chromatin states
    rate_per_chromatin_state(muts_input_file, ax = ax2, names='chr,start,Annotation'.split(','), usecols = [0,1,9],
                             d = 300, chr_in=0, window_size=1000000)
    
    #f.subplots_adjust(hspace=0)
    
    sns.despine(right=True, top=True, bottom=True, left=False)
    plt.savefig("../analysis/Fig1max1.svg", bbox_inches='tight')
    plt.savefig("../analysis/Fig1max1.pdf", bbox_inches='tight')
    plt.close()
    
if __name__ == '__main__':
    draw_plot()
    