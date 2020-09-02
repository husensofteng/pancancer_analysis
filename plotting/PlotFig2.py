'''
Created on May 29, 2017

@author: husensofteng
'''
import matplotlib
#matplotlib.use('TkAgg')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
#plt.style.use('seaborn-ticks')
sns.set_style('white', {'text.color': '.15'})
from matplotlib.pyplot import tight_layout, ylabel
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, ticker
import pandas as pd
from operator import itemgetter
import os,sys
import argparse

from utils import *

def draw_motifalignment(ax, fig):
    
    all_scores_2 = [[('C', 0.02247014831444764), ('T', 0.057903843733384308), ('A', 0.30370837683591219),('G', 0.44803586793255664)],
              [('T', 0.046608227674354567), ('G', 0.048827667087419063), ('A', 0.084338697696451109), ('C', 0.92994511407402669)],
              [('G', 0.0), ('T', 0.011098351287382456), ('A', 0.022196702574764911), ('C', 0.9164301607015951)],
              [('C', 0.020803153636453006), ('T', 0.078011826136698756), ('G', 0.11268374886412044),('A', 0.65529933954826969)],
              ]
    all_scores_1 = [[('C', 0.02247014831444764), ('T', 0.057903843733384308), ('A', 0.30370837683591219),('T', 0.44803586793255664)],
              [('T', 0.046608227674354567), ('C', 0.048827667087419063), ('G', 0.084338697696451109), ('A', 0.92994511407402669)],
              [('C', 0.0), ('T', 0.011098351287382456), ('G', 0.022196702574764911), ('T', 0.9164301607015951)],
              [('C', 0.020803153636453006), ('T', 0.078011826136698756), ('A', 0.11268374886412044),('G', 0.65529933954826969)],
              ]
    xs = -7
    ys = 0
    ax.plot([i for i in np.arange(xs,50)], [ys for i in np.arange(xs,50)], color='grey', linewidth=1.0)
    
    #scattered motifs
    for xi in [xs+0, xs+20, xs+40]:
        draw_motif(ax, fig, motif_pwm=all_scores_1, size=10, x_shift=xi, ypos=ys)
        
    for xi in [xs+10, xs+30, xs+50]:
        draw_motif(ax, fig, motif_pwm=all_scores_2, size=10, x_shift=xi, ypos=ys)
    #aligned
    for yi in [-15, -20, -25]:
            draw_motif(ax, fig, motif_pwm=all_scores_1, size=10, x_shift=10, ypos=yi)
        
    for yi in [-15, -20, -25]:
        draw_motif(ax, fig, motif_pwm=all_scores_2, size=10, x_shift=30, ypos=yi)
    
    #dashed lines to connect
    line_length = 10
    x_range = 20
    x_shift_dashed_line = 5
    ax.plot([i for i in np.arange(2, 2+line_length)], [i*-1 for i in np.arange(2,2+line_length)], linestyle='--', color='grey', linewidth=0.5)
    ax.plot([(12+(i/16.0))-x_shift_dashed_line for i in np.arange(11, 11+line_length)], [(i*-1) for i in np.arange(2,2+line_length)], linestyle='--', color='grey', linewidth=0.5)
    ax.plot([(22-((i-22)/2.0))-x_shift_dashed_line for i in np.arange(22, 22+line_length)], [(i*-1) for i in np.arange(2,2+line_length)], linestyle='--', color='grey', linewidth=0.5)
    ax.plot([(32-(i-32))-x_shift_dashed_line for i in np.arange(32, 32+line_length)], [(i*-1) for i in np.arange(2,2+line_length)], linestyle='--', color='grey', linewidth=0.5)
    ax.plot([(42-(i-42))-x_shift_dashed_line for i in np.arange(42, 42+line_length)], [(i*-1) for i in np.arange(2,2+line_length)], linestyle='--', color='grey', linewidth=0.5)
    ax.plot([(52-(i-52))-x_shift_dashed_line for i in np.arange(52, 52+line_length)], [(i*-1) for i in np.arange(2,2+line_length)], linestyle='--', color='grey', linewidth=0.5)
    
    #text boxes
    draw_text(ax, x=-18, y=0, text="Mutated motifs\nacross the genome", horizontalalignment='left', fontsize=10)
    draw_text(ax, x=-18, y=-18, text="Alignment of\nmotifs per TF", horizontalalignment='left', fontsize=10)
    draw_text(ax, x=-18, y=-30, text="Significance test", horizontalalignment='left', fontsize=10)
    draw_marker(ax, x=11.5, y=-30, marker='*', color='grey', markersize=10)
    
    ax.set_xlim(-15,57)
    ax.set_ylim(-35,10)
    ax.set_frame_on(False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    return

def draw_sigTFs(ax, sig_input_file, cohort, ylabel):
    
    df = pd.read_table(sig_input_file, sep='\t')
    df = df[df['Cohorts']==cohort]
    df = df[df['Mean #Mutated Motifs in Simulated Sets']>=100.0]
    df['FDR'] = np.where(df['FDR']==0.0, 1e-300, df['FDR'])#to replace zero with a number that can be converted to log10
    df['FDR'] = df['FDR'].apply(lambda x: np.log10(x)*-1)
    df['Enrichment'] = df['#Mutated Motifs']/df['Mean #Mutated Motifs in Simulated Sets']
    
    ax.scatter(df['Enrichment'], df['FDR'], color='grey')
    labels_to_plot_df = df[df['FDR']>100]
    for i, r in labels_to_plot_df.iterrows():
        rotation = 0
        y_shift = 30
        x_shift = 0
        if 'TF Positions' in r.keys():
            label = r['TF Positions'].split('_')[0] + '-'+ r['TF Positions'].split('#')[-1]
            rotation = -45
            y_shift = 60
            x_shift = 0.5
        else:
            label = r['TFs'].split('_')[0]
        
        draw_text(ax, x=r['Enrichment']+x_shift, y=r['FDR']-y_shift, text=label, horizontalalignment='center', color='red', fontsize=8, rotation=rotation)
    #draw_text(ax, x, y, text, color, fontsize, horizontalalignment, rotation, verticalalignment)
    
    ax.minorticks_off()
    step_fdr =100
    if df['FDR'].max()<100:
        step_fdr=10
    step_enrich =1
    if df['Enrichment'].max()>5:
        step_enrich=2
    ax.set_xticks(np.arange(0,df['Enrichment'].max()+1,step_enrich))
    ax.set_yticks(np.arange(0,df['FDR'].max()+step_fdr, step_fdr))
    ax.set_xlim(0,df['Enrichment'].max()+0.5)
    ax.set_ylim(0,df['FDR'].max()+10)
    ax.set_xlabel("Encirhment over background")
    if ylabel:
        ax.set_ylabel('Adjusted P-values (-log10)')
    
    return

def draw_tfs_per_cohort(ax, fig, sig_input_file, tfs=['CEBPB_MA0466.2', 'CTCF_MA0139.1', 'CEBPG_MA0838.1', 
                                                      ]):#'Tcf12_MA0521.1', 'Gabpa_MA0062.2', 'ELK4_MA0076.2', 'ETV1_MA0761.1']):#'ZNF740_MA0753.1', 
    df = pd.read_table(sig_input_file, sep='\t')
    df = df[df['TFs'].isin(tfs)]
    
    #df = df[df['Mean #Mutated Motifs in Simulated Sets']>=100.0]
    df['TFs'] = df['TFs'].apply(lambda x: x.split('_')[0].upper())
    df['Cohorts'] = df['Cohorts'].apply(lambda  
                                        x: x.replace('All-tumors-without-Lymphatic-system-Skin-Melanoma', 'ATELM').replace('-system-tumors', ''))
    df['FDR'] = np.where(df['FDR']==0.0, 1e-300, df['FDR'])#to replace zero with a number that can be converted to log10
    df['FDR'] = df['FDR'].apply(lambda x: np.log10(x)*-1)
    df['FDR'] = np.where(df['FDR']<0.0, 0.0, df['FDR'])
    #df['#Mutated Motifs'] = np.where(df['#Mutated Motifs']<100.0, 0.0, df['#Mutated Motifs'])
    df['Enrichment'] = df['#Mutated Motifs']/df['Mean #Mutated Motifs in Simulated Sets']
    #conditions to keep
    df = df[df['#Mutated Motifs']>=100]
    #df = df[df['FDR']>=2]
    df = df[df['Enrichment']>=1]
    print (df['Enrichment'].min())
    print (len(df))
    df_pivot = df.pivot(index='Cohorts', columns='TFs', values='Enrichment')
    ax.set_frame_on(False)
    cbar_ax = fig.add_axes([.94, 0.05, .015, .15])
    plt_res = sns.heatmap(df_pivot, ax=ax, robust=True, center=1, square=True, 
                          cbar_ax=cbar_ax, cbar=True, cbar_kws={'ticks': [0,1.5,3]}, cmap = "YlGnBu")
    #plt.xticks(rotation=90)
    plt_res.set_yticklabels(plt_res.get_yticklabels(), rotation=0)
    plt_res.set_xticklabels(plt_res.get_xticklabels(), rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt_res.set_xlabel(plt_res.get_xlabel(), None)
    plt_res.set_ylabel(plt_res.get_ylabel(), None)
    ax.minorticks_off()
    
    return

def draw_motiflog(ax, fig, motif_pwm, fontsize=20):
    draw_motif(ax, fig, size=fontsize, motif_pwm=motif_pwm, x_shift=0, ypos=0, add_axis=True)

def draw_bar_chart(ax, sig_input_file, tf='CEBPB_MA0466.2', cohort='All-tumors-without-Lymphatic-system-Skin-Melanoma', motif_length=10):
    df = pd.read_table(sig_input_file, sep='\t')
    tf_pos = [tf+'#'+str(i) for i in range(1,motif_length+1)]
    df =  df[df['Cohorts']==cohort]
    df = df[df['TF Positions'].isin(tf_pos)]
    df['TF Positions'] = df['TF Positions'].apply(lambda x: int(x.split('#')[-1]))
    df.sort_values(by='TF Positions', inplace=True)
    plt_results = sns.barplot(x='TF Positions', y='#Mutated Motifs', ax=ax, data=df, ci=None, estimator=sum, color='grey', edgecolor = "none")#, palette="PRGn")
    sns.pointplot(data=df, x='TF Positions', y='Mean #Mutated Motifs in Simulated Sets', ci=None, ax=ax, linewidth=1.0, color='red', markers='')
    if df['#Mutated Motifs'].max()>500:# and df['#Mutated Motifs'].max() % 500 > 100:
        plt_results.set_yticks(np.arange(0,df['#Mutated Motifs'].max()+500,500))
    else:
        plt_results.set_yticks(np.arange(0,df['#Mutated Motifs'].max()+10,100))
    plt_results.axes.set_ylabel("Mutation frequency")
    plt_results.set_xticklabels(plt_results.get_xticklabels(), rotation=0)
    plt_results.set_xlabel(tf.split('_')[0]+' Motif')
    
    
def draw_fig2(sig_input_file, sig_input_file_pos, pwm_file, output_dir):
    fig = plt.figure(figsize=(12,8), linewidth=1.0)#design a figure with the given size
    gs = gridspec.GridSpec(4, 7, height_ratios=[4,4,2,4], width_ratios=[4,4,0.5,4,3,4,4], wspace=0.0, hspace=0.0)#create 4 rows and three columns with the given ratio for each
    #sketch
    ax0 = fig.add_subplot(gs[0, 0:5])
    #scatter plots
    ax1 = fig.add_subplot(gs[1, 0:3])
    ax2 = fig.add_subplot(gs[1, 3:5])
    #heatmap
    ax3 = fig.add_subplot(gs[:, 5:])
    #motif logos
    ax4 = fig.add_subplot(gs[2, 0:2])
    ax5 = fig.add_subplot(gs[2, 2:5])
    #ax6 = fig.add_subplot(gs[2, 4:6])
    #bar chars
    ax7 = fig.add_subplot(gs[3, 0:2])
    ax8 = fig.add_subplot(gs[3, 2:5])
    #ax9 = fig.add_subplot(gs[3, 4:6])
    
    gs.tight_layout(fig, pad=2, h_pad=0.0, w_pad=0.0)
    
    sns.despine(top=True, right= True)
    #draw alignment schema
    #sig_input_file = '/home/huum/projs/regMotifs/analysis_exclVEP/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_sigTFs_0.05.tsv'
    #sig_input_file_pos = '/home/huum/projs/regMotifs/analysis_exclVEP/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_sigTFpos_0.05.tsv'
    #sig_input_file = '/Users/karolinasg/Documents/pcawg/NEW_RESULTS_removig_VEP_23_october/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_sigTFs_0.05.tsv'
    #sig_input_file_pos = '/Users/karolinasg/Documents/pcawg/NEW_RESULTS_removig_VEP_23_october/merged200bp_extended200bp_nofullexon_pancan/combined_rand103setsTFsigQval0.05_sigTFpos_0.05.tsv'
    
    
    draw_motifalignment(ax0,fig)
    draw_sigTFs(ax1, sig_input_file, cohort='All-tumors-without-Lymphatic-system-Skin-Melanoma', ylabel='Adjusted P-values (-log10)')
    draw_sigTFs(ax2, sig_input_file_pos, cohort='All-tumors-without-Lymphatic-system-Skin-Melanoma', ylabel=False)
    draw_tfs_per_cohort(ax3, fig, sig_input_file)
    
    #pwm_file = '/Users/karolinasg/Documents/pcawg/analysis/JASPAR_CORE_2016_vertebrates.meme'
    #pwm_file = '/home/huum/projs/regMotifs/datafiles/Motifs/JASPAR_CORE_2016_vertebrates.meme'

    tf_pwms = get_freq_per_motif(pwm_file)
    
    motif_name = 'CEBPB_MA0466.2'
    draw_bar_chart(ax7, sig_input_file_pos, tf=motif_name, motif_length=len(tf_pwms[motif_name]))
    draw_motiflog(ax4, fig, tf_pwms[motif_name], fontsize=35)
    
    motif_name = 'CTCF_MA0139.1'
    draw_bar_chart(ax8, sig_input_file_pos, tf=motif_name, motif_length=len(tf_pwms[motif_name]))
    draw_motiflog(ax5, fig, tf_pwms[motif_name], fontsize=25)
    
    '''motif_name = 'ZNF740_MA0753.1'
    draw_bar_chart(ax9, sig_input_file_pos, tf=motif_name, motif_length=len(tf_pwms[motif_name]))
    draw_motiflog(ax6, fig, tf_pwms[motif_name])
    '''
    #plt.savefig("/Users/karolinasg/Documents/pcawg/NEW_RESULTS_removig_VEP_23_october/plots/Fig2.pdf")#, bbox_inches='tight')
    plot_dir = output_dir + "/Fig2.pdf"
    plt.savefig(plot_dir)#, bbox_inches='tight')
    plt.close()
    
    return




def parse_args():
    '''Parse command line arguments'''
    
    parser = argparse.ArgumentParser(description='Plot Figure 2')
    parser.add_argument('--sig_input_file', default='', help='')
    parser.add_argument('--sig_input_file_pos', default='', help='')
    parser.add_argument('--pwm_file', default='', help='')
    parser.add_argument('--output_dir', default='', help='')
    

    
    return parser.parse_args(sys.argv[1:])


if __name__ == '__main__':
    
    args = parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    
    print("Generate Fig 2")

    draw_fig2(args.sig_input_file, args.sig_input_file_pos, args.pwm_file, args.output_dir)
