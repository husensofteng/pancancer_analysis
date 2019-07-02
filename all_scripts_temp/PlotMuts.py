'''
Created on May 23, 2017

@author: husensofteng
'''
import pandas as pd
import matplotlib.pyplot as plt
import math
from collections import Counter
import numpy as np
import seaborn as sns

sns.set(style="ticks")
#plt.style.use('ggplot')
sns.set_style("whitegrid")
sns.set_context("paper")#talk

def get_xaxis(n):
    return round(((int(math.ceil(n/50000)))*1.0)/6000.0, 6)
    
def get_yaxis(n):
    return n-int(math.ceil(n/10000000))*10000000

def get_unique_state(s):
    return Counter(s.split(',')).most_common(1)[0][0]

def getchrnum(chr):
    return int(chr.replace('X', '23').replace('Y', '24').replace('M', '25').replace('chr', ''))

def plot_muts(muts_input_file, names='chr,start,DNase1,TFBinding'.split(','), usecols = [0,1,10,11],
              motifs=True):
    
    df = pd.read_table(muts_input_file, sep='\t', skiprows=0, header=None, usecols=usecols, names=names)
    #plt.figure(figsize=(12, 12))
    df['chr'] = df['chr'].apply(getchrnum)
    df['x'] = df['chr'] + df['start'].apply(float).apply(get_xaxis)
    if motifs:
        df['Active'] = np.where((df['DNase1']>1e-300) | (df['TFBinding']>1e-300), 'Yes', 'No')#col='ChromatinState', col_wrap=5, 
    else:
        df['Active'] = np.where(('TFBinding' in df['Annotation']) | ('DNase1' in df['Annotation']), 'Yes', 'No')
    #df['x'] = df.start.apply(get_xaxis, args=(int(df['chr'].replace('X','23').replace('Y','24').replace('M','25').replace('chr','')),))
    #df['y'] = df.start.apply(get_yaxis)
    #df['ChromatinState'] = df['ChromatinState'].apply(get_unique_state)
    
    active_windows = []
    active_window_counts = []
    unactive_windows = []
    unactive_window_counts = []
    chromosomes = []
    for active_label, df_activity in df.groupby('Active'):
        for chr_label, df_chr in df_activity.groupby('chr'):
            chromosomes.append(chr_label)
            windows = []
            window_counts = []
            for window_label, df_window in df_chr.groupby('x'):
                if len(df_window)==0:
                    windows.append('nan')
                    window_counts.append('nan')
                windows.append(window_label)
                window_counts.append(len(df_window))
                
                if len(df_window)==0:
                    windows.append('nan')
                    window_counts.append('nan')
                    
            if active_label=='Yes':
                active_windows.extend(windows)
                active_window_counts.extend(window_counts)
                active_windows.append('nan')
                active_window_counts.append('nan')
            else:
                unactive_windows.extend(windows)
                unactive_window_counts.extend(window_counts)
                unactive_windows.append('nan')
                unactive_window_counts.append('nan')
    return active_windows, active_window_counts, unactive_windows, unactive_window_counts, chromosomes
    
    #dfg.vals.plot(kind="kde", ax=ax, label=label)
    #plt.legend()    
    '''number_muts_per_mb = {}
    for i,window in enumerate(df['x']):
        try :
            number_muts_per_mb[window][df['Active'][i]]+=1
        except KeyError:    
            try:
                number_muts_per_mb[window][df['Active'][i]] = 0
            except KeyError:
                number_muts_per_mb[window] = {df['Active'][i]:0}
    print number_muts_per_mb
    #df_total_numbers_per_mb_per_state = pd.DataFrame(names=['MB', 'Active', 'Number'])
    #print(pd.DataFrame.from_dict(number_muts_per_mb, orient='index').head(10))
    active_items =[]
    unactive_items =[]
    
    for k in sorted(number_muts_per_mb.keys()):
        for v in number_muts_per_mb[k]:
            if v == "Yes":
                active_items.append([k,v,number_muts_per_mb[k][v]])
            else:
                unactive_items.append([k,v,number_muts_per_mb[k][v]])
    
        #active_items.append([k,v,'nan'])
        #unactive_items.append([k,v,'nan'])
               
    df_av = pd.DataFrame(active_items, columns=['Window', 'Active', 'Number'])
    df_uv = pd.DataFrame(unactive_items, columns=['Window', 'Active', 'Number'])
    #df_v['Number'] = df_v['Number'].apply(math.log)
    
    sns.set_style("white")
    plt.plot(df_uv['Window'], df_uv['Number'], 'grey', df_av['Window'], df_av['Number'], 'r', linewidth=0.5)
    sns.despine(right=True, top=True)
    plt.savefig(fig_out+"_avgperWindow" +".pdf")
    plt.close()
    
    '''
    
    '''plt_results = sns.lmplot(x='Window', y='Number', data=df_v, hue='Active', palette=dict(Yes='red', No='grey'), #col='ChromatinState', col_wrap=5, size=10,
                             fit_reg=False, legend=True, legend_out=True, scatter=True, markers='o', scatter_kws={"s": 1})
    plt_results.set_axis_labels("1Mb Window (hg19)", "#Mutations (log)").set(ylim=[0,df_v['Number'].max()], xlim=[1, df_v['Window'].max()])
    plt_results.savefig(fig_out+"_avgperWindow" +".pdf")
    '''
    
    '''plt_results = sns.lmplot(x='x', y='y', data=df, hue='Active', palette=dict(Yes='red', No='grey'), col='ChromatinState', col_wrap=5, size=10,
                             fit_reg=False, legend=True, legend_out=True, scatter=True, markers='o', scatter_kws={"s": 1})
    plt_results.set_axis_labels("10Kb Window (chr5)", "Mutations").set(ylim=[0,df['y'].max()], xlim=[0, df['x'].max()+1])
    #sns.despine(right=True, top=True)
    plt_results.savefig(fig_out+"_ActivePerChromatinState" +".pdf")
    plt.close()
    #df.plot.scatter(x='x', y='y')
    #plt.show()
    '''
    #print df.head(2)

def plot_lines(windows, chromosomes, fig_out):
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(14,6))
    plt.plot(windows[0], windows[1], windows[2], 
             windows[3], windows[4], windows[5],
             windows[6], windows[7], windows[8],
             windows[9], windows[10], windows[11],)
    plt.xticks(chromosomes)
    plt.savefig(fig_out+"_perchr" +".pdf")
    
if __name__ == '__main__':  
    muts_input_file = 'analysis/observed_agreement_22May2017_annotated.bed10'
    motif_muts_input_file = 'analysis/motifmuts_all.bed12'
    fig_out='analysis/mutsplots_all'
    active_windows, active_window_counts, unactive_windows, unactive_window_counts, chromosomes = plot_muts(muts_input_file, names='chr,start,Annotation'.split(','),
              usecols = [0,1,9], motifs=False)
    print len(active_windows)
    windows = [active_windows, active_window_counts, '#e6e6e6', unactive_windows, unactive_window_counts, '#ccccff']
    active_windows, active_window_counts, unactive_windows, unactive_window_counts, chromosomes = plot_muts(motif_muts_input_file)
    print len(active_windows)
    windows.extend([active_windows, active_window_counts, '#999999', unactive_windows, unactive_window_counts, '#ff3333'])
    plot_lines(windows, chromosomes, fig_out)
    