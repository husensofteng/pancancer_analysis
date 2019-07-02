'''
Created on Dec 7, 2016

@author: husensofteng
'''
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.backends.backend_pdf

matplotlib.style.use('ggplot')
    
def plot_dataframe(results_df, plot_per_group=False, group_by_cols = ['name'], output_fig_name="DBOut.png"):
    
    plt.figure();
    print results_df.head()
    print results_df.tail()
    fig = ""
    if plot_per_group and len(group_by_cols)>0:
        ax = results_df.boxplot(by=group_by_cols)
        fig = ax[0][0].get_figure()
    else:
        ax = results_df.plot.box()
        fig = ax.get_figure()
    
    fig.savefig(output_fig_name)
    
    return

if __name__ == '__main__':
    results_df = ""#get it elsewhere
    plot_dataframe(results_df, plot_per_group=False, group_by_cols = ['name'], output_fig_name="DBOut.png")
    