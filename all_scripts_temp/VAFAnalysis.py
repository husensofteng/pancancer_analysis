'''
Created on 25 Jul 2017

@author: husensofteng
'''
import matplotlib
matplotlib.use('TkAgg')
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
import sys

def plot_vafs(input_file):
    
    #sns.boxplot(x, y, hue, data, order, hue_order, orient, color, palette, saturation, width, fliersize, linewidth, whis, notch, ax)
    df = pd.read_table(input_file, sep='\t', index_col=None)
    print df.head(2)
    print type(df['VAF'])
    sns.boxplot(x=df["Type"], y=df["VAF"])
    plt.savefig(input_file+".pdf")
    return

if __name__ == '__main__':
    plot_vafs(sys.argv[1])
