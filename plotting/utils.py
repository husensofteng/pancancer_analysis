'''
Created on Jun 4, 2017

@author: husensofteng
'''
import matplotlib
from numpy.lib.function_base import average
from collections import Counter
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
import pandas as pd
from operator import itemgetter
import sys
rcParams['svg.fonttype'] = 'none'
    
def plot_sin(ax, shift_x=10, y_shift=0, R=4, B=4, A=3, color='blue', linestyle='-'):
    x = np.arange(0, R, 0.01)
    y = (A*np.sin((np.pi*x)/B))#change 10 to larger value to get a larger curve
    ax.plot(x+shift_x, y+y_shift, color=color, linestyle=linestyle, linewidth=1)
    
def stem_line(ax, x, y, marker, markerfacecolor, markeredgecolor, stemline_color):
    #muts without motif or signal or not effective
    markerline, stemlines, baseline = ax.stem(x, y)
    plt.setp(markerline, marker=marker, markersize=10.0 , markerfacecolor=markerfacecolor, markeredgewidth=2.0, markeredgecolor= markeredgecolor)
    plt.setp(stemlines, color= stemline_color, linestyle='--', linewidth= 2)
    plt.setp(baseline, color='w', linestyle='-', linewidth= 2)
    
def draw_marker(ax, x,y, marker="*", color='green', markersize=10.0, markeredgecolor='black', markeredgewidth=0):
    ax.plot(x,y,color=color, marker=marker, markersize=markersize, markeredgecolor=markeredgecolor, markeredgewidth=markeredgewidth)
    
def draw_text(ax, x,y, text="", color='0.15', fontsize=10, horizontalalignment='center', rotation=0, verticalalignment='center'):
    ax.text(x,y,text, color=color, fontsize=fontsize, horizontalalignment=horizontalalignment, rotation=rotation, verticalalignment=verticalalignment)

def get_xaxis(n, window_size = 50000, d=6000.0):
    return round(((int(math.ceil(n/window_size)))*1.0)/(d*1.0), 6)
    
def get_unique_state(s):
    return Counter(s.split(',')).most_common(1)[0][0]

def getchrnum(chr):
    return int(chr.replace('X', '23').replace('Y', '24').replace('M', '25').replace('chr', ''))

def replace_state(x):
    if 'Tss' in x or 'BivFlnk' in x:
        return 'Tss'
    elif 'Tx' in x:
        return 'Tx'
    elif 'Enh' in x:
        return 'Enh'
    elif 'Repr' in x:
        return 'Repr'
    elif 'Quies' in x or 'Rpts' in x or 'Het' in x:
        return 'Quies'
    else:
        return 'NA'
    
def get_state(s):
    try:
        if 'ChromHMM' in s:
            states = [x.split(':')[1].split('_')[-1].replace('','').replace('','').replace('','').replace('','').replace('','').replace('','') for x in s.split('|') if 'ChromHMM:' in x]
            x = Counter(states).most_common(1)[0][0]
            return replace_state(x) 
        else:
            return 'NA'
    except TypeError:
        return 'NA'
    
class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def draw_motif(ax, fig, fontfamily='Arial', size=80, motif_pwm=[], x_shift = 0, ypos=0, add_axis=False): 
    COLOR_SCHEME = {'G': 'orange', 
                'A': 'green', 
                'C': 'blue', 
                'T': 'red'}
    
    BASES = list(COLOR_SCHEME.keys())
    #ax = fig.add_subplot(111, aspect='equal')#)
    
    font = FontProperties()
    font.set_size(size)
    font.set_weight('light')
    if add_axis:
        ax.set_xticks(range(1,len(motif_pwm)+1))
        #ax.set_yticks(range(ypos,ypos+3))
        #ax.set_xticklabels([''])#range(1,len(all_scores)+1), rotation=0)
        #ax.set_yticklabels(np.arange(ypos,ypos+3,1))    
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
    sns.despine(ax=ax, trim=True, left=True,bottom=True)
    
    trans_offset = transforms.offset_copy(ax.transData, 
                                          fig=fig, 
                                          x=1, 
                                          y=0, 
                                          units='dots')
    
    for index, scores in enumerate(motif_pwm):
        yshift_window = 0
        for base, score in sorted(dict(scores).items(), key=lambda z: z[1], reverse=False): 
            #print 'index,base,score, yshift_window, ypos'
            #print index,base,score, yshift_window, ypos
            txt = ax.text(index+x_shift+1, 
                          ypos, 
                          base, 
                          transform=trans_offset,
                          fontsize=size, 
                          color=COLOR_SCHEME[base],
                          ha='center',
                          fontproperties=font,
                         )
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift_window = (window_ext.height*score)
            trans_offset = transforms.offset_copy(txt._transform, 
                                                  fig=fig,
                                                  y=yshift_window,
                                                  units='points')
        trans_offset = transforms.offset_copy(ax.transData, 
                                              fig=fig, 
                                              x=1, 
                                              y=0, 
                                              units='points')    
    
    #ax.yaxis.set_visible(False)
    #ax.xaxis.set_visible(False)
    #ax.set_frame_on(False)

def get_freq_per_motif(motif_PFM_input_file):
    "given a PFM file return a dict, a key for each tf and the freq as value"
    PFM_motifs_lines = [""]
    with open(motif_PFM_input_file, 'r') as PFM_motifs_infile:
        PFM_motifs_lines = PFM_motifs_infile.readlines()
    
    PFM_motifs_dict = {}
    nucleotides  = ['A', 'C', 'G', 'T']#default nucleotides
    motif_info_sep = ' '
    motif_name = ""
    if 'MEME' in PFM_motifs_lines[0]: 
        motif_info_sep = ' '
        for line in PFM_motifs_lines:
            if 'ALPHABET=' in line:
                nucleotides = []
                ALPHABETS = line.split('=')[1].strip()
                for alph in ALPHABETS:
                    nucleotides.append(alph.upper())
            
            if line.strip()!="" and not line.startswith('letter') and not line.startswith('URL'):
                if line.startswith('MOTIF'):
                    motif_name = line.strip().split(motif_info_sep)[2]+'_'+line.strip().split(motif_info_sep)[1]
                else:
                    if motif_name!="":#if it has been initialized
                        if motif_name not in PFM_motifs_dict:
                            PFM_motifs_dict[motif_name] = []
                        split_line = line.strip().split()
                        freq_per_allele = []
                        for s in split_line:
                            try:
                                freq_per_allele.append(float(s.strip()))
                            except ValueError:
                                continue
                        if len(freq_per_allele)==len(nucleotides): #freq of the 4 nucleotides
                            nucl_weigts = {}
                            #calcualte information content for each position
                            ic_pos_i = 2
                            for i,nucl in enumerate(nucleotides):
                                try:
                                    ic_pos_i += float(freq_per_allele[i])* math.log(float(freq_per_allele[i]), 2)
                                except ValueError:
                                    ic_pos_i+=0.0
                            
                            #multiply probability of each nucl by the information content at that position
                            for i,nucl in enumerate(nucleotides):
                                nucl_weigts[nucl] = float(freq_per_allele[i])*ic_pos_i
                            
                            PFM_motifs_dict[motif_name].append(nucl_weigts)#, C: float(split_line[1]), G: float(split_line[2]), T: float(split_line[3])})
    return PFM_motifs_dict

