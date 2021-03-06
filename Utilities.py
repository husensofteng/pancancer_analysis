'''
Created on Dec 9, 2016

@author: husensofteng
'''
import sys, os
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
import pickle
import json
from pybedtools import BedTool, set_tempdir, cleanup
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import binom, hypergeom
from scipy.stats import binom
import shutil
from shutil import copyfile
from multiprocessing import Pool
from itertools import product
import csv
import math
import glob
#from score_motifs_tissuepertable import open_connection, close_connection

#import matplotlib.backends.backend_pdf

#matplotlib.style.use('ggplot')
#temp_dir = 'tmp_pybedtoos'
#if not os.path.exists(temp_dir):
#    os.mkdir(temp_dir)  

#use scratch
#temp_dir = tmp_dir


def unify_muts(annotated_mutations_input_file, output_extension, filter_mut_motifs=True, filter_cond = "", operation_on_unify='mean'):
    
    annotated_mutations_grouped_file=annotated_mutations_input_file + output_extension + "_groupedbymut"
    annotated_mutations_grouped_output_file=annotated_mutations_grouped_file+"withmotifinfo"
    fsep = '\t'
    vsep = '#'
    if not os.path.exists(annotated_mutations_grouped_output_file):
        awk_stmt = """awk 'BEGIN{{FS=OFS="{fsep}"}}{{{filter_cond} {{gsub("X","23",$1); gsub("Y","24",$1);gsub("MT","25",$1);gsub("M","25",$1);gsub("chr","",$1); $2=($2/1)*1; $3=($3/1)*1;   
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10+$11,$10+$11"{vsep}"$10"{vsep}"$11"{vsep}"$12"{vsep}"$13"{vsep}"$14"{vsep}"$15"{vsep}"$16"{vsep}"$17"{vsep}"$18"{vsep}"$19"{vsep}"$20"{vsep}"$21"{vsep}"$22"{vsep}"$23"{vsep}"$24"{vsep}"$25"{vsep}"$26"{vsep}"$27"{vsep}"$28"{vsep}"$29"{vsep}"$30"{vsep}"$31"{vsep}"$32}}}}' {infile} | 
        sort -k1,1n -k2,2n -k3,3n -k4 -k5 -k6 -k7 -k8 -k9 | 
        groupBy -g 1-9 -c 10,11 -o {operation_on_unify},collapse > {grouped_file}""".format(
                                            fsep=fsep, vsep=vsep, filter_cond=filter_cond, infile=annotated_mutations_input_file, operation_on_unify=operation_on_unify, grouped_file=annotated_mutations_grouped_file)
        #print(awk_stmt)
        os.system(awk_stmt)
        #print("groupBy mutation is done")
    
    return annotated_mutations_grouped_file
    
def get_max_motif_in_grouped_muts(annotated_mutations_grouped_file):
    
    annotated_mutations_grouped_output_file=annotated_mutations_grouped_file+"withmotifinfo"
    fsep = '\t'
    vsep = '#'
    if os.path.exists(annotated_mutations_grouped_output_file):
        return annotated_mutations_grouped_output_file
    #print("Reading the grouped mutations file")
    score_index_in_grouped_file = 9
    motif_index_in_grouped_file = 10
    score_index_in_motif_info = 0
    ref_alt_index_in_motif_info = 4
    mutpos_index_in_motif_info = 5
    motifname_index_in_motif_info = 9
    with open(annotated_mutations_grouped_file, 'r') as grouped_file, open(annotated_mutations_grouped_output_file, 'w') as annotated_mutations_grouped_outfile_bed12:
        l = grouped_file.readline()
        while l:
            sl = l.strip().split(fsep)
            max_mut_score = float(sl[score_index_in_grouped_file]) 
            motifs_info = [s.strip().split(vsep) for s in sl[motif_index_in_grouped_file].split(',')]
            for motif_info in motifs_info:
                mut_motif_score = float(motif_info[score_index_in_motif_info])
                if mut_motif_score>=max_mut_score:#append the max score to the chromatin status of the overlapping motif
                    annotated_mutations_grouped_outfile_bed12.write(
                        l.strip() + '\t' + motif_info[ref_alt_index_in_motif_info]+ '\t' + 
                        motif_info[mutpos_index_in_motif_info]+ '\t' + 
                        motif_info[motifname_index_in_motif_info] + '\t' + 
                        vsep.join(motif_info) + '\n')
                    break
            l = grouped_file.readline()
            
    return annotated_mutations_grouped_output_file

def get_scores(annotated_mutations_grouped_file):
    
    fsep = '\t'
    vsep = '#'
    
    overall_scores_file = annotated_mutations_grouped_file+"_overallscores"
    overall_scores_chromatin_status_file = annotated_mutations_grouped_file+"_chromatinscores"
    overall_scores_phenotype_file = annotated_mutations_grouped_file+"_phenotypescores"
    mut_scores_overall_combined_file = annotated_mutations_grouped_file+"_combined"
    if (os.path.exists(overall_scores_file) and os.path.exists(overall_scores_chromatin_status_file) and os.path.exists(overall_scores_phenotype_file) 
        and os.path.exists(mut_scores_overall_combined_file)):
        print('Loading scores')
        with open(overall_scores_file, 'r') as gp:
            mut_scores_overall = pickle.load(gp)
        with open(overall_scores_chromatin_status_file, 'r') as gp:
            mut_scores_overall_per_chromatin_status = pickle.load(gp)
        with open(overall_scores_phenotype_file, 'r') as gp:
            mut_scores_overall_per_phenotype = pickle.load(gp)
        with open(mut_scores_overall_combined_file, 'r') as gp:
            mut_scores_overall_combined = pickle.load(gp)
        
        print(mut_scores_overall_per_chromatin_status.keys())
        print(mut_scores_overall_per_phenotype.keys())
        print(mut_scores_overall_combined.keys())
        
        return mut_scores_overall_combined, mut_scores_overall, mut_scores_overall_per_chromatin_status, mut_scores_overall_per_phenotype
    print("Reading the grouped file")
    motif_index_in_grouped_file = 11
    chromatin_status_index_in_motif_info = 13
    phenotype_index_in_grouped_file = 5
    score_index_in_motif_info = 0
    mut_scores_overall = []
    mut_scores_overall_combined = {'scores':[], 'chromatin_status': [], 'cancer_types':[]}
    mut_scores_overall_per_chromatin_status = {}
    mut_scores_overall_per_phenotype = {}
    with open(annotated_mutations_grouped_file) as grouped_file:
        l = grouped_file.readline()
        while l:
            sl = l.strip().split(fsep)
            motif_info = sl[motif_index_in_grouped_file].split(vsep)
            max_mut_score = float(motif_info[score_index_in_motif_info])
            mut_scores_overall.append(max_mut_score)
            
            phenotype = sl[phenotype_index_in_grouped_file]
            try:
                mut_scores_overall_per_phenotype[phenotype].append(max_mut_score)
            except KeyError:
                mut_scores_overall_per_phenotype[phenotype] = [max_mut_score]
            
            chromatin_status = motif_info[chromatin_status_index_in_motif_info]
            try:
                mut_scores_overall_per_chromatin_status[chromatin_status].append(max_mut_score)
            except KeyError:
                mut_scores_overall_per_chromatin_status[chromatin_status] = [max_mut_score]
            
            mut_scores_overall_combined['scores'].append(max_mut_score)
            mut_scores_overall_combined['chromatin_status'].append(chromatin_status)
            mut_scores_overall_combined['cancer_types'].append(phenotype)
            
            l = grouped_file.readline()
            
    #print(mut_scores_overall)
    with open(overall_scores_file, 'w') as gp:
        pickle.dump(mut_scores_overall, gp)
    with open(overall_scores_chromatin_status_file, 'w') as gp:
        pickle.dump(mut_scores_overall_per_chromatin_status, gp)
    with open(overall_scores_phenotype_file, 'w') as gp:
        pickle.dump(mut_scores_overall_per_phenotype, gp)
    with open(mut_scores_overall_combined_file, 'w') as gp:
        pickle.dump(mut_scores_overall_combined, gp)
    
    return mut_scores_overall_combined, mut_scores_overall, mut_scores_overall_per_chromatin_status, mut_scores_overall_per_phenotype

def get_mean_and_sd_from_file(simulated_input_file, scores_index = 9, index_mut_type = 6, report_overlall_score = False):
    
    if os.path.exists(simulated_input_file+"_statsdict"):
        with open(simulated_input_file+"_statsdict", 'r') as simulated_infile:
            d = simulated_infile.readline().strip()
            return json.loads(d)
    
    scores_SNPs = []
    scores_MNPs = []
    scores = []
    with open(simulated_input_file, 'r') as simulated_infile:
        l = simulated_infile.readline()
        while l:
            sl = l.strip().split('\t')
            scores.append(float(sl[scores_index]))
            if not report_overlall_score:
                if sl[index_mut_type]=="SNP":
                    scores_SNPs.append(float(sl[scores_index]))
                else:
                    scores_MNPs.append(float(sl[scores_index]))
            l = simulated_infile.readline()
    
    stats_dict = {'scores':scores,'avg': np.mean(scores), 'std': np.std(scores)}
    if not report_overlall_score:
        if len(scores_SNPs)>0:
            stats_dict['avgSNPs'] = np.mean(scores_SNPs)
            stats_dict['stdSNPs'] = np.std(scores_SNPs) 
        if len(scores_MNPs)>0:
            stats_dict['avgMNPs'] = np.mean(scores_MNPs)
            stats_dict['stdMNPs'] = np.std(scores_MNPs)
    
    with open(simulated_input_file+"_statsdict", 'w') as simulated_outfile:
        json.dump(stats_dict, simulated_outfile)
    
    return stats_dict



def assess_stat_muts(muts_input_file, simulated_input_file, observed_output_file, observed_onlysig_output_file, score_index_observed_elements=9, score_index_sim_elements=9, index_mut_type = 6, mut_sig_threshold=0.05, stats_dict = {}):
    
    if os.path.exists(observed_output_file) and os.path.exists(observed_onlysig_output_file):
        return observed_output_file, observed_onlysig_output_file, 'NA'
    
    if len(stats_dict.keys())==0:
        stats_dict = get_mean_and_sd_from_file(simulated_input_file, scores_index=score_index_sim_elements, index_mut_type = index_mut_type)
    print("getting pval for muts in {} using {}: {}".format(muts_input_file, simulated_input_file, stats_dict))
    
    n_sig = 0
    with open(muts_input_file, 'r') as observed_infile, open(observed_output_file, 'w') as observed_outfile, open(observed_onlysig_output_file, 'w') as observed_onlysig_outfile:
        l = observed_infile.readline()
        p_values = []
        lines_to_write = []
        while l:
            sl = l.strip().split('\t')
            p_value = 0.0
            if 'avgSNPs' in stats_dict.keys() and 'avgMNPs' in stats_dict.keys(): 
                if sl[index_mut_type]=='SNP':
                    p_value = get_pval(float(sl[score_index_observed_elements]), stats_dict['avgSNPs'], stats_dict['stdSNPs'])
                else:
                    p_value = get_pval(float(sl[score_index_observed_elements]), stats_dict['avgMNPs'], stats_dict['stdMNPs'])
            else:
                p_value = get_pval(float(sl[score_index_observed_elements]), stats_dict['avg'], stats_dict['std'])
            
            p_values.append(p_value)
            lines_to_write.append(l.strip())
            
            l = observed_infile.readline()
        
        p_values_adjusted = adjust_pvales(p_values)
        for i,lo in enumerate(lines_to_write):
            observed_outfile.write(lo + '\t' + str(p_values[i]) + '\t' + str(p_values_adjusted[i]) + '\n')
            if p_values[i]<mut_sig_threshold:
                observed_onlysig_outfile.write(lo + '\t' + str(p_values[i]) + '\t' + str(p_values_adjusted[i]) + '\n')
                
            #observed_onlysig_outfile.write(l.strip() + '\t' + str(p_value) + '\n')
        return observed_output_file, observed_onlysig_output_file, n_sig

def merge_muts(muts_input_file, merged_muts_output_ext, filter_mut_motifs=False, filter_col_index=15, filter_value=0.05, mut_score_index=9, motifs_col_index =10, ref_alt_col_index=11, mutpos_col_index=12, motifname_col_index=13, motif_col_index=14, distance_to_merge=20):
    
    
    merged_muts_output_file = muts_input_file + merged_muts_output_ext
    
    fsep = '\t'
    vsep = '#'
    msep = '*'
    
    MotifInfo, matching_motifs_sep, MaxMotif_sep = 'MotifInfo', 'MatchingMotifs', 'MaxMotif'
    
    if os.path.exists(merged_muts_output_file):
        return merged_muts_output_file
    #selected_regions_file = 'analysis/encode_merge_cdsAndSpliceSitesSubtracted.bed3'
    cond = ""
    if filter_mut_motifs:
        cond = 'if({}<{})'.format(filter_col_index, filter_value)
    ''''cols info:
    awk print: 1 chr, 2 start, 3 end, 4 score, 5 phenotype, 6 ref_alt, 7 mutpos, 8 motifname, 
               10 donorID, 11 mutinfo*motifsinfo*maxmotifinfo
    
    merge (-c -o) mean(score), count(mutinfo), count(donorID), distinct: phenotype, 
    collapse: ref_alt, mutpos, motifname, score, phenotype, distinct(donorID), donorID, mut_motifs_info'''
    
    awk_stmt = """awk 'BEGIN{{FS=OFS="{fsep}"}}{{{cond}{{ gsub(",", "MotifInfo"); print $1,$2,$3,$10,$6,$12,$13,$14,$9,$1"{vsep}"$2"{vsep}"$3"{vsep}"$4"{vsep}"$5"{vsep}"$6"{vsep}"$7"{vsep}"$8"{vsep}"$9"{vsep}"$10"{matching_motifs_sep}"$11"{MaxMotif_sep}"$15}}}}' {muts_input_file} | sort -k1,1n -k2,2n -k3,3n | mergeBed -i stdin -d {distance_to_merge} -c 4,10,9,5,6,7,8,4,5,9,9,10 -o mean,count_distinct,count_distinct,distinct,collapse,collapse,collapse,collapse,collapse,distinct,collapse,collapse > {merged_muts_output_file}""".format(**locals()) # | awk 'BEGIN{{FS=OFS="\t"}}{{if($2!=$3){{$2=$2+1; $3=$3-1}}; if($2>$3){{$2=$2-1; $3=$3+1}}; print}}
    #awk_stmt = """awk 'BEGIN{{FS=OFS="{fsep}"}}{{{cond}{{ gsub(",", "MotifInfo"); print $1,$2,$3,$10,$6,$12,$13,$14,$9,$1"{vsep}"$2"{vsep}"$3"{vsep}"$4"{vsep}"$5"{vsep}"$6"{vsep}"$7"{vsep}"$8"{vsep}"$9"{vsep}"$10"{matching_motifs_sep}"$11"{MaxMotif_sep}"$15}}}}' {muts_input_file} | sort -k1,1n -k2,2n -k3,3n | intersectBed -split -wo -a stdin -b {selected_regions_file} | sort -k11,11n -k12,12n -k13,13n | groupBy -g 11,12,13 -c 4,10,9,5,6,7,8,4,5,9,9,10 -o sum,count_distinct,count_distinct,distinct,collapse,collapse,collapse,collapse,collapse,distinct,collapse,collapse > {merged_muts_output_file}""".format(**locals())
    #print(awk_stmt)
    os.system(awk_stmt)
    #print("merge mutations process is done")
    
    return merged_muts_output_file
    
def get_chr_lengths(chr_lengths_file):
    chr_lengths = {}
    with open(chr_lengths_file, 'r') as chr_min_max_ifile:
        lines = chr_min_max_ifile.readlines()
        chr_names = range(1,26)
        for l in lines:
            if not l.startswith('chr') and l!="" and not l.startswith('//') and not l.startswith('#'):
                sl = l.strip().split('\t')
                chr_name = sl[0].replace('X', '23').replace('Y','24').replace('MT','25').replace('M','25')
                if int(chr_name) in chr_names: 
                    chr_lengths[int(chr_name)] = int(sl[1])
    return chr_lengths

def sum_fscore_motif_breaking_score(feature,fscore_index, motif_breaking_score_index):
    if(feature[fscore_index] != '.'):
        sums = float(feature[fscore_index]) + float(feature[motif_breaking_score_index])
        feature[fscore_index] = str(sums)
    return feature

def empirical_pval(sl, stats_dict_scores):
    p_values=[]
    for score in sl:
        scores_higher_than_observed = [i for i in stats_dict_scores if i >= score]
        p_value= len(scores_higher_than_observed)/(len(stats_dict_scores))
        if p_value==0.0:
            p_value=1/103
        p_values.append(p_value)
    
    return p_values

def empirical_pval_global(dict_lines_observed_split, stats_dict_scores, pval_file):
    with open(pval_file, 'a') as pval_ifile:
        for index in dict_lines_observed_split:
            score = dict_lines_observed_split[index]
            scores_higher_than_observed = [i for i in stats_dict_scores if i >= score]
            p_value= len(scores_higher_than_observed)/(len(stats_dict_scores))
            if p_value==0.0:
                p_value=1/103
            pval_ifile.write(str(index) + '\t' + str(p_value)+'\n')
    return pval_file

def empirical_pval_local_window(dict_lines_observed_split, pval_file):
   # dict_p_values={}
    with open(pval_file, 'a') as pval_ifile:
        for index in dict_lines_observed_split:
    
            simulated_score_vec=dict_lines_observed_split[index][1]
            scores_len=len(simulated_score_vec)
            element_score = float(dict_lines_observed_split[index][0])
            scores_higher_than_observed = [i for i in simulated_score_vec if i >= element_score]
            p_value= len(scores_higher_than_observed)/scores_len
            if p_value==0.0:
                p_value=1/103
            pval_ifile.write(str(index) + '\t' + str(p_value)+ '\n')
    return pval_file


def split_dict_equally(input_dict, chunks=2):
    # prep with empty dicts
    return_list = [dict() for idx in range(chunks)]
    idx = 0
    for k,v in input_dict.items():
        return_list[idx][k] = v
        if idx < chunks-1:  # indexes start at 0
            idx += 1
        else:
            idx = 0
    return return_list

def assess_stat_elements_local_domain(observed_input_file, simulated_input_file, merged_elements_statspvalues, merged_elements_statspvaluesonlysig, 
                                      chr_lengths_file, local_domain_window=25000, 
                                      merged_mut_sig_threshold = 0.05, score_index_observed_elements=4, score_index_sim_elements=4, p_value_on_score=False):
    
    if os.path.exists(merged_elements_statspvalues) and os.path.exists(merged_elements_statspvaluesonlysig):
        return merged_elements_statspvalues, merged_elements_statspvaluesonlysig, 'NA'
    
    simulated_input_file_sort=simulated_input_file+'_sort'
         
         
         
    if not os.path.exists(simulated_input_file_sort):
            os.system("""sort -k1,1n -k2,2n {} > {}""".format(simulated_input_file,simulated_input_file_sort))
    
    #extend elements size
    "replace chr, X, Y, add Line Number to use as window ID and sort by chr,start"
    observed_input_file_temp_file = observed_input_file+"_temp" 
    cmd = """awk 'BEGIN{{OFS="\t"}}{{gsub("chr","",$1); gsub("X", 23, $1); gsub("Y", 24, $1); print $1,$2,$3,NR,$4}}' {} | bedtools slop -g /proj/snic2020-16-50/nobackup/pancananalysis/pancan12Feb2020/cancer_datafiles/chr_order_hg19.txt -b {}| sort -k1,1n -k2,2n > {}""".format(
        observed_input_file, local_domain_window,observed_input_file_temp_file)
    os.system(cmd)
    
    os.system("""awk '{{print $0>"{}""_"$1".bed"}}' {}""".format(
       observed_input_file_temp_file, observed_input_file_temp_file))
    print(glob.glob(observed_input_file_temp_file+'*.bed'))
    pval_files=[]
    for observed_input_file_temp_file_per_chr in glob.glob(observed_input_file_temp_file+'_*.bed'):
        print(observed_input_file_temp_file_per_chr)
    #chr_lengths = get_chr_lengths(chr_lengths_file)
    
        dict_lines_observed = {}
        line_number = 1
        with open(observed_input_file_temp_file_per_chr, 'r') as observed_infile: #, open(observed_input_file_temp_file, 'w') as observed_input_file_temp_ofile:
            l = observed_infile.readline().strip().split('\t')
            while l and len(l)>2:
                dict_lines_observed[int(float(l[3]))] = [l[4],[]]
        #             extended_element_start = (int(l[1])-local_domain_window)
        #             extended_element_end = int(l[2])+local_domain_window
        #             if extended_element_start<0:
        #                 extended_element_start = 0
        #                 extended_element_end += 0 - (int(l[1])-local_domain_window)
        #             if extended_element_end>chr_lengths[int(l[0])]:
        #                 extended_element_end = chr_lengths[int(l[0])]
        #                 extended_element_start -= (int(l[2])+local_domain_window) - chr_lengths[int(l[0])]
        #                 
        #             observed_input_file_temp_ofile.write(l[0] + '\t' + str(extended_element_start) + '\t' + str(extended_element_end) + '\t' + str(line_number) + '\n')
                line_number+=1
                l = observed_infile.readline().strip().split('\t')
        print('observed_input_file: ', observed_input_file_temp_file_per_chr)
        observed_input_file_temp_file_per_chr_sort=observed_input_file_temp_file_per_chr+'_sort'

        if not os.path.exists(observed_input_file_temp_file_per_chr_sort):
            os.system("""sort -k1,1n -k2,2n {} > {}""".format(observed_input_file_temp_file_per_chr, observed_input_file_temp_file_per_chr_sort))
        observed_input_file_obj = BedTool(observed_input_file_temp_file_per_chr_sort)
        
        

            
            
        simulated_input_file_temp = simulated_input_file+"_temp"
        observed_input_file_obj.map(BedTool(simulated_input_file_sort), c=4, o=['collapse'], g='/proj/snic2020-16-50/nobackup/pancananalysis/pancan12Feb2020/cancer_datafiles/chr_order_hg19.txt').saveas(simulated_input_file_temp)
       
       
    
       
        
            

    

        with open(simulated_input_file_temp, 'r') as simulated_input_file_temp_ifile:
            l = simulated_input_file_temp_ifile.readline().strip().split('\t')
            
            while l and len(l)>1:
                
                sim_scores = []
                for x in l[5].split(','):
                    try:
                        sim_scores.append(float(x))
                    except ValueError:
                        sim_scores.append(0.0)
                dict_lines_observed[int(float(l[3]))][1].extend(sim_scores)
                l = simulated_input_file_temp_ifile.readline().strip().split('\t')

        #split dictionery into chunks
        dict_lines_observed_chunks=split_dict_equally(dict_lines_observed, 100)
        
        pval_file=observed_input_file_temp_file_per_chr+'_elem_pval_local'+str(local_domain_window)
    
        if os.path.exists(pval_file):
            os.remove(pval_file)
    

    
        #print('p-value on score local')
        pm = Pool(15)
        pm.starmap(empirical_pval_local_window, product(dict_lines_observed_chunks, [pval_file]))
        pm.close()
        pm.join()
        
        
   
        os.remove(simulated_input_file_temp)
        os.remove(observed_input_file_temp_file_per_chr_sort)
        pval_files.append(pval_file)
        

    combined_pval_file=observed_input_file+'_elem_pval_local'+str(local_domain_window)
    
    
    dict_pvals={}
    p_values=[]
    if not os.path.exists(combined_pval_file):
        with open(combined_pval_file, 'w') as combined_pval_outfile:
            for pval_file in pval_files:
                with open(pval_file, 'r') as pval_ifile:
                    combined_pval_outfile.write(pval_ifile.read())
                    os.remove(pval_file)
                    
    with open(combined_pval_file, 'r') as combined_pval_ifile:                
        l = combined_pval_ifile.readline().strip().split('\t')
    # 
        while l and len(l)>1:
            p_values.append(float(l[1]))
            dict_pvals[int(float(l[0]))]=float(l[1])
            l = combined_pval_ifile.readline().strip().split('\t')
                        

  
    
    
    #l=1


#     p_values = []
#     pvalues_adjusted = []
#     n_sig = 0
#     lines = []
#     dict_pvals={}
#     p_values=pval_df[1]
#     with open(pval_file, 'r') as pval_ifile:
#         l = pval_ifile.readline().strip().split('\t')
# 
#         while l and len(l)>1:
#             p_values.append(float(l[1]))
#             dict_pvals[int(float(l[0]))]=float(l[1])
#             l = pval_ifile.readline().strip().split('\t')
    pvalues_adjusted = p_values

    lambda_factor=np.median(stats.chi2.isf(p_values,1))/stats.chi2.ppf(0.5, 1)
    lambda_values_file=observed_input_file+"_lambda_values_local_window_"+str(local_domain_window)+'.txt'
    lambda_file=open(lambda_values_file, "w")
    lambda_file.write(observed_input_file.split('/')[-1] + '\t' + str(lambda_factor)+'\n')
    lambda_file.close()
    
    #if os.path.exists(pval_file):
    #    os.remove(pval_file)
        
           
    
    
    
    
    #l=1 
    n_sig = 0
    with open(observed_input_file, 'r') as observed_infile, open(merged_elements_statspvalues, 'w') as merged_elements_statspvalues_outfile, open(merged_elements_statspvaluesonlysig, 'w') as merged_elements_statspvaluesonlysig_outfile:
        l = observed_infile.readline()
        l_number=1
        while l:
            sl = l.strip().split('\t')
            #print(dict_pvals[l_number])
            #print(l.strip())
            merged_elements_statspvalues_outfile.write(l.strip() + '\t' + str(dict_pvals[l_number]) + '\t' + str(dict_pvals[l_number]) + '\n')
            if dict_pvals[l_number]<merged_mut_sig_threshold:
                n_sig+=1
                merged_elements_statspvaluesonlysig_outfile.write(l.strip() + '\t' + str(dict_pvals[l_number]) + '\t' + str(dict_pvals[l_number]) + '\n')
            l = observed_infile.readline()
            l_number+=1 
    
    cleanup()
    return merged_elements_statspvalues, merged_elements_statspvaluesonlysig, n_sig




def assess_stat_elements(observed_input_file, simulated_input_file, 
                         merged_elements_statspvalues, 
                         merged_elements_statspvaluesonlysig, 
                         merged_mut_sig_threshold = 0.05, 
                         score_index_observed_elements=4, 
                         score_index_sim_elements=4, p_value_on_score=False):
    
    if os.path.exists(merged_elements_statspvalues) and os.path.exists(merged_elements_statspvaluesonlysig):
        return merged_elements_statspvalues, merged_elements_statspvaluesonlysig, 'NA'
    
    stats_dict = get_mean_and_sd_from_file(simulated_input_file, scores_index=score_index_sim_elements, report_overlall_score=True)
    #print("getting pval for elements in {} using {}: {}".format(observed_input_file, simulated_input_file))
    
    #extend each merged region by wbp and combine all columns into one; intersect the entire list with the combine simulated list of elemenets, group by the column ID, take average and std from the grouping;  
    
    
    if p_value_on_score:
        
        dict_lines_observed_score = {}
        line_number = 1
        with open(observed_input_file, 'r') as observed_infile: #, open(observed_input_file_temp_file, 'w') as observed_input_file_temp_ofile:
            l = observed_infile.readline().strip().split('\t')
            while l and len(l)>3:
                dict_lines_observed_score[line_number] = float(l[3])
                line_number+=1
                l = observed_infile.readline().strip().split('\t')
        dict_lines_observed_chunks=split_dict_equally(dict_lines_observed_score, 100)
        
        pval_file=observed_input_file+'_elem_pval'
        if os.path.exists(pval_file):
            os.remove(pval_file)
        pm = Pool(15)
        pm.starmap(empirical_pval_global, product(dict_lines_observed_chunks,  [stats_dict['scores']], [pval_file]))
        pm.close()
        pm.join()
        
        dict_pvals={}
        p_values=[]
        with open(pval_file, 'r') as pval_ifile:
            l = pval_ifile.readline().strip().split('\t')
        
            while l and len(l)>1:
                p_values.append(float(l[1]))
                dict_pvals[int(float(l[0]))]=float(l[1])
                l = pval_ifile.readline().strip().split('\t')
        pvalues_adjusted = p_values
                
        lambda_factor=np.median(stats.chi2.isf(p_values,1))/stats.chi2.ppf(0.5, 1)
        #print(observed_input_file+ "_lambda_values_whole_genome.txt")
        lambda_file=open(observed_input_file+ "_lambda_values_whole_genome.txt", "w")
        lambda_file.write(observed_input_file.split('/')[-1] + '\t' + str(lambda_factor)+'\n')
        lambda_file.close()    
        
        
        lines = []
        n_sig = 0
        with open(observed_input_file, 'r') as observed_infile, open(merged_elements_statspvalues, 'w') as merged_elements_statspvalues_outfile, open(merged_elements_statspvaluesonlysig, 'w') as merged_elements_statspvaluesonlysig_outfile:
            l = observed_infile.readline()
            l_number=1
            while l:
                sl = l.strip().split('\t')
                #print(dict_pvals[l_number])
                #print(l.strip())
                merged_elements_statspvalues_outfile.write(l.strip() + '\t' + str(dict_pvals[l_number]) + '\t' + str(dict_pvals[l_number]) + '\n')
                if dict_pvals[l_number]<merged_mut_sig_threshold:
                    n_sig+=1
                    merged_elements_statspvaluesonlysig_outfile.write(l.strip() + '\t' + str(dict_pvals[l_number]) + '\t' + str(dict_pvals[l_number]) + '\n')
                l = observed_infile.readline()
                l_number+=1
        
        if os.path.exists(pval_file):
            os.remove(pval_file)
            
    else:
        p_values = []
        lines = []
        pvalues_adjusted = []
        n_sig = 0
        with open(observed_input_file, 'r') as observed_infile, open(merged_elements_statspvalues, 'w') as merged_elements_statspvalues_outfile, open(merged_elements_statspvaluesonlysig, 'w') as merged_elements_statspvaluesonlysig_outfile:
            l = observed_infile.readline()
            while l:
                sl = l.strip().split('\t')
                #get avg and std from simulated merged elements located within w bps of this region
                p_value = get_pval(float(sl[score_index_observed_elements]), stats_dict['avg'], stats_dict['std'])
                p_values.append(p_value)
                lines.append(l.strip())
                l = observed_infile.readline()
            if len(p_values)>0:
                pvalues_adjusted = adjust_pvales(p_values)
            
            for i,l in enumerate(lines):
                merged_elements_statspvalues_outfile.write(l.strip() + '\t' + str(p_values[i]) + '\t' + str(pvalues_adjusted[i]) + '\n')
                if pvalues_adjusted[i]<merged_mut_sig_threshold:
                    n_sig+=1
                    merged_elements_statspvaluesonlysig_outfile.write(l.strip() + '\t' + str(p_values[i]) + '\t' + str(pvalues_adjusted[i]) + '\n')
    
    
    return merged_elements_statspvalues, merged_elements_statspvaluesonlysig, n_sig




def get_tf_pval(cohort, sig_muts_per_tf_mutation_input_files, p_value_on_score, motif_name_index, 
                f_score_index, motif_breaking_score_index, 
                filter_cond, fsep, sig_tfs_file, sig_tfpos_file,
                filter_on_signal = True, dnase_index = 24, fantom_index = 25, 
                num_other_tfs_index = 27):
    
    print('sig_muts_per_tf_mutation_input_files: ', sig_muts_per_tf_mutation_input_files)
    observed_mut_motifs = sig_muts_per_tf_mutation_input_files[0]
    if os.path.isfile(sig_tfs_file) and os.path.isfile(sig_tfpos_file):
        return sig_tfs_file, sig_tfpos_file
    
    observed_mut_motifs_temp = observed_mut_motifs+'_sigTFs_mintfscores_temp'
    print('Calculating pval for TFs in ', cohort)
    '''filter the mutations by motif-breaking score and gene-expression as given in filter_cond
        the mutations in the input file are already checked for signicance (or TF binding>0)'''
    
    os.system("""awk 'BEGIN{{FS=OFS="{fsep}"}}{{{filter_cond}{{print ${motif_name_index},${f_score_index}+${mut_break_score_index}}}}}' {observed_mut_motifs} | sort -k1 | groupBy -g 1 -c 2 -o min,mean,stdev,median > {observed_mut_motifs_temp}""".format(
                        filter_cond=filter_cond, observed_mut_motifs=observed_mut_motifs, 
                        motif_name_index=motif_name_index+1, f_score_index=f_score_index+1, 
                        mut_break_score_index=motif_breaking_score_index+1, 
                        observed_mut_motifs_temp=observed_mut_motifs_temp, fsep=fsep))
    
    tf_min_scores_in_sig_obs_motifs = {}
    with open(observed_mut_motifs_temp, 'r') as i_observed_mut_motifs_temp:
        lines = i_observed_mut_motifs_temp.readlines()
        for l in lines:
            #min
            
            #Min
            #tf_min_scores_in_sig_obs_motifs[l.strip().split('\t')[0]] = float(l.strip().split('\t')[1])
            #mean
            tf_min_scores_in_sig_obs_motifs[l.strip().split('\t')[0]] = float(float(l.strip().split('\t')[2])-float(l.strip().split('\t')[3]))
    print('tf_min_scores_in_sig_obs_motifs: ', tf_min_scores_in_sig_obs_motifs) 
    gene_expression_index = 31
    tf_binding_index = 30
    mut_motif_pos_index = 13
    breaking_score_threshold = 0.3
    tf_counts_in_sim_sets = {}
    tfpos_counts_in_sim_sets = {}
    
    tf_counts_in_this_sim_set = {}
    tfpos_counts_in_this_sim_set = {}
    sim_file=sig_muts_per_tf_mutation_input_files[0]
    with open(sim_file) as i_sim_file:
         l = i_sim_file.readline().strip().split('\t')
         while l and len(l)>gene_expression_index:
             if ((l[gene_expression_index]=='nan' or float(l[gene_expression_index])>0) and 
                 float(l[motif_breaking_score_index])>=breaking_score_threshold):
                 try:
                     tf_counts_in_this_sim_set[l[motif_name_index]] +=1
                 except KeyError:
                     tf_counts_in_this_sim_set[l[motif_name_index]] = 1
                  
                 try:
                     tfpos_counts_in_this_sim_set[l[motif_name_index]+"#"+l[mut_motif_pos_index]] +=1
                 except KeyError:
                     tfpos_counts_in_this_sim_set[l[motif_name_index]+"#"+l[mut_motif_pos_index]] = 1
             l = i_sim_file.readline().strip().split('\t')
    for tf in tf_counts_in_this_sim_set.keys():
             try:
                 tf_counts_in_sim_sets[tf].append(tf_counts_in_this_sim_set[tf])
             except KeyError:
                 tf_counts_in_sim_sets[tf] = [tf_counts_in_this_sim_set[tf]]
          
    for tf in tfpos_counts_in_this_sim_set.keys():
         try:
             tfpos_counts_in_sim_sets[tf].append(tfpos_counts_in_this_sim_set[tf])
         except KeyError:
             tfpos_counts_in_sim_sets[tf] = [tfpos_counts_in_this_sim_set[tf]]


    for sim_file in sig_muts_per_tf_mutation_input_files[1:]: #count for all files incl. observed
        tf_counts_in_this_sim_set = {}
        tfpos_counts_in_this_sim_set = {}
        with open(sim_file) as i_sim_file:
            l = i_sim_file.readline().strip().split('\t')
            while l and len(l)>gene_expression_index:
                if ((l[gene_expression_index]=='nan' or float(l[gene_expression_index])>0) and 
                    float(l[motif_breaking_score_index])>=breaking_score_threshold):

                    
                    '''means no fdr and signal check has been applied therefore only keep motifs that have:
                        (f_score+breaking_score> minimum score obtained for the same motif in the obs set
                        (instead of FDR calculations to ensure only reasonable mutations are counted)
                        and there is dnase signal 
                        or there is tf binindg signal 
                        this is usually the case for sim sets. (the idea is to get mut similar to those in obs set)'''
                    try:
                        min_obs_score_this_motif = tf_min_scores_in_sig_obs_motifs[l[motif_name_index]]
                    except KeyError:
                        min_obs_score_this_motif = None
                    
                    if min_obs_score_this_motif:
                        if ((float(l[f_score_index])) >= min_obs_score_this_motif) :
                              if((float(l[dnase_index])>0.0) or# or float(l[fantom_index])>0.0 or float(l[num_other_tfs_index])>0.0 
                             (float(l[tf_binding_index])>0 and l[tf_binding_index]!="nan")):
                                try:
                                    tf_counts_in_this_sim_set[l[motif_name_index]] +=1
                                except KeyError:
                                    tf_counts_in_this_sim_set[l[motif_name_index]] = 1
                                
                                try:
                                    tfpos_counts_in_this_sim_set[l[motif_name_index]+"#"+l[mut_motif_pos_index]] +=1
                                except KeyError:
                                    tfpos_counts_in_this_sim_set[l[motif_name_index]+"#"+l[mut_motif_pos_index]] = 1
                    
                l = i_sim_file.readline().strip().split('\t')
        for tf in tf_counts_in_this_sim_set.keys():
            try:
                tf_counts_in_sim_sets[tf].append(tf_counts_in_this_sim_set[tf])
            except KeyError:
                tf_counts_in_sim_sets[tf] = [tf_counts_in_this_sim_set[tf]]
        
        for tf in tfpos_counts_in_this_sim_set.keys():
            try:
                tfpos_counts_in_sim_sets[tf].append(tfpos_counts_in_this_sim_set[tf])
            except KeyError:
                tfpos_counts_in_sim_sets[tf] = [tfpos_counts_in_this_sim_set[tf]]
    tf_p_values = []
    tf_names = []
    for tf in sorted(tf_counts_in_sim_sets.keys()):
        if len(tf_counts_in_sim_sets[tf])<len(sig_muts_per_tf_mutation_input_files):
            for i in range(len(tf_counts_in_sim_sets[tf]), len(sig_muts_per_tf_mutation_input_files)):
                tf_counts_in_sim_sets[tf].append(0)
        num_tf_obs =  tf_counts_in_sim_sets[tf][0]
        num_tf_sim = tf_counts_in_sim_sets[tf][1:]
        tf_names.append(tf)
        num_tf_sim_mean = np.mean(num_tf_sim)
        num_tf_sim_sd = np.std(num_tf_sim)
        if num_tf_sim_sd==0.0:
            num_tf_sim_sd = 1.0
        if p_value_on_score:
            tf_p_values.append(empirical_pval((num_tf_obs,), num_tf_sim))
            #adjusted_tf_p_values = []
            #adjusted_tf_p_values = tf_p_values
        else:
            tf_p_values.append(get_pval(num_tf_obs, num_tf_sim_mean, num_tf_sim_sd))
            adjusted_tf_p_values = []
            if len(tf_p_values)>0:
                adjusted_tf_p_values = adjust_pvales(tf_p_values)
            else:
                print('tf_p_values nothing:', tf_p_values)
                
    tf_p_values_vec = [j for i in tf_p_values for j in i]

    with open(sig_tfs_file, 'w') as ofile:
        for i,tf in enumerate(tf_names):
            ofile.write(tf + '\t' + str(tf_p_values_vec[i]) + '\t' + str(tf_p_values_vec[i]) + '\t' + str(tf_counts_in_sim_sets[tf][0]) + '\t' + str(np.mean(tf_counts_in_sim_sets[tf][1:])) + '\t' + ','.join([str(x) for x in tf_counts_in_sim_sets[tf][1:]])+ '\n')
    
    tfpos_p_values = []
    tfpos_names = []
    for tfpos in sorted(tfpos_counts_in_sim_sets.keys()):
        if len(tfpos_counts_in_sim_sets[tfpos])<len(sig_muts_per_tf_mutation_input_files):
            for i in range(len(tfpos_counts_in_sim_sets[tfpos]), len(sig_muts_per_tf_mutation_input_files)):
                tfpos_counts_in_sim_sets[tfpos].append(0)
        
        num_tfpos_obs =  tfpos_counts_in_sim_sets[tfpos][0]
        num_tfpos_sim = tfpos_counts_in_sim_sets[tfpos][1:]
        tfpos_names.append(tfpos)
        num_tfpos_sim_mean = np.mean(num_tfpos_sim)
        num_tfpos_sim_sd = np.std(num_tfpos_sim)
        if num_tfpos_sim_sd==0.0:
            num_tfpos_sim_sd = 1.0
        if p_value_on_score:
            tfpos_p_values.append(empirical_pval((num_tfpos_obs,), num_tfpos_sim))
        else:
            tfpos_p_values.append(get_pval(num_tfpos_obs, num_tfpos_sim_mean, num_tfpos_sim_sd))
    
            adjusted_tfpos_p_values = []
            if len(tfpos_p_values)>0:
                adjusted_tfpos_p_values = adjust_pvales(tfpos_p_values)
                
    tfpos_p_values_vec = [j for i in tfpos_p_values for j in i]
    with open(sig_tfpos_file, 'w') as ofile:
        for i,tfpos in enumerate(tfpos_names):
            ofile.write(tfpos + '\t' + str(tfpos_p_values_vec[i]) + '\t' + str(tfpos_p_values_vec[i]) + '\t' + str(tfpos_counts_in_sim_sets[tfpos][0]) + '\t' + str(np.mean(tfpos_counts_in_sim_sets[tfpos][1:])) + '\t' + ','.join([str(x) for x in tfpos_counts_in_sim_sets[tfpos][1:]])+ '\n')
    if os.path.exists(observed_mut_motifs_temp):
        os.remove(observed_mut_motifs_temp)
    return sig_tfs_file, sig_tfpos_file


def get_pval(element_score, avg, sd):
    try:
        z_score = (element_score - avg)/sd
    except ZeroDivisionError: #in case sd is zero
        z_score = (element_score - avg)
    p_value = stats.norm.sf(z_score)
    return p_value

def adjust_pvales(pvalues):
    significant_bool_report, corrected_p_values_array, alphacSidak, alphacBonf = multipletests(pvalues, alpha=0.05, method='fdr_bh', returnsorted=False) #returns 4 things: a boolean array contains True or False for each value meaning wether the value after correction compared to the given alpha is significant or not, an array of the values after correction, a single for corrected alpha for Sidak method, a single value for corrected alpha for Bonferroni method 
    return corrected_p_values_array.tolist()
    
def process_input_file(observed_input_file, simulated_input_files, 
                       combined_simulated_muts_output_file, 
                       combined_simulated_muts_merged_output_file, 
                       output_extension, distance_to_merge, 
                       filter_cond, mut_sig_threshold):
    
    if os.path.exists(combined_simulated_muts_merged_output_file) and os.path.exists(combined_simulated_muts_output_file):
        pass
    else:
        with open(combined_simulated_muts_output_file, 'w') as combined_simulated_muts_outfile, open(combined_simulated_muts_merged_output_file, 'w') as combined_simulated_muts_merged_oufile:
            for simulated_input_file in simulated_input_files:
                
                unified_muts_file = simulated_input_file + output_extension + "_groupedbymut" 
                print(unified_muts_file)
                unified_muts_file = unify_muts(simulated_input_file, unified_muts_file, filter_mut_motifs=True, filter_cond=filter_cond)
                with open(unified_muts_file, 'r') as unified_muts_readfile:
                    combined_simulated_muts_outfile.write(unified_muts_readfile.read())
                
                unified_muts_file_wihtmotifinfo = unified_muts_file+"withmotifinfo"
                print(unified_muts_file_wihtmotifinfo)
                unified_muts_file_wihtmotifinfo = get_max_motif_in_grouped_muts(annotated_mutations_grouped_file=unified_muts_file, annotated_mutations_grouped_output_file=unified_muts_file_wihtmotifinfo)
                
                calculated_pvalues_unified_muts_file_wihtmotifinfo = unified_muts_file_wihtmotifinfo+"_statmuts"
                sig_calculated_pvalues_unified_muts_file_wihtmotifinfo = unified_muts_file_wihtmotifinfo+"_statmutsonlysig"
                calculated_pvalues_unified_muts_file_wihtmotifinfo, sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, n_sig_muts = assess_stat_muts(muts_input_file=unified_muts_file_wihtmotifinfo, simulated_input_file=unified_muts_file_wihtmotifinfo, observed_output_file=calculated_pvalues_unified_muts_file_wihtmotifinfo, observed_onlysig_output_file=sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, score_index_observed_elements=9, score_index_sim_elements=9, mut_sig_threshold=mut_sig_threshold)
                
                merged_muts_output_file = sig_calculated_pvalues_unified_muts_file_wihtmotifinfo+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)
                merged_muts_output_file = merge_muts(muts_input_file=sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, merged_muts_output_file=merged_muts_output_file, filter_mut_motifs=False, filter_col_index=15, filter_value=0.05, mut_score_index=9, motifs_col_index =10, ref_alt_col_index=11, mutpos_col_index=12, motifname_col_index=13, motif_col_index=14, distance_to_merge=20)
                with open(merged_muts_output_file, 'r') as merged_muts_read_file:
                    combined_simulated_muts_merged_oufile.write(merged_muts_read_file.read())
                
    print(combined_simulated_muts_output_file)
    print(combined_simulated_muts_merged_output_file)
    
    unified_muts_file = observed_input_file + output_extension + "_groupedbymut" 
    unified_muts_file = unify_muts(observed_input_file, unified_muts_file, filter_mut_motifs=True, filter_cond=filter_cond)
    
    unified_muts_file_wihtmotifinfo = unified_muts_file+"withmotifinfo"
    unified_muts_file_wihtmotifinfo = get_max_motif_in_grouped_muts(annotated_mutations_grouped_file=unified_muts_file, annotated_mutations_grouped_output_file=unified_muts_file_wihtmotifinfo)
    
    calculated_pvalues_unified_muts_file_wihtmotifinfo = unified_muts_file_wihtmotifinfo+"_statmuts"
    sig_calculated_pvalues_unified_muts_file_wihtmotifinfo = unified_muts_file_wihtmotifinfo+"_statmutsonlysig"
    calculated_pvalues_unified_muts_file_wihtmotifinfo, sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, n_sig_muts = assess_stat_muts(muts_input_file=unified_muts_file_wihtmotifinfo, simulated_input_file=combined_simulated_muts_output_file, observed_output_file=calculated_pvalues_unified_muts_file_wihtmotifinfo, observed_onlysig_output_file=sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, score_index_observed_elements=9, score_index_sim_elements=9, mut_sig_threshold=mut_sig_threshold)
    merged_muts_output_file = sig_calculated_pvalues_unified_muts_file_wihtmotifinfo+"_mergedmuts{distance_to_merge}bp".format(distance_to_merge=distance_to_merge)
    merged_muts_output_file = merge_muts(muts_input_file=sig_calculated_pvalues_unified_muts_file_wihtmotifinfo, merged_muts_output_file=merged_muts_output_file, filter_mut_motifs=False, filter_col_index=15, filter_value=0.05, mut_score_index=9, motifs_col_index =10, ref_alt_col_index=11, mutpos_col_index=12, motifname_col_index=13, motif_col_index=14, distance_to_merge=20)
    merged_elements_statspvalues = merged_muts_output_file+"_statspvalues"
    merged_elements_statspvaluesonlysig = merged_muts_output_file+"_statspvaluesonlysig"
    merged_elements_statspvalues, merged_elements_statspvaluesonlysig, n_sig = assess_stat_elements(observed_input_file=merged_muts_output_file, simulated_input_file=combined_simulated_muts_merged_output_file, merged_elements_statspvalues=merged_elements_statspvalues, merged_elements_statspvaluesonlysig=merged_elements_statspvaluesonlysig, score_index_observed_elements=3, score_index_sim_elements=3)
    print('Number of Sig elements: ', n_sig, 'initial #sig muts ', n_sig_muts)
    
    return merged_elements_statspvaluesonlysig

def get_scores_per_window(observed_input_files_objs, observed_input_file, tmp_dir, window_size, ext, simulated_input_file):
    
    
    simulated_input_file_tmp_overallTFs_local = tmp_dir +'/'+ observed_input_file.split('/')[-1] + '_' + simulated_input_file.split('/')[-1] + ext
    
    if os.path.exists(simulated_input_file_tmp_overallTFs_local):
        return  simulated_input_file_tmp_overallTFs_local
    
    
    
    #"remove chr, X>23,Y>24 and print in string format, check position if number"
    #simulated_input_file_fixed_sorted = tmp_dir + '/' + observed_input_file.split('/')[-1] + simulated_input_file.split('/')[-1]   
        
    #awk_stmt = r"""awk 'BEGIN{{FS=OFS="\t"}}{{gsub("X","23", $1); gsub("Y","24", $1); gsub("chr","", $1);if($10==".") print $1, $2*1, $3*1, $11, $18; else print $1, $2*1, $3*1, $10+$11, $18}}' {simulated_file} | sort -k1,1n -k2,2n -k3,3n > {simulated_outfile_temp}""".format(simulated_file = simulated_input_file, simulated_outfile_temp = simulated_input_file_fixed_sorted)

    #os.system(awk_stmt)
    
    #sim_chrs_dir = simulated_input_file_fixed_sorted+'_sim_chrs/'
    #if not os.path.isdir(sim_chrs_dir):
    #    os.makedirs(sim_chrs_dir)
    
    #os.system("""awk '{{print $0>>"{}"$1".bed"}}' {}""".format(
    #    sim_chrs_dir, simulated_input_file_fixed_sorted))
    
    simulated_input_file_tmp_overallTFs_local_temp = simulated_input_file_tmp_overallTFs_local + '_temp'
        
    sim_chr = simulated_input_file.split('/')[-1].split('.')[0]
    #print(sim_chr)
    #print(observed_input_files_objs[sim_chr])
    "if observed mutation file for any chromosome doesn't exist, return en ampty file"
    try: 
        obs_chr_obj =  BedTool(observed_input_files_objs[sim_chr])
    except KeyError:
        open(simulated_input_file_tmp_overallTFs_local, 'a').close()
        return(simulated_input_file_tmp_overallTFs_local)
        
    sim_chr_obj = BedTool(simulated_input_file)
    print("Intersecting ", simulated_input_file)
    sim_chr_file_intersected = simulated_input_file+ '_intersected'
    obs_chr_obj.map(sim_chr_obj, c=4, o=['mean', 'stdev', 'count']).saveas(sim_chr_file_intersected)
    #obs_chr_obj.map(sim_chr_obj, c=4, o=['collapse']).saveas(sim_chr_file_intersected)
    
    #window_id_fscroe_file = """awk 'BEGIN{{FS=OFS="\t"}}{{if($5!=".") print $4,$5}}' {sim_intersected} >> {sim_scores_combined}""".format(
    #            sim_intersected=sim_chr_file_intersected, sim_scores_combined=simulated_input_file_tmp_overallTFs_local_temp)
    window_id_fscroe_file = """awk 'BEGIN{{FS=OFS="\t"}}{{if($5!=".") print $4,$5,$6,$7, $8}}' {sim_intersected} >> {sim_scores_combined}""".format(
                sim_intersected=sim_chr_file_intersected, sim_scores_combined=simulated_input_file_tmp_overallTFs_local_temp)
    os.system(window_id_fscroe_file)
    
    
#     for chr_file in os.listdir(sim_chrs_dir):
#         if chr_file.endswith('.bed'):
#             sim_chr_file = sim_chrs_dir+chr_file
#             obs_chr_obj = observed_input_files_objs[chr_file.replace('.bed', '')]
#             sim_chr_obj = BedTool(sim_chr_file)
#             
#             print("Intersecting ", sim_chr_file)
#             sim_chr_file_intersected = sim_chr_file+'_intersected'
#             obs_chr_obj.map(sim_chr_obj, c=4, o=['mean', 'stdev', 'count']).saveas(sim_chr_file_intersected)
#             #obs_chr_obj.window(sim_chr_obj, w = window_size).saveas(sim_chr_file_intersected)
#             
#             #col 4: windowID; col18: tf-binding score; col9:fscore 
#             window_id_fscroe_file = """awk 'BEGIN{{FS=OFS="\t"}}{{if($7!=0) print $4,$5,$6,$7}}' {sim_intersected} >> {sim_scores_combined}""".format(
#                 sim_intersected=sim_chr_file_intersected, sim_scores_combined=simulated_input_file_tmp_overallTFs_local_temp)
#             os.system(window_id_fscroe_file)
#             #os.remove(sim_chr_file_intersected)
#             #os.remove(sim_chr_file)
    
    if os.path.isfile(simulated_input_file_tmp_overallTFs_local_temp):
        shutil.move(simulated_input_file_tmp_overallTFs_local_temp, simulated_input_file_tmp_overallTFs_local)
    
    #os.remove(simulated_input_file_tmp_overallTFs_local_temp)
    #shutil.rmtree(sim_chrs_dir)
    
    
    #print('cleanup')
    cleanup()
    return simulated_input_file_tmp_overallTFs_local

# def groupby_per_mut(score_input_files_objs, tmp_dir):
#     
#     score_group_file = 
#     d = score_input_files_objs.groupby(g=[1], c=2, o=['mean', 'stdev','count'])
# 
#     
#     
#     return


def get_simulated_mean_sd_per_TF_motif_background_window(cohort_full_name, annotated_input_file, 
                                                         simulated_annotated_input_files, 
                                       mutations_cohorts_dir,
                                       cohort_mean_sd_per_tf_overall_output_dict_file, 
                                       chr_lengths_file,
                                       background_window_size = 50000,
                                       motif_name_index = 17, f_score_index = 9, 
                                       motif_breaking_score_index = 10, chromatin_cat_index=22, tmp_dir = '$SNIC_TMP', n_cores_fscore=10):
    
    if os.path.exists(cohort_mean_sd_per_tf_overall_output_dict_file):
        with open(cohort_mean_sd_per_tf_overall_output_dict_file, 'r') as dict_simulated_mean_sd_per_TF_motif_ifile:
            dict_type_mean_std_scores = json.loads(dict_simulated_mean_sd_per_TF_motif_ifile.readline())
        return  dict_type_mean_std_scores
    
    print("Extracting avg and std per TF and overall from the simulation sets... onto: ", cohort_mean_sd_per_tf_overall_output_dict_file)
    cohort = cohort_full_name.split('/')[-1]
    
    "replace chr, X, Y, add Line Number to use as window ID and sort by chr,start"
    observed_input_file_sorted = tmp_dir+'/'+ annotated_input_file.split('/')[-1] + '_fixed_sorted'
    cmd = """awk 'BEGIN{{OFS="\t"}}{{gsub("chr","",$1); gsub("X", 23, $1); gsub("Y", 24, $1); print $1,$2,$3,NR}}' {} | bedtools slop -g /proj/snic2020-16-50/nobackup/pancananalysis/pancan12Feb2020/cancer_datafiles/chr_order_hg19.txt -b {}| sort -k1,1n -k2,2n > {}""".format(
        annotated_input_file, background_window_size, observed_input_file_sorted)
    os.system(cmd)
    
    obs_chrs_dir = tmp_dir+'/'+cohort + '_chrs/'
    if not os.path.isdir(obs_chrs_dir):
        os.makedirs(obs_chrs_dir)
    
    observed_input_files_objs = {}    
    os.system("""awk '{{print $0>>"{}"$1".bed"}}' {}""".format(
        obs_chrs_dir, observed_input_file_sorted))
    
    for chr_file in os.listdir(obs_chrs_dir):
        
        if chr_file.endswith('.bed'):
            #chr_file_window=obs_chrs_dir+chr_file+'_window'
            #chr_file_window_sorted = chr_file_window +'_sorted'
            #BedTool(obs_chrs_dir+chr_file).slop(b=background_window_size,genome='hg19').saveas(chr_file_window)
            #os.system("""sort -k1,1n -k2,2n {} > {}""".format(
        #chr_file_window, chr_file_window_sorted))
            observed_input_files_objs[chr_file.replace('.bed', '')] = obs_chrs_dir+chr_file
    
    sim_chrs_dir = tmp_dir + '/'+ cohort + '_sim/'
    if not os.path.isdir(sim_chrs_dir):
        os.makedirs(sim_chrs_dir)
    
    "combine simulated files with the same chr"
    for sim_file in simulated_annotated_input_files:
        sim_file_tmp = sim_file +'_tmp'
        os.system("""awk 'BEGIN{{FS=OFS="\t"}}{{gsub("X","23", $1); gsub("Y","24", $1); gsub("chr","", $1);if($10==".") print $1, $2*1, $3*1, $11, $18; else print $1, $2*1, $3*1, $10+$11, $18}}' {} > {}""".format(
            sim_file, sim_file_tmp)) 
        os.system("""awk '{{print $0>>"{}"$1".bed"}}' {}""".format(
        sim_chrs_dir, sim_file_tmp))
        os.remove(sim_file_tmp)
        
    sim_input_files =[]
    for sim_file in os.listdir(sim_chrs_dir):
        sim_file_sorted = sim_chrs_dir +'/' + sim_file +'sorted'
        os.system("""sort -k1,1n -k2,2n -k3,3n {} > {}""".format(
        sim_chrs_dir+ sim_file, sim_file_sorted))
        sim_input_files.append(sim_file_sorted)
        #split files by 1 000 000 bps
        #os.system("""awk 'BEGIN{{n=1}}{{x=$3;if(x>n*10000000){{++n}}{{print > "{sim_file_sorted}""_split_"n}}}}' {sim_file_sorted}""".format(
        #     sim_file_sorted=sim_file_sorted))

    #sim_input_files = [sim_chrs_dir+'/'+x for x in os.listdir(sim_chrs_dir) if 'sorted_split_' in x]

    #sim_input_files.append(sim_file_sorted)
    #print(sim_input_files)
        
          #awk_stmt = r"""awk 'BEGIN{{FS=OFS="\t"}}{{gsub("X","23", $1); gsub("Y","24", $1); gsub("chr","", $1);if($10==".") print $1, $2*1, $3*1, $11, $18; else print $1, $2*1, $3*1, $10+$11, $18}}' {simulated_file} | sort -k1,1n -k2,2n -k3,3n > {simulated_outfile_temp}""".format(simulated_file = simulated_input_file, simulated_outfile_temp = simulated_input_file_fixed_sorted)

    #awk_stmt = r"""awk 'BEGIN{{FS=OFS="\t"}}{{gsub("X","23", $1); gsub("Y","24", $1); gsub("chr","", $1); print $1, $2*1, $3*1, $10, $11, $18}}' {simulated_file} | sort -k1,1n -k2,2n > {simulated_outfile_temp}""".format(simulated_file = simulated_input_file, simulated_outfile_temp = simulated_input_file_fixed_sorted)
    #awk_stmt = r"""awk 'BEGIN{{FS=OFS="\t"}}{{gsub("X","23", $1); gsub("Y","24", $1); gsub("chr","", $1); printf ("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2*1, $3*1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32)}}' {simulated_file} > {simulated_outfile_temp}""".format(simulated_file = simulated_input_file, simulated_outfile_temp = simulated_input_file_fixed)
    #os.system(awk_stmt)  
        
    ext = '_scoresPerWindow'
    obs_scores_files = []
    if n_cores_fscore>1:
        p = Pool(n_cores_fscore)
        obs_scores_files = p.starmap(get_scores_per_window, product(
            [observed_input_files_objs],[observed_input_file_sorted], [tmp_dir], [background_window_size],[ext], 
            sim_input_files))
        p.close()
        p.join()
    else:
        for simulated_annotated_input_file in sim_input_files:
            obs_scores_files.append(get_scores_per_window(observed_input_files_objs, observed_input_file_sorted, tmp_dir, 
                                  background_window_size, ext, simulated_annotated_input_file))
    print(obs_scores_files)
    
    print('Combining the scores')
    


    simulated_mean_sd_outfiles = tmp_dir + '/' + cohort + '_allscores'
     
    if not os.path.isfile(simulated_mean_sd_outfiles):
        #merge files from the same category, sort by the line number and group by position, TF motif, chromatin cat. and line number
        with open(simulated_mean_sd_outfiles, 'w') as sim_fn:
            for obs_scores_file in obs_scores_files:
                with open(obs_scores_file, 'r') as score_fn:
                    sim_fn.write(score_fn.read())
     
    "create a dictionery for mean, std scores for all categories"         
    dict_simulated_mean_sd = {} 
    with open(simulated_mean_sd_outfiles, 'r') as simulated_mean_sd_ifile:        
        l = simulated_mean_sd_ifile.readline().strip().split('\t')            
        while l and len(l)>1:    
            #fscore = map(float, l[1].split(','))          
            #fscore = [float(x) for x in l[1].split(',')]
            #dict_simulated_mean_sd[l[0]] = {'mean': np.mean(fscore),          
            #                                       "std": np.std(fscore),                                                        
            #                                       "nummotifs": len(fscore)}
            dict_simulated_mean_sd[l[0]] = {'mean': l[1],          
                                                  "std": l[2],                                                        
                                                   "nummotifs": l[3]}
            l = simulated_mean_sd_ifile.readline().strip().split('\t')
     
    #save the dictionery per category
    dict_type_mean_std_scores = {}
    dict_type_mean_std_scores['overallTFs'] = dict_simulated_mean_sd
    
    with open(cohort_mean_sd_per_tf_overall_output_dict_file, 'w') as dict_simulated_mean_sd_per_TF_motif_outfile:
            json.dump(dict_type_mean_std_scores, dict_simulated_mean_sd_per_TF_motif_outfile)
    
    shutil.rmtree(sim_chrs_dir)
    shutil.rmtree(obs_chrs_dir)   
    cleanup() 
    
    return  dict_type_mean_std_scores


def get_simulated_mean_sd_per_TF_motif(simulated_annotated_input_files, 
                                       cohort_mean_sd_per_tf_overall_output_dict_file, 
                                       obs_chrs_dir, observed_input_files_objs,
                                       motif_name_index = 17, f_score_index = 9, 
                                       motif_breaking_score_index = 10, chromatin_cat_index=22):
    
    if os.path.exists(cohort_mean_sd_per_tf_overall_output_dict_file):
        with open(cohort_mean_sd_per_tf_overall_output_dict_file, 'r') as dict_simulated_mean_sd_per_TF_motif_ifile:
            dict_simulated_mean_sd_per_TF_motif = json.loads(dict_simulated_mean_sd_per_TF_motif_ifile.readline())
        return dict_simulated_mean_sd_per_TF_motif
    print("Extracting avg and std per TF and overall from the simulation sets... onto: ", cohort_mean_sd_per_tf_overall_output_dict_file)
    overall_score_values = []
    
    dict_simulated_score_per_TF_motif = {}
    dict_simulated_mean_sd_per_TF_motif = {}
    
    dict_simulated_score_per_chromatin_cat = {}
    dict_simulated_mean_sd_per_chromatin_cat = {}
    
    dict_simulated_score_per_TF_motif_per_chromatin_cat = {}
    dict_simulated_mean_sd_per_TF_motif_per_chromatin_cat = {}
    
    for simulated_annotated_input_file in simulated_annotated_input_files:
        with open(simulated_annotated_input_file, 'r') as simulated_annotated_ifile:
            l = simulated_annotated_ifile.readline().strip().split('\t')
            while l and len(l)>motif_name_index:
                try:
                    dict_simulated_score_per_TF_motif[l[motif_name_index]].append(float(l[f_score_index])+float(l[motif_breaking_score_index]))
                except KeyError:
                    dict_simulated_score_per_TF_motif[l[motif_name_index]] = [float(l[f_score_index])+float(l[motif_breaking_score_index])]
                
                try:
                    dict_simulated_score_per_chromatin_cat[l[chromatin_cat_index]].append(float(l[f_score_index])+float(l[motif_breaking_score_index]))
                except KeyError:
                    dict_simulated_score_per_chromatin_cat[l[chromatin_cat_index]] = [float(l[f_score_index])+float(l[motif_breaking_score_index])]
                
                try:
                    dict_simulated_score_per_TF_motif_per_chromatin_cat[l[motif_name_index]][l[chromatin_cat_index]].append(float(l[f_score_index])+float(l[motif_breaking_score_index]))
                except KeyError:
                    try:
                        dict_simulated_score_per_TF_motif_per_chromatin_cat[l[motif_name_index]][l[chromatin_cat_index]] = [float(l[f_score_index])+float(l[motif_breaking_score_index])]
                    except KeyError:
                        dict_simulated_score_per_TF_motif_per_chromatin_cat[l[motif_name_index]] = {l[chromatin_cat_index] : [float(l[f_score_index])+float(l[motif_breaking_score_index])]}
                        
                overall_score_values.append(float(l[f_score_index])+float(l[motif_breaking_score_index]))
                l = simulated_annotated_ifile.readline().strip().split('\t')
    
    #get mean and std of scores per TF_motif
    for tf in dict_simulated_score_per_TF_motif.keys():
        tf_mean = np.mean(dict_simulated_score_per_TF_motif[tf])
        tf_std = np.std(dict_simulated_score_per_TF_motif[tf])
        num_motifs = len(dict_simulated_score_per_TF_motif[tf])
        dict_simulated_mean_sd_per_TF_motif[tf] = {'mean': tf_mean, 
                                                   "std": tf_std, 
                                                   "nummotifs": num_motifs}
    
    #get mean and std of scores per chromatin category
    for chromatin_cat in dict_simulated_score_per_chromatin_cat.keys():
        chromatin_cat_mean = np.mean(dict_simulated_score_per_chromatin_cat[chromatin_cat])
        chromatin_cat_std = np.std(dict_simulated_score_per_chromatin_cat[chromatin_cat])
        num_motifs = len(dict_simulated_score_per_chromatin_cat[chromatin_cat])
        dict_simulated_mean_sd_per_chromatin_cat[chromatin_cat] = {'mean': chromatin_cat_mean, 
                                                                   "std": chromatin_cat_std, 
                                                                   "nummotifs": num_motifs} 
    
    #get mean and std of scores per TF_motif per chromatin category
    for tf in dict_simulated_score_per_TF_motif_per_chromatin_cat.keys():
        for chromatin_cat in dict_simulated_score_per_TF_motif_per_chromatin_cat[tf].keys():
            tf_chromatin_cat_mean = np.mean(dict_simulated_score_per_TF_motif_per_chromatin_cat[tf][chromatin_cat])
            tf_chromatin_cat_std = np.std(dict_simulated_score_per_TF_motif_per_chromatin_cat[tf][chromatin_cat])
            num_motifs = len(dict_simulated_score_per_TF_motif_per_chromatin_cat[tf][chromatin_cat])
            try:
                dict_simulated_mean_sd_per_TF_motif_per_chromatin_cat[tf][chromatin_cat] = {'mean': tf_chromatin_cat_mean, 
                                                   "std": tf_chromatin_cat_std, 
                                                   "nummotifs": num_motifs}
            except KeyError:
                dict_simulated_mean_sd_per_TF_motif_per_chromatin_cat[tf] = {chromatin_cat: {'mean': tf_chromatin_cat_mean, 
                                                   "std": tf_chromatin_cat_std, 
                                                   "nummotifs": num_motifs}}
    
    #get mean and std of scores across the genome regardless motif or chromatin category
    overall_num_motifs = len(overall_score_values)
    number_subsets = 1
    subset_size = 10000000
    if overall_num_motifs > subset_size: 
        number_subsets = int(overall_num_motifs / subset_size) 
    
    subsets = np.array_split(overall_score_values, number_subsets)
    means = []
    stds = []
    for subset in subsets:
        means.append(np.mean(subset))
        stds.append(np.std(subset))
    overallTFs_mean = np.mean(means)
    overallTFs_std = np.mean(stds)
    #overallTFs_mean = np.mean(overall_score_values)
    #overallTFs_std = np.std(overall_score_values)
    #overall_num_motifs = len(overall_score_values)
    num_tfs = len(dict_simulated_score_per_TF_motif.keys())
    dict_type_mean_std_scores = {}
    dict_type_mean_std_scores['overallTFs'] = {'mean': overallTFs_mean, 'std': overallTFs_std, 'nummotifs_total': overall_num_motifs, 'nummotifs_avg': overall_num_motifs/num_tfs}
    dict_type_mean_std_scores['perTF'] = dict_simulated_mean_sd_per_TF_motif
    dict_type_mean_std_scores['perChromatinCat'] = dict_simulated_mean_sd_per_chromatin_cat
    dict_type_mean_std_scores['perTF_perChromatinCat'] = dict_simulated_mean_sd_per_TF_motif_per_chromatin_cat
    with open(cohort_mean_sd_per_tf_overall_output_dict_file, 'w') as dict_simulated_mean_sd_per_TF_motif_outfile:
        json.dump(dict_type_mean_std_scores, dict_simulated_mean_sd_per_TF_motif_outfile)
    
    return dict_type_mean_std_scores

def get_muts_sig_per_TF(annoted_input_file, dict_type_mean_std_scores, 
                        annoted_output_file_extension, annoted_output_file_extension_onlysig, 
                        background_window = False,
                        motif_name_index = 17, f_score_index = 9, chromatin_index = 22,
                        motif_breaking_score_index = 10,
                        filter_on_qval=True, sig_cat='overallTFs',
                        sig_thresh=0.05,
                        filter_on_signal = True, dnase_index = 24, fantom_index = 25, 
                        num_other_tfs_index = 27, tf_binding_index=30):
    
    if sig_thresh>=1.0:
        return annoted_input_file
    annoted_output_file = annoted_input_file + annoted_output_file_extension
    annoted_output_file_onlysig = annoted_input_file + annoted_output_file_extension_onlysig
    
    if os.path.exists(annoted_output_file_onlysig):
        return annoted_output_file_onlysig
    
    if filter_on_qval:
        print("Calculating Q-values using {} for set: {}".format(sig_cat, annoted_input_file))
    else:
        print("Calculating P-values using {} for set: {}".format(sig_cat, annoted_input_file))
    
    if os.path.exists(annoted_output_file):
        with open(annoted_output_file, 'r') as annoted_output_ifile, open(annoted_output_file_onlysig, 'w') as annoted_input_ofile_onlysig:
            l = annoted_output_ifile.readline()
            while l:
                sl = l.strip().split('\t')
                if len(sl) >= motif_breaking_score_index:
                    
                    if sl[tf_binding_index]!="nan":
                        if float(sl[tf_binding_index]) > 0:
                            annoted_input_ofile_onlysig.write(l)
                            l = annoted_output_ifile.readline()
                            continue
                    pvalues = dict(json.loads(sl[motif_breaking_score_index+1].split('@')[0]).replace(';', ','))
                    adj_pvalues = dict(json.loads(sl[motif_breaking_score_index+1].split('@')[1].replace(';', ',')))
                    sig_level = 1.0
                    if background_window:
                        if filter_on_qval:
                            sig_level = float(adj_pvalues)
                        else:
                           sig_level = float(pvalues)
                    else:
                        if filter_on_qval:
                            sig_level = float(adj_pvalues[sig_cat])
                        else:
                            sig_level = float(pvalues[sig_cat])
                    
                        
                    if filter_on_signal:
                        if (sig_level < sig_thresh and 
                            (float(sl[dnase_index])>0.0)):# or float(sl[fantom_index])>0.0 or float(sl[num_other_tfs_index])>0.0
                            annoted_input_ofile_onlysig.write(l)
                    else:
                        if sig_level < sig_thresh:
                            annoted_input_ofile_onlysig.write(l)
                    
                l = annoted_output_ifile.readline()
        return annoted_output_file_onlysig
    
    dict_pvals = [] #store p-values for seach line per category in the same order as the  lines
    dict_line_indices = {} #store index of the lines from the input files to keep track of the categories
    lines = []
    
    
    if background_window:
        with open(annoted_input_file, 'r') as observed_annoted_input_ifile:
            lines = observed_annoted_input_ifile.readlines()
            print("Computing P-values (TF motifs) for: ", annoted_input_file)
            for line_index, line in enumerate(lines):
                l = line.strip().split('\t')
                if len(l)<motif_name_index:
                    continue
                p_value = 1.0
                try: 
                    #check if the background exist
                    avg = float(dict_type_mean_std_scores[sig_cat][str(line_index)]['mean'])
                    sd = float(dict_type_mean_std_scores[sig_cat][str(line_index)]['std'])
                    p_value = get_pval(float(l[f_score_index]) + float(l[motif_breaking_score_index]), 
                                         avg=avg, 
                                         sd=sd)
                except KeyError:
                    #no simulated mutations in the background to compare; set p-value as 1
                    p_value = 0.0
                try:
                    dict_pvals.append(p_value)
                    dict_line_indices = (line_index)
                except KeyError:
                    dict_pvals= [float(p_value)]
                    dict_line_indices = [line_index]
            
            print("Computing adjusted P-values for {}".format(annoted_input_file))
            adjusted_dict_pvals = {} 
            adjusted_dict_pvals = adjust_pvales(dict_pvals)
            
            with open(annoted_output_file, 'w') as annoted_input_ofile, open(annoted_output_file_onlysig, 'w') as annoted_input_ofile_onlysig:
            
                for line_index, line in enumerate(lines):
                    sl = line.strip().split('\t')    
                    #if dict_p_values_per_tf[tf][i]<0.05 or dict_p_values_overall_per_tf[tf][i]<0.05:
                    #iterate through the pvalues list and report them all for this line
                    pvals = {}
                    adjust_pvals = {}
                    pvals = dict_pvals[line_index]
                    adjust_pvals = adjusted_dict_pvals[line_index]
                    sl[motif_breaking_score_index+1] = (json.dumps({sig_cat: pvals})+'@'+json.dumps({sig_cat: adjust_pvals})).replace(',', ';')
                    
                    annoted_input_ofile.write('\t'.join(sl) + '\n')
                    sig_level = 1.0
                    if filter_on_qval:
                        sig_level = adjust_pvals
                    else:
                        sig_level = pvals
                    #select the motif if its ChiP-seq signal larger than 0
                    if sl[tf_binding_index]!="nan":
                        if float(sl[tf_binding_index]) > 0.0:
                            annoted_input_ofile_onlysig.write('\t'.join(sl) + '\n')
                            continue
                    
                    if filter_on_signal:
                        if (sig_level<sig_thresh and 
                            float(sl[dnase_index])>0.0):# or float(tf_motif[fantom_index])>0.0 or float(tf_motif[num_other_tfs_index])>0.0
                            annoted_input_ofile_onlysig.write('\t'.join(sl) + '\n')
                    else:
                        if sig_level<sig_thresh: #or adjusted_dict_p_values_per_tf[tf][i]<0.05:
                            annoted_input_ofile_onlysig.write('\t'.join(sl) + '\n')
                
    
    else:
        dict_pvals = {} #store p-values for seach line per category in the same order as the  lines
        dict_line_indices = {} #store index of the lines from the input files to keep track of the categories
        lines = []
        with open(annoted_input_file, 'r') as observed_annoted_input_ifile:
            lines = observed_annoted_input_ifile.readlines()
            print("Computing P-values (TF motifs) for: ", annoted_input_file)
            for line_index, line in enumerate(lines):
                l = line.strip().split('\t')
                if len(l)<motif_name_index:
                    continue
                p_value = 1.0
                pval_type = "overallTFs"
                try:
                    p_value = get_pval(float(l[f_score_index]) + float(l[motif_breaking_score_index]), 
                                     avg=dict_type_mean_std_scores[pval_type]['mean'], 
                                     sd=dict_type_mean_std_scores[pval_type]['std'])
                except KeyError:
                    p_value = 0.0
    
                try:
                    dict_pvals[pval_type].append(p_value)
                    dict_line_indices[pval_type].append(line_index)
                except KeyError:
                    dict_pvals[pval_type] = [p_value]
                    dict_line_indices[pval_type] = [line_index]
    
                pval_type = "perTF"
                try:
                    p_value = get_pval(float(l[f_score_index]) + float(l[motif_breaking_score_index]), 
                                 avg=dict_type_mean_std_scores[pval_type][l[motif_name_index]]['mean'], 
                                 sd=dict_type_mean_std_scores[pval_type][l[motif_name_index]]['std'])
                except KeyError:
                    p_value = 0.0
    
                try:
                    dict_pvals[pval_type][l[motif_name_index]].append(p_value)
                    dict_line_indices[pval_type][l[motif_name_index]].append(line_index)
                except KeyError:
                    try:
                        dict_pvals[pval_type][l[motif_name_index]] = [p_value]
                        dict_line_indices[pval_type][l[motif_name_index]] = [line_index]
                    except KeyError:
                        dict_pvals[pval_type] = {l[motif_name_index] : [p_value]}
                        dict_line_indices[pval_type] = {l[motif_name_index] :[line_index]}
    
    
                pval_type = "perChromatinCat"
                try:
                    p_value = get_pval(float(l[f_score_index]) + float(l[motif_breaking_score_index]), 
                                 avg=dict_type_mean_std_scores[pval_type][l[chromatin_index]]['mean'], 
                                 sd=dict_type_mean_std_scores[pval_type][l[chromatin_index]]['std'])
                except KeyError:
                    p_value = 0.0
    
                try:
                    dict_pvals[pval_type][l[chromatin_index]].append(p_value)
                    dict_line_indices[pval_type][l[chromatin_index]].append(line_index)
                except KeyError:
                    try:
                        dict_pvals[pval_type][l[chromatin_index]] = [p_value]
                        dict_line_indices[pval_type][l[chromatin_index]] = [line_index]
                    except KeyError:
                        dict_pvals[pval_type] = {l[chromatin_index] : [p_value]}
                        dict_line_indices[pval_type] = {l[chromatin_index] :[line_index]}
    
    
                pval_type = "perTF_perChromatinCat"
                try:
                    p_value = get_pval(float(l[f_score_index]) + float(l[motif_breaking_score_index]), 
                                 avg=dict_type_mean_std_scores[pval_type][l[motif_name_index]][l[chromatin_index]]['mean'], 
                                 sd=dict_type_mean_std_scores[pval_type][l[motif_name_index]][l[chromatin_index]]['std'])
                except KeyError:
                    p_value = 0.0
    
                try:
                    dict_pvals[pval_type][l[motif_name_index]][l[chromatin_index]].append(p_value)
                    dict_line_indices[pval_type][l[motif_name_index]][l[chromatin_index]].append(line_index)
                except KeyError:
                    try:
                        dict_pvals[pval_type][l[motif_name_index]][l[chromatin_index]] = [p_value]
                        dict_line_indices[pval_type][l[motif_name_index]][l[chromatin_index]] = [line_index]
                    except KeyError:
                        try:
                            dict_pvals[pval_type][l[motif_name_index]] = {l[chromatin_index] : [p_value]}
                            dict_line_indices[pval_type][l[motif_name_index]] = {l[chromatin_index] :[line_index]}
                        except KeyError:
                            dict_pvals[pval_type] = {l[motif_name_index]: {l[chromatin_index] : [p_value]}}
                            dict_line_indices[pval_type] = {l[motif_name_index]: {l[chromatin_index] :[line_index]}}
    
            print("Computing adjusted P-values for {}".format(annoted_input_file))
            adjusted_dict_pvals = {} 
       
    
            adjusted_dict_pvals["overallTFs"] = adjust_pvales(dict_pvals["overallTFs"])
            pval_type = "perTF"
            adjusted_dict_pvals[pval_type] = {}
            for tf in dict_pvals[pval_type].keys():
                adjusted_dict_pvals[pval_type][tf] = adjust_pvales(dict_pvals[pval_type][tf])
            pval_type = "perChromatinCat"
            adjusted_dict_pvals[pval_type] = {}
            for chrom_cat in dict_pvals[pval_type].keys():
                adjusted_dict_pvals[pval_type][chrom_cat] = adjust_pvales(dict_pvals[pval_type][chrom_cat])
            pval_type = "perTF_perChromatinCat"
            adjusted_dict_pvals[pval_type] = {}
            for tf in dict_pvals[pval_type].keys():
                adjusted_dict_pvals[pval_type][tf] = {}
                for chrom_cat in dict_pvals[pval_type][tf].keys():
                    adjusted_dict_pvals[pval_type][tf][chrom_cat] = adjust_pvales(dict_pvals[pval_type][tf][chrom_cat])
            
            with open(annoted_output_file, 'w') as annoted_input_ofile, open(annoted_output_file_onlysig, 'w') as annoted_input_ofile_onlysig:
                
                for line_index, line in enumerate(lines):
                    sl = line.strip().split('\t')    
                    #if dict_p_values_per_tf[tf][i]<0.05 or dict_p_values_overall_per_tf[tf][i]<0.05:
                    #iterate through the pvalues list and report them all for this line
                    pvals = {}
                    adjust_pvals = {}
        
                    pval_type = "overallTFs"
                    pvals[pval_type] = dict_pvals[pval_type][dict_line_indices[pval_type].index(line_index)]
                    adjust_pvals[pval_type] = adjusted_dict_pvals[pval_type][dict_line_indices[pval_type].index(line_index)]
                
                    pval_type = "perTF"
                    sub_type = sl[motif_name_index]
                    pvals[pval_type] = dict_pvals[pval_type][sub_type][dict_line_indices[pval_type][sub_type].index(line_index)]
                    adjust_pvals[pval_type] = adjusted_dict_pvals[pval_type][sub_type][dict_line_indices[pval_type][sub_type].index(line_index)]
                
                    pval_type = "perChromatinCat"
                    sub_type = sl[chromatin_index]
                    pvals[pval_type] = dict_pvals[pval_type][sub_type][dict_line_indices[pval_type][sub_type].index(line_index)]
                    adjust_pvals[pval_type] = adjusted_dict_pvals[pval_type][sub_type][dict_line_indices[pval_type][sub_type].index(line_index)]
                    
                    pval_type = "perTF_perChromatinCat"
                    sub_type = sl[motif_name_index]
                    sub_sub_type = sl[chromatin_index]
                    pvals[pval_type] = dict_pvals[pval_type][sub_type][sub_sub_type][dict_line_indices[pval_type][sub_type][sub_sub_type].index(line_index)]
                    adjust_pvals[pval_type] = adjusted_dict_pvals[pval_type][sub_type][sub_sub_type][dict_line_indices[pval_type][sub_type][sub_sub_type].index(line_index)]
                    
                        
                    sl[motif_breaking_score_index+1] = (json.dumps(pvals)+'@'+json.dumps(adjust_pvals)).replace(',', ';')
                    
                    annoted_input_ofile.write('\t'.join(sl) + '\n')
                    sig_level = 1.0
                    if filter_on_qval:
                        sig_level = adjust_pvals[sig_cat]
                    else:
                        sig_level = pvals[sig_cat]
                    #select the motif if its ChiP-seq signal larger than 0
                    if sl[tf_binding_index]!="nan":
                        if float(sl[tf_binding_index]) > 0.0:
                            annoted_input_ofile_onlysig.write('\t'.join(sl) + '\n')
                            continue
                    
                    if filter_on_signal:
                        if (sig_level<sig_thresh and 
                            float(sl[dnase_index])>0.0):# or float(tf_motif[fantom_index])>0.0 or float(tf_motif[num_other_tfs_index])>0.0
                            annoted_input_ofile_onlysig.write('\t'.join(sl) + '\n')
                    else:
                        if sig_level<sig_thresh: #or adjusted_dict_p_values_per_tf[tf][i]<0.05:
                            annoted_input_ofile_onlysig.write('\t'.join(sl) + '\n')
        
    
    return annoted_output_file_onlysig
 
def calculate_p_value_motifregions(mutated_regions_list, num_muts_per_sample_dict, 
                                   total_number_of_regions_tested, 
                                   index_mutation_frequency=12, index_sample_ids=-1, 
                                   index_elment_start_coordinate=1, 
                                   index_elment_stop_coordinate=2, 
                                   genome_size=3000000000.0):
    
    reported_p_values = []#this holds p-values of all the regions, it will be used for p-value correction and after correction it is written to the output file as an additional column after the calculated p-value, the full list of p-values is need to make the correction test
    for element in mutated_regions_list:
        mutation_frequency = int(element[index_mutation_frequency])
        sample_ids = element[index_sample_ids].split(',')
        avg_proportion_of_mutations_in_the_samples_of_this_region = 0.0
        samples_counted = []
        
        for sample_id in sample_ids:
            if sample_id not in samples_counted:
                samples_counted.append(sample_id)
                if sample_id in num_muts_per_sample_dict.keys():
                    avg_proportion_of_mutations_in_the_samples_of_this_region += ((float(num_muts_per_sample_dict[sample_id]))/genome_size)
        p = avg_proportion_of_mutations_in_the_samples_of_this_region#/(len(samples_counted)*1.0)
        n = (int(element[index_elment_stop_coordinate]) - int(element[index_elment_start_coordinate])) #* len(sample_id_and_number_of_mutations_per_sample_dict.keys()) #region length (end-start) multiplied by the total number of tested samples
        k = mutation_frequency
        p_val_of_this_region = 1 - (binom.cdf(k, n, p))
        reported_p_values.append(p_val_of_this_region)
    '''Extend the number of tested elements according to the given total_number_of_regions_tested; set pval of this regions not selected as 1'''
    n_elements = len(reported_p_values)
    for i in range(n_elements, total_number_of_regions_tested):
        reported_p_values.append(1) 
    
    print("correcting p-values for multiple testing")
    if len(reported_p_values)>0:
        significant_bool_report, corrected_p_values_array, alphacSidak, alphacBonf = multipletests(reported_p_values, alpha=0.05, method='fdr_bh', returnsorted=False) #returns 4 things: a boolean array contains True or False for each value meaning wether the value after correction compared to the given alpha is significant or not, an array of the values after correction, a single for corrected alpha for Sidak method, a single value for corrected alpha for Bonferroni method 
        corrected_p_values_list = corrected_p_values_array.tolist()
        for l in range(0, len(mutated_regions_list)):
            mutated_regions_list[l].append(str(reported_p_values[l]))
            mutated_regions_list[l].append(str(corrected_p_values_list[l]))
    return mutated_regions_list
    
def get_number_of_mutations_per_sample_list_and_write_to_file(mutations_file, numberofmutationspersample_output_file, index_sample_ids=8):
    
    num_muts_per_sample_dict = {}  
    if not os.path.exists(numberofmutationspersample_output_file):
        print("Counting number of mutations per sample from the initial mutation file")
        with open(mutations_file, "r") as mutations_infile:
            mutations_line = mutations_infile.readline().strip().split('\t')
            while len(mutations_line)>index_sample_ids:
                sample_id_of_this_mutation = mutations_line[index_sample_ids].strip()
                try:
                    num_muts_per_sample_dict[sample_id_of_this_mutation] +=1
                except KeyError:
                    num_muts_per_sample_dict[sample_id_of_this_mutation] = 1
                mutations_line = mutations_infile.readline().strip().split('\t')
                
        with open(numberofmutationspersample_output_file, 'w') as numberofmutationspersample_writefile:
            for sample_id in num_muts_per_sample_dict.keys():#write the sample and its number of mutations to the output file
                numberofmutationspersample_writefile.write(sample_id + "\t" + str(num_muts_per_sample_dict[sample_id]) +"\n")
    else:
        #in case the number of mutations per sample was already available then just read them and insert them to the returned list
        with open(numberofmutationspersample_output_file, 'r') as numberofmutationspersample_readfile:
            numberofmutationspersample_lines = numberofmutationspersample_readfile.readlines()
            for line in numberofmutationspersample_lines:
                num_muts_per_sample_dict[line.split('\t')[0].strip()] = int(line.split('\t')[1].strip()) # append sample id # append number of mutations in this sample id
    return num_muts_per_sample_dict


def get_unique_muts_from_collection(mutations_file, unique_muts_file, sample_id_index=8, mut_type_index=6, prioritize_SNP_over_indel=True):
    
    muts_per_position_per_sample = {}
    with open(mutations_file, 'r') as ifile:
        l= ifile.readline().strip().split('\t')
        while len(l)>sample_id_index:
            k = '::'.join([l[0], l[1], l[sample_id_index]])
            try:
                muts_per_position_per_sample[k].append("duplicate")
                muts_per_position_per_sample[k].append(l)
            except KeyError:
                muts_per_position_per_sample[k] = l
            l= ifile.readline().strip().split('\t')
    
    with open(unique_muts_file, 'w') as ofile:
        for k in muts_per_position_per_sample.keys():
            if 'duplicate' not in muts_per_position_per_sample[k]:
                ofile.write('\t'.join(muts_per_position_per_sample[k]) + '\n')
                
    return unique_muts_file

def get_unique_mutsmotifs_from_collection(mutations_file, unique_muts_file, sample_id_index=8, mut_type_index=6, prioritize_SNP_over_indel=True, 
                                          motif_start_index=15, motif_end_index=16, motif_name_index=17):
    
    muts_per_position_per_sample = {}
    with open(mutations_file, 'r') as ifile:
        l= ifile.readline().strip().split('\t')
        while len(l)>motif_name_index:
            k = '::'.join([l[0], l[1], l[sample_id_index], l[motif_start_index], l[motif_end_index], l[motif_name_index]])
            try:
                muts_per_position_per_sample[k].append("duplicate")
                print(muts_per_position_per_sample[k])
                muts_per_position_per_sample[k].append(l)
                print(muts_per_position_per_sample[k])
                sys.exit(0)
            except KeyError:
                muts_per_position_per_sample[k] = l
            l= ifile.readline().strip().split('\t')
    
    with open(unique_muts_file, 'w') as ofile:
        for k in muts_per_position_per_sample.keys():
            if 'duplicate' not in muts_per_position_per_sample[k]:
                ofile.write('\t'.join(muts_per_position_per_sample[k]) + '\n')
                
    return unique_muts_file

def calculate_pval_for_genesets(geneset_enrichement_results_input_file, index_total_number_genes_per_set=2, index_number_enriched_genes=3, total_number_of_genes_in_the_universe=27000, total_number_of_genes_tried_in_the_search=3135, header_line = True, number_of_tried_gene_sets=24, keywords_to_filter_out_with=[]):#although not all are recognized in the pathways)
    
    infile = open(geneset_enrichement_results_input_file, 'r')
    calculated_p_value_out_file  = '.'.join(geneset_enrichement_results_input_file.split('.')[0:-1]) + "_calculated_pval." + geneset_enrichement_results_input_file.split('.')[-1]
    calculated_p_value_sig_out_file  = '.'.join(geneset_enrichement_results_input_file.split('.')[0:-1]) + "_calculated_pvalsig." + geneset_enrichement_results_input_file.split('.')[-1]
    calculated_p_value_sig_out_file_keywords  = '.'.join(geneset_enrichement_results_input_file.split('.')[0:-1]) + "_calculated_pval_sig_kw." + geneset_enrichement_results_input_file.split('.')[-1]
    
    outfile = open(calculated_p_value_out_file, "w")
    outfile_sig = open(calculated_p_value_sig_out_file, "w")
    outfile_sig_keywords = open(calculated_p_value_sig_out_file_keywords, "w")
    
    M = total_number_of_genes_in_the_universe
    N = total_number_of_genes_tried_in_the_search
    sep = '\t'
    if header_line:
        header = infile.readline()
        outfile.write(header.strip() +  sep + "p-value" + sep + "q-value" + "\n")
        outfile_sig.write(header.strip() + sep + "p-value" + sep + "q-value"  + "\n")
        outfile_sig_keywords.write(header.strip() + sep + "p-value" + sep + "q-value" + "\n")
    inlines = infile.readlines()
    calculated_pvalues = []
    for line in inlines:
        split_line = line.strip().split(sep)
        if split_line[index_number_enriched_genes].isdigit() and split_line[index_total_number_genes_per_set].isdigit():
            n = int(split_line[index_total_number_genes_per_set])
            x = int(split_line[index_number_enriched_genes])
            pval = hypergeom.sf(x,M,n,N) #previous 1-h.cdf (but gives negative values
            calculated_pvalues.append(pval)
            
    corrected_p_values_list = []
    if len(calculated_pvalues)>1: #if only one geneset is tried then there is no need for multiple test correction
        significant_bool_report, corrected_p_values_array, alphacSidak, alphacBonf = multipletests(calculated_pvalues, alpha=0.05, method='fdr_bh', returnsorted=False) #returns 4 things: a boolean array contains True or False for each value meaning wether the value after correction compared to the given alpha is significant or not, an array of the values after correction, a single for corrected alpha for Sidak method, a single value for corrected alpha for Bonferroni method 
        corrected_p_values_list = corrected_p_values_array.tolist()
    elif len(calculated_pvalues)==1:
        corrected_p_values_list = calculated_pvalues
    else:
        print("No genesets are reported, check the filters, params and gene names")
    number_of_sig_enriched_genesets = 0
    for l in range(0, len(inlines)):
        outfile.write(sep.join(inlines[l].strip().split(sep)) + sep + str(calculated_pvalues[l]) + sep + str(corrected_p_values_list[l]) +"\n")
        #write only the significant ones and satisfying the given condition
        if calculated_pvalues[l]<0.05:
            number_of_sig_enriched_genesets+=1 
            outfile_sig.write(sep.join(inlines[l].strip().split(sep)) + sep + str(calculated_pvalues[l]) + sep + str(corrected_p_values_list[l]) + "\n")
            if len(keywords_to_filter_out_with)>0:
                for keyword in keywords_to_filter_out_with:
                    if keyword.upper() in inlines[l].split(sep)[0].upper() or keyword.upper() in inlines[l].split(sep)[1].upper(): #if the given keyword(s) was found in the gene set name or discreption then report 
                        outfile_sig_keywords.write(sep.join(inlines[l].strip().split()) + sep + str(calculated_pvalues[l]) + sep + str(corrected_p_values_list[l]) + "\n")
                        break
    print("Number of significantly enriched genesets: "  + str(number_of_sig_enriched_genesets))
    return calculated_p_value_out_file, calculated_p_value_sig_out_file, calculated_p_value_sig_out_file_keywords
    
def find_overlap_genesets_genelist(geneset_input_file, genelist_input_file, 
                                   enriched_genes_output_file, 
                                   total_number_of_genes_in_the_universe=27000, 
                                   min_number_of_genes_be_enriched_for_geneset_to_be_reported = 10, 
                                   index_gene_name=0, index_gene_names_start=3, 
                                   keywords_to_filter_out_with=[], 
                                   only_keep_the_sig_file = True, 
                                   min_number_of_genes_in_geneset_to_consider_the_geneset = 10, 
                                   header_line = False,
                                   sample_ids_given=False):
    
    with open(geneset_input_file, 'r') as ifile:
        genesets_lines = ifile.readlines()
    with open(genelist_input_file, 'r') as ifile:
        genelist_lines = [x.strip().split('\t') for x in ifile.readlines()]
    enriched_genes_outfile = open(enriched_genes_output_file, 'w')
    
    #read the gene names from the gene input list
    print("Number of genes provided for search: " + str(len(genelist_lines)))
    print("Total number of genesets to try: " + str(len(genesets_lines)))
    
    number_of_tried_genesets = 0
    if header_line:
        enriched_genes_outfile.write('ID' + "\t" + 'description' + "\t"+ 'total_number_of_genes_in_this_geneset' + "\t" + 'number_of_enriched_genes' + "\t"+ 'enriched_genes' + "\n")

    for geneset in genesets_lines:
        split_geneset_info = geneset.split('\t')
        total_number_of_genes_in_this_geneset = 0
        if index_gene_names_start>2:
            total_number_of_genes_in_this_geneset = int(split_geneset_info[index_gene_names_start-1])#if there were three columns before the gene names start, it means the total number of genes are already given per gene set in addition to the pathway info and ID. This is especially important for KEGG pathways converted from KEGG IDs since the gene names are wtitten in many variations so it is not possible to caount them to get the actual number of genes in the set
        else:
            total_number_of_genes_in_this_geneset = len(set(split_geneset_info[index_gene_names_start::]))
        if total_number_of_genes_in_this_geneset > min_number_of_genes_in_geneset_to_consider_the_geneset:
            number_of_tried_genesets +=1
            #only search in those sets that fulfill the criteria
            enriched_genes = []
            regmut_samples = []
            mut_samples = []
            for gene_line in genelist_lines:
                if gene_line[index_gene_name].strip() in split_geneset_info[index_gene_names_start::]:
                    enriched_genes.append(gene_line[index_gene_name])
                    if sample_ids_given:
                        mut_samples.extend(gene_line[-1].split(','))
                        regmut_samples.extend(gene_line[-2].split(','))
                    
            if len(set(enriched_genes)) >= min_number_of_genes_be_enriched_for_geneset_to_be_reported:
                enriched_genes_outfile.write(split_geneset_info[0] + "\t" + split_geneset_info[1].replace(' - Homo sapiens (human)','') + "\t"+ str(total_number_of_genes_in_this_geneset) + "\t" + str(len(set(enriched_genes))) + '\t' + str(len(set(regmut_samples))) + "\t" + str(len(set(mut_samples))) + '\t' +','.join(set(regmut_samples)) + "\t" + ','.join(set(mut_samples)) + "\t"+ ','.join(set(enriched_genes)) +"\n")
            
    enriched_genes_outfile.close()
    #calculate p-values for each gene set/pathway
    calculated_p_value_out_file, calculated_p_value_sig_out_file, calculated_p_value_sig_out_file_keywords  = calculate_pval_for_genesets(enriched_genes_output_file, index_total_number_genes_per_set=2, index_number_enriched_genes=3, total_number_of_genes_in_the_universe=total_number_of_genes_in_the_universe, total_number_of_genes_tried_in_the_search=len(genelist_lines), header_line = header_line, number_of_tried_gene_sets=number_of_tried_genesets, keywords_to_filter_out_with=keywords_to_filter_out_with)
    
    if len(keywords_to_filter_out_with)<1:
        os.remove(calculated_p_value_sig_out_file_keywords)
    del genesets_lines
    del genelist_lines
    
    if only_keep_the_sig_file:
        os.remove(calculated_p_value_out_file)
        os.remove(enriched_genes_output_file)    

    return calculated_p_value_sig_out_file

def get_simulated_mean_sd_per_TF_motif_background_window_correction(cohort_full_name, annotated_input_file, simulated_annotated_input_files, 
                                       mutations_cohorts_dir,
                                       cohort_mean_sd_per_tf_overall_output_dict_file, 
                                       chr_lengths_file,
                                       background_window_size = 50000,
                                       motif_name_index = 17, f_score_index = 9, 
                                       motif_breaking_score_index = 10, chromatin_cat_index=22, tmp_dir = '$SNIC_TMP'):
    
    
    copyfile(cohort_mean_sd_per_tf_overall_output_dict_file, cohort_mean_sd_per_tf_overall_output_dict_file + '_copy.dict') 
    
    print("Extracting avg and std per TF and overall from the simulation sets... onto: ", cohort_mean_sd_per_tf_overall_output_dict_file)
    cohort = cohort_full_name.split('/')[-1]
    tmp_dir_intersect = mutations_cohorts_dir + '/' + cohort + '_tmp_pybedtoos/'
    if not os.path.exists(tmp_dir_intersect):
        os.mkdir(tmp_dir_intersect) 
    
    chr_lengths = get_chr_lengths(chr_lengths_file)
    #divided observed mutations files into subfiles. Extend mutations with the backgroud window
    splited_file_name = tmp_dir  + '/' + cohort + '_splited'
    splited_file_name_sorted = splited_file_name + '_sorted'
    print(splited_file_name)

    splited_file_name_local = tmp_dir_intersect  + '/' + cohort + '_splited'
    
    #lines_per_file = 10000
    if not os.path.exists(splited_file_name_local):
        line_number = 0
        with open(annotated_input_file, 'r') as observed_infile, open(splited_file_name, "w") as splited_ifile:
            l = observed_infile.readline().strip().split('\t')
            while l and len(l)>3:
                motif_start = (int(l[1])-background_window_size)
                motif_end = int(l[2])+background_window_size
                motif_names = l[motif_name_index]
                chrom_cat = l[chromatin_cat_index]
                chr_name = l[0].replace('chr','')
                chr_name2 = chr_name.replace('X', '23').replace('Y','24').replace('MT','25').replace('M','25')

                if motif_start<0:
                    motif_start = 0
                    motif_end += 0 - (int(l[1])-background_window_size)
                if motif_end>chr_lengths[int(chr_name2)]:
                    motif_end = chr_lengths[int(chr_name2)]
                    motif_start -= (int(l[2])+background_window_size) - chr_lengths[int(chr_name2)]
                #if line_number % lines_per_file == 0:
                #    if splited_file:
                #        splited_file.close()
                #    splited_file_name = splited_files_name + '_{}'.format(line_number)
                #    splited_files_list.append(splited_file_name)
                #    splited_file = open(splited_file_name, "w")
                # save background window, motif name, chromatin cat, and number of line
                splited_ifile.write(chr_name + '\t' + str(motif_start) + '\t' +   str(motif_end) + '\t' + str( motif_names) + '\t' +chrom_cat + '\t' + str(line_number) + '\n')
                line_number+=1
                l = observed_infile.readline().strip().split('\t')
            #if splited_file:
            #    splited_file.close()
        #copy file from scratch to project folder
        awk_stmt_split_sort = """grep -E '^X|^Y|^M' {splited_file_name} |sort -k1,1n -k2,2n > {splited_file_name_sorted}""".format(splited_file_name = splited_file_name, splited_file_name_sorted = splited_file_name_sorted )
        os.system(awk_stmt_split_sort)
        
        copyfile(splited_file_name_sorted, splited_file_name_local)      
        os.remove(splited_file_name)

    observed_input_file_obj = BedTool(splited_file_name_sorted)
    #define motif breaking score and fscore for the intersected files
    new_motif_breaking_score_index = motif_breaking_score_index + 6
    new_fscore_index = f_score_index + 6
    #define extensions for the merged files for all categories
    simulated_input_file_tmp_overallTFs_extension ="_tmp_overallTFs"
    #simulated_input_file_tmp_perTF_extension = "_tmp_perTF"
    #simulated_input_file_tmp_perChromatinCat_extension = "_tmp_perChromatinCat"
    #simulated_input_file_tmp_perTF_perChromatinCat_extension = "_tmp_perTF_perChromatinCat"
    simulated_files_temp = []
    # intersection the observed mutation file with the simulated file to find the background
    # group by: 1 - position for all TFs within the window
    #2 - positions of the same TF motifs within the window in the simulated sets

    for simulated_input_file in simulated_annotated_input_files:
            simulated_input_file_name = simulated_input_file.split('/')[-1]
            
            
            simulated_input_file_tmp_overallTFs = tmp_dir +'/' + simulated_input_file_name + '_' + splited_file_name.split('_')[-1] + simulated_input_file_tmp_overallTFs_extension
            simulated_input_file_tmp_overallTFs_local = tmp_dir_intersect + simulated_input_file_name + '_' + splited_file_name.split('_')[-1] + simulated_input_file_tmp_overallTFs_extension

            #simulated_input_file_tmp_TFs = tmp_dir +'/' + simulated_input_file_name + '_' + splited_file_name.split('_')[-1] + simulated_input_file_tmp_perTF_extension
            #simulated_input_file_tmp_chromatin = tmp_dir +'/' + simulated_input_file_name + '_' + splited_file_name.split('_')[-1] + simulated_input_file_tmp_perChromatinCat_extension
            #simulated_input_file_tmp_TFs_chromatin = tmp_dir +'/' + simulated_input_file_name + '_' + splited_file_name.split('_')[-1] + simulated_input_file_tmp_perTF_perChromatinCat_extension
            
            
            #check if mutation position is string and convert to intiger
            #remove from the simulation file rows where mut positions are string and compare number of lines
            simulated_input_file_position = tmp_dir + simulated_input_file_name + '_pos'
            simulated_ifile_pos_temp = simulated_input_file_position + '_tmp'
            awk_stmt_sim  = """awk 'BEGIN{{FS=OFS="\t"}} {{if ( $2 ~ "^[0-9][0-9]*$" && $3 ~ "^[0-9][0-9]*$") print $0}} ' {sim_ifile} > {sim_ofile}""".format(sim_ifile = simulated_input_file, sim_ofile =simulated_input_file_position )
            os.system(awk_stmt_sim)
            count = len(open(simulated_input_file).readlines(  ))
            count2 = len(open(simulated_input_file_position).readlines(  ))
            if(count != count2):
                awk_tmp =r"""awk 'BEGIN{{FS=OFS="\t"}}{{ printf ("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32)}}' {sim_ifile} > {sim_ofile} """.format(sim_ifile = simulated_input_file, sim_ofile = simulated_ifile_pos_temp)
                os.system(awk_tmp)
                simulated_input_file = simulated_ifile_pos_temp 
            #os.remove(simulated_input_file_position)
            print(simulated_input_file)
            if not os.path.exists(simulated_input_file_tmp_overallTFs):
                #check if 'chr' is present
                with open(simulated_input_file, 'r') as simulated_ifile:
                    line = simulated_ifile.readline()
                    print(line)
                    if line[0:3] == 'chr':
                        simulated_ifile_temp = tmp_dir + simulated_input_file_name + '_tmp'
                        #simulated_ifile_temp = simulated_input_file + '_tmp'
                        awk_stmt = """cat {simulated_file} | sed 's/^...//'  > {simulated_outfile_temp}""".format(simulated_file = simulated_input_file, simulated_outfile_temp = simulated_ifile_temp)
                        os.system(awk_stmt)
                        simulated_files_temp = simulated_ifile_temp
                        simulated_input_file = simulated_ifile_temp
                print(simulated_input_file)
                simulated_input_file_sorted = tmp_dir + simulated_input_file_name + '_sorted'
                awk_stmt_sort = """grep -E '^X|^Y|^M' {simulated_input_file} | sort -k1,1n -k2,2n  > {simulated_input_file_sorted}""".format(simulated_input_file = simulated_input_file,simulated_input_file_sorted = simulated_input_file_sorted )
                os.system(awk_stmt_sort)
                simulated_input_file_obj = BedTool(simulated_input_file_sorted)                
                #intersect the simulated file with the observed mutation file. Provide a sum of f_score and motif breaking score
                observed_input_file_obj_inter = observed_input_file_obj.intersect(simulated_input_file_obj, wo = True).each(sum_fscore_motif_breaking_score, new_fscore_index, new_motif_breaking_score_index).saveas()
                #group files to obtain the mean and stdev for the functional score
                #os.remove(simulated_input_file + '_sorted')
                try: 
                    observed_input_file_obj_inter.groupby(g=[1,2,3,4,5,6], c=16, o=['mean', 'stdev', 'count']).saveas(simulated_input_file_tmp_overallTFs)
                    #observed_input_file_obj_inter.filter(lambda x: str(x[3]) == str(x[23])).groupby(g=[1,2,3,4,5,6], c=16, o=['mean', 'stdev', 'count']).saveas(simulated_input_file_tmp_TFs)
                    #observed_input_file_obj_inter.filter(lambda x: str(x[4]) == str(x[28])).groupby(g=[1,2,3,4,5,6], c=(new_fscore_index+1), o=['mean', 'stdev', 'count']).saveas(simulated_input_file_tmp_chromatin)
                    #observed_input_file_obj_inter.filter(lambda x: (str(x[4]) == str(x[28])) & (str(x[3]) == str(x[23]))).groupby(g=[1,2,3,4,5,6], c=(new_fscore_index+1), o=['mean', 'stdev', 'count']).saveas(simulated_input_file_tmp_TFs_chromatin)
                except KeyError:
                    open(simulated_input_file_tmp_overallTFs, 'a').close()
                    #open(simulated_input_file_tmp_perTF, 'a').close()
                    #open(simulated_input_file_tmp_chromatin, 'a').close()
                    #open(simulated_input_file_tmp_perTF_perChromatinCat_extension, 'a').close()
                #print(os.listdir(tmp_dir_intersect))
                copyfile(simulated_input_file_tmp_overallTFs, simulated_input_file_tmp_overallTFs_local)         
                #if "_tmp" in simulated_input_file:
                #    os.remove(simulated_input_file)
                #os.remove(simulated_input_file_sorted)
            cleanup()   
    #list of categories for simulated_mean_sd_files
    simulated_mean_sd_cat = ["overallTFs"]
    #simulated_mean_sd_cat = ["overallTFs", "perTF", "perChromatinCat", "perTF_perChromatinCat"]
    #print(simulated_mean_sd_cat)
    
    #create a dictionery for mean, std scores for all categories
    with open(cohort_mean_sd_per_tf_overall_output_dict_file, 'r') as dict_simulated_mean_sd_per_TF_motif_ifile:
        dict_type_mean_std_scores = json.loads(dict_simulated_mean_sd_per_TF_motif_ifile.readline())
    
    #dict_type_mean_std_scores = {}
    
    mean_f_score_index = 5
    for cat_type in simulated_mean_sd_cat:
        print('Combining dictionery')
        simulated_mean_sd_files = tmp_dir + '/' +'*_tmp_' + cat_type
        simulated_mean_sd_outfiles = tmp_dir + '/' + cohort + '_'+ cat_type
        simulated_mean_sd_outfiles_local = mutations_cohorts_dir + '/' + cohort + '_'+ cat_type

        #merge files from the same category, sort by the line number and group by position, TF motif, chromatin cat. and line number
        awk_comm = """cat {files} | 
        sort -k6 -V | 
        groupBy -g 1-6 -c 7,8,9 -o mean,mean,sum > {file_out} """.format(files = simulated_mean_sd_files, file_out = simulated_mean_sd_outfiles)
        os.system(awk_comm)
        #save the scores per line number
        dict_simulated_mean_sd = {}


        with open(simulated_mean_sd_outfiles, 'r') as simulated_mean_sd_ifile:
            l = simulated_mean_sd_ifile.readline().strip().split('\t')
            while l and len(l)>3:  
                dict_simulated_mean_sd[l[mean_f_score_index]] = {'mean': l[mean_f_score_index +1 ], 
                                                       "std": l[mean_f_score_index +2 ], 
                                                       "nummotifs": l[mean_f_score_index +3 ]}
                l = simulated_mean_sd_ifile.readline().strip().split('\t')
       #     #save the dictionery per category
        dict_type_mean_std_scores[cat_type].update(dict_simulated_mean_sd)
    #print(dict_type_mean_std_scores)

    with open(cohort_mean_sd_per_tf_overall_output_dict_file, 'w') as dict_simulated_mean_sd_per_TF_motif_outfile:
            json.dump(dict_type_mean_std_scores, dict_simulated_mean_sd_per_TF_motif_outfile)
    
    #if os.path.exists(tmp_dir_intersect):
    #    shutil.rmtree(tmp_dir_intersect)
            
    return  dict_type_mean_std_scores

def get_sig_merged_elements_oncodrive(unified_mutation_input_files, mutation_input_files, cohort_full_name, 
                            output_extension, distance_to_merge, 
                            merged_mut_sig_threshold, local_domain_window, 
                            chr_lengths_file, sig_elements_output_file, 
                            sig_thresh):
    
    
    if os.path.exists(sig_elements_output_file):
        return sig_elements_output_file
    
    "Merge the observed mutations into elements"
    unified_observed_muts_file_wihtmotifinfo = unified_mutation_input_files
    merged_muts_output_file = unified_observed_muts_file_wihtmotifinfo+"_mergedmuts{distance_to_merge}bp".format(
                                            distance_to_merge=distance_to_merge)
    Utilities.merge_muts(muts_input_file=unified_observed_muts_file_wihtmotifinfo, 
                         merged_muts_output_file=merged_muts_output_file, 
                         filter_mut_motifs=False, filter_col_index=15, 
                         filter_value=merged_mut_sig_threshold, mut_score_index=9, 
                         motifs_col_index =10, ref_alt_col_index=11, mutpos_col_index=12, 
                         motifname_col_index=13, motif_col_index=14, 
                         distance_to_merge=distance_to_merge)
    
    
    '''Prepare mutation file'''
    print('Calcuate pval for each element using OncodriveFML')
    mutation_file_oncodrive = mutation_input_files + '_oncodrive'
    fsep = '\t'
    awk_stmt_mut = """awk 'BEGIN{{FS=OFS="{fsep}"}}{{print $1,$2,$4,$5,$8,$6}}' {infile} | sort -k1,1n -k2,2n | uniq -u | awk 'BEGIN{{FS=OFS="\t"}}{{gsub("23","X", $1); gsub("24","Y", $1); gsub("chr","", $1); print $0}}' > {mutation_file}""".format(
                                                fsep=fsep,  infile=mutation_input_files, mutation_file=mutation_file_oncodrive+'_header')
    os.system(awk_stmt_mut)
    
    #add header
    awk_stmt_mut2 = """echo "CHROMOSOME\tPOSITION\tREF\tALT\tSAMPLE\tCANCER_TYPE" | cat - {infile} > {mutation_file}""".format(
                                                 infile=mutation_file_oncodrive+'_header', mutation_file=mutation_file_oncodrive)
    os.system(awk_stmt_mut2)
    
    os.remove(mutation_file_oncodrive +'_header')
    
    '''Prepare elements file'''
    element_file_oncodrive = merged_muts_output_file + '_oncodrive'
    
    filter_cond = 'if($6>=1)' #remove elements with one mutation
    awk_stmt_elem = """awk 'BEGIN{{FS=OFS="{fsep}"}}{{{filter_cond} {{print $1,$2,$3,$15}}}}' {infile} | sort -k1,1n -k2,3n | uniq -u | awk 'BEGIN{{FS=OFS="\t"}}{{gsub("23","X", $1); gsub("24","Y", $1); print $0}}' > {element_file}""".format(
                                                fsep=fsep, filter_cond= filter_cond, infile=merged_muts_output_file, element_file=element_file_oncodrive+'_header')
    os.system(awk_stmt_elem)
    
    filesize = os.path.getsize(element_file_oncodrive+'_header')
    #check if elements exist
    if filesize == 0:
        os.remove(element_file_oncodrive +'_header')
        with open(sig_elements_output_file, 'w') as fp: 
            pass
    else:
        
        #add header
        awk_stmt_elem2 = """echo "CHROMOSOME\tSTART\tEND\tELEMENT" | cat - {infile} > {element_file}""".format(
                                                     infile=element_file_oncodrive+'_header', element_file=element_file_oncodrive)
        os.system(awk_stmt_elem2)
        
        os.remove(element_file_oncodrive +'_header')
        
        '''Calcuate pval for each element using oncodrivefml'''
        
        tmp_dir = mutation_file_oncodrive + '_tmp'
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir) 
            
        awk_stm_oncodrive ="""oncodrivefml -i {mutation_file} -e {element_file} -s wgs -c /proj/snic2020-16-50/nobackup/pancananalysis/pancan12Feb2020/cancer_datafiles/oncodrivefml_v2.conf -o {oncodrive_dir}""".format(mutation_file = mutation_file_oncodrive,
                                                                                                        element_file = element_file_oncodrive, oncodrive_dir = tmp_dir)
        
        #print(awk_stm_oncodrive)
        os.system(awk_stm_oncodrive)
        
        #oncodrive result: tsv file
        oncodrive_out_file = [tmp_dir+'/'+x for x in os.listdir(tmp_dir) if '.tsv' in x]
        #merge elements with oncodrive results
        merged_elements_statspvalues = merged_muts_output_file+"_statspvalues"    
        
        element=pd.read_csv(merged_muts_output_file, sep="\t",  header=None)
        oncodrive_element=pd.read_csv(oncodrive_out_file[0], sep="\t")
        
        merged_element = element.merge(oncodrive_element, left_on=14, right_on='GENE_ID')
        #remove unnecessary columns
        merged_element_removed_columns = merged_element.drop(['GENE_ID','MUTS', 'MUTS_RECURRENCE', 'SAMPLES','SNP', 'MNP','INDELS', 'SYMBOL','P_VALUE_NEG', 'Q_VALUE_NEG'], axis=1)
        merged_element_removed_columns.to_csv(merged_elements_statspvalues, index=False, sep='\t', header =False)
        #find significant elements in oncodrive results
        awk_stm_sig_elem = """awk 'BEGIN{{FS=OFS="{fsep}"}}{{if ($17<= {sig_thresh} && $17 != "") print $0,$16,$17; else if ($16<= {sig_thresh} && $17 == "") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$16,$16,$16}}' {infile} > {merged_elements_statspvaluesonlysig}""".format(
        fsep=fsep, sig_thresh=sig_thresh,infile = merged_elements_statspvalues,  merged_elements_statspvaluesonlysig=sig_elements_output_file)
        os.system(awk_stm_sig_elem)
       
    os.remove(element_file_oncodrive)   
    os.remove(mutation_file_oncodrive)
    return sig_elements_output_file


def get_sig_muts(elements_input_file, mutations_input_file, sig_muts_file, motif_breaking_score_index, tf_binding_index, dnase_index):
    
    
    if os.path.exists(sig_muts_file):
        return sig_muts_file
    
    breaking_score_threshold = 0.3
    
    elements_input_file_tmp=elements_input_file+'_tmp'
    os.system("""awk 'BEGIN{{FS=OFS="\t"}}{{gsub("23","X", $1); gsub("24","Y", $1); print "chr"$0}}' {}>{}""".format(elements_input_file,elements_input_file_tmp ))

    sig_muts_file_tmp=sig_muts_file+"_tmp"
    sig_muts_file_tmp2=sig_muts_file_tmp+"2"
         
    BedTool(elements_input_file_tmp).intersect(BedTool(mutations_input_file), wb=True).saveas(sig_muts_file_tmp)
    os.system("""awk 'BEGIN{{FS=OFS="\t"}}{{print $22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53}}' {} | sort -k1,1 -k2,2n -V | uniq  >{}""".format(sig_muts_file_tmp,sig_muts_file_tmp2))
    print('YES')
    with open(sig_muts_file_tmp2, 'r') as annoted_output_ifile, open(sig_muts_file, 'w') as annoted_input_ofile_onlysig:
            l = annoted_output_ifile.readline()
            while l:
                sl = l.strip().split('\t')
                if (len(sl) >= 3 and (float(sl[motif_breaking_score_index])>=breaking_score_threshold)):
                        if sl[tf_binding_index]!="nan":
                            if float(sl[tf_binding_index]) > 0:
                                annoted_input_ofile_onlysig.write(l)
                                l = annoted_output_ifile.readline()
                                continue
                                
                            
                        if (float(sl[dnase_index])>0.0):# or float(sl[fantom_index])>0.0 or float(sl[num_other_tfs_index])>0.0
                            annoted_input_ofile_onlysig.write(l)
                l = annoted_output_ifile.readline()
    
    #os.remove(elements_input_file_tmp)
    #os.remove(sig_muts_file_tmp)        
    return sig_muts_file