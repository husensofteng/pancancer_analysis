'''
Created on Mar 15, 2017

@author: husensofteng
'''
import os,sys
import numpy as np
import pandas as pd
from collections import Counter
from pybedtools import BedTool
from ProcessCohorts import process_cohorts
from decimal import Decimal
from AnnotateMutations import get_annotated_muts
from Utilities import get_number_of_mutations_per_sample_list_and_write_to_file, calculate_p_value_motifregions, find_overlap_genesets_genelist 
def get_chr_lengths(chr_lengths_file):
    chr_lengths = {}
    with open(chr_lengths_file, 'r') as chr_min_max_ifile:
        lines = chr_min_max_ifile.readlines()
        for l in lines:
            if not l.startswith('chr') and l!="" and not l.startswith('//') and not l.startswith('#'):
                sl = l.strip().split('\t')
                chr_name = "chr"+sl[0]
                chr_lengths[chr_name] = int(sl[1])
    return chr_lengths

def generate_extended_regions(regions, extended_output_file, chr_lengths, window):
    
    with open(extended_output_file, 'w') as ofile:
        for region in regions:
            extended_element_start = (int(region[1])-window)
            extended_element_end = int(region[2])+window
            if extended_element_start<0:
                extended_element_start = 0
                extended_element_end += 0 - (int(region[1])-window)
            if extended_element_end>chr_lengths[region[0]]:
                extended_element_end = chr_lengths[region[0]]
                extended_element_start -= (int(region[2])+window) - chr_lengths[region[0]]
            ofile.write(region[0] + '\t' + str(extended_element_start) + '\t' + str(extended_element_end) + '\t' + str(region[3]) + '\n')
    return extended_output_file

def get_nearby_genes(regions_input_file, genes_input_file,
                     gene_types_to_consider = ['protein_coding', 'RNA'], gene_status_to_consider = ['KNOWN'],
                     n=3, upstream=True, downstream=True, overlapping = True, max_dist = 100000,
                     genesets_to_check_for_overlap_input_files=[]):
    
    regions_input_file_intersect_genes = regions_input_file+"intersect_genes"
    BedTool(regions_input_file).intersect(BedTool(genes_input_file), wo=True).saveas(regions_input_file_intersect_genes)#groupby(g=4, c=[], o=[])
    
    regions_ugenes_dict = {}
    regions_dgenes_dict = {}#contains gene info of each unique region - based on index
    regions_ogenes_dict = {}
    with open(regions_input_file_intersect_genes, 'r') as ifile:
        l = ifile.readline().strip().split('\t')
        while l:
            try:
                if l[11] in gene_status_to_consider and l[12] in gene_types_to_consider:
                    region_start = int(l[3].split(':')[1].split('-')[0])
                    region_end = int(l[3].split(':')[1].split('-')[1])
                    gene_start = int(l[5])
                    gene_end = int(l[6])
                    
                    if ((region_start >= gene_start and region_end <= gene_end) or
                        (region_start>=gene_start and region_start<=gene_end and region_end>=gene_end) or
                        (region_start<=gene_start and region_end<=gene_end and region_end>=gene_start)):
                        if overlapping:
                            d = region_start-gene_start
                            if l[8]=='-':
                                d = region_end-gene_end
                            gene_info = l[10]+"::"+l[9]
                            try:
                                if len(regions_ogenes_dict[l[3]])<=n:#no limit on overlapping genes
                                    regions_ogenes_dict[l[3]].append([gene_info, d])
                            except KeyError:
                                regions_ogenes_dict[l[3]] = [[gene_info, d]]
                    
                    elif (region_start > gene_end):#genes before the region start (upstream genes) 
                        if upstream:
                            d = region_start-gene_start
                            if l[8] == '-':
                                d = region_start-gene_end
                            gene_info = l[10]+"::"+l[9]
                            if d<max_dist:
                                try:
                                    if len(regions_ugenes_dict[l[3]])<n:
                                        regions_ugenes_dict[l[3]].append([gene_info, d])
                                    else:#there are n genes already, remove a gene with a larger distance
                                        regions_ugenes_dict[l[3]] = sorted(regions_ugenes_dict[l[3]],key=lambda l:l[1], reverse=False)
                                        if regions_ugenes_dict[l[3]][n-1][1] > d:
                                            regions_ugenes_dict[l[3]][n-1] = [gene_info, d]
                                except KeyError:
                                    regions_ugenes_dict[l[3]] = [[gene_info, d]]
                        
                    elif (region_start < gene_end):
                        if downstream:
                            d = gene_start-region_end#region_start
                            if l[8] == '-':
                                d = gene_end-region_end#region_start
                            gene_info = l[10]+"::"+l[9]
                            if d<max_dist:
                                try:
                                    if len(regions_dgenes_dict[l[3]])<n:
                                        regions_dgenes_dict[l[3]].append([gene_info, d])
                                    else:#there are n genes already, remove a gene with a larger distance
                                        regions_dgenes_dict[l[3]] = sorted(regions_dgenes_dict[l[3]],key=lambda l:l[1], reverse=False)
                                        if regions_dgenes_dict[l[3]][n-1][1] > d:
                                            regions_dgenes_dict[l[3]][n-1] = [gene_info, d]
                                except KeyError:
                                    regions_dgenes_dict[l[3]] = [[gene_info, d]]
                                
            except IndexError:
                l = ifile.readline().strip().split('\t')
                break
            l = ifile.readline().strip().split('\t')
    os.remove(regions_input_file_intersect_genes)
    
    regions_genes_dict = {}
    for reg in regions_ogenes_dict:
        genes = regions_ogenes_dict[reg]
        for i,g in enumerate(genes):
            try:
                regions_genes_dict[reg].append(g[0]+"::{i}O{dist}".format(i=i, dist=g[1]))
            except KeyError:
                regions_genes_dict[reg] = [g[0]+"::{i}O{dist}".format(i=i, dist=g[1])]
    
    for reg in regions_ugenes_dict:
        if reg in regions_ogenes_dict.keys():
            continue
        genes = regions_ugenes_dict[reg]
        for i, g in enumerate(genes):
            try:
                regions_genes_dict[reg].append(g[0]+"::{i}U{dist}".format(i=i, dist=g[1]))
            except KeyError:
                regions_genes_dict[reg] = [g[0]+"::{i}U{dist}".format(i=i, dist=g[1])]
    
    for reg in regions_dgenes_dict:
        if reg in regions_ogenes_dict.keys():# or reg in regions_ugenes_dict.keys():
            continue
        genes = regions_dgenes_dict[reg]
        for i,g in enumerate(genes):
            try:
                regions_genes_dict[reg].append(g[0]+"::{i}D{dist}".format(i=i, dist=g[1]))
            except KeyError:
                regions_genes_dict[reg] = [g[0]+"::{i}D{dist}".format(i=i, dist=g[1])]
    return regions_genes_dict

def aggregate_results(regions_input_file):
    #generated_sig_merged_element_files = process_cohorts()
    """
    Desired output format:
    Position
    #cohorts    cohorts(unique)
    Information from functional mutations (aggregate across all merged cohorts) (RegMuts):
        - Summary score: sum of all unique mutations
        - Summary FDR: min FDR from all merged elements
        - Number of unique mutations in the element
        - Number of unique mutations per cancer type: vector of Cancer_type:#muts
        - Number of unique samples in the element
        - Number of unique samples per cancer type: vector of Cancer_type:#samples with mut
        *- Number of mutations per TF-motif (vector)
        *- Number of mutations per chromatin state (vector)
        *- Number of unique samples per TF-motif (vector)
        *- Number of unique samples per chromatin state (vector)  
    Info from all mutations (AllMuts):
        - Number of unique mutations in the element
        - Number of unique mutations per cancer type: vector of Cancer_type:#muts
        - Number of unique samples in the element
        - Number of unique samples per cancer type: vector of Cancer_type:#samples with mut
        *- Number of unique mutations per chromatin state
        *- Number of unique samples per chromatin state
    *Infor from all mutated motifs (sig and not sig) (MotifMuts):
        - Number of unique mutations per TF-motif
        - Number of unique samples per TF-motif
        - Number of unique mutations per chromatin state
        - Number of unique samples per chromatin state
    Mutations info: a unique list of all functional mutations in the element
        - list mutated motifs per mutation to be supplied for plotting:
            - heatmap of chromHMM per cancer type
            - heatmap of TF motif per cancer type
            - heatmap of TF motif per chrom HMM
            - DNase/Other enrichment per cancer type
    """
    
    #to keep the order of the columns define a list with the same dict key names
    aggregated_lines = []
    summaries_dict = {'#Elements':0, '#RegMuts':0, '#Samples(RegMuts)':0, '#Muts':0, '#Samples':0}
    
    with open(regions_input_file, 'r') as final_input_ifile:
        lines = [l.strip().split('\t') for l in final_input_ifile.readlines()]
        for l in lines:
            cols_dict = {'chr':'', 'start':'', 'end': '', 'Position':'', 
                 'Cohorts':[], '#Cohorts':0, 'Score':0.0, 'FDR':0.0, 
                     '#RegMuts':0, '#Samples(RegMuts)':0, 'Cancer-Types:#RegMuts':{}, 'Cancer-Types:#Samples(RegMuts)':{},
                     '#Muts':0, '#Samples':0, 'Cancer-Types:#Muts':[], 'Cancer-Types:#Samples':[],
                     'StatsMuts':[],'StatsSamples':[],
                     'RegMuts':[],'Muts':[], 'Mutated-Moitfs':[],'Mutated-Motifs-All':[], 'Max-RegMotif': [],
                     'SamplesMuts':[]}
            
            
            cols_dict['chr'] = l[0]
            #cols_dict['start'] = l[1]
            #cols_dict['end'] = l[2]
            #cols_dict['Position'] = '{}:{}-{}'.format(l[0],l[1],l[2])
            
            cols_dict['Cohorts'] = set(l[4].split(','))
            cols_dict['#Cohorts'] = len(set(l[4].split(',')))
            
            #Functional mutations from sig elements across the cohorts
            mutations_in_cohorts = []
            cohorts_info = [x.strip().split('~') for x in l[10].split(',')]
            
            cols_dict['FDR'] = float(cohorts_info[0][18])
            for cohort_info in cohorts_info:
                if float(cohort_info[18]) > cols_dict['FDR']:
                    cols_dict['FDR'] = float(cohort_info[18])
                
                mutations_in_cohorts.extend(cohort_info[14].split('|'))
            muts_per_cancer_type = {}
            motifs_per_cancer_type = {}
            
            samples_regmuts_per_chromhmm = {}
            samples_regmuts_per_TFmotif = {}
            samples_regmuts_per_TFmotif_position = {}
            samples_regmuts_per_TFmotif_allele = {}
            samples_regmuts_per_TFmotif_position_allle = {}
            mutated_motifs = []
            
            samples_in_this_element = []
            regmuts_already_added = []
            for mut in mutations_in_cohorts:
                regmut_info = mut.split('MatchingMotifs')[0].split('#')
                motifs_info = [x.split('#') for x in mut.split('MatchingMotifs')[1].split('MaxMotif')[0].split('MotifInfo')]
                max_motifs_info = mut.split('MatchingMotifs')[1].split('MaxMotif')[1].split('#')
                '''add one mut and one max motif for each unique mutation (there can be duplicates since one mut can be present in more than one cohort)'''
                
                if regmut_info not in regmuts_already_added:
                    if regmut_info[8] not in samples_in_this_element:
                        samples_in_this_element.append(regmut_info[8])
                        cols_dict['#Samples(RegMuts)']+=1
                    max_motifs_info[3] = regmut_info[5]+"#"+regmut_info[8]
                    cols_dict['Max-RegMotif'].append('#'.join(max_motifs_info))
                    cols_dict['Score']+=float(regmut_info[9])
                    cols_dict['#RegMuts']+=1
                    cols_dict['RegMuts'].append('#'.join(regmut_info))
                    regmuts_already_added.append(regmut_info)
                    
                    '''Append sample ID to each ChromHMM state in order to get enrichment of mutations per chromatin states per element'''
                    '''get chromHMM state from the max motif of the regMut. All motifs (approximately) have the same state'''
                    try:
                        samples_regmuts_per_chromhmm[max_motifs_info[14]].append(regmut_info[8])
                    except KeyError:
                        samples_regmuts_per_chromhmm[max_motifs_info[14]]= [regmut_info[8]]
                    try:
                        muts_per_cancer_type[regmut_info[5]].append(regmut_info)
                    except KeyError:
                        muts_per_cancer_type[regmut_info[5]] = [regmut_info]
                    
                regmotifs_already_added = []
                for motif_info in motifs_info:
                    '''when more than one instance of a single motif is overlapping the same mutation just count it once.'''
                    if motif_info[9] in regmotifs_already_added:
                        continue
                    regmotifs_already_added.append(motif_info[9])
                    
                    '''mutations_in_cohorts list contains duplicate mutated motifs since it collects muts from all cohorts
                    and a single motif may be mutated in several cohorts but they all have the same info
                    except col4 which contains p-values and that is cohort specific. Therefore, col4 is replaced with cancer type and sample ID'''
                    motif_info[3] = regmut_info[5]+"#"+regmut_info[8]
                    if motif_info in mutated_motifs:
                        continue
                    '''
                    Append sample ID to each TF_motif in order to get the list of samples that have a particular
                    TF-Motif mutated. 
                    '''
                    try:
                        samples_regmuts_per_TFmotif[motif_info[9]].append(regmut_info[8])
                    except KeyError:
                        samples_regmuts_per_TFmotif[motif_info[9]]=[regmut_info[8]]
                    try:
                        samples_regmuts_per_TFmotif_allele[motif_info[9]+'#'+motif_info[4]].append(regmut_info[8])
                    except KeyError:
                        samples_regmuts_per_TFmotif_allele[motif_info[9]+'#'+motif_info[4]]=[regmut_info[8]]
                    try:
                        samples_regmuts_per_TFmotif_position[motif_info[9]+'#'+motif_info[5]].append(regmut_info[8])
                    except KeyError:
                        samples_regmuts_per_TFmotif_position[motif_info[9]+'#'+motif_info[5]]=[regmut_info[8]]
                    try:
                        samples_regmuts_per_TFmotif_position_allle[motif_info[9]+'#'+motif_info[5]+'#'+motif_info[4]].append(regmut_info[8])
                    except KeyError:
                        samples_regmuts_per_TFmotif_position_allle[motif_info[9]+'#'+motif_info[5]+'#'+motif_info[4]]=[regmut_info[8]]
                    
                    '''Append sample ID to each ChromHMM state in order to get enrichment of mutations per chromatin states per element'''
                    try:
                        if motif_info not in motifs_per_cancer_type[regmut_info[5]]:
                            motifs_per_cancer_type[regmut_info[5]].append(motif_info)
                            mutated_motifs.append(motif_info)
                    except KeyError:
                        motifs_per_cancer_type[regmut_info[5]] = [motif_info]
                        mutated_motifs.append(motif_info)
                            
                        
            for cancer_type in muts_per_cancer_type:
                cols_dict['Cancer-Types:#RegMuts'][cancer_type] = 0
                samples_per_cancer_type = [] 
                for mut in  muts_per_cancer_type[cancer_type]:
                    cols_dict['Cancer-Types:#RegMuts'][cancer_type]+=1
                    if mut[8] not in samples_per_cancer_type:
                        samples_per_cancer_type.append(mut[8])
                    
                cols_dict['Cancer-Types:#Samples(RegMuts)'][cancer_type]=len(samples_per_cancer_type)
            
            cols_dict['Mutated-Moitfs'] = ['#'.join(x) for x in mutated_motifs]
            '''Process All-Mutations in the element (functional and others)'''
            mutations_info = [x.strip().split('#') for x in l[11].split(',')]
            cancer_types = []
            cancer_types_samples = {}
            cancer_types_samples_count = {}
            
            mut_start_positions = []
            mut_end_positions = []
            
            samples_per_chromhmm = {}
            samples_per_allele_mut = {}
            
            for mut_info in mutations_info:
                try:
                    mut_start_positions.append(int(mut_info[0].split(':')[1].split('-')[0]))
                    mut_end_positions.append(int(mut_info[0].split(':')[1].split('-')[1]))
                except ValueError:
                    pass
                if mut_info[5] not in cols_dict['SamplesMuts']:
                    cols_dict['SamplesMuts'].append(mut_info[5])
                cancer_types.append(mut_info[2])
                try:
                    if mut_info[5] not in cancer_types_samples[mut_info[2]]:
                        cancer_types_samples[mut_info[2]].append(mut_info[5])
                        cancer_types_samples_count[mut_info[2]]+=1
                except KeyError:
                    cancer_types_samples[mut_info[2]] = [mut_info[5]]
                    cancer_types_samples_count[mut_info[2]] = 1
                
                '''the last column in the mutation info contains annotations'''
                state = "NA"
                if 'ChromHMM:' in mut_info[-1]:
                    for x in mut_info[-1].split('|'):
                        if x.startswith("ChromHMM:"):
                            state = x 
                try:
                    samples_per_chromhmm[state].append(mut_info[5])
                except KeyError:
                    samples_per_chromhmm[state] = [mut_info[5]]
                '''Get samples per allele: e.g C>A:sample1,sample2...'''
                try:
                    samples_per_allele_mut[mut_info[1]].append(mut_info[5])
                except:
                    samples_per_allele_mut[mut_info[1]] = [mut_info[5]]
            try:
                cols_dict['start'] = str(sorted(mut_start_positions)[0])#l[1]
                cols_dict['end'] = str(sorted(mut_end_positions)[-1])#l[2]
                cols_dict['Position'] = '{}:{}-{}'.format(l[0],cols_dict['start'],cols_dict['end'])
            except IndexError:
                #if no mutation was found in the region, skip the line and move on to the next region
                continue
            cols_dict['Muts'] = ['#'.join(x) for x in mutations_info]
            cols_dict['#Muts'] = len(mutations_info)
            cols_dict['#Samples'] = len(cols_dict['SamplesMuts'])
            c = Counter(cancer_types)
            cols_dict['Cancer-Types:#Muts'] = [x+":"+str(c[x]) for x in c]
            cols_dict['Cancer-Types:#Samples'] = [x+":"+str(cancer_types_samples_count[x]) for x in cancer_types_samples_count]
            
            '''Process Mutated-Motifs All (not only those that are significant based on regulatory mutations)'''
            mutated_motifsall = [x.strip().split('#') for x in l[12].split(',')]
            samples_per_tf_motifsall = {}
            samples_per_tf_position_motifsall = {}
            samples_per_tf_allele_motifsall = {}
            samples_per_tf_per_position_per_allele_motifsall = {}
            samples_per_chromhmm_motifsall = {}
            motifs_already_added = []
            muts_already_added = []
            for motif_info in mutated_motifsall:
                '''avoid recounting the same mutation in the same motif. 
                This may happen due to multiple instances of the same motif in the mutation position'''
                if '#'.join(motif_info[14:18])+"#"+motif_info[8] in motifs_already_added:
                    continue
                motifs_already_added.append('#'.join(motif_info[14:18])+"#"+motif_info[8])
                '''Get samples per TF-motif'''
                try:
                    samples_per_tf_motifsall[motif_info[17]].append(motif_info[8])
                except KeyError:
                    samples_per_tf_motifsall[motif_info[17]] = [motif_info[8]]
                '''Get samples per TF Motif position; to identify the number of muts per motif position'''
                try:
                    samples_per_tf_position_motifsall[motif_info[17]+"#"+motif_info[13]].append(motif_info[8])
                except KeyError:
                    samples_per_tf_position_motifsall[motif_info[17]+"#"+motif_info[13]] = [motif_info[8]]
                '''Get samples per allele; to identify the number of mutation types per allele such as CTCF#A>C:3'''
                try:
                    samples_per_tf_allele_motifsall[motif_info[17]+"#"+motif_info[12]].append(motif_info[8])
                except KeyError:
                    samples_per_tf_allele_motifsall[motif_info[17]+"#"+motif_info[12]] = [motif_info[8]]
                try:
                    samples_per_tf_per_position_per_allele_motifsall[motif_info[17]+"#"+motif_info[13]+"#"+motif_info[12]].append(motif_info[8])
                except KeyError:
                    samples_per_tf_per_position_per_allele_motifsall[motif_info[17]+"#"+motif_info[13]+"#"+motif_info[12]] = [motif_info[8]]
            
                '''Almost all motifs that overlap with single mutation are expected to have the same chromHMM state.
                Therefore take state of the first motif'''
                if '#'.join(motif_info[0:3])+"#"+motif_info[8] in muts_already_added:
                    continue
                muts_already_added.append('#'.join(motif_info[0:3])+"#"+motif_info[8])
                try:
                    samples_per_chromhmm_motifsall[motif_info[22]].append(motif_info[8])
                except KeyError:
                    samples_per_chromhmm_motifsall[motif_info[22]] = [motif_info[8]]
            
            '''Get counts from dicts'''
            samples_dicts = {'RegMuts-ChromHMM':samples_regmuts_per_chromhmm,#Regmuts
                             'RegMuts-Motifs':samples_regmuts_per_TFmotif, 'RegMuts-MotifPositions':samples_regmuts_per_TFmotif_position, 
                             'RegMuts-MotifAllele':samples_regmuts_per_TFmotif_allele, 'RegMuts-MotifPositions-Allele':samples_regmuts_per_TFmotif_position_allle, 
                             'Muts-ChromHMM':samples_per_chromhmm, 'Muts-Allele':samples_per_allele_mut,#All muts 
                             'Muts-Motifs':samples_per_tf_motifsall, 'Muts-MotifPositions':samples_per_tf_position_motifsall, 'Muts-MotifAllele':samples_per_tf_allele_motifsall,#MutMotifsAll 
                             'Muts-MotifPositions-Allele':samples_per_tf_per_position_per_allele_motifsall, 'Muts-Motifs-ChromHMM':samples_per_chromhmm_motifsall}
            for k in sorted(samples_dicts.keys()):
                mutcounts = {}
                samplecounts = {}
                for x in samples_dicts[k]:
                    mutcounts[x]=len(samples_dicts[k][x])
                    samplecounts[x]=len(set(samples_dicts[k][x]))
                
                c = ','.join([x.replace('ChromHMM:','')+":"+str(mutcounts[x]) for x in sorted(mutcounts, key=mutcounts.get, reverse=True)]) 
                cols_dict['StatsMuts'].append(k+'['+c+']')
                
                c_unique = ','.join([x.replace('ChromHMM:','')+":"+str(samplecounts[x]) for x in sorted(samplecounts, key=samplecounts.get, reverse=True)]) 
                cols_dict['StatsSamples'].append(k+'['+c_unique+']')
            
            '''Report CancerType:ChromatinType:number_of_times'''
                    
            #instead of writing the dict put them in a list to keep the columns order 
            cols_to_write = [cols_dict['chr'], cols_dict['start'], cols_dict['end'], cols_dict['Position'], 
                             ','.join(cols_dict['Cohorts']), cols_dict['#Cohorts'], cols_dict['Score'], cols_dict['FDR'], 
                             cols_dict['#RegMuts'], cols_dict['#Samples(RegMuts)'], 
                             ','.join([x+":"+str(cols_dict['Cancer-Types:#RegMuts'][x]) for x in cols_dict['Cancer-Types:#RegMuts']]), 
                             ','.join([x+":"+str(cols_dict['Cancer-Types:#Samples(RegMuts)'][x]) for x in cols_dict['Cancer-Types:#Samples(RegMuts)']]),
                             cols_dict['#Muts'], cols_dict['#Samples'], ','.join(cols_dict['Cancer-Types:#Muts']), ','.join(cols_dict['Cancer-Types:#Samples']),
                             #'Nearby-Genes(Downstream/Upstream:Distance;COSMIC;KEGG;PCAWG)'
                             ','.join(cols_dict['StatsMuts']),','.join(cols_dict['StatsSamples']),
                             ','.join(cols_dict['RegMuts']), ','.join(cols_dict['Muts']), ','.join(cols_dict['Mutated-Moitfs']), ','.join(cols_dict['Max-RegMotif']),
                             ','.join(cols_dict['SamplesMuts'])
                             ]
            aggregated_lines.append(cols_to_write)
            
            summaries_dict['#Elements'] +=1 
            summaries_dict['#RegMuts'] += cols_dict['#RegMuts']
            summaries_dict['#Samples(RegMuts)'] += cols_dict['#Samples(RegMuts)']
            summaries_dict['#Muts'] += cols_dict['#Muts']
            summaries_dict['#Samples'] += cols_dict['#Samples']
            
    return aggregated_lines, summaries_dict

def generate_genesets_genes_dict(geneset_files):
    
    cosmic_genes_file = "analysis/cancer_gene_census.csv"
    kegg_pathways_file = "analysis/kegg_pathways_fromdb_madeAgenesetPerPathway.gmt"
    pcawg_drivers_file = "analysis/PCAWG_cancer_drivers_fulllist.txt"
    
    genesets_genes_dict = {'KCP':[], 'COSMIC': [], 'PCD': []}
    with open(cosmic_genes_file, 'r') as cosmic_genes_ifile:
        lines = cosmic_genes_ifile.readlines()
        for l in lines[1::]:
            sl = l.strip().split(',')
            if sl[0]!="":
                genesets_genes_dict['COSMIC'].append(sl[0])
            for s in sl:
                if s.startswith('ENSG'):
                    genesets_genes_dict['COSMIC'].append(s)
    
    with open(kegg_pathways_file, 'r') as kegg_pathways_ifile:
        lines = kegg_pathways_ifile.readlines()
        for l in lines:
            sl = l.split('\t')
            if sl[0]=="path:hsa05200":
                genesets_genes_dict['KCP'] = sl[3::]
    
    with open(pcawg_drivers_file, 'r') as pcawg_drivers_ifile:
        lines = pcawg_drivers_ifile.readlines()
        for l in lines:
            sl = l.strip().split('\t')
            if sl[4]!="":
                genesets_genes_dict['PCD'].extend(sl[4].split(';'))
            for s in sl[1].split('::'):
                if s.startswith('ENSG'):
                    genesets_genes_dict['PCD'].append(s)
    
    return genesets_genes_dict

def get_enriched_gene_geneset(regions_genes_dict, genesets_genes_dict):
    enriched_genesets_dict = {}
    enriched_genesets_dict_overall = {}
    genes_all = {}
    genes_all_per_side = {}
    for region in regions_genes_dict.keys():
        for i, gene_info in enumerate(regions_genes_dict[region]):
            gene_name = gene_info.split('::')[0]
            gene_id = gene_info.split('::')[1]
            gene_place_dist = gene_info.split('::')[2].split('O')[0].split('D')[0].split('U')[0] + gene_info.split('::')[2][len(gene_info.split('::')[2].split('O')[0].split('D')[0].split('U')[0])]
            
            try:
                if gene_name not in genes_all[gene_place_dist[-1]]:
                    genes_all[gene_place_dist[-1]].append(gene_name)
            except KeyError:
                genes_all[gene_place_dist[-1]] = [gene_name]
            
            try:
                if gene_name not in genes_all_per_side[gene_place_dist]:
                    genes_all_per_side[gene_place_dist].append(gene_name)
            except KeyError:
                genes_all_per_side[gene_place_dist] = [gene_name]
                
            for geneset in sorted(genesets_genes_dict.keys()):
                if (gene_name in genesets_genes_dict[geneset] or 
                    gene_id in genesets_genes_dict[geneset] or
                    gene_id.split('.')[0] in genesets_genes_dict[geneset]):
                    regions_genes_dict[region][i]+="::"+geneset
                    
                    try:
                        if gene_name not in enriched_genesets_dict[geneset+gene_place_dist]:
                            enriched_genesets_dict[geneset+gene_place_dist].append(gene_name)
                    except KeyError:
                        enriched_genesets_dict[geneset+gene_place_dist] = [gene_name]
                    
                    try:
                        if gene_name not in enriched_genesets_dict_overall[geneset]:
                            enriched_genesets_dict_overall[geneset].append(gene_name)
                    except KeyError:
                        enriched_genesets_dict_overall[geneset] = [gene_name]
    
    print '\t'.join([r + ':' + str(len(genes_all[r])) for r in genes_all.keys()])
    print '\t'.join([r + ':' + str(len(genes_all_per_side[r])) for r in genes_all_per_side.keys()])
    print '\t'.join([r + ':' + str(len(enriched_genesets_dict_overall[r])) for r in enriched_genesets_dict_overall.keys()])
    print '\t'.join([r + ':' + str(len(enriched_genesets_dict[r])) for r in enriched_genesets_dict.keys()]) 
    return regions_genes_dict, genes_all, genes_all_per_side, enriched_genesets_dict_overall, enriched_genesets_dict

def write_aggregated_lines_to_outfile(aggregated_lines, cols_to_write, 
                                      regions_genes_dict, summary_info_to_write, summary_dicts_to_write,
                                      region_types_dict,
                                      aggregated_output_file):
    with open(aggregated_output_file, 'w') as output_ofile:
        output_ofile.write("Summary Info\n")
        
        for dict_name in sorted(summary_info_to_write.keys()):
            output_ofile.write(dict_name + '\t' + '\t'.join([r + ':' + str((summary_info_to_write[dict_name][r])) for r in summary_info_to_write[dict_name].keys()]) + '\n')
        
        for dict_name in sorted(summary_dicts_to_write.keys()):
            output_ofile.write(dict_name + '\t' + '\t'.join([r + ':' + str(len(summary_dicts_to_write[dict_name][r])) for r in summary_dicts_to_write[dict_name].keys()]) + '\n')
        
        output_ofile.write('\t'.join(cols_to_write) + '\t' + 'Feature_types'+'\n')
        
        for l in aggregated_lines:
            genes = "None"
            try:
                genes = ','.join(regions_genes_dict[l[3]])
            except KeyError:
                pass
            
            
            feature_type = 'intergenic'
            feature_types  = ['intergenic']
            try:
                feature_to_select_index = 0
                n = 0
                feature_types = [x[0]+":"+str(x[1]) for x in region_types_dict[l[3]]]
                for i, feature in enumerate(region_types_dict[l[3]]):#find index of the feature that has the largest num of muts
                    if feature[0] != 'gene':
                        if feature[1]>n:
                            feature_to_select_index = i
                            n = feature[1]
                        
                feature_type = region_types_dict[l[3]][feature_to_select_index][0]
                if feature_type=='gene':
                    feature_type = 'intronic'
            except KeyError:
                pass
            
            output_ofile.write('\t'.join([str(x) for x in l]) + '\t' + genes + '\t' + feature_type + '\t' + ','.join(feature_types)+'\n')
            
    return aggregated_output_file

def overlaps(ms, me, rs, re):
    if ms>=rs and ms<=re:#if m starts in r then it overlaps
        return True
    elif me>=rs and me<=re:#if m ends in r then it overlaps
        return True
    elif ms<=rs and me >= re:#r is within m
        return True
    return False


def get_region_type(aggregated_lines, genes_segments_input_file, gene_types_to_consider, gene_status_to_consider, feature_types_to_consider, muts_col_index=19):
    region_temp_file = 'aggregated_lines_temp_regions'
    with open(region_temp_file, 'w') as ofile:
        for reg in aggregated_lines:
            ofile.write('\t'.join(reg[0:4]) + '\t' + reg[muts_col_index] + '\n')
    BedTool(genes_segments_input_file).intersect(BedTool(region_temp_file), wo=True).saveas(region_temp_file+'genesfeatures')
    
    region_types_dict = {}#store feature type and number of muts in the feature type
    with open(region_temp_file+'genesfeatures', 'r') as ifile:
        l = ifile.readline().strip().split('\t')
        while len(l)>3:
            if l[8] in gene_types_to_consider and l[9] in gene_status_to_consider and l[3] in feature_types_to_consider:
                #check the number of muts that have start larger than the start of this feature
                #keep the feature that have the largest number of muts
                #for each positionID append its overlapping feature types
                num_muts_in_this_feature = 0
                for mut in l[14].split(','):
                    mut_start = int(mut.split('#')[0].split(':')[1].split('-')[0])
                    mut_end = int(mut.split('#')[0].split(':')[1].split('-')[1]) 
                    if overlaps(mut_start, mut_end, int(l[1]), int(l[2])):
                    #if mut_start>int(l[11]) or (mut_start<int(l[11]) and mut_end>int(l[11])):#the mut is within the region
                        num_muts_in_this_feature+=1
                #append feature type and the number of muts in it for each regionID
                try:
                    region_types_dict[l[13]].append([l[3], num_muts_in_this_feature])
                except KeyError:
                    region_types_dict[l[13]] = [[l[3], num_muts_in_this_feature]]
                
            l = ifile.readline().strip().split('\t')
        
    #os.remove(region_temp_file)
    #os.remove(region_temp_file+'genesfeatures')
    
    return region_types_dict

def get_features_from_gencode(gencode_input_file, gencode_output_file):
    
    with open(gencode_input_file, 'r') as ifile, open(gencode_output_file, 'w') as ofile:
            l = ifile.readline()
            while l:
                if l.startswith('#'):
                    print 'Skipping: ', l
                    l = ifile.readline()
                    continue
                sl = l.strip().split('\t')
                if len(sl)>8:
                    info_dict = {}
                    for info in sl[8].split(';'):
                        if '=' in info:
                            info_dict[info.split('=')[0]] = info.split('=')[1]
                            
                    #if info_dict['gene_status'] == 'KNOWN':
                    try:
                        ol = [sl[0], sl[3], sl[4], sl[2], sl[1], sl[6], info_dict['ID'], info_dict['gene_id']+"::"+info_dict['gene_name'], info_dict['gene_type'], info_dict['gene_status']]                            
                        ofile.write('\t'.join(ol) + '\n')
                        #if it was a gene then write a promoter region to the file too
                        if sl[2]=='gene':
                            promoter_start = int(sl[3])-2000
                            promoter_end = int(sl[3])-1
                            if sl[6]=='-':
                                promoter_start = int(sl[4])+1
                                promoter_end = int(sl[4])+2000
                            if promoter_start<1:
                                promoter_start = 1
                            if promoter_end<1:
                                promoter_end = 1
                            ol = [sl[0], str(promoter_start), str(promoter_end), 'proximal_promoter', sl[1], sl[6], 'proximal_promoter:'+info_dict['ID'], info_dict['gene_id']+"::"+info_dict['gene_name'], info_dict['gene_type'], info_dict['gene_status']]                            
                            ofile.write('\t'.join(ol) + '\n')
                        
                    except KeyError:
                        print "Key not found: ", sl, info_dict
                else:
                    print "Length<8: ", l, sl            
                l = ifile.readline()
                if l=="":
                    break
                
    return gencode_output_file
    
def getSigElements(generated_sig_merged_element_files, n, max_dist, window):
    
    upstream=True
    downstream=True
    overlapping = True
    
    ext = ""
    try:
        ext = generated_sig_merged_element_files[0].split('/')[-1].split('.bed9')[1].replace('groupedbymutwithmotifinfo_','').replace('_statspvalues', '')
    except IndexError:
        print "error: ", generated_sig_merged_element_files
        sys.exit()

    aggregated_output_file = 'analysis_exclVEP/combined{ext}_merged_intersectedmuts_grouped_aggregated{n}{up}{dw}maxdist{max_dist}kb_within{window}kb.tsv'.format(ext=ext, n=n, up="Up", dw="Dw", max_dist=max_dist/1000, window=window/1000)
    if os.path.exists(aggregated_output_file):
        return aggregated_output_file
    
    print generated_sig_merged_element_files
    annotated_motifs = 'mutations_files/obsann22May2017_exclVEP.bed9'
    tracks_dir = 'datafiles/chromatin_marks_all_cells_onlynarrowpeaks'
    observed_mutations_all = 'mutations_files/obsagr22May2017_exclVEP.bed9'#_notInExonsProteinCodingProcessedTranscriptIG.bed9'#'mutations_files/observedunique.bed9'
    #regions_input_file = 'analysis/combined_onlysig_merged_intersectedmuts_grouped_recurrent.col12'
    combined_mut_grouped_file = 'analysis_exclVEP/combined{ext}_merged_intersectedmuts_grouped_recurrent.col12'.format(ext=ext)
    if not os.path.exists(combined_mut_grouped_file):
        combined_file_all = combined_mut_grouped_file+'_temp'
        with open(combined_file_all, 'w') as regions_input_ofile:
            for cohort_sigregions_file in generated_sig_merged_element_files:
                cohort_name = cohort_sigregions_file.split('/')[-1].split('_')[0]
                with open(cohort_sigregions_file, 'r') as cohort_sigregions_ifile:
                    l = cohort_sigregions_ifile.readline().strip().split('\t')
                    while l and len(l)>10:
                        regions_input_ofile.write('\t'.join(l[0:3]) + '\t' + cohort_name + '\t' + '~'.join([x.replace(',', '|') for x in l]) + '\n')
                        l = cohort_sigregions_ifile.readline().strip().split('\t')
        nbp_to_extend = 200
        combined_file_all_merged = combined_mut_grouped_file+'_temp_merged'
        awk_stmt = ("""awk 'BEGIN{{FS=OFS="\t"}}{{if(($3-$2)<{nbp_to_extend}){{s={nbp_to_extend}-($3-$2); $2=$2-int(s/2); $3=$3+int(s/2);}}; print $0}}' {combined_file_all} | sort -k1,1n -k2,2n -k3,3n | mergeBed -i stdin -c 4,4,5 -o count_distinct,collapse,collapse | awk 'BEGIN{{FS=OFS="\t"}}{{gsub("23","X", $1); gsub("24","Y", $1); print "chr"$0}}' > {combined_file_all_merged} 
                    """).format(combined_file_all=combined_file_all, nbp_to_extend = nbp_to_extend, combined_file_all_merged=combined_file_all_merged)
        '''awk_stmt = ("""sort -k1,1n -k2,2n -k3,3n -i {combined_file_all} | mergeBed -i stdin -c 4,4,5 -o count_distinct,collapse,collapse |
                    awk 'BEGIN{{FS=OFS="\t"}}{{gsub("23","X", $1); gsub("24","Y", $1); if(($3-$2)<{nbp_to_extend}){{s={nbp_to_extend}-($3-$2); $2=$2-int(s/2); $3=$3+int(s/2);}}; print "chr"$0}}' > {combined_file_all_merged} 
                    """).format(combined_file_all=combined_file_all, nbp_to_extend = nbp_to_extend, combined_file_all_merged=combined_file_all_merged)
        '''
        #print awk_stmt
        os.system(awk_stmt)
        muts_overlapping_combined_file_all = combined_file_all_merged+"_muts"
        BedTool(observed_mutations_all).intersect(BedTool(combined_file_all_merged)).saveas(muts_overlapping_combined_file_all)
        #annotat the overlapping muts
        muts_overlapping_combined_file_all_annotated = muts_overlapping_combined_file_all+'_annotated'
        muts_overlapping_combined_file_all_annotated = get_annotated_muts(muts_input_file=muts_overlapping_combined_file_all, tracks_dir=tracks_dir, muts_out=muts_overlapping_combined_file_all_annotated, filter_on_dnase1_or_tf_peak=False)
        #sort and merge the file
        combined_mut_grouped_file_with_annotated_muts = combined_mut_grouped_file + "_withannotatedmuts"
        print "Combining results"
        awk_stmt = ("""intersectBed -wo -loj -a {combined_file_all_merged} -b {observed_mutations_all} | 
                    awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$3,$4,$5,$10">"$11,$12,$15,$6,$7":"$8"-"$9"#"$10">"$11"#"$12"#"$13"#"$14"#"$15"#"$16}}' | 
                    groupBy -g 1-5 -c 8,8,8,6,7,9,10 -o count,count_distinct,collapse,collapse,collapse,distinct,collapse > {combined_mut_grouped_file_with_annotated_muts}""".format(
                    combined_file_all_merged=combined_file_all_merged, observed_mutations_all=muts_overlapping_combined_file_all_annotated, combined_mut_grouped_file_with_annotated_muts=combined_mut_grouped_file_with_annotated_muts))
        os.system(awk_stmt)#awk '$7>1'
        
        #get all mutated motifs in the extended element
        #combined_mut_grouped_file_with_annotated_muts_with_motifs = combined_mut_grouped_file + "_withannotatedmuts_motifs"
        awk_stmt = ("""intersectBed -wo -loj -a {combined_mut_grouped_file_with_annotated_muts} -b {annotated_motifs} | 
                    awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,{motif_cols}}}' | 
                    groupBy -g 1-12 -c 13 -o collapse > {combined_mut_grouped_file_with_annotated_muts_with_motifs}""".format(
                    combined_mut_grouped_file_with_annotated_muts=combined_mut_grouped_file_with_annotated_muts, annotated_motifs=annotated_motifs, 
                    combined_mut_grouped_file_with_annotated_muts_with_motifs=combined_mut_grouped_file,
                    motif_cols = '"#"'.join(["$"+str(x) for x in range(13,45)])))#motif cols are starting from col12 and end in col44
        os.system(awk_stmt)#awk '$7>1'
        
        '''
        os.remove(combined_file_all)
        os.remove(combined_file_all_merged)
        os.remove(muts_overlapping_combined_file_all)
        os.remove(muts_overlapping_combined_file_all_annotated)
        '''
        
        #os.remove(regions_input_file+'_temp')
    print "Aggregating for final results"
    aggregated_lines, summaries_dict = aggregate_results(combined_mut_grouped_file)
    
    num_muts_per_sample_dict = get_number_of_mutations_per_sample_list_and_write_to_file(mutations_file=observed_mutations_all, numberofmutationspersample_output_file=observed_mutations_all+"numbermutspersample.txt", index_sample_ids=8)
    aggregated_lines = calculate_p_value_motifregions(aggregated_lines, num_muts_per_sample_dict, index_mutation_frequency=12, index_sample_ids=-1, index_elment_start_coordinate=1, index_elment_stop_coordinate=2, genome_size=3000000000.0, total_number_of_regions_tested=833999)#len(aggregated_lines)) 
    
    #aggregated_lines = compute_fdr_per_element(aggregated_lines, mutations_input_file=observed_mutations_all, sample_ids_index_muts_file=8, num_muts_index = 12, index_sample_ids=-1, index_elment_start_coordinate=1, index_elment_stop_coordinate=2, genome_size=3000000000.0)
    '''
    #Get pvalue for each Element
    if not os.path.exists(annotated_mutations_statcalc_output_file):
        print "Calculating p-values"
        print "getting mut frequency per sample"
        sample_id_and_number_of_mutations_per_sample_dict = get_number_of_mutations_per_sample_list_and_write_to_file(Mutations_dir_list, "", index_sample_ids=8)
        number_of_elements_tested = file_len(annotated_mutations_final_output_file_scored_merged)
        if header:
            number_of_elements_tested-=1
            calculate_p_value_motifregions(annotated_mutations_final_output_file_scored_merged, sample_id_and_number_of_mutations_per_sample_dict, mutated_regions_pval_outfile=annotated_mutations_statcalc_output_file, index_mutation_frequency=5, index_sample_ids=4, index_elment_start_coordinate=1, index_elment_stop_coordinate=2, genome_size=3100000000.0, total_number_tested_regions=number_of_elements_tested)
    '''
    cols_to_write = ['chr', 'start', 'end', 'Position', 'Cohorts', '#Cohorts', 'Score', 'FDR', 
                     '#RegMuts', '#Samples(RegMuts)', 'Cancer-Types:#RegMuts', 'Cancer-Types:#Samples(RegMuts)',
                     '#Muts', '#Samples', 'Cancer-Types:#Muts', 'Cancer-Types:#Samples','StatsMuts', 'StatsSamples',
                     'RegMuts','Muts', 'Mutated-Moitfs', 'Max-RegMotif', 
                     'SamplesMuts', 'ElementPval', 'ELementFDR', 
                     'Nearby-Genes(Name::ID::O|U|Ddistance::COSMIC|KCP|PCD)',
                     'Feature_type'
                     ]
    
    chr_lengths_file = 'datafiles/chr_lengths_hg19.txt'
    genes_input_file = "datafiles/GeneExp/gencode.v19.annotation.gff3_onlygenes.bed"
    genocode_genes_segments_input_file = "datafiles/GeneExp/gencode.v19.annotation.gff3"
    
    #gene_types_to_consider = ['protein_coding', 'lincRNA', 'miRNA', 'snRNA', 'snoRNA', 'rRNA', 'Mt_tRNA', 'Mt_rRNA', 'antisense', 'sense_intronic', 'sense_overlapping', '3prime_overlapping_ncrna']
    gene_types_to_consider = ['protein_coding', 'lincRNA',
                              'IG_V_gene', 'IG_C_gene', 'IG_J_gene', 'IG_D_gene', 
                              'TR_V_gene', 'TR_C_gene', 'TR_J_gene', 'TR_D_gene', 
                              'processed_transcript']
    gene_status_to_consider = ['KNOWN']
    
    #aggregated_output_file = 'analysis/combined{ext}_merged_intersectedmuts_grouped_aggregated{n}{up}{dw}maxdist{max_dist}kb_within{window}kb.tsv'.format(ext=ext, n=n, up="Up", dw="Dw", max_dist=max_dist/1000, window=window/1000)
    
    
    chr_lengths = get_chr_lengths(chr_lengths_file)
    
    extended_output_file = aggregated_output_file+"_extendedtemp"
    extended_output_file = generate_extended_regions(regions=aggregated_lines, extended_output_file=extended_output_file, chr_lengths=chr_lengths, window=window)
    
    regions_genes_dict = get_nearby_genes(regions_input_file=extended_output_file, 
                        genes_input_file = genes_input_file, 
                        gene_types_to_consider = gene_types_to_consider, gene_status_to_consider = gene_status_to_consider,
                        n=n, upstream=upstream, downstream=downstream, overlapping = overlapping, max_dist = max_dist,
                        genesets_to_check_for_overlap_input_files=[])
    os.remove(extended_output_file)
    
    geneset_files = ['']
    genesets_genes_dict = generate_genesets_genes_dict(geneset_files)
    enrichment_regions_genes_dict, genes_all, genes_all_per_side, enriched_genesets_dict_overall, enriched_genesets_dict = get_enriched_gene_geneset(
                                                                                                                        regions_genes_dict, genesets_genes_dict)
    summary_dicts_to_write = {"All genes:": genes_all, "All genes per dir:": genes_all_per_side ,"Enriched genes:": enriched_genesets_dict_overall, "Enriched genes per dir:": enriched_genesets_dict}
    summary_info_to_write = {'Element Info': summaries_dict}
    
    gencode_output_file="datafiles/GeneExp/gencode.v19.annotation.gff3_extractedinfo"
    if not os.path.exists(gencode_output_file):
        get_features_from_gencode(gencode_input_file="datafiles/GeneExp/gencode.v19.annotation.gff3", gencode_output_file=gencode_output_file)
    gene_types_to_consider = ['protein_coding',
                              'IG_V_gene', 'IG_C_gene', 'IG_J_gene', 'IG_D_gene', 
                              'TR_V_gene', 'TR_C_gene', 'TR_J_gene', 'TR_D_gene', 
                              'processed_transcript']
    
    region_types_dict = get_region_type(aggregated_lines=aggregated_lines, genes_segments_input_file=gencode_output_file, 
                                        gene_types_to_consider=gene_types_to_consider, gene_status_to_consider=gene_status_to_consider,
                                        feature_types_to_consider=['CDS', 'UTR','proximal_promoter', 'gene','start_codon', 'stop_codon'])
    
    write_aggregated_lines_to_outfile(aggregated_lines, cols_to_write, 
                                      enrichment_regions_genes_dict, summary_info_to_write, summary_dicts_to_write, 
                                      region_types_dict,
                                      aggregated_output_file)
    
    return aggregated_output_file

def combine_sig_TFs(sig_tfs_files, tf_label='TFs'):
    ext = ""
    try:
        ext = sig_tfs_files[0].split('/')[-1].split('.bed9')[1]
    except IndexError:
        print "error: ", generated_sig_merged_element_files
        sys.exit()
    
    aggregated_output_file = 'analysis_exclVEP/combined{ext}.tsv'.format(ext=ext)
    if os.path.exists(aggregated_output_file):
        return aggregated_output_file
    header_cols = [tf_label, 'P-Val', 'FDR', '#Mutated Motifs', 'Mean #Mutated Motifs in Simulated Sets','#Mutated Motifs in Simulated Sets', 'Cohorts']
    with open(aggregated_output_file, 'w') as ofile:
        ofile.write('\t'.join(header_cols) + '\n')
        for sig_tf_file in sig_tfs_files:
            cohort_name = sig_tf_file.split('/')[-1].split('_')[0]
            with open(sig_tf_file, 'r') as ifile:
                l = ifile.readline().strip().split('\t')
                while l and len(l)>4:
                    l = [l[0],l[1],l[2],l[3], str(np.mean([int(x) for x in l[5].split(',')])), l[5]]
                    ofile.write('\t'.join(l) + '\t'+ cohort_name +'\n')
                    l = ifile.readline().strip().split('\t')
    print 'results are in : ', aggregated_output_file           
    return aggregated_output_file

def get_gene_enrichments(elements_input_file, elements_output_file, header_lines_to_skip=6, skip_exon_elements=True):
    elements_input = pd.read_table(elements_input_file, sep='\t', skiprows=header_lines_to_skip, header=0)
    elements_input_filtered = elements_input[(elements_input['#Samples(RegMuts)']>1)]
    genes_dict = {}
    for i,element in elements_input_filtered.iterrows():
        if skip_exon_elements:
            if element['Feature_type']=='CDS':
                continue
        regmuts = []
        muts = []
        for regmut in element['RegMuts'].split(','):
            regmuts.append(regmut.split('#')[8])
        for mut in element['Muts'].split(','):
            muts.append(mut.split('#')[5])
        for g in element['Nearby-Genes(Name::ID::O|U|Ddistance::COSMIC|KCP|PCD)'].split(','):
            element_info = element['Position']+'::'+str(element['FDR'])+'::'+element['Feature_type'].replace(',','|')+'::'+g
            gene_name = g.split('::')[0]
            gene_id = "None"
            if len(g.split('::'))>2:
                gene_id = g.split('::')[1]
            try:
                genes_dict[gene_name]['#RegMuts']+=len(regmuts)
                genes_dict[gene_name]['#Muts']+=len(muts)
                genes_dict[gene_name]['RegMutSamples'].update(regmuts)
                genes_dict[gene_name]['MutSamples'].update(muts)
                genes_dict[gene_name]['Elements'].append(element_info)
            except KeyError:
                genes_dict[gene_name] = {'Gene_ID':gene_id, '#RegMuts': len(regmuts), '#Muts': len(set(muts)), 'RegMutSamples':set(regmuts), 'MutSamples':set(muts), 'Elements':[element_info]}
    
    genes_samples = []
    samples = []
    with open(elements_output_file, 'w') as ofile:
        for g in sorted(genes_dict.keys()):
            gene = genes_dict[g]
            gene['#RegMutSamples']=len(gene['RegMutSamples'])
            gene['#MutSamples']=len(gene['MutSamples'])
            
            ofile.write('\t'.join([g, gene['Gene_ID'], str(gene['#RegMuts']), str(gene['#RegMutSamples']), str(gene['#Muts']), str(gene['#MutSamples']),
                         str(len(gene['Elements'])), ','.join(gene['Elements']), ','.join(gene['RegMutSamples']), ','.join(gene['MutSamples'])]) + '\n')
            
            '''
            if gene['#RegMutSamples']>=10 and gene['#MutSamples']>=20:
                for sample in gene['MutSamples']:
                    samples.append(sample)
                    if sample in gene['RegMutSamples']:
                        genes_samples.append([sample, g, 'AMP', 'CNA', len(gene['RegMutSamples']), len(gene['MutSamples'])])
                        #samples_per_gene.write('\t'.join([sample, g, 'AMP', 'CNA', str(len(gene['RegMutSamples'])), str(len(gene['MutSamples']))])+'\n')
                    else:
                        genes_samples.append([sample, g, 'mut', 'TRUNC', len(gene['RegMutSamples']), len(gene['MutSamples'])])
                        #samples_per_gene.write('\t'.join([sample, g, 'mut', 'TRUNC', str(len(gene['RegMutSamples'])), str(len(gene['MutSamples']))])+'\n')
                for sample in gene['RegMutSamples']:
                    if sample not in gene['MutSamples']:
                        genes_samples.append([sample, g, 'AMP', 'CNA', len(gene['RegMutSamples']), len(gene['MutSamples'])])
                        samples.append(sample)
                        #samples_per_gene.write('\t'.join([sample, g, 'AMP', 'CNA', str(len(gene['RegMutSamples'])), str(len(gene['MutSamples']))])+'\n')
    genes_samples.sort(key=lambda x: x[5], reverse=True)
    total_numnber_samples = 2520
    
    for s in range(len(set(samples)), total_numnber_samples):
        genes_samples.append(['Sample'+str(s)])
    
    with open(elements_output_file+"_samples_per_gene", 'w') as samples_per_gene:
        for sample_gene in genes_samples:
            if len(sample_gene)==6:
                samples_per_gene.write('\t'.join([sample_gene[0], sample_gene[1]+':'+str(sample_gene[4])+':'+str(sample_gene[5]), sample_gene[2],sample_gene[3]]) + '\n')
            else:
                samples_per_gene.write(sample_gene[0] + '\n')
    '''
    return elements_output_file

def get_sample_pathways(calculated_p_value_sig_out_file, output_file, total_number_samples=2520):
    
    with open(calculated_p_value_sig_out_file, 'r') as ifile, open(output_file, 'w') as ofile:
        pathways = [x.strip().split('\t') for x in ifile.readlines()]
        all_samples = []
        pathways.sort(key=lambda x: int(x[5]), reverse=True)
        for pathway in pathways:
            if int(pathway[4])>200 and int(pathway[5])>200 and float(pathway[-1])<0.01:
                pathway[1] = pathway[1].replace(' signaling pathway','').replace(' ','_').replace('-','_') + ':' + ':'.join([pathway[4], pathway[5], '%.1E' % Decimal(pathway[-1])])
                for sample in pathway[7].split(','):
                    all_samples.append(sample)
                    if sample in pathway[6].split(','):
                        ofile.write('\t'.join([sample, pathway[1], 'AMP', 'CNA']) + '\n')
                    else:
                        ofile.write('\t'.join([sample, pathway[1], 'mut', 'TRUNC']) + '\n')
                
                for sample in pathway[6].split(','):
                    if sample not in pathway[7].split(','):
                        all_samples.append(sample)
                        ofile.write('\t'.join([sample, pathway[1], 'AMP', 'CNA']) + '\n')
                
        for s in range(len(set(all_samples)), total_number_samples):
            ofile.write('Sample'+str(s)+'\n')
        
    return


if __name__ == '__main__':
    
    cohort_names_input = sys.argv[1]#'meta_tumor_cohorts_v2_22May2017/cohorts_to_run_definedPCAWG'
    #generated_sig_merged_element_files, sig_tfs_files, sig_tfpos_files = process_cohorts(cohort_names_input)
    mutation_input_dir = 'mutations_cohorts_output_exclVEP'
    generated_sig_merged_element_files =  [mutation_input_dir+'/'+x for x in os.listdir(mutation_input_dir) if 'statspvalueslocalw25000onlysig0.05' in x]
    sig_tfs_files = [mutation_input_dir+'/'+x for x in os.listdir(mutation_input_dir) if 'sigTFs_0.05' in x]

    sig_tfpos_files = [mutation_input_dir+'/'+x for x in os.listdir(mutation_input_dir) if 'sigTFpos_0.05' in x]    
    print sig_tfpos_files
    n=int(sys.argv[2])#3
    max_dist = int(sys.argv[3])#500000
    window = int(sys.argv[4])#500000
    
    aggregated_output_file = getSigElements(generated_sig_merged_element_files, n, max_dist, window)
    combine_sig_TFs(sig_tfs_files)
    combine_sig_TFs(sig_tfpos_files, tf_label='TF Positions')
    
    geneset_input_file = 'analysis/kegg_pathways_fromdb_madeAgenesetPerPathway.gmt'
    elements_output_file = get_gene_enrichments(elements_input_file=aggregated_output_file, elements_output_file=aggregated_output_file+"_GenesInclCDS.tsv", skip_exon_elements=False)
    calculated_p_value_sig_out_file = find_overlap_genesets_genelist(geneset_input_file, elements_output_file, elements_output_file+'_pathways.tsv', total_number_of_genes_in_the_universe=20278, 
                                   min_number_of_genes_be_enriched_for_geneset_to_be_reported = 10, index_gene_name=0, index_gene_names_start=3, 
                                   keywords_to_filter_out_with=[], only_keep_the_sig_file = False, min_number_of_genes_in_geneset_to_consider_the_geneset = 10, header_line = False,
                                   sample_ids_given=True)
    #produce a list of all genes including exon elements
    #get_gene_enrichments(elements_input_file=aggregated_output_file, elements_output_file=aggregated_output_file+"_GenesInclExons.tsv", skip_exon_elements=False)
    
    #get_sample_pathways(calculated_p_value_sig_out_file, elements_output_file+'_pathwaySamples.txt')
    
    
    
    
    
