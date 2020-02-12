'''
Created on Dec 11, 2016

@author: husensofteng
'''
import sys, os
import pandas as pd
import numpy as np
import string
import psycopg2
import json
from multiprocessing import Pool
from psycopg2.extras import DictCursor
import time

from score_motifs_tissuepertable import open_connection, close_connection, get_col_names_from_table
from WeightFeatures import *

'''
Get a list of mutations (bed9 format) and:
    - filter mutations affecting motifs with TF expression=0 or TFbinding=0 and breaking score<0.3
        and groupby mutation to get a single line for each mutation [use: unify_muts(muts_file, filter_mut_motifs=True)]
    - for each mutation add a column to hold information about the motif with the max score [use get_max_motif_in_grouped_muts(annotated_mutations_grouped_file=grouped_muts_file)]
    - for each mutation in the observed set get a pvalue using scores of annotated mutations in the simulated set
    - for each mutation in the simulated set get a pvalue using scores of annotated mutations in the simulated set
    - filter insignifant mutations
    - merge nearby mutations and get a final score for each element.
    - compare the obtained score of each observed mutation to the distribution of scores of merged simulated elemens
    - correct the p-values
    
'''

def updateColNames(col_names):
    for i in range(0,len(col_names)):
        col_names[i] = '_'.join((col_names[i].replace('(','').replace(')','').replace('-','__')).split()).lower()
    return col_names

def updateColName(col_name):
    return '_'.join((col_name.replace('(','').replace(')','').replace('-','__')).split()).lower()
    
def get_cond_stmt(conds = []):
    if len(conds)>0:
        return ' where {} '.format(' and '.join(conds))
    else:
        return ''
    
def get_limit_smt(number_rows_select='all'):
    limit_number_rows_select_stmt = ""
    if number_rows_select!="all":
        if number_rows_select>0:
            limit_number_rows_select_stmt = ' limit {}'.format(str(number_rows_select))
    return limit_number_rows_select_stmt

def run_query(col_names_to_retrieve, cond_statement, limit_number_rows_select_stmt, cell_table, conn, n):
    curs = conn.cursor(name = "countcurs"+n, cursor_factory=DictCursor)
    #print "Query statement: ", 'select {} from {} {} {}'.format(col_names_to_retrieve, cell_table, cond_statement, limit_number_rows_select_stmt)
    #return pd.read_sql_query('select {} from {}{} {}'.format(col_names_to_retrieve, cell_table, cond_statement, limit_number_rows_select_stmt), conn)#, index_col=['chr', 'start', 'end', 'name'])
    curs.execute('select {} from {}{} {}'.format(col_names_to_retrieve, cell_table, cond_statement, limit_number_rows_select_stmt))#.fetchall()
    if curs is not None:
        return curs.fetchall()
        curs.close()
    else:
        curs.close()
        return []

def get_col_names_from_cells_assays(col_names, cells, assays, col_names_from_db):

    cells_to_report = []
    if len(cells)>0:
        if cells[0]=='all':
            for c in col_names_from_db:
                if '___' in c:
                    cells_to_report.append(c.split('___')[0])
        else:
            cells_to_report = cells
    for c in cells_to_report:
        for a in assays:
            if a == 'all':
                for ca in col_names_from_db:
                    if c in ca:#if the cell name was in the col_names
                        if ca not in col_names:
                            col_names.append(ca)
                break
            elif c+"___"+a in col_names_from_db:
                if c+"___"+a not in col_names:
                    col_names.append(c+"___"+a)
    return ','.join(col_names)
    
def get_cell_info_for_mutations(mutations_input_file, mutations_output_file, 
                                db_name, 
                                phenotype_collections, motif_info_cols,
                                cols_indices_to_report_from_file, cols_names_to_report_from_file,
                                mut_ref_index = 3, mut_alt_index = 4, phenotype_index_infile = 5, 
                                number_rows_select='all', sep='\t', region_name_index = 6,
                                TF_motif_weights_dict={}, annotation_weights={}):
    
    if os.path.exists(mutations_output_file):
        return mutations_output_file
    print "Annotating {} ...".format(mutations_input_file)
    conn = open_connection(db_name)
    #curs = conn.cursor()
    limit_stmt = get_limit_smt(number_rows_select=number_rows_select) 
    t =  time.time()
    number_lines_processed = 0
    with open(mutations_input_file, 'r') as mutations_infile, open(mutations_output_file, 'w') as outfile:
        line = mutations_infile.readline()
        while line:
            sline = line.strip().split(sep)
            if ((line.startswith('#') or line.startswith('//') or len(sline)<3) or 
                ( ((int(float(sline[2]))-int(float(sline[1]))) + 1 != len(sline[mut_ref_index]) and sline[mut_ref_index]!='-' and sline[mut_alt_index]!='-'))):#skip mis appropriate lines
                print 'Warning -- skipped line: ', line
                line = mutations_infile.readline()
                continue
            updated_chr = sline[0].replace('X', '23').replace('Y', '24').replace('MT','25').replace('M','25')
            chr_table = updated_chr+'motifs'
            if not updated_chr.startswith('chr'):
                chr_table = 'chr'+updated_chr+'motifs'
            try:
                phenotype_table = phenotype_collections[sline[phenotype_index_infile]][0]
                phenotype_cols = phenotype_collections[sline[phenotype_index_infile]][1:]
            except KeyError:
                print "WARNING: Skipped line, No collection is available for: ", sline[phenotype_index_infile]
                line = mutations_infile.readline()
                continue
            conds = []
            conds.append('(' + ("posrange && int4range({0},{1},'[]')".format(int(float(sline[1])), int(float(sline[2]))) + ')'))
            conds.append(' {}.mid={}.mid '.format(chr_table, phenotype_table))
            cond_stmt = get_cond_stmt(conds)
            
            col_names_to_retrieve = motif_info_cols[:]
            col_names_to_retrieve.extend(phenotype_cols)
            
            rows = run_query(','.join(col_names_to_retrieve), cond_stmt, limit_stmt, chr_table+','+phenotype_table, conn, str(number_lines_processed))
            cols_to_report_from_file = []
            for c in cols_indices_to_report_from_file:
                cols_to_report_from_file.append(sline[c])
            
            for row in rows:
                #compute a final score from the retrieved results
                lrow = list(row)
                motif_score = compute_motif_score(row, annotation_weights=annotation_weights)
                
                #get motif-breaking score for the each mutated_motif
                breaking_score, breaking_score_cumulative, mut_sig, motif_mut_pos = get_motif_breaking_score(TF_motif_weights_dict, row['name'], row['strand'], row['motifstart'], row['motifend'], 
                                                                  int(float(sline[1])), int(float(sline[2])), sline[mut_ref_index].upper(), sline[mut_alt_index].upper())
                #breaking_score, mut_sig, motif_mut_pos = get_motif_breaking_score(TF_motif_weights_dict, row[3], row[6], row[1], row[2], int(sline[1]), int(sline[2]), sline[mut_ref_index].upper(), sline[mut_alt_index].upper())
                
                outfile.write('\t'.join([str(x) for x in cols_to_report_from_file]) + '\t' + str(motif_score) + '\t' + 
                              str(breaking_score) + '\t'+ str(breaking_score_cumulative) + '\t' + mut_sig + '\t' + motif_mut_pos + '\t' +
                              '\t'.join([str(x) for x in lrow]) + '\n')
            line = mutations_infile.readline()
            
            number_lines_processed+=1
            if number_lines_processed % 100000 == 0:
                print '{} Lines are processed from {}'.format(number_lines_processed, mutations_input_file)
                print time.time()-t
                t = time.time()
                close_connection(conn)
                conn = open_connection(db_name)
            
    close_connection(conn)
    print "Finished: ", mutations_output_file
    return mutations_output_file
   
def get_collection_per_phenotype(db_name, phenotype_collection_matchings_input):
    """for each line in the phenotype_collection_matchings_input generate a list of cols
        phenotype_collection_matchings_input should either be a list or a system file (one element/line per phenotype).
    """
    
    phenotype_collection_matchings_input_lines = []
    phenotype_collection_cols = {}
    try:
        with open(phenotype_collection_matchings_input, 'r') as phenotype_collection_matchings_input_file:
            phenotype_collection_matchings_input_lines = phenotype_collection_matchings_input_file.readlines()
    except (IOError, TypeError):
        phenotype_collection_matchings_input_lines = phenotype_collection_matchings_input
        
    conn = open_connection(db_name)
    for l in phenotype_collection_matchings_input_lines:
        if l.startswith('//') or l.startswith('#') or l=="":
            continue
        phenotype_input = l.strip().split('=')[0]
        if '=' not in l:
            phenotype = updateColName(phenotype_input)
            available_cols = get_col_names_from_table(phenotype, conn)
            phenotype_collection_cols[phenotype_input] = [phenotype]#index 0 is for the table name and the remaining elements refer to the col names
            for c in available_cols:
                if phenotype+'.'+c+' as '+c not in phenotype_collection_cols[phenotype_input]:
                    phenotype_collection_cols[phenotype_input].append(phenotype+'.'+c+' as '+c)
        else:
            matchings = l.strip().split('=')[1]
            for matching in matchings.split(';'):
                matching_collection = updateColName(matching.split(':')[0])
                available_cols = get_col_names_from_table(matching_collection, conn)
                phenotype_collection_cols[phenotype_input] = [matching_collection]#index 0 is for the table name and the remaining elements refer to the col names
                if ':' in matching:
                    matching_collection_assays = updateColNames(matching.split(':')[1].split(','))
                    for matching_collection_assay in matching_collection_assays:
                        if matching_collection_assay in available_cols:
                            if matching_collection+'.'+matching_collection_assay+' as '+matching_collection_assay not in phenotype_collection_cols[phenotype_input]:
                                phenotype_collection_cols[phenotype_input].append(matching_collection+'.'+matching_collection_assay+' as '+matching_collection_assay)
                else:
                    for c in available_cols:
                        if matching_collection+'.'+c+' as '+c not in phenotype_collection_cols[phenotype_input]:
                            phenotype_collection_cols[phenotype_input].append(matching_collection+'.'+c+' as '+c)
        if len(phenotype_collection_cols[phenotype_input])<=1:
            print "Nothing to be done for: ", phenotype_input
            phenotype_collection_cols.pop(phenotype_input)
    return phenotype_collection_cols


def get_motif_breaking_score(TF_motif_weights_dict, motif_name, motif_strand, motif_start, motif_end, mut_start, mut_end, ref_allele, alt_allele):
    
    if motif_strand=='-':
        ref_allele = ref_allele.translate(string.maketrans('ACGT','TGCA'))
        alt_allele = alt_allele.translate(string.maketrans('ACGT','TGCA'))
    
    breaking_score = 0.0
    breaking_score_cumulative = 0.0
    mut_sig = ""
    motif_mut_pos_start = 0
    motif_mut_pos_end = 0
    motif_length = motif_end-motif_start
    if mut_start >= motif_start and mut_end <=motif_end:#motif contains the mutation
        if motif_strand=='+':
            motif_mut_pos_start = mut_start-motif_start
            motif_mut_pos_end = mut_end-motif_start
        else:
            motif_mut_pos_start = motif_end-mut_end
            motif_mut_pos_end = motif_end-mut_start
    elif mut_start < motif_start and (mut_end >=motif_start and mut_end <=motif_end):#mut stretches to the left of the motif
        bp_to_strip = motif_start-mut_start
        if motif_strand == '+':
            motif_mut_pos_start = 0
            motif_mut_pos_end = mut_end-motif_start
        else:
            motif_mut_pos_start = motif_end-mut_end
            motif_mut_pos_end = motif_length
        
        if not ref_allele == '-':#if it is not insertion
            ref_allele = ref_allele[bp_to_strip:]
        if not alt_allele == '-' and not ref_allele == '-':#if it is not deletion nor insertion (don't touch insertions)
            alt_allele = alt_allele[bp_to_strip:]
            
            
    elif (mut_start >= motif_start and mut_start <= motif_end) and mut_end >motif_end:#mut stretches to the right of the motif
        if not ref_allele == '-':
            bp_to_strip = len(ref_allele)-(mut_end-motif_end)
            ref_allele = ref_allele[:bp_to_strip]
        if not alt_allele == '-' and not ref_allele == '-':
            #bp_to_strip = len(ref_allele)-(mut_end-motif_end)
            alt_allele = alt_allele[:bp_to_strip]
        
        if motif_strand=='+':
            motif_mut_pos_start = mut_start-motif_start
            motif_mut_pos_end = motif_length
        else:
            motif_mut_pos_start = 0
            motif_mut_pos_end = motif_end-mut_start
        
    elif mut_start < motif_start and mut_end > motif_end:#motif contains the mutation
        motif_mut_pos_start = 0
        motif_mut_pos_end = motif_length
        bp_to_strip = motif_start-mut_start
        if not ref_allele == '-':
            ref_allele = ref_allele[bp_to_strip:bp_to_strip+motif_length+1]
        if not alt_allele == '-' and not ref_allele=='-':
            alt_allele = alt_allele[bp_to_strip:bp_to_strip+motif_length+1]
    
    '''print TF_motif_weights_dict[motif_name]
    print motif_name, len(TF_motif_weights_dict[motif_name]), motif_strand
    print motif_start, motif_end
    print mut_start, mut_end
    print motif_mut_pos_start, motif_mut_pos_end
    print breaking_score
    print ref_allele, '>', alt_allele
    print mut_sig
    '''
    if ref_allele == '-' or alt_allele == '-':
        breaking_score = 1.0
        breaking_score_cumulative = (motif_mut_pos_end-motif_mut_pos_start)+1#number of deleted or inserted bps
        if breaking_score_cumulative>motif_length+1:#in cases where an insertion streches more than the motif then just count the number of bases in the motif and ignore the rest
            breaking_score_cumulative = motif_length+1
    else:
        for i, mut_pos in enumerate(range(motif_mut_pos_start, motif_mut_pos_end+1)):
            try:
                breaking_score += abs(TF_motif_weights_dict[motif_name][mut_pos][ref_allele[i]] - TF_motif_weights_dict[motif_name][mut_pos][alt_allele[i]])
                breaking_score_cumulative+=breaking_score
            except KeyError:
                continue
            #keep the breaking_score at max 1
        if breaking_score>=1.0:
            breaking_score = 1.0
    mut_sig = ref_allele+">"+alt_allele        
    
    
    motif_mut_pos = str(motif_mut_pos_start+1) + '-' + str(motif_mut_pos_end+1)
    if motif_mut_pos_start==motif_mut_pos_end:
        motif_mut_pos = str(motif_mut_pos_start+1)
    
    return breaking_score, breaking_score_cumulative, mut_sig, motif_mut_pos

'''
desc: reports motif-breaking info for mutations in each motif sites. It check the difference between nucleotide frequency of the ref allele and the mutated-to allele of each mutation in the PWM of the anchor motif.
in: motifs_PFM matrix (get it from ENCODE, for instance), mutated motifs infput file (contains mutations info), output file name
out: all mutated motifs with adding breaking info to the mutations, all mutations at the motifs with adding breaking info  
out: only mutated motifs have motif-breaking mutation (difference between ref and mut allele > given_Threshold) and the fourth returned file is the list of motif breaking mutations (diff>Threshold)
'''
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
                            for i,nucl in enumerate(nucleotides):
                                nucl_weigts[nucl] = float(freq_per_allele[i])
                            PFM_motifs_dict[motif_name].append(nucl_weigts)#, C: float(split_line[1]), G: float(split_line[2]), T: float(split_line[3])})
    return PFM_motifs_dict


def compute_motif_score(row, annotation_weights={}):
    motif_score = 0.0
    for c in row.keys():
        if row[c]=='nan':
            continue
        wc = 0
        cv = 0
        try:
            if c=='dnase__seq' or c=='tfbinding' or c=='fantom':
                cv = make_binary(row[c])
                if annotation_weights[c]>0:
                    wc = cv*annotation_weights[c]
            elif c=='numothertfbinding':
                cv = row[c]
                if cv>3:    
                    cv = 3 #limit max number of other TFs to 10 to limit the score range (10 is already very high)
                if annotation_weights[c]>0:
                    wc = cv*annotation_weights[c]
            elif c=='replidomain':# or c=='cellname' or c=='name':
                cv = row[c]
                if annotation_weights[cv]>0:
                    wc = annotation_weights[cv]
            elif c == 'chromhmm':
                cv = combine_lables(row[c])
                if annotation_weights[cv]>0:
                    wc = annotation_weights[cv]
            elif c=='tfexpr': 
                cv = float(row[c])
                if annotation_weights[c]>0:
                    if cv>0:
                        cv = np.log2(float(row[c]))
                        wc = cv*annotation_weights[c]
            if wc > 0:
                motif_score+=wc
        except KeyError:
            continue
    return motif_score#np.exp(motif_score)


def get_annotation_weights(db_name, motif_info_col_names, col_names_to_weight_param, weights_params_dict_file):
    
    if os.path.exists(weights_params_dict_file):
        with open(weights_params_dict_file, 'r') as weights_params_dict_readfile:
            l = weights_params_dict_readfile.readline()
            logit_params = json.loads(l)
            return logit_params
        
    datafiles_motifs_dir = 'datafiles/Motifs/motifs_split_chr'
    
    datafiles_HepG2_geneexpr_dir = 'datafiles/GeneExp/ENCODEGeneExpr/HepG2/HepG2.bed' 
    datafiles_K562_geneexpr_dir = 'datafiles/GeneExp/ENCODEGeneExpr/K562/K562.bed'
    datafiles_GM12878_geneexpr_dir = 'datafiles/GeneExp/ENCODEGeneExpr/GM12878/GM12878.bed'
    datafiles_MCF7_geneexpr_dir = 'datafiles/GeneExp/ENCODEGeneExpr/MCF-7/MCF-7.bed'
    
    training_dir_results = 'datafiles/TrainingSets/Weight_features_analysis'
    training_dir_Ernst = 'datafiles/TrainingSets/Ernst_NatGen_2016'
    training_dir_Tewhey = 'datafiles/TrainingSets/Tewhey_Cell2016' 
    training_dir_Vockley = 'datafiles/TrainingSets/Vockley_Cell_2016'
    
    logit_params = get_param_weights(col_names_to_weight_param, db_name, motif_info_col_names, datafiles_motifs_dir, 
                      training_dir_results, training_dir_Ernst, training_dir_Tewhey, training_dir_Vockley,
                      datafiles_HepG2_geneexpr_dir, datafiles_K562_geneexpr_dir, datafiles_GM12878_geneexpr_dir, datafiles_MCF7_geneexpr_dir)
    
    with open(weights_params_dict_file, 'w') as weights_params_dict_outfile:
        json.dump(logit_params.params.to_dict(), weights_params_dict_outfile)
    #print logit_params.summary()
    return logit_params.params.to_dict()#log values    


if __name__ == '__main__':
    db_name = 'regmotifs'
    phenotype_collection_matchings_input = 'PhenotypeCollectionMatchings' #'Biliary-AdenoCA=Biliary-AdenoCA:replidomain,tfbinding;Bladder-TCC:dnase-seq;Bone-Epith#Liver-HCC#dfdfd'.split('#')
    weights_params_dict_file = "weightsdict_used.txt"
    #generate phenotype_collection_cols (contains the list of columns to be used for each named phenotype) from phenotype_collection_matchings and available_cols
    #available_cols = ['chromhmm', 'contactingdomain', 'dnase__seq', 'fantom', 'loopdomain', 'numothertfbinding', 'othertfbinding', 'replidomain', 'tfbinding', 'tfexpr']
    phenotype_collections = get_collection_per_phenotype(db_name, phenotype_collection_matchings_input)
    #expected output: phenotype_collection_cols = {'Ovary-AdenoCA': ['ovary__adenoca___chromhmm', 'ovary__adenoca___contactingdomain', 'ovary__adenoca___dnase__seq'], 'Liver-HCC':['liver__hcc___chromhmm', 'liver__hcc___contactingdomain', 'liver__hcc___dnase__seq', 'liver__hcc___fantom']}
    
    motif_info_cols = ['chr', 'motifstart', 'motifend', 'name', 'score', 'pval', 'strand']
    cols_indices_to_report_from_file=[0,1,2,3,4,5,6,7,8]
    cols_names_to_report_from_file=['Mut_chr', 'Mut_start', 'Mut_end', 'Mut_ref', 'Mut_alt','Cancer_type', 'Mut_type', 'Sample_ID', 'Donor_id'] 
    
    motif_PFM_input_file = "datafiles/Motifs/JASPAR_CORE_2016_vertebrates.meme"
    TF_motif_weights_dict = get_freq_per_motif(motif_PFM_input_file)
    
    col_names_to_weight_param = ['ChromHMM'.lower(), 'DNase__seq'.lower(), 'FANTOM'.lower(), 'NumOtherTFBinding'.lower(), 'RepliDomain'.lower(), 'TFBinding'.lower(), 'TFExpr'.lower(), 'score'.lower()]#sys.argv[1].split(',')#
    annotation_weights = get_annotation_weights(db_name, motif_info_cols, col_names_to_weight_param, weights_params_dict_file)
    print annotation_weights
    mutations_input_file = sys.argv[1]
    mutations_output_file = sys.argv[2]
    
    nprocess = 0
    num_cores = 0
    print sys.argv
    if len(sys.argv)>4:
        num_cores = int(sys.argv[3])
    if num_cores>1 and os.path.isdir(mutations_input_file) and os.path.isdir(mutations_output_file) and len(sys.argv)>=6:
        start_index= int(sys.argv[4])
        end_index = int(sys.argv[5])
        p = Pool(num_cores)
        
        '''other_mut_files = ['mutations_files/observed.bed9', 'mutations_files/simulation_broad.bed9', 'mutations_files/simulation_Sangerneutral.bed9', 'mutations_files/simulation_dkfz.bed9']
        other_mut_outfiles = ['mutations_files/observed_nocur_annotated.bed9', 'mutations_files/simulation_broad_nocur_annotated.bed9', 'mutations_files/simulation_Sangerneutral_nocur_annotated.bed9', 'mutations_files/simulation_dkfz_nocur_annotated.bed9']
        for i,f in enumerate(other_mut_files):
            p.apply_async(get_cell_info_for_mutations, args=(other_mut_files[i], 
                                                           other_mut_outfiles[i], 
                                                           db_name,
                                                           phenotype_collections, motif_info_cols,
                                                           cols_indices_to_report_from_file, cols_names_to_report_from_file,
                                                           3, 4, 5,
                                                           'all', '\t', 6,
                                                           TF_motif_weights_dict, annotation_weights))
            nprocess+=1
        '''
        for i in range(start_index, end_index+1):
            p.apply_async(
                          get_cell_info_for_mutations, args=(mutations_input_file+'/randomised_' + str(i) + '_f100k_w50kb_nonparametric.bed9', 
                                                               mutations_output_file+'/randomised_' + str(i) + '_f100k_w50kb_nonparametric_nocur_annotated.bed9', 
                                                               db_name,
                                                               phenotype_collections, motif_info_cols,
                                                               cols_indices_to_report_from_file, cols_names_to_report_from_file,
                                                               3, 4, 5,
                                                               'all', '\t', 6,
                                                               TF_motif_weights_dict, annotation_weights))
            nprocess+=1
            if nprocess%num_cores==0:
                p.close()
                p.join()
                print 'Number of sets finished{}'.format(nprocess)
                p = Pool(num_cores)
        p.close()
        p.join()
        '''split_mutations_input_file = []
        temp_split_dir = mutations_output_file+"_dir"
        if not os.path.exists(temp_split_dir):
            os.mkdir(temp_split_dir)
            cmd = """awk '{{print $0 >> "{}/"$1}}' {}""".format(temp_split_dir, mutations_input_file)
            os.system(cmd)
        split_mutations_input_file = [temp_split_dir + "/" + f for f in os.listdir(temp_split_dir) if '_annotated' not in f]
        split_mutations_output_file = [temp_split_dir + "/" + f + "_annotated" for f in os.listdir(temp_split_dir) if '_annotated' not in f]
        print split_mutations_input_file
        print split_mutations_output_file
        
        p = Pool(num_cores)
        for i in range(0, len(split_mutations_input_file)):
            p.apply_async(get_cell_info_for_mutations, args=(split_mutations_input_file[i], split_mutations_output_file[i], db_name,
                                phenotype_collections, motif_info_cols,
                                cols_indices_to_report_from_file, cols_names_to_report_from_file,
                                3, 4, 5, 
                                'all', '\t', 6,
                                TF_motif_weights_dict, annotation_weights))
        p.close()
        p.join()
        
        with open(mutations_output_file, 'w') as mutations_outfile: 
            for i in range(0, len(split_mutations_output_file)):
                with open(split_mutations_output_file[i], 'r') as split_mutations_outfile:
                    mutations_outfile.write(split_mutations_outfile.read())
        #os.removedirs(temp_split_dir)
    '''
    else:
        annotated_mutations_output_file = get_cell_info_for_mutations(mutations_input_file=mutations_input_file, mutations_output_file=mutations_output_file, 
                                    db_name=db_name,
                                    phenotype_collections=phenotype_collections, motif_info_cols=motif_info_cols,
                                    cols_indices_to_report_from_file=cols_indices_to_report_from_file, cols_names_to_report_from_file=cols_names_to_report_from_file,
                                    mut_ref_index = 3, mut_alt_index = 4, phenotype_index_infile = 5,
                                    TF_motif_weights_dict = TF_motif_weights_dict, annotation_weights=annotation_weights)
    
    
