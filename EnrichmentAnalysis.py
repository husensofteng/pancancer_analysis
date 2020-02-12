'''
Created on Aug 12, 2016

@author: Husen M. Umer
'''
import os,sys
from scipy.stats import binom, hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
#preprocess the input regions file 
#get closest genes according to the options (all genes/only certain biotypes), distance, max_number_genes. Return a list of gene names
def calculate_pval_for_genesets(geneset_enrichement_results_input_file, index_total_number_genes_per_set=2, index_number_enriched_genes=3, total_number_of_genes_in_the_universe=27000, total_number_of_genes_tried_in_the_search=3135, header_line = True, number_of_tried_gene_sets=24, keywords_to_filter_out_with=[]):#although not all are recognized in the pathways)
    
    infile = open(geneset_enrichement_results_input_file, 'r')
    calculated_p_value_out_file  = '.'.join(geneset_enrichement_results_input_file.split('.')[0:-1]) + "_calculated_pval." + geneset_enrichement_results_input_file.split('.')[-1]
    calculated_p_value_sig_out_file  = '.'.join(geneset_enrichement_results_input_file.split('.')[0:-1]) + "_calculated_pval_sig." + geneset_enrichement_results_input_file.split('.')[-1]
    calculated_p_value_sig_out_file_keywords  = '.'.join(geneset_enrichement_results_input_file.split('.')[0:-1]) + "_calculated_pval_sig_keywords." + geneset_enrichement_results_input_file.split('.')[-1]
    
    outfile = open(calculated_p_value_out_file, "w")
    outfile_sig = open(calculated_p_value_sig_out_file, "w")
    outfile_sig_keywords = open(calculated_p_value_sig_out_file_keywords, "w")
    
    M = total_number_of_genes_in_the_universe
    N = total_number_of_genes_tried_in_the_search
    sep = '\t'
    if header_line:
        header = infile.readline()
        outfile.write(sep.join(header.strip().split(sep)[0:-1]) + sep + "p-value" + sep + "q-value"  + sep + header.strip().split(sep)[-1] + "\n")
        outfile_sig.write(sep.join(header.strip().split(sep)[0:-1]) + sep + "p-value" + sep + "q-value"  + sep + header.strip().split(sep)[-1] + "\n")
        outfile_sig_keywords.write(sep.join(header.strip().split(sep)[0:-1]) + sep + "p-value" + sep + "q-value"  + sep + header.strip().split(sep)[-1] + "\n")
    inlines = infile.readlines()
    calculated_pvalues = []
    for line in inlines:
        split_line = line.split(sep)
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
        print "No genesets are reported, check the filters, params and gene names"
    number_of_sig_enriched_genesets = 0
    for l in range(0, len(inlines)):
        outfile.write(sep.join(inlines[l].strip().split(sep)[0:-1]) + sep + str(calculated_pvalues[l]) + sep + str(corrected_p_values_list[l]) + sep + inlines[l].strip().split(sep)[-1] +"\n")
        #write only the significant ones and satisfying the given condition
        if calculated_pvalues[l]<0.05:
            number_of_sig_enriched_genesets+=1 
            outfile_sig.write(sep.join(inlines[l].strip().split(sep)[0:-1]) + sep + str(calculated_pvalues[l]) + sep + str(corrected_p_values_list[l]) + sep + inlines[l].strip().split(sep)[-1] +"\n")
            if len(keywords_to_filter_out_with)>0:
                for keyword in keywords_to_filter_out_with:
                    if keyword.upper() in inlines[l].split(sep)[0].upper() or keyword.upper() in inlines[l].split(sep)[1].upper(): #if the given keyword(s) was found in the gene set name or discreption then report 
                        outfile_sig_keywords.write(sep.join(inlines[l].strip().split(sep)[0:-1]) + sep + str(calculated_pvalues[l]) + sep + str(corrected_p_values_list[l]) + sep + inlines[l].strip().split(sep)[-1] + "\n")
                        break
    print "Number of significantly enriched genesets: "  + str(number_of_sig_enriched_genesets)
    return calculated_p_value_out_file, calculated_p_value_sig_out_file, calculated_p_value_sig_out_file_keywords
    
def find_overlap_genesets_genelist(geneset_input_file, genelist_input_file, enriched_genes_output_file, total_number_of_genes_in_the_universe=27000, 
                                   min_number_of_genes_be_enriched_for_geneset_to_be_reported = 10, index_gene_name=0, index_gene_names_start=2, 
                                   keywords_to_filter_out_with=[], only_keep_the_sig_file = True, min_number_of_genes_in_geneset_to_consider_the_geneset = 10, header_line = False):
    genesets_infile = open(geneset_input_file, 'r')
    genelist_infile = open(genelist_input_file, 'r')
    enriched_genes_outfile = open(enriched_genes_output_file, 'w')
    genesets_lines = genesets_infile.readlines()
    genelist_lines = genelist_infile.readlines()
    
    #read the gene names from the gene input list
    genes_in_genelist = []
    for gene_line in genelist_lines:
        if gene_line.split("\t")[index_gene_name].strip() not in genes_in_genelist:
            genes_in_genelist.append(gene_line.split("\t")[index_gene_name].strip())
    print "Number of genes provided for search: " + str(len(genes_in_genelist))
    print "Total number of genesets tried: " + str(len(genesets_lines))
    
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
            for gene_line in genes_in_genelist:
                if gene_line.split('\t')[index_gene_name].strip() in split_geneset_info[index_gene_names_start::]:
                    enriched_genes.append(gene_line.strip().split('\t')[index_gene_name])
            if len(set(enriched_genes)) >= min_number_of_genes_be_enriched_for_geneset_to_be_reported:
                enriched_genes_outfile.write(split_geneset_info[0] + "\t" + split_geneset_info[1] + "\t"+ str(total_number_of_genes_in_this_geneset) + "\t" + str(len(set(enriched_genes))) + "\t"+ ','.join(set(enriched_genes)) +"\n")
            
    genesets_infile.close()
    genelist_infile.close()
    enriched_genes_outfile.close()
    #calculate p-values for each gene set/pathway
    calculated_p_value_out_file, calculated_p_value_sig_out_file, calculated_p_value_sig_out_file_keywords  = calculate_pval_for_genesets(enriched_genes_output_file, index_total_number_genes_per_set=2, index_number_enriched_genes=3, total_number_of_genes_in_the_universe=total_number_of_genes_in_the_universe, total_number_of_genes_tried_in_the_search=len(set(genes_in_genelist)), header_line = header_line, number_of_tried_gene_sets=number_of_tried_genesets, keywords_to_filter_out_with=keywords_to_filter_out_with)
    
    if len(keywords_to_filter_out_with)<1:
        os.remove(calculated_p_value_sig_out_file_keywords)
    del genesets_lines
    del genelist_lines
    
    if only_keep_the_sig_file:
        os.remove(calculated_p_value_out_file)
        os.remove(enriched_genes_output_file)    

    return calculated_p_value_sig_out_file

if __name__ == '__main__':
    
    if len(sys.argv)<3:
        print "Usage: python EnrichmentAnalysis.py genesets_file gene_list_file [index_gene_names_in_gene_list_file] [index_gene_names_in_genesets_file] [total_number_of_genes_in_the_genome_used]"
    
    genesets_file = sys.argv[1]#e.g kegg pathways (format: PathwayID, Pathway description, total_number_of_genes_in_the_pathways, gene_symbols(comma separated) 
    gene_list_input_files = []
    if os.path.isdir(sys.argv[2]):
        for x in os.listdir(sys.argv[2]):
            gene_list_input_files.append(sys.argv[2] + "/" + x)
    else:
        gene_list_input_files = [sys.argv[2]]
    
    index_gene_names_in_gene_list_file = 0
    index_gene_names_start_in_geneset = 3
    total_number_of_genes_in_the_universe = 27000
    if len(sys.argv)==4:
        index_gene_names_in_gene_list_file = int(sys.argv[4])
    if len(sys.argv)==5:
        index_gene_names_start_in_geneset = int(sys.argv[5])
    if len(sys.argv)==5:
        total_number_of_genes_in_the_universe = int(sys.argv[6])
    
    for gene_list_file in gene_list_input_files:
        calculated_p_value_sig_out_file = find_overlap_genesets_genelist(geneset_input_file=genesets_file, genelist_input_file=gene_list_file, enriched_genes_output_file=gene_list_file+"_enrichimentanalysis", total_number_of_genes_in_the_universe=total_number_of_genes_in_the_universe, min_number_of_genes_be_enriched_for_geneset_to_be_reported = 1, index_gene_name=index_gene_names_in_gene_list_file, index_gene_names_start=index_gene_names_start_in_geneset, keywords_to_filter_out_with=[], only_keep_the_sig_file = False, min_number_of_genes_in_geneset_to_consider_the_geneset=10, header_line = True) 
        