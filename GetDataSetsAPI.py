'''
Created on Dec 15, 2016

@author: husensofteng
'''
import sys, os
import requests 
import urllib
import numpy as np
import json

def get_cellnames_from_cellinfodict(cellinfodict_inputfile, cell_names_start_with="#"):
    cellinfo_lines = []
    with open(cellinfodict_inputfile, 'r') as cellinfodict_infile:
        cellinfo_lines = cellinfodict_infile.readlines()
    
    cell_lines = []
    cell_name = ""
    for line in cellinfo_lines:
        if line.startswith('//'):#skip lines start with // researved for comments in the config file
            continue
        elif line.startswith('***'):
            break
        elif line.startswith(cell_names_start_with):
            cell_name=line.strip().strip(cell_names_start_with)
            if cell_name not in cell_lines:
                cell_lines.append(cell_name)
    return cell_lines


def get_data_API(biosamples_out_dir="./", biosample_term_names_to_get=[]):
    HEADERS = {'accept': 'application/json'}
    biosample_term_names_to_get_str = '&biosample_term_name='+ '&biosample_term_name='.join(biosample_term_names_to_get)
    #Changed RNA-seq to total+RNA-seq in the URL
    URL = "https://www.encodeproject.org/search/?type=Experiment&assay_slims=Transcription&assay_title=total+RNA-seq&status=released&assembly=hg19&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=tsv&files.analysis_step_version.analysis_step.pipelines.title=RNA-seq+of+long+RNAs+%28paired-end%2C+stranded%29&frame=object{}".format(biosample_term_names_to_get_str)
    response_json_dict = requests.get(URL, headers=HEADERS).json()
    #print json.dumps(response_json_dict, indent=4, separators=(',', ': '))
    bio_sample_files = {}
    '''
    iterate through each file in @graph list, for each file 'files' key send a request to get the file's info in json
    check file_type, assembly,... and the get link to the file in 'href'
    '''
    for cell in range(0, len(response_json_dict['@graph'])):#it contains a list, one element for each experiment
        bio_sample = response_json_dict['@graph'][cell]['biosample_term_name']
        if 'hg19' not in response_json_dict['@graph'][cell]['assembly'] or response_json_dict['@graph'][cell]['status']!='released':
            continue
        if bio_sample not in bio_sample_files.keys():
            bio_sample_files[bio_sample] = []
        for f in response_json_dict['@graph'][cell]['files']:
            f_info = requests.get("https://www.encodeproject.org/{}/".format(f), headers=HEADERS).json()
            if 'assembly' in f_info.keys() and f_info['file_type']=='tsv' and f_info['output_type']=='gene quantifications':
                if f_info['assembly'] == 'hg19':
                    if not os.path.exists(biosamples_out_dir+bio_sample):
                        os.mkdir(biosamples_out_dir+bio_sample)
                    print biosamples_out_dir+bio_sample + ':' + f_info['href']
                    bio_sample_files[bio_sample].append(f_info['href'])
                    if not os.path.exists(biosamples_out_dir+bio_sample+'/'+f_info['href'].split('/')[-1]):
                        urllib.urlretrieve('https://www.encodeproject.org/{}'.format(f_info['href']), biosamples_out_dir+bio_sample+'/'+f_info['href'].split('/')[-1])
    
    print "Finished downloading files for: ", bio_sample_files.keys()
    #now the content of bio_sample_files can be further processed to combine files of each bio_sample
    
def process_RNA_seq_datafolder(input_dir_path, num_header_lines=1, gene_id_index=0, gene_value_index=5,sep='\t'):
    bio_samples = {}
    bio_samples_output_dict = "bio_samples_dict.json"
    if os.path.exists(bio_samples_output_dict):
        with open(bio_samples_output_dict, 'r') as bio_samples_infile_dict:
            bio_samples = json.load(bio_samples_infile_dict)
        return bio_samples
    for bio_sample_input_dir in os.listdir(input_dir_path):
        if os.path.exists(input_dir_path+'/'+bio_sample_input_dir+'/'+bio_sample_input_dir+".bed"):
            continue
        bio_samples[bio_sample_input_dir] = {}
        for bio_sample_input_file in os.listdir(input_dir_path+'/'+bio_sample_input_dir):
            if bio_sample_input_file.split('.')[-1] != "tsv":
                continue
            with open(input_dir_path+'/'+bio_sample_input_dir+'/'+bio_sample_input_file, 'r') as bio_sample_infile:
                for i in range(0, num_header_lines):
                    bio_sample_infile.readline()#skip the header lines
                genes = bio_sample_infile.readlines()
                for g in genes:
                    if g=="" or g.startswith('//'):
                        continue
                    gene_id = g.strip().split(sep)[gene_id_index]
                    gene_value = float(g.strip().split(sep)[gene_value_index])
                    if gene_id not in bio_samples[bio_sample_input_dir].keys():
                        bio_samples[bio_sample_input_dir][gene_id] = []
                    bio_samples[bio_sample_input_dir][gene_id].append(gene_value)
    with open(bio_samples_output_dict, 'w') as bio_samples_output_dict_outfile:
        json.dump(bio_samples, bio_samples_output_dict_outfile)
    return bio_samples

def get_expr_per_bio_sample(bio_samples, gencode_id_info_dict, output_dir_path):
    
    genes_with_gencode_info_dict = {}
    for bio_sample in bio_samples:
        genes_with_gencode_info_dict[bio_sample] = {}

        if os.path.exists(output_dir_path+'/'+bio_sample + '/' + bio_sample+".bed"):
            return genes_with_gencode_info_dict
        
        with open(output_dir_path+'/'+bio_sample + '/' + bio_sample+".bed", 'w') as bio_sample_outfile:
            for gene_id in bio_samples[bio_sample]:
                if gene_id in gencode_id_info_dict.keys():
                    genes_with_gencode_info_dict[bio_sample][gene_id] = np.mean(bio_samples[bio_sample][gene_id])
                    bio_sample_outfile.write('\t'.join(gencode_id_info_dict[gene_id]) + '\t' + str(np.mean(bio_samples[bio_sample][gene_id])) + '\n')
    return genes_with_gencode_info_dict
    #combine results of all bio samples to one file (one col per bio sample) --- extend the the GTEx file
    
def get_gene_names_and_ids_from_genecode(genecode_genes_input_file, genecode_genes_only_genes_bed_output_file, num_header_lines=0, sep='\t'):
    gencode_id_info_dict = {}
    gene_id_name_dict = {}
    gencode_id_info_dict_savefile = genecode_genes_only_genes_bed_output_file+"saveddict"
    if os.path.exists(gencode_id_info_dict_savefile):
        with open(gencode_id_info_dict_savefile, 'r') as gencode_id_info_dict_file:
            return json.load(gencode_id_info_dict_file)
    with open(genecode_genes_input_file, 'r') as genecode_genes_infile, open(genecode_genes_only_genes_bed_output_file, 'w') as genecode_genes_only_genes_bed_outfile:
        for i in range(0, num_header_lines):
            genecode_genes_infile.readline()
            
        l = genecode_genes_infile.readline()
        while l:
            if l.startswith('#') or l=="":
                l = genecode_genes_infile.readline()
                continue
            sl = l.strip().split(sep)
            if sl[2]=="gene":
                sl_info_dict = dict([x.split('=') for x in sl[8].split(';')])
                gencode_id_info_dict[sl_info_dict['gene_id']] = [sl[0], sl[3], sl[4], sl[5], sl[6], sl_info_dict['gene_id'], 
                                                              sl_info_dict['gene_name'], sl_info_dict['gene_status'], sl_info_dict['gene_type']]
                gene_id_name_dict[sl_info_dict['gene_id']] = sl_info_dict['gene_name']
                genecode_genes_only_genes_bed_outfile.write(sep.join(gencode_id_info_dict[sl_info_dict['gene_id']]) + '\n')
            l = genecode_genes_infile.readline()
    
    with open(gencode_id_info_dict_savefile, 'w') as gencode_id_info_dict_savefile_outfile:
        json.dump(gencode_id_info_dict, gencode_id_info_dict_savefile_outfile)            
    return gencode_id_info_dict#, gene_id_name_dict

if __name__ == '__main__':
    
    biosamples_dir_path=sys.argv[1]
    if not os.path.exists(biosamples_dir_path):
        os.mkdir(biosamples_dir_path)
    cellinfodict_inputfile = sys.argv[2] 
    biosample_term_names_to_get = get_cellnames_from_cellinfodict(cellinfodict_inputfile)
    #biosample_term_names_to_get = ['A549','GM12878','Gastric','HCT116','HEK293','HeLa-S3','HepG2','IMR-90','Ishikawa','K562','Kidney_curated','LNCaP clone FGC','MCF 10A','MCF-7','Ovary','Panc1','Pancreas_curated','Prostate_curated','SK-N-SH','T47D','U2OS','astrocyte','epithelial cell of esophagus','keratinocyte','osteoblast','urothelium cell line']
    genecode_genes_input_file=sys.argv[3]
    genecode_genes_only_genes_bed_output_file=genecode_genes_input_file+"_onlygenes.bed"
    
    get_data_API(biosamples_dir_path+'/', biosample_term_names_to_get)
    bio_samples = process_RNA_seq_datafolder(biosamples_dir_path)#give a dir that contains all bio_samples that have previously been downloaded
    gencode_id_info_dict = get_gene_names_and_ids_from_genecode(genecode_genes_input_file=genecode_genes_input_file, genecode_genes_only_genes_bed_output_file=genecode_genes_only_genes_bed_output_file)
    genes_with_gencode_info_dict = get_expr_per_bio_sample(bio_samples, gencode_id_info_dict, output_dir_path=biosamples_dir_path)
    
    