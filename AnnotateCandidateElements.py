import os,sys
from collections import Counter
from AnnotateMutations import get_annotated_muts
from pybedtools import BedTool

def load_element_coordinates(combined_bed12_elements_input_files):
    
    element_id_index_bed12 = 3
    elements_bed12_dict = {}
    for f in combined_bed12_elements_input_files:
        with open(f, 'r') as ifile:
            l = ifile.readline().strip().split('\t')
            while len(l)> element_id_index_bed12:
                elements_bed12_dict[l[element_id_index_bed12]] = l
                l = ifile.readline().strip().split('\t')
    
    return elements_bed12_dict

def get_element_coordinates_bed12(candidate_elements_input_file, elements_bed12_dict):
    
    index_element_id = 1
    with open(candidate_elements_input_file, 'r') as candidate_elements_ifile, open(candidate_elements_input_file+'.bed12', 'w') as candidate_elements_ofile:
        l = candidate_elements_ifile.readline().strip().split('\t')
        while len(l)>index_element_id:
            new_element_id = ';'.join(l)
            element_bed12_info = elements_bed12_dict[l[index_element_id]]
            element_bed12_info[3] = new_element_id
            candidate_elements_ofile.write('\t'.join(element_bed12_info) + '\n')
            l = candidate_elements_ifile.readline().strip().split('\t')
            
    return candidate_elements_input_file+'.bed12'
    
def get_cohort_cancer_type(cohorts_file, cohort_files):
    
    cohort_cancer_type_dict = {}
    with open(cohorts_file, 'r') as cohorts_ifile:
        lines = cohorts_ifile.readlines()
        for l in lines:
            if l!="" and not l.startswith('//'):
                if '=' not in l and '::' not in l and "*" not in l:
                    cohort_cancer_type_dict[l.strip().upper()] = l.strip()
                elif '=' in l:
                    cohort_name = l.split('=')[1].split('_')[0] + '_' + l.split('=')[1].split('/')[1].split('.')[0].replace('_tumors', '')
                    if cohort_name == 'meta_Female_reproductive_system':
                        cohort_name = 'meta_Female_reproductive_tract'
                    elif cohort_name == 'meta_Hematopoietic':
                        cohort_name = 'meta_Hematopoietic_system'
                    elif cohort_name == 'meta_Lymph':
                        cohort_name = 'meta_Lymphatic_system'
                    
                    cohort_cancer_type_dict[cohort_name.upper()] = []
                    file_to_read_from = l.strip().split('=')[1]
                    with open(file_to_read_from, 'r') as ifile:
                        l = ifile.readline()
                        while len(l)>1:
                            cancer_type = l.strip().split('\t')[1]
                            if cancer_type not in cohort_cancer_type_dict[cohort_name.upper()]:
                                cohort_cancer_type_dict[cohort_name.upper()].append(cancer_type)
                            l = ifile.readline()
    
    #add all cancer types except lymph and skin to Pancan-no-skin-melanoma-lymph cohort 
    cohort_cancer_type_dict['Pancan-no-skin-melanoma-lymph'.upper()] = ['Biliary-AdenoCA' ,'Bladder-TCC' ,'Bone-Cart' ,'Bone-Epith' ,'Bone-Leiomyo' ,'Bone-Osteosarc' ,'Breast-AdenoCa' ,'Breast-DCIS' ,'Breast-LobularCa' ,'Cervix-AdenoCA' ,'Cervix-SCC' ,'CNS-GBM' ,'CNS-Medullo' ,'CNS-Oligo' ,'CNS-PiloAstro' ,'ColoRect-AdenoCA' ,'Eso-AdenoCa' ,'Head-SCC' ,'Kidney-ChRCC' ,'Kidney-RCC' ,'Liver-HCC' ,'Lung-AdenoCA' ,'Lung-SCC' ,'Myeloid-AML' ,'Myeloid-MDS' ,'Myeloid-MPN' ,'Ovary-AdenoCA' ,'Panc-AdenoCA' ,'Panc-Endocrine' ,'Prost-AdenoCA' ,'Stomach-AdenoCA' ,'Thy-AdenoCA' ,'Uterus-AdenoCA']
    
    return cohort_cancer_type_dict

def get_muts_from_corresponding_cancertypes(muts_elements_input_file, cohort_cancer_type_dict):
    
    index_cancer_type = 5
    index_cohort = 12
    index_cohort_name_in_element_id = 3
    muts_elements_output_file = muts_elements_input_file+'onlycorrespondig_cancertypeMuts.bed9' 
    
    with open(muts_elements_input_file, 'r') as ifile, open(muts_elements_output_file, 'w') as ofile:
        l = ifile.readline().strip().split('\t')
        while len(l)>index_cohort:
            cohort_name = l[index_cohort].split(';')[index_cohort_name_in_element_id].upper()
            if l[index_cancer_type] in cohort_cancer_type_dict[cohort_name]:
                l[8] = l[index_cohort]
                ofile.write('\t'.join(l[0:9]) + '\n')
                
            l = ifile.readline().strip().split('\t')
    
    return muts_elements_output_file

def get_unique_annotations_per_element(muts_elements_output_file_annotated, element_id_index=8, annotations_index=9):
    
    muts_elements_output_file_annotated_grouped = muts_elements_output_file_annotated+"grouped.tsv"
    elements_inf_dict = {}
    elements_inf_dict_cells = {}
    with open(muts_elements_output_file_annotated, 'r') as muts_elements_output_file_annotated_ifile:
        l = muts_elements_output_file_annotated_ifile.readline().strip().split('\t')
        while len(l)>annotations_index:
            element_id = l[element_id_index]
            if element_id not in elements_inf_dict.keys():
                elements_inf_dict[element_id] = {'ChromHMM':[], 'DNase-seq': [], 'TFBinding':[], 'RepliDomain':[], 'FANTOM':[], 'ContactingDomain':[], 'LoopDomain':[]}
                elements_inf_dict_cells[element_id] = []
            
            if l[annotations_index]=='NaN':
                l = muts_elements_output_file_annotated_ifile.readline().strip().split('\t')
                continue
            elements_inf_dict_cells[element_id].extend(l[annotations_index+1].split('|'))
            
            for annotation in l[annotations_index].split('|'):
                if annotation.split(':')[0] in elements_inf_dict[element_id].keys():
                    if annotation.split(':')[0]=='TFBinding':
                        elements_inf_dict[element_id][annotation.split(':')[0]].extend(annotation.split(':')[1].split(';'))
                    else:
                        try:
                            elements_inf_dict[element_id][annotation.split(':')[0]].append(float(annotation.split(':')[1]))
                        except (ValueError, TypeError):
                            elements_inf_dict[element_id][annotation.split(':')[0]].append(annotation.split(':')[1])
            l = muts_elements_output_file_annotated_ifile.readline().strip().split('\t')
    print len(elements_inf_dict)
    
    with open(muts_elements_output_file_annotated_grouped, 'w') as muts_elements_output_file_annotated_grouped_ofile:
        element = 'ID;Name;target;tissue;Gene'
        element_annos = sorted(['ChromHMM', 'ContactingDomain', 'DNase-seq', 'FANTOM', 'RepliDomain', 'TFBinding', 'LoopDomain'])#sorted(['ChromHMM','DNase1', 'FANTOM-Peak', 'HIC-Contacting-Domain', 'Replication-Domain', 'TF-Peaks'])
        cell_names = 'Cells/Tissues of Evidence'
        muts_elements_output_file_annotated_grouped_ofile.write('\t'.join(element.split(';')) + '\t' + '\t'.join(element_annos) + '\t' + cell_names + '\n')
        for element in elements_inf_dict.keys():
            element_annos = []
            for anno in sorted(elements_inf_dict[element].keys()):
                annos = elements_inf_dict[element][anno]
                if len(annos)==0:
                    element_annos.append('NA')
                else:
                    try:
                        if sum(annos)/(len(annos)*1.0)>0.0:
                            element_annos.append('1')
                        else:
                            element_annos.append('NA')
                    except TypeError:
                        element_annos.append(','.join(k+":"+str(v) for k,v in dict.items(Counter(annos))))
                        #element_annos.append(Counter(annos))
            cell_names = ','.join(set(elements_inf_dict_cells[element]))
            if len(elements_inf_dict_cells[element])==0:
                cell_names = 'NA'
            if len(element.split(';'))>5:
                muts_elements_output_file_annotated_grouped_ofile.write('\t'.join(element.split(';')[0:4]) + '\t' + ';'.join(element.split(';')[4::]) + '\t' + '\t'.join(element_annos) + '\t' + cell_names + '\n')
            elif len(element.split(';'))>4:
                muts_elements_output_file_annotated_grouped_ofile.write('\t'.join(element.split(';')) + '\t' + '\t'.join(element_annos) + '\t' + cell_names + '\n')
            else:
                muts_elements_output_file_annotated_grouped_ofile.write('\t'.join(element.split(';')) + '\t\t' + '\t'.join(element_annos) + '\t' + cell_names + '\n')
            
    return muts_elements_output_file_annotated_grouped

if __name__ == '__main__':
    
    dir_in = 'analysis/candidate_drivers_pcawg_10Apr2017/'
    candidate_elements_input_file = dir_in+"driver_candidates_10Apr2917.tsv"
    elements_dir = dir_in+'ElementSetsAll/'
    mutations_input_file = 'mutations_files/observedunique.bed9'
    combined_bed12_elements_input_files = [elements_dir + x for x in os.listdir(elements_dir) if x.endswith('.bed')]
    
    elements_bed12_dict = load_element_coordinates(combined_bed12_elements_input_files=combined_bed12_elements_input_files)
    
    candidate_elements_bed12 = get_element_coordinates_bed12(candidate_elements_input_file, elements_bed12_dict)
    #intersect observed.bed9 with the candidate_elements_bed12 file
    muts_candidate_elements_bed12 = candidate_elements_bed12+'.muts'
    if not os.path.exists(muts_candidate_elements_bed12):
        BedTool(mutations_input_file).intersect(BedTool(candidate_elements_bed12), split = True, wo=True).saveas(muts_candidate_elements_bed12)
    #muts_elements_input_file = '/Users/husensofteng/Documents/workspace/ActiveMotifs/analysis/driver_candidates_17Mar2917.tsv.bed12.muts'
    
    cohorts_file = 'meta_tumor_cohorts_Oct2016/cohorts_to_run'
    cohort_files='meta_tumor_cohorts_Oct2016'
    cohort_cancer_type_dict = get_cohort_cancer_type(cohorts_file=cohorts_file, cohort_files=cohort_files)
    
    muts_elements_corresponding_cancertypes = get_muts_from_corresponding_cancertypes(muts_candidate_elements_bed12, cohort_cancer_type_dict)

    #annotate the muts_elements_output_file in AnnotateMutations module
    tracks_dir = 'datafiles/chromatin_marks_all_cells_onlynarrowpeaks'
    muts_elements_corresponding_cancertypes_annotated = muts_elements_corresponding_cancertypes+'_annotated'
    get_annotated_muts(muts_input_file=muts_elements_corresponding_cancertypes, tracks_dir=tracks_dir, muts_out=muts_elements_corresponding_cancertypes_annotated)
    #muts_elements_output_file_annotated  = '/Users/husensofteng/Documents/workspace/ActiveMotifs/analysis/candidate_drivers_annotated.out'
    get_unique_annotations_per_element(muts_elements_corresponding_cancertypes_annotated, element_id_index=8, annotations_index=9)
    
    