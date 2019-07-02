'''
Created on May 27, 2017

@author: husensofteng
'''
import glob, os
from pybedtools import BedTool

def extract_cohort_info(source_dir, target_dir, tfbs_file, 
                        cohorts_to_extract=[]):
    for cohort in cohorts_to_extract:
        for f in glob.glob(source_dir+'/'+cohort+'*.bed9'):
            print f.split('/')[-1]
            BedTool(f).intersect(tfbs_file, split=True, f = 1.0).saveas(target_dir+'/'+f.split('/')[-1])
    

if __name__ == '__main__':

    cohorts_to_extract = ['Carcinoma-tumors_',
                          'All-tumors-without-Lymphatic-Melanoma-Liver-Esophagus-tumors_',
                          'All-tumors-without-Lymphatic-system-Skin-Melanoma_']
    extract_cohort_info(source_dir='mutations_cohorts_output', target_dir='mutations_cohorts_output_tfbs', 
                        tfbs_file = 'analysis/encode_merge_cdsAndSpliceSitesSubtracted.bed', 
                        cohorts_to_extract=cohorts_to_extract)