# pancancer_analysis
Scripts used to analyse the PCAWG datasets

## 1. Motif detection and annotation
GenerateMotifsFIMO.py based on datafiles/Motifs/

## 2. Download datasets

- ENCODE_searchTerms_29Nov.txt: shows steps for generating metadata for downlaoding 

- ParseCellInfo.py: Downloads and processes datasets based on datafiles/CellInfoDict

- GetDataSetsAPI.py: dowanlods RNA-seq from ENCODE and reports expression level per gene per sample 


## 3. Generate database tables and retrieve data

- score_motifs_tissuepertable.py: extract motifs per TF and annotates them based on data obtained from (1) and (2). Reads parameters from score_motifs_param_regulatorymotifs.conf

- GetDBData_TrainingSets.py: get motifs from MPRAs and inactive genes

- WeightFeatures.py: Get weights per feature for motif scoring based on active and inactive motifs

- RetrieveDB.py: get motif info from the motif datbase

4. Mutation annotation


## 5. Element annotation
...

## 6. Downstream analysis
