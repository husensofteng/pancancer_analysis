library(ActiveDriverWGS)
library(tidyverse)
library(data.table)
#args
df_elem = args[1]
df_mut = args[2]
n_mut_ele = args[3]

#elements
data_my <- read.csv(df_elem, sep = '\t', header = TRUE, skip = 6)


#at least 2 mutations in an element
data_my_2mut <- data_my[which(data_my$X.Muts >= n_mut_ele),] 

#elements to test
data_elem <- data_my_2mut[,c(1:4)] 
colnames(data_elem) <- c('chr',
                         'start',
                         'end',
                         'id')
data_elem$chr <- gsub('23','X',data_elem$chr)
data_elem$chr <- gsub('24','Y',data_elem$chr)  
data_elem$start <- as.numeric(data_elem$start)
data_elem$end <- as.numeric(data_elem$end)
data_elem$id <- as.character(data_elem$id)


#mutations to test
data_mut <- fread(df_mut)
mut <- data_mut
colnames(mut) <- c('chr', 'pos1', 'pos2','ref', 'alt', 'cohort','var_type', 'sample_id','patient', 'pval','fsore',paste0('V', seq(from = 12, to = 32)))
mut_in <- mut[,c('chr', 'pos1', 'pos2','ref', 'alt', 'patient')] %>%  unique
mut_in$chr <- gsub('23','X',mut_in$chr)
mut_in$chr <- gsub('24','Y',mut_in$chr)
mut_in$pos1 <- as.numeric(mut_in$pos1)
mut_in$pos2 <- as.numeric(mut_in$pos2)
result = ActiveDriverWGS(mutations = mut_in,
                         elements = data_elem)

saveRDS(results, file= paste0('ActiveDriverWGS_results',n_mut_ele,'.RDS'))