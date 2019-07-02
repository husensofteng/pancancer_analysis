#!/bin/bash
f=merged200bp_extended200bp_nofullexon_ATELM/combined_rand103setsTFsigQval0.05_meanTFExprMotifBreaking03Filters_mergedmuts200bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_aggregated0UpDwmaxdist2kb_within500kb.tsv
echo "Number of elements:"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | wc -l

echo "Number of RegMuts:"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f19 | awk '{gsub("\t","\n"); print}' | awk '{gsub(",","\n"); print}' | awk -F# '{if($0~":"){split($1,s,":"); split(s[2],ss,"-"); gsub("chrX","23",s[1]); gsub("chrY","24",s[1]); gsub("chrM","25",s[1]); gsub("chr","",s[1]); print s[1]"\t"ss[1]"\t"ss[2]"\t"$3"\t"$6} else print $1"\t"$2"\t"$3"\t"$6"\t"$9}' | wc -l

echo "number of samples that had 6 or fewer regmuts"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f19 | awk '{gsub(",","\n"); print}' | cut -f9 -d '#' | sort | uniq -c | awk '$1<=6' | wc -l

echo "number of samples that had 20 or more regmuts"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f19 | awk '{gsub(",","\n"); print}' | cut -f9 -d '#' | sort | uniq -c | awk '$1>=20' | wc -l

echo "Number of Samples with RegMuts"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f19 | awk '{gsub("\t","\n"); print}' | awk '{gsub(",","\n"); print}' | awk -F# '{if($0~":"){split($1,s,":"); split(s[2],ss,"-"); gsub("chrX","23",s[1]); gsub("chrY","24",s[1]); gsub("chrM","25",s[1]); gsub("chr","",s[1]); print s[1]"\t"ss[1]"\t"ss[2]"\t"$3"\t"$6} else print $1"\t"$2"\t"$3"\t"$6"\t"$9}' | sort -k1n -k2n -k3n -k4 -k5 | uniq | cut -f5 | sort | uniq | wc -l

echo "Number of Muts"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f20 | awk '{gsub("\t","\n"); print}' | awk '{gsub(",","\n"); print}' | awk -F# '{if($0~":"){split($1,s,":"); split(s[2],ss,"-"); gsub("chrX","23",s[1]); gsub("chrY","24",s[1]); gsub("chrM","25",s[1]); gsub("chr","",s[1]); print s[1]"\t"ss[1]"\t"ss[2]"\t"$3"\t"$6} else print $1"\t"$2"\t"$3"\t"$6"\t"$9}'| wc -l

echo "Number of samples with Muts"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f20 | awk '{gsub("\t","\n"); print}' | awk '{gsub(",","\n"); print}' | awk -F# '{if($0~":"){split($1,s,":"); split(s[2],ss,"-"); gsub("chrX","23",s[1]); gsub("chrY","24",s[1]); gsub("chrM","25",s[1]); gsub("chr","",s[1]); print s[1]"\t"ss[1]"\t"ss[2]"\t"$3"\t"$6} else print $1"\t"$2"\t"$3"\t"$6"\t"$9}' | sort -k1n -k2n -k3n -k4 -k5 | uniq | cut -f5 | sort | uniq | wc -l

echo "Number of RegMuts per ChromHMM state:"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f17 | awk '{gsub("],","\n"); print}' | awk '{split($0,s,"["); if(s[1]=="RegMuts-ChromHMM") print s[2]}' | awk '{gsub(",","\n"); print}' | awk '{gsub(":","\t"); print}' | sort -k1 | groupBy -g 1 -c 2 -o sum | sort -k2nr

echo "Number of muts per motif in TSS-like states"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f21 | awk '{gsub(",","\n"); print}' | awk '{gsub(",","\n"); print}' | awk '$0~"Tss" || $0~"BivFlnk"' | cut -f11 -d '#' | sort | uniq -c | sort -k1nr | head

echo "Number of muts per motif in Quies-like states"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f21 | awk '{gsub(",","\n"); print}' | awk '{gsub(",","\n"); print}' | awk '$0~"Quies" || $0~"Rpts" || $0~"Het"' | cut -f11 -d '#' | sort | uniq -c | sort -k1nr | head

echo "number of samples with regmuts in in Tss-like states"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f21 | awk '{gsub(",","\n"); print}' | awk '{gsub(",","\n"); print}' | awk '($0~"Tss" || $0~"BivFlnk")' | cut -f5 -d '#' | sort | uniq | wc -l

echo "number of samples with regmuts in SP1 or SP2 in Tss-like states"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f21 | awk '{gsub(",","\n"); print}' | awk '{gsub(",","\n"); print}' | awk '($0~"Tss" || $0~"BivFlnk") && ($0~"SP2_" || $0~"SP1_")' | cut -f5 -d '#' | sort | uniq | wc -l

echo "number of samples with regmuts in EGR1 in Tss-like states"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f21 | awk '{gsub(",","\n"); print}' | awk '{gsub(",","\n"); print}' | awk '($0~"Tss" || $0~"BivFlnk") && ($0~"EGR1_")' | cut -f5 -d '#' | sort | uniq | wc -l

echo "number of samples with regmuts in CTCF in Quies-like states"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f21 | awk '{gsub(",","\n"); print}' | awk '{gsub(",","\n"); print}' | awk '($0~"Quies" || $0~"Rpts" || $0~"Het") && ($0~"CTCF_")' | cut -f5 -d '#' | sort | uniq | wc -l

echo "number of samples with regmuts in CEBPB in Quies-like states"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f21 | awk '{gsub(",","\n"); print}' | awk '{gsub(",","\n"); print}' | awk '($0~"Quies" || $0~"Rpts" || $0~"Het") && ($0~"CEBPB_")' | cut -f5 -d '#' | sort | uniq | wc -l

echo "Number of intergenic,intronic,... regions in elements"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f27 | sort | uniq -c | sort -k1nr

echo "Number of Muts per ChromHMM state:"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f17 | awk '{gsub("],","\n"); print}' | awk '{split($0,s,"["); if(s[1]=="Muts-ChromHMM") print s[2]}' | awk '{gsub(",","\n"); print}' | awk '{gsub(":","\t"); print}' | sort -k1 | groupBy -g 1 -c 2 -o sum | sort -k2nr

echo "Number of mutations per TF peak:"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f20 | awk '{gsub("\t","\n"); print}' | awk '{gsub(",","\n"); print}' | awk -F# '{if($0~":"){if($7~"TFBinding" && $7~"ChromHMM:") print $7}}' | awk '{gsub("\\|","\n"); print}' | awk '$1~"TFBin"' | cut -f2 -d ':' | awk '{gsub(";","\n"); print}' | sort | uniq -c | sort -k1nr | head


echo "Number of CTCF-peak mutations that also have RAD21 peaks:"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f20 | awk '{gsub("\t","\n"); print}' | awk '{gsub(",","\n"); print}' | awk -F# '{if($0~":"){if($7~"TFBinding" && $7~"ChromHMM:") print $7}}' | awk '$1~"CTCF" && ($1~"RAD21")' | wc -l

echo "Number of CTCF-peak mutations that also have RAD21 peaks and SMC3 peaks"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f20 | awk '{gsub("\t","\n"); print}' | awk '{gsub(",","\n"); print}' | awk -F# '{if($0~":"){if($7~"TFBinding" && $7~"ChromHMM:") print $7}}' | awk '$1~"CTCF" && ($1~"RAD21" && $1~"SMC3")' | wc -l

echo "Number of CTCF-peak mutations per Chromatin state:"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | cut -f20 | awk '{gsub("\t","\n"); print}' | awk '{gsub(",","\n"); print}' | awk -F# '{if($0~":"){if($7~"TFBinding" && $7~"ChromHMM:") print $7}}' | awk '$1~"CTCF"' | cut -f1 -d '|' | cut -f2 -d ':' | sort | uniq -c | sort -k1nr 

echo "Total number of mutations in all elements that are associated to each gene:"
awk 'NR>7 && $10>1 && $8<0.05 && $24<0.05' $f | sort -k14nr | awk 'BEGIN{FS=OFS="\t"}{split($26,s, ","); for(i=1;i<=length(s); i++){split(s[i],ss,":"); print ss[1],$10,$14,$27,$8}}' | sort -k1 | groupBy -g 1 -c 2,3,3,4,5 -o sum,sum,collapse,collapse,collapse | sort -k3nr | awk 'NR>1' | awk '$2>=10' | head

echo "Single elements that are highly recurrent (#muts>=10):"
awk 'NR>7 && $10>1 && $14>=10' $f | sort -k14nr | cut -f10,8,9,13,14,26-  | head

echo "Single elemnets that have 10 or more CFRMs and are highly recurrent (#Muts>=10):"
awk 'NR>7 && $10>=10 && $14>=10' $f | sort -k14nr | cut -f10,8,9,13,14,26- 

echo "#Genes harbor 10 or more CFRMs:"
awk '$4>=10 && $6>=10' $f""_Genes.tsv | wc -l

echo "COSMIC or KCP Genes harbor 10 or more CFRMs:"
awk '$4>=10 && $6>=10' $f""_Genes.tsv | awk '$0~"COSMIC" || $0~"KCP"' | sort -k6nr | cut -f1,4,6
