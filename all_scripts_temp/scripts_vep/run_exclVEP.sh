#!/bin/sh
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 99:00:00
date
#cat /home/huum/projs/regMotifs/mutations_simulations_files_103sets/*_nonparametric_annotated.bed9 | awk 'BEGIN {FS=OFS="\t"}{if($4=="-") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $1 "\t" ($3) "\t" $2 "\t" $4 "/" $5 "\t" $1 "_" ($3) "_" $4 "/" $5; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "/" $5 "\t" $1 "_" ($2) "_" $4 "/" $5 }'  | sort -u -T /home/huum/projs/  > /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_id


#cat /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_paa.vep.txt2 /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_pab.vep.txt2 /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_pac.vep.txt2 /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_pad.vep.txt2 /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_pae.vep.txt2 /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_paf.vep.txt2 mutations_simulations_nonparametric_annotated.bed9.out2_pag.vep.txt2 /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_pah.vep.txt2 /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_pai.vep.txt2 /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_paj.vep.txt2 | sort -u -T /home/huum/projs/regMotifs/ > /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_all.vep.txt2

#cat /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_all.vep.txt2 |  awk '!/#/{print }'| awk '/IMPACT=HIGH/{print $0}'  > /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_all.vep.txt_HIGH

#awk   'NR==FNR{a[$1];next} (($10) in a)' /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_all.vep.txt_HIGH /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_id_v2  > /home/huum/projs/regMotifs/mutations_cohorts_output_mut/mutations_simulations_nonparametric_annotated.bed9.out2_id_HIGH

#FILES="/home/huum/projs/regMotifs/mutations_cohorts_output_mut/simulation_broad_annotated.bed9"
#/home/huum/projs/regMotifs/mutations_cohorts_output_mut/simulation_Sangerneutral_annotated.bed9
#/home/huum/projs/regMotifs/mutations_cohorts_output_mut/simulation_dkfz_annotated.bed9"
#FILES="nonparametric_annotated.bed9"
#var=$(echo $FILES | grep -oE "[^/]+$" /home/huum/projs/regMotifs/mutations_simulations_files_103sets/)


#for f in $FILES
#do	
#	cat $f"".out2.vep.txt2 | awk '!/#/{print }'| awk '/IMPACT=HIGH/{print $0}'  >  $f"".out2.vep.txt_HIGH
#	awk 'NR==FNR{a[$1];next} (($10) in a)' $f"".out2.vep.txt_HIGH $f"".out2_id > $f"".out2_id_HIGH
#var=$(echo $f | grep -oE "[^/]+$")
#	awk 'NR==FNR{a[$1$2$3$4$5];next} !(($1$2$3$4$5) in a)' $f"".out2_id_HIGH /home/huum/projs/regMotifs/mutations_simulations_files_103sets/""$var > /home/huum/projs/regMotifs/mutations_simulations_files_103sets_exclVEP/""$var""exclVEP
#done

awk 'NR==FNR{a[$1$2$3$4$5];next} !(($1$2$3$4$5) in a)' /home/huum/projs/regMotifs/mutations_cohorts_output_mut/observed_annotated_agreement_22May2017.bed9.out2_id_HIGH /home/huum/projs/regMotifs/mutations_files/observed_annotated_agreement_22May2017.bed9 > /home/huum/projs/regMotifs/mutations_files/observed_annotated_agreement_22May2017_exclVEP.bed9
#for d in mutations_simulations_files_103sets/*; do
# echo $d
# tar -czvf $d".tar.gz" $d
#done

#tar -czvf mutations_cohorts_output_randomised_files.tar.gz mutations_cohorts_output_randomised_files
#python ProcessSigElements.py meta_tumor_cohorts_v2_22May2017/cohorts_to_run_definedPCAWG 0 2000 500000
#python plotting/GeneExprAnalysis.py 
#python ProcessCohorts.py
#python TFBSAnalysisPCAWG.py
#python PlotMuts.py
#cat datafiles/Motifs/motifs_split_chr/* | intersectBed -wo -a analysis/onservedmuts_annotated.bed10 -b stdin > analysis/onservedmuts_annotated.bed10_motifs
#python AnnotateMutations.py mutations_files/observed_agreement_22May2017.bed9 datafiles/chromatin_marks_all_cells_onlynarrowpeaks analysis/observed_agreement_22May2017_annotated2.bed10
#python ProcessSigElements.py meta_tumor_cohorts_v2_22May2017/cohorts_to_run_definedPCAWG 0 2000 500000
#python convert_randomcol6tobd9.py randomized_sets50kb/ randomized_sets50kb_bed9/
#python AnnotateCandidateElements.py
#pg_ctl restart -l /home/huum/tools/local_installations/postgres9.6.1/data/logfile -D /home/huum/tools/local_installations/postgres9.6.1/data/
#sleep 30
#python RetrieveDB.py randomized_sets50kb_bed9 randomized_sets50kb_bed9 8 1 47
#pg_restore -j 8 -d regmotifs /home/huum/projs/dbdump/regmotifsdb -v 
#psql -d regmotifs -c '\d'
#wait
#python ProcessCohorts.py meta_tumor_cohorts_v2_22May2017/cohorts_to_run_definedPCAWG_ATELM

#for f in All-tumors_*_annotated.bed9; do awk 'BEGIN{FS=OFS="\t"}{if($1!~"chr"){gsub("23","X",$1); gsub("24","Y", $1); gsub("25", "M", $1); print "chr"$0} else {print;}}' "$f" | intersectBed -v -a stdin -b ../datafiles/GeneExp/gencode.v19.annotation.gtf_onlyENSEBLProteinCodingExons.bed9 > ../mutations_cohorts_output_Alltumors_noexon/"$f"; done

#cat ../regMotifs_22Mar/datafiles/chromatin_marks_all_cells_onlynarrowpeaks/*.bed | intersectBed -wo -a mutations_files/observed.bed9 -b stdin > mutations_files/observed_intersect_tracks.bed13

#python AnnotateMutations.py analysis/combined_rand23setsTFsigQval0.2_meanTFExprMotifBreaking03Filters_mergedmuts1000bpSimSig1.0localw25000onlysig0.05_merged_intersectedmuts_grouped_recurrent_aggregated0UpDwmaxdist100kb_within500kb_muts.bed9 datafiles/chromatin_marks_all_cells_onlynarrowpeaks tmp/obs.out

#python AnnotateMutations.py analysis/driver_candidates_17Mar2917.tsv.bed12.mutsonlycorrespondig_cancertypeMuts.bed9 datafiles/chromatin_marks_all_cells_onlynarrowpeaks analysis/candidate_drivers_annotated.out

echo "The process is done"
date
