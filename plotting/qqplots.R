library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)


gg_qqplot <- function(ps, ci = 0.95) {
  
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df, aes(x= expected, y= observed, colour= observed > 1.30103)) +
    scale_color_manual(name = 'PC1 > 0', values = setNames(c('#56B4E9','#999999'),c(T, F)))+
    
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point() +
    #geom_point(aes(x= expected, y= observed, fill= observed > 1.30103), shape = 1, size = 3) +
    #scale_fill_manual(name = 'PC1 > 0', values = setNames(c('#56B4E9','#999999'),c(T, F)))+
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}


#directory to the cohort file
dir_elem<- '/Users/karsm933/Documents/regMotif/2020_MS_pancacer/lambda/'

data_frame_elem=c(ATELM ='All-tumors-without-Lymphatic-system-Skin-Melanoma_observed_annotated_agreement_22May2017.bed9_meanTFExprMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts200bp_statspvaluesSimSig1.0_statspvalueslocalw25000',
                  Adenocarcinoma_tumors='Adenocarcinoma-tumors_observed_annotated_agreement_22May2017.bed9_meanTFExprMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts200bp_statspvaluesSimSig1.0_statspvalueslocalw25000',
                  Carcinoma_tumors='Carcinoma-tumors_observed_annotated_agreement_22May2017.bed9_meanTFExprMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts200bp_statspvaluesSimSig1.0_statspvalueslocalw25000',
                  Digestive_tract_tumors='Digestive-tract-tumors_observed_annotated_agreement_22May2017.bed9_meanTFExprMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts200bp_statspvaluesSimSig1.0_statspvalueslocalw25000',
                  Skin_Melanoma='Skin-Melanoma_observed_annotated_agreement_22May2017.bed9_meanTFExprMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts200bp_statspvaluesSimSig1.0_statspvalueslocalw25000',
                  Liver_HCC='Liver-HCC_observed_annotated_agreement_22May2017.bed9_meanTFExprMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts200bp_statspvaluesSimSig1.0_statspvalueslocalw25000')
#data_frame_elem='CNS-PiloAstro_observed_annotated_agreement_22May2017.bed9_meanTFExprMotifBreaking03Filters_groupedbymutwithmotifinfo_mergedmuts200bp_statspvaluesSimSig1.0_statspvalueslocalw25000'
#data_frame_elem_names<-c(Pilo='CNS-PiloAstro')
data_frame_elem_names<-c( ATELM ='ATELM',
                          Adenocarcinoma_tumors='Adenocarcinoma-tumors',
                          Carcinoma_tumors='Carcinoma-tumors',
                          Digestive_tract_tumors='Digestive-tract-tumors',
                          Skin_Melanoma='Skin-Melanoma',
                          Liver_HCC='Liver-HCC')
plot_vec<-NULL
for(i in 1:2){#length(data_frame_elem)){
  print(data_frame_elem[i])
  ps <- NULL
  ps <- fread(paste0(dir_elem,data_frame_elem[i]), select = c(18), stringsAsFactors = FALSE, colClasses = 'numeric')
  ps <- as.vector(as.matrix(ps))
  plot_vec[[i]] <- gg_qqplot(ps) +
    theme_bw() +
    labs(title=data_frame_elem_names[i])+
    annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      hjust = -0.15,
      vjust = 1 + 0.15 * 3,
      label = sprintf("λ = %.2f", inflation(ps)),
      size = 6
    ) +
    theme(
      axis.ticks = element_line(size = 0.5),
      panel.grid = element_blank(),
      legend.position =  "none"
      # panel.grid = element_line(size = 0.5, color = "grey80")
    )
}





png(file=paste0(dir_elem,'plot_newwww.png'), width = 8, height = 9,units = "in", res=300)
#pdf(paste0(dir_elem,'plot.pdf'), paper="a4")
plot_grid(
  plot_vec[[1]], plot_vec[[2]],
  plot_vec[[3]], plot_vec[[4]],
  plot_vec[[5]], plot_vec[[6]],
  labels =c('(a)','(b)'), '(c)','(d)','(e)','(f)'),
            ncol = 2)
dev.off() 

png(file=paste0(dir_elem,'plot.png'), width = 8, height = 9,units = "in", res=300)
#pdf(paste0(dir_elem,'plot.pdf'), paper="a4")
plot_grid(
  plot_vec[[1]], plot_vec[[2]],
  plot_vec[[3]], plot_vec[[4]],
  plot_vec[[5]], plot_vec[[6]],
  labels =c('(a)','(b)','(c)','(d)','(e)','(f)'), ncol = 2
)
dev.off()


