#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(data.table)
require(cowplot)

plot_enrich_compare <- function(data, graphname){
  c <- signif(cor(log10(data$rep1), log10(data$rep2)),2)
  print (graphname)
  print (paste("R:", c, sep=''))
  textsize   <- 7
  p <-  ggplot(data,aes(x=log10(rep1), y=log10(rep2))) +
            geom_point(size=0.05) +
            theme_cowplot(12) +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                  text = element_text(size=textsize),
                  axis.text=element_text(size=textsize,face="bold",colour='black'),
                  axis.text.x=element_text(angle = 0, hjust = 0.5, size=textsize, vjust=0.5, face="bold"),
                  axis.title=element_text(size=textsize,face="bold",colour='black')) +
            ylab("Replicate 2") +
            xlab("Replicate 1")
  ggsave(filename=graphname,p,height=1.5,width=1.5,dpi=300)
  }

filename <- 'result/Pan99_mut_fit.tsv'
mydata <- read_tsv(filename)
mydata_noAb <- mydata %>%
		    rename(rep1=fit_noAb_rep1) %>%
		    rename(rep2=fit_noAb_rep2)
mydata_2G12 <- mydata %>%
		    rename(rep1=fit_2G12_rep1) %>%
		    rename(rep2=fit_2G12_rep2)
plot_enrich_compare(mydata_noAb, 'graph/compare_rep_Pan99_noAb.png')
plot_enrich_compare(mydata_2G12, 'graph/compare_rep_Pan99_2G12.png')
