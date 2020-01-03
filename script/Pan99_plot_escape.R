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

plot_escape <- function(data, graphname){
  c <- signif(cor(log10(data$fit_noAb), log10(data$fit_2G12)),2)
  print (graphname)
  print (paste("R:", c, sep=''))
  textsize   <- 7
  p <-  ggplot(data,aes(x=log10(fit_noAb), y=log10(fit_2G12), color=color)) +
            geom_point(size=0.05) +
            scale_color_manual(values=c('red','blue','black')) +
            theme_cowplot(12) +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                  text = element_text(size=textsize),
                  legend.title=element_blank(),
                  legend.position="none",
                  axis.text=element_text(size=textsize,face="bold",colour='black'),
                  axis.text.x=element_text(angle = 0, hjust = 0.5, size=textsize, vjust=0.5, face="bold"),
                  axis.title=element_text(size=textsize,face="bold",colour='black')) +
            ylab("with 2G12") +
            xlab("without 2G12")
  ggsave(filename=graphname,p,height=1.5,width=1.5,dpi=300)
  }

coloring <- function(fit_noAb, fit_2G12){
  if (fit_noAb > 1 & fit_noAb < 10 & fit_2G12 > 10){return ('esp')}
  if (fit_noAb > 10){return ('fit')}
  else {return ('no')}
  }

filename <- paste('result/Pan99_mut_fit.tsv', sep='')
mydata <- read_tsv(filename) %>%
            mutate(color=mapply(coloring, fit_noAb, fit_2G12))
write.table(select(filter(mydata,color=='esp'),mut,name,num), file='result/Pan99_esp.tsv', row.names=F, col.names=F, quote=F)
write.table(select(filter(mydata,color=='fit'),mut,name,num), file='result/Pan99_fit.tsv', row.names=F, col.names=F, quote=F)
plot_escape(mydata, paste('graph/compare_condition_Pan99.png', sep=''))
