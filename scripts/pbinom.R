library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(grid)
library(gridExtra)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <common> <rare> <nucleotide_counts> <cds windows bed file> <non-coding windows bed file> <path to output file> <window_size> <color> <prefix for output plot and table> ', call.=FALSE)
} 

#open input files, add column names and merge common and rare count tables
common<- read.table(args[1], sep = "\t", header = F)
rare<-read.table(args[2], sep = "\t", header = F)
nuc<-read.table(args[3], sep = "\t", header = T)
cds<- read.table(args[4], sep = "\t", header = F)
nc<- read.table(args[5], sep = "\t", header = F)
colnames(cds)<-c("chrom", "start", "end")
colnames(nc)<-c("chrom", "start", "end")
cds$type<- "cds"
nc$type<-"non_coding"
anno<-rbind(cds,nc)
f = args[1]
code = strsplit(basename(f), ".", fixed = T) %>% sapply(extract2, 1)
colnames(common)<-c("chrom", "start", "end","commonCounts")
colnames(rare)<-c("chrom", "start", "end","rareCounts")
all<- merge(rare, common, by = c("chrom", "start", "end"))
allType<-merge(all, anno, by = c("chrom", "start", "end")) %>% distinct()

# merge rare/common table and table with nucleotides numbers
final<-allType %>% unite(chrpos, c("chrom", "start"), sep = ":", remove = F) %>% unite(Gene, c("chrpos", "end"), sep = "-", remove = F) %>% merge(nuc, by = "Gene") %>% select(-chrpos)
final %>% write.table(args[6], sep = "\t", quote = F, col.names = T, row.names = F)
### calculate rare and common frequencies
rareProb<-final %>% mutate(med = (rareCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(rareProb)
commonProb<-final %>% mutate(med = (commonCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(commonProb)

#### rare/common frequencies in CDS only
rareProbC<-final %>% filter(type == "cds") %>% mutate(med = (rareCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(rareProbC)
commonProbC<-final %>% filter(type == "cds") %>% mutate(med = (commonCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(commonProbC)

#### rare/common frequencies in non-coding regions only
rareProbN<-final %>% filter(type == "non_coding") %>% mutate(med = (rareCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(rareProbN)
commonProbN<-final %>% filter(type == "non_coding") %>% mutate(med = (commonCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(commonProbN)

# calculate probabilities of finding windows with excess of rare variants and plot all probabilities
window_size = args[7]
mycol = args[8]
final %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProb))) %>% filter(pbinom != 1) %>% ggplot(aes(start, -log10(1-pbinom))) +geom_point(aes(alpha = 0.3 ), color = mycol ) + theme_bw() + scale_alpha(guide = 'none') + ylim(0,6) + ggtitle(paste(code, window_size, sep = "-"))
ggsave(paste(args[8], "pbinom.png", sep = "."), width = 14, heigh = 4) 
#create table containing windows with excess of rare variants
final %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProb))) %>% mutate(log1pbinom = -log10(1-pbinom)) %>% filter(pbinom != 1) %>% filter(log1pbinom > 4) %>% write.table(paste(args[8], "tsv", sep = "."), sep = "\t",quote = F, col.names=T, row.names=F)

# calculate probabilities of finding windows with excess of rare variants and plot all probabilities (in CDS )
pcds<-final %>% filter(type == "cds") %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProbC))) %>% filter(pbinom != 1) %>% ggplot(aes(start, -log10(1-pbinom))) +geom_point(aes(alpha = 0.3 ), color = mycol ) + theme_bw() + scale_alpha(guide = 'none') + ylim(0,6) + ggtitle(paste(code, window_size, "CDS", sep = "-"))

# create table with significant windows
sigCds<-final %>% filter(type == "cds") %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProbC))) %>% mutate(log1pbinom = -log10(1-pbinom)) %>% filter(pbinom != 1) %>% filter(log1pbinom > 4)

# calculate probabilities of finding windows with excess of rare variants and plot all probabilities (in non coding regions)
pnc<-final %>% filter(type == "non_coding") %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProbN))) %>% filter(pbinom != 1) %>% ggplot(aes(start, -log10(1-pbinom))) +geom_point(aes(alpha = 0.3 ), color = mycol ) + theme_bw() + scale_alpha(guide = 'none') + ylim(0,6) + ggtitle(paste(code, window_size, "Non coding" ,sep = "-"))

sigNc<-final %>% filter(type == "non_coding") %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProbN))) %>% mutate(log1pbinom = -log10(1-pbinom)) %>% filter(pbinom != 1) %>% filter(log1pbinom > 4)

rbind(sigCds, sigNc) %>% write.table(paste(args[8], "cds_nc.tsv", sep = "."), sep = "\t",quote = F, col.names=T, row.names=F)

# make panel with cds and non coding plots

myplot<-grid.arrange(pcds, pnc, nrow = 2)
ggsave(paste(args[8], "cds_nc_panel.png", sep = "."), plot = myplot, width = 14, heigh = 9)

