library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(magrittr)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <common> <rare> <nucleotide_counts> <path to output file> <window_size> <color> <prefix for output plot and table> ', call.=FALSE)
} 

#open input files, add column names and merge common and rare count tables
common<- read.table(args[1], sep = "\t", header = F)
rare<-read.table(args[2], sep = "\t", header = F)
nuc<-read.table(args[3], sep = "\t", header = T)
f = args[1]
code = strsplit(basename(f), ".", fixed = T) %>% sapply(extract2, 1)
colnames(common)<-c("chrom", "start", "end","commonCounts")
colnames(rare)<-c("chrom", "start", "end","rareCounts")
all<- merge(rare, common, by = c("chrom", "start", "end"))
#print(all)

# merge rare/common table and table with nucleotides numbers
final<-all %>% unite(chrpos, c("chrom", "start"), sep = ":", remove = F) %>% unite(Gene, c("chrpos", "end"), sep = "-", remove = F) %>% merge(nuc, by = "Gene") %>% select(-chrpos)
final %>% write.table(args[4], sep = "\t", quote = F, col.names = T, row.names = F)
### calculate rare and common frequencies
rareProb<-final %>% mutate(med = (rareCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(rareProb)
commonProb<-final %>% mutate(med = (commonCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(commonProb)

# calculate probabilities of finding windows with excess of rare variants and plot all probabilities
window_size = args[5]
mycol = args[6]
all %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProb))) %>% filter(pbinom != 1) %>% ggplot(aes(start, -log10(1-pbinom))) +geom_point(aes(alpha = 0.3 ), color = mycol ) + theme_bw() + scale_alpha(guide = 'none') + ylim(0,6) + ggtitle(paste(code, window_size, sep = "-"))
ggsave(paste(args[7], "pbinom.png", sep = "."), width = 14, heigh = 4) 
#create table containing windows with excess of rare variants
all %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProb))) %>% mutate(log1pbinom = -log10(1-pbinom)) %>% filter(pbinom != 1) %>% filter(log1pbinom > 4) %>% write.table(paste(args[7], "tsv", sep = "."), sep = "\t",quote = F, col.names=T, row.names=F)



