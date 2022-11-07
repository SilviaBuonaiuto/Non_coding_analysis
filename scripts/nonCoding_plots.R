library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(DescTools)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <common> <rare> <nucleotide_counts> <window_size> <color> <prefix for output plot and table> ', call.=FALSE)
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
print(head(final))
### calculate rare and common frequencies
rareProb<-final %>% mutate(med = (rareCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(rareProb)
commonProb<-final %>% mutate(med = (commonCounts/(rareCounts+commonCounts))) %>% filter(rareCounts+commonCounts >0) %>% summarise(median(med))
print(commonProb)

# calculate probabilities of finding windows with excess of rare variants and plot all probabilities
window_size = args[4]
mycol = args[5]
all %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProb))) %>% filter(pbinom != 1) %>% ggplot(aes(start, -log10(1-pbinom))) +geom_point(aes(alpha = 0.3 ), color = mycol ) + theme_bw() + scale_alpha(guide = 'none') + ylim(0,6) + ggtitle(paste(code, window_size, sep = "-"))
ggsave(paste(args[6], "pbinom.png", sep = "."), width = 14, heigh = 4) 
#create table containing windows with excess of rare variants
all %>% mutate(pbinom = pbinom(rareCounts, rareCounts + commonCounts, as.numeric(rareProb))) %>% mutate(log1pbinom = -log10(1-pbinom)) %>% filter(pbinom != 1) %>% filter(log1pbinom > 4) %>% write.table(paste(args[6], "tsv", sep = "."), sep = "\t",quote = F, col.names=T, row.names=F)

#Linear regression of observed data
lmObs<-lm(rareCounts~ A+C+G+CG , data = final)
modObs = lmObs$coefficients[1]+lmObs$coefficients[2]*final$A+lmObs$coefficients[3]*final$C+lmObs$coefficients[4]*final$G+lmObs$coefficients[5]*final$CG

#Simulation of expected data and linear regression
numA = rbinom(nrow(final), as.numeric(window_size), .2)
numC = rbinom(nrow(final), as.numeric(window_size), .3)
numG = rbinom(nrow(final), as.numeric(window_size), .3)
numT = rbinom(nrow(final), as.numeric(window_size), .2)
numCG = rbinom(nrow(final), as.numeric(window_size), .1)
Y2 = rbinom(nrow(final), as.numeric(window_size), (.3*numA+.25*numC+.24*numG+.3*numT+.25*numCG)/1000)
lmExp = lm(Y2 ~ numA+numC+numG+numCG)
modExp = lmExp$coefficients[1]+lmExp$coefficients[2]*numA+lmExp$coefficients[3]*numC+lmExp$coefficients[4]*numG+lmExp$coefficients[5]*numCG

#Gtest to compare observed and expected and plot
final %>% mutate(normO = (modObs-min(modObs))/(max(modObs)-min(modObs)), normE = (modExp-min(modExp))/(max(modExp)-min(modExp))) %>% rowwise() %>% mutate(gtest = GTest(cbind(normO, normE))$p.value) %>% ggplot(aes(start, -log10(gtest))) + geom_point() + theme_bw() + ylab("-log10(p-value)")
ggsave(paste(args[6], "gtest.png", sep = "."), width = 10, heigh = 4)





