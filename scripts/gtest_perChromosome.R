library(dplyr)
library(tidyr)
library(ggplot2)
library(DescTools)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <common counts file> <rare counts file> <nucleotides file> <prefix for output tables and plots> <window_size> <plot color> ', call.=FALSE)
} 

### open tables containing rare variants, common variants, number of nucleotides and CPGs and define window size
common<- read.table(args[1], sep = "\t", header = F)
rare<-read.table(args[2], sep = "\t", header = F)
nuc<-read.table(args[3], sep = "\t", header = T)
#f = args[1]
#code = strsplit(basename(f), ".", fixed = T) %>% sapply(extract2, 1)
colnames(common)<-c("chrom", "start", "end","commonCounts")
colnames(rare)<-c("chrom", "start", "end","rareCounts")
all<- merge(rare, common, by = c("chrom", "start", "end"))

# merge rare/common table and table with nucleotides numbers
final<-all %>% unite(chrpos, c("chrom", "start"), sep = ":", remove = F) %>% unite(Gene, c("chrpos", "end"), sep = "-", remove = F) %>% merge(nuc, by = "Gene") %>% select(-chrpos)
final %>% write.table(paste(args[4], "completeTable.tsv", sep = "."), sep = "\t", quote = F, col.names = T, row.names = F)

###### define windown size
window_size = args[5]
window_size = as.numeric(window_size)

### calculate expected number of common and rare variants using linear model

lmExpR = lm(final$rareCounts ~ final$A+final$C+final$G+final$CG)
lmExpC = lm(final$commonCounts ~ final$A+final$C+final$G+final$CG)

modExpR = lmExpR$coefficients[1]+lmExpR$coefficients[2]*final$A+lmExpR$coefficients[3]*final$C+lmExpR$coefficients[4]*final$G+lmExpR$coefficients[5]*final$CG

modExpC = lmExpC$coefficients[1]+lmExpC$coefficients[2]*final$A+lmExpC$coefficients[3]*final$C+lmExpC$coefficients[4]*final$G+lmExpC$coefficients[5]*final$CG

#### Gtest to compare observed and expected variants and scatter plot representing chromosome and -log10(pvalues)

mycol = args[6]
chromosome = final %>% select(chrom) %>% distinct()
final %>% mutate(expR = modExpR, expC = modExpC) %>% rowwise() %>% mutate(gtest = GTest(x = c(rareCounts, commonCounts), p = c(expR, expC)/(expR+expC))$p.value, stat = -log10(gtest)) %>% ggplot(aes(start, stat)) + geom_point(color = mycol) + theme_bw() + ylab("-log10(p-value)") + ggtitle(paste(chromosome, window_size, sep = "-"))
ggsave(paste(args[4], "gtest.png", sep = "."), width = 10, heigh = 4)

#final %>% mutate(expR = modExpR, expC = modExpC) %>% rowwise() %>% mutate(gtest = GTest(x = c(rareCounts, commonCounts), p = c(expR, expC)/(expR+expC))$p.value, stat = -log10(gtest)) %>% filter(stat >20) %>% write.table(paste(args[4], "gtest_significant.tsv", sep = "."), sep = "\t", quote = F, col.names = T, row.names = F)

myd<-final %>% mutate(expR = modExpR, expC = modExpC) %>% rowwise() %>% mutate(gtest = GTest(x = c(rareCounts, commonCounts), p = c(expR, expC)/(expR+expC))$p.value, stat = -log10(gtest))

##### QQplot
# Uniform distribution
library(gap)
#qqunif(myd$gtest)
lambda=(gcontrol2(myd$gtest))$lambda
png(paste(args[4], "qqPlot.png", sep = "."))
gcontrol2(myd$gtest, col = "black")
text(0.5,11, expression(paste(lambda, "=", sep = " ")))
text(0.8,11, paste(round(lambda,2), sep = ""))
title(paste(chromosome, "QQ plot", sep = " "))
dev.off()
