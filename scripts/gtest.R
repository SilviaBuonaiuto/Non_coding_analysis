library(dplyr)
library(tidyr)
library(ggplot2)
library(DescTools)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <input file> <window_size> <color> <prefix for output plot and table> ', call.=FALSE)
} 

### open table containing rare and common variants counts and nucleotides and CPGs counts and define window size
final<-read.table(args[1], sep = "\t", header = T)
window_size = args[2]
window_size = as.numeric(window_size)

### calculate expected number of common and rare variants using linear model

lmExpR = lm(final$rareCounts ~ final$A+final$C+final$G+final$CG)
lmExpC = lm(final$commonCounts ~ final$A+final$C+final$G+final$CG)

modExpR = lmExpR$coefficients[1]+lmExpR$coefficients[2]*final$A+lmExpR$coefficients[3]*final$C+lmExpR$coefficients[4]*final$G+lmExpR$coefficients[5]*final$CG

modExpC = lmExpC$coefficients[1]+lmExpC$coefficients[2]*final$A+lmExpC$coefficients[3]*final$C+lmExpC$coefficients[4]*final$G+lmExpC$coefficients[5]*final$CG

#### Gtest to compare observed and expected and plot

mycol = args[3]
chromosome = final %>% select(chrom) %>% distinct()
final %>% mutate(expR = modExpR, expC = modExpC) %>% rowwise() %>% mutate(gtest = GTest(x = c(rareCounts, commonCounts), p = c(expR, expC)/(expR+expC))$p.value, stat = -log10(gtest)) %>% ggplot(aes(start, stat)) + geom_point(color = mycol) + theme_bw() + ylab("-log10(p-value)") + ggtitle(paste(chromosome, window_size, sep = "-"))
ggsave(paste(args[4], "gtest.png", sep = "."), width = 10, heigh = 4)

df %>% filter(stat > 20) %>% write.table(paste(args[4], "gtest_significant.tsv", sep = "."), sep = "\t", quote = F, col.names = T, row.names = F)

##### QQplot