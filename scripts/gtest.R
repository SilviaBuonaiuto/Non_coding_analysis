library(dplyr)
library(tidyr)
library(ggplot2)
library(DescTools)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <input file> <window_size> <color> <prefix for output plot and table> ', call.=FALSE)
} 

final<-read.table(args[1], sep = "\t", header = T)
window_size = args[2]
window_size = as.numeric(window_size)

#Y2 = rbinom(nrow(final), window_size, (final$A/window_size+final$C/window_size+final$T/window_size+final$G/window_size+final$CG/window_size)/100)

lmExpR = lm(final$rareCounts ~ final$A+final$C+final$G+final$CG)
lmExpC = lm(final$commonCounts ~ final$A+final$C+final$G+final$CG)

modExpR = lmExpR$coefficients[1]+lmExpR$coefficients[2]*final$A+lmExpR$coefficients[3]*final$C+lmExpR$coefficients[4]*final$G+lmExpR$coefficients[5]*final$CG

modExpC = lmExpC$coefficients[1]+lmExpC$coefficients[2]*final$A+lmExpC$coefficients[3]*final$C+lmExpC$coefficients[4]*final$G+lmExpC$coefficients[5]*final$CG

#Gtest to compare observed and expected and plot (sommare statistiche per expected e rare capire come e come ottnere pvalue da rappresentare)
#trovare modo di sostituire Inf nel df per poi rappresentare bene nel plot
mycol = args[3]
chromosome = final %>% select(chrom) %>% distinct()
df<-final %>% mutate(expR = modExpR, expC = modExpC) %>% rowwise() %>% mutate(gtest = GTest(x = c(rareCounts, commonCounts), p = c(expR/window_size, expC/window_size))$p.value, stat = -log10(gtest)) 
df[sapply(df, is.infinite)] <- 350
df %>% ggplot(aes(start, -log10(gtest))) + geom_point(color = mycol) + theme_bw() + ylab("-log10(p-value)") + ggtitle(paste(chromosome, window_size, sep = "-"))
ggsave(paste(args[4], "gtest.png", sep = "."), width = 10, heigh = 4)

