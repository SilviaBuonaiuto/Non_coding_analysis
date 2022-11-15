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
Y2 = rbinom(nrow(final), window_size, (final$A/window_size+final$C/window_size+final$T/window_size+final$G/window_size+final$CG/window_size)/100)
lmExp = lm(Y2 ~ final$A+final$C+final$T+final$G+final$CG)
modExp = lmExp$coefficients[1]+lmExp$coefficients[2]*final$A+lmExp$coefficients[3]*final$C+lmExp$coefficients[4]*final$G+lmExp$coefficients[5]*final$CG

#Gtest to compare observed and expected and plot
#final %>% mutate(normO = (modObs-min(modObs))/(max(modObs)-min(modObs)), normE = (modExp-min(modExp))/(max(modExp)-min(modExp))) %>% rowwise() %>% mutate(gtest = GTest(cbind(normO, normE))$p.value) %>% ggplot(aes(start, -log10(gtest))) + geom_point() + theme_bw() + ylab("-log10(p-value)")
mycol = args[3]
final %>% mutate(Exp = modExp) %>% rowwise() %>% mutate(gtest = GTest(cbind(rareCounts, Exp))$p.value) %>% ggplot(aes(start, -log10(gtest))) + geom_point(color = mycol) + theme_bw() + ylab("-log10(p-value)")
ggsave(paste(args[4], "gtest.png", sep = "."), width = 10, heigh = 4)

