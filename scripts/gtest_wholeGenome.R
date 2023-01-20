library(tidyr)
library(dplyr)
library(ggplot2)
library(bracer)
library(DescTools)
library(gap)
library(ggtext)
library(quantsmooth)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
stop('specify arguments: <input files directory> <prefix for output table and plots>', call.=FALSE)
}

# open single chromosomes input files and merge all chromosome in one df
#final<-read.table(args[1], sep = "\t", header = T)
final = data.frame()
fileList=glob(args[1], engine = "r")
for (f in fileList){
df<-read.table(f, sep = "\t", header = T)
final <-rbind(final,df)
}

#define window_size
window_size = 400
#predict number of nucleotides using linear model
lmExpR = lm(final$rareCounts ~ final$A+final$C+final$G+final$CG)
lmExpC = lm(final$commonCounts ~ final$A+final$C+final$G+final$CG)
#calculate variants based on previous linear model
modExpR = lmExpR$coefficients[1]+lmExpR$coefficients[2]*final$A+lmExpR$coefficients[3]*final$C+lmExpR$coefficients[4]*final$G+lmExpR$coefficients[5]*final$CG
modExpC = lmExpC$coefficients[1]+lmExpC$coefficients[2]*final$A+lmExpC$coefficients[3]*final$C+lmExpC$coefficients[4]*final$G+lmExpC$coefficients[5]*final$CG
#perform g-test and obtain p-values
chromosome = final %>% select(chrom) %>% distinct()
toRemove<-final %>% mutate(expR = modExpR, expC = modExpC) %>% filter(rareCounts == 0 & commonCounts == 0)
myd<-final %>% mutate(expR = modExpR, expC = modExpC) %>% anti_join(toRemove) %>% rowwise() %>% mutate(gtest = GTest(x = c(rareCounts, commonCounts), p = c(expR, expC)/(expR+expC))$p.value, stat = -log10(gtest))
write.table(myd, paste(args[2], "finalTable_gtest.tsv", sep = "."), sep = "\t", quote = F, col.names = T, row.names = F)

#qq plot
lambda=(gcontrol2(myd$gtest))$lambda
png(paste(args[2], "qqplot_wholeGenome.png", sep = "."))
gcontrol2(myd$gtest, col = "black")
text(2,30, expression(paste(lambda, "=", sep = " ")))
text(2.4,30, paste(round(lambda,2), sep = ""))
title("Q-Q Plot")
dev.off()

# Manhattan plot
myd$chromosome<-numericCHR(myd$chrom, prefix="chr")
df<-myd %>% select(Gene, chromosome, start, gtest)
sig_data<- df %>% subset(gtest < 0.05)
ns_data<- df %>% subset(gtest >= 0.05) %>% group_by(chromosome) %>% sample_frac(0.1)
data<- bind_rows(sig_data, ns_data)
data_cum <- data %>% group_by(chromosome) %>% summarise(max_bp = max(start)) %>% mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% select(chromosome, bp_add)
data<- data %>% inner_join(data_cum, by = "chromosome") %>% mutate(bp_cum = start + bp_add)
axis_set <- data %>% group_by(chromosome) %>% summarize(center = mean(bp_cum))
ylim<- data %>% filter(gtest == min(gtest)) %>% mutate(ylim = abs(floor(log10(gtest))) + 10) %>% pull(ylim)
sig<- 5e-8
manhplot<-ggplot(data, aes(x = bp_cum, y = -log10(gtest), color = as.factor(chromosome), size = -log10(gtest))) + geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + geom_point(alpha = 0.75) + scale_x_continuous(label = axis_set$chromosome, breaks = axis_set$center) + scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) + scale_size_continuous(range = c(0.5,3)) + labs(x = NULL, y = "-log10(p-Value)") + theme_minimal() + theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5), axis.text.y = element_text(size = 12))
ggsave(paste(args[2], "manhattanPlot.png", sep = "."), plot = manhplot, width = 15, heigh = 8)
