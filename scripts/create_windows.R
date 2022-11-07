library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <accessible_regions_bedFile> <window_size> <chromosome> <outputFile> ', call.=FALSE)
}

#open bed file containing accessible regions coordinates and rename columns
mask<-read.table(args[1], sep = "\t", header = F)
colnames(mask)<-c("chrom", "start", "end", "comment")
#assign windows size and chromosome name
dimension = args[2]
chromosome = args[3]
#evaluate regions length, select the ones with length > than window size and put them in vectors
regions<-mask %>% mutate(length = end-start) %>% filter(length >= as.numeric(dimension)) %>% rowwise() %>% mutate(my_seq = map2(start, end, seq, by=as.numeric(dimension))) %>% pull(my_seq)
#create dataframe with windows and write in a tab separated file (bed file)
data.frame(unlist(regions)) %>% rename(start = "unlist.regions.") %>% mutate(CHROM = paste(chromosome), end = start+as.numeric(dimension)) %>% select(CHROM, start, end) %>% write.table(args[4], sep = "\t", quote = F, col.names = F, row.names = F)
