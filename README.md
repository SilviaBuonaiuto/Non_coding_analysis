# Non_coding_analysis
#### Required tools
- bcftools (http://samtools.github.io/bcftools/bcftools.html)
- bedtools (https://bedtools.readthedocs.io/en/latest/)

#### 1. Split chromosomes (only accessible regions) 
```
for c in $(seq 1 22); do \
cat 20160622.allChr.mask.bed | \
grep "chr$c" > chr$c.accessible.bed; \
done 
```
#### 2. Create bed file with windows of the desired length (https://github.com/SilviaBuonaiuto/Non_coding_analysis/tree/main/scripts/create_windows.R)
```
for c in $(seq 1 22); do \
Rscript create_windows.R chr$c.accessible.bed \
window_size \
chr$c \
output_file.chr$c.tsv ; \
done
```
#### 3. Count rare and common variants in windows
- Extract info from vcf in tsv format
```
mkdir common
mkdir rare

#RARE

for c in $(seq 1 22); do \
bcftools view \
-f PASS \
-v snps \
-Q 0.001 \
chr$c.BRAVO_TOPMed_Freeze_8.vcf.gz | \
bcftools query -f'[%CHROM\t%POS\t%POS\t%REF\t%ALT\t%VRT\t%AN\t%AC\t%AF\t%Het\t%Hom\n]'\
-o rare/chr$c.rare.info.tsv; \
done

#COMMON
for c in $(seq 1 22); do \
bcftools view \
-f PASS \
-v snps \
-q 0.001 \
chr$c.BRAVO_TOPMed_Freeze_8.vcf.gz | \
bcftools query -f'[%CHROM\t%POS\t%POS\t%REF\t%ALT\t%VRT\t%AN\t%AC\t%AF\t%Het\t%Hom\n]'\
-o common/chr$c.common.info.tsv; \
done
```

- Create bed file with rare/common variants positions
```
for c in $(seq 1 22); do \
for freq in rare common; do \
cut -f 1,2,3 chr$c.$freq.info.tsv > \
$freq/chr$c.$freq.bed; \
done; \
done
```
- Intersect bed file containing windows with the one containing rare/common variants

```
for c in $(seq 1 22); do \
for freq in rare common; do \
bedtools intersect \
-a chr$c.accessible.bed \
-b chr$c.$freq.bed
-c > \
$freq/chr$c.$freq.count; \
done; \
done
```
#### 4. Use reference fasta file to count the number of nucleotides per window
- Download reference genome
- Split fasta file for each chromosome

```
for c in $(seq 1 22); do \
bedtools getfasta \
-fi reference.fa \
-bed chr$c.accessible.bed \
-fo hg38.chr$c.accessible.window_size.fa ; \
done
```
- Count nucleotides in each window and write table (https://github.com/SilviaBuonaiuto/Non_coding_analysis/tree/main/scripts/fastaAnalysis.py)

```
for c in $(seq 1 22); do \
python3 fastaAnalysis.py \
-i hg38.chr$c.accessible.window_size.fa
-o hg38.chr$c.nucleotides.tsv ; \
done 
```

#### 5. Per chromosome analysis : calculate probability of finding windows with excess of rare variants with G-test and compare observed and expected models (https://github.com/SilviaBuonaiuto/Non_coding_analysis/tree/main/scripts/gtest_perChromosome.R)

```
paste -d'\n' chromosomes.txt colors.txt | \
while read f1 && read f2; do \
Rscript gtest_perChromosome.R \
$f1.common.count \
$f1.rare.count \
hg38.$f1.nucleotides.tsv \
$f1 \ (prefix for output tables and plots)
window_size (number) \
$f2 \
```

#### 6. Whole-genome analysis : calculate probability of finding windows with excess of rare variants with G test on observed/expected data (https://github.com/SilviaBuonaiuto/Non_coding_analysis/blob/main/scripts/gtest_wholeGenome.R)

```
Rscript gtest_wholeGenome.R \
input files path \
prefix for output table and plots
```
