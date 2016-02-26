## Get coverage data for the different alignments
## Christian P
## 26 February 2016

### Get the genes file, just have to do this once
#grep -P '\tgene\t' /proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7.gff > genes.txt
#cat genes.txt | cut -f 1,4,5 > genes.bed # 9 is the metadata

## Get the coverage data
bedtools coverage -a genes.bed -b ../aln/F4R4F.realn.bam -d > F.res
echo "F4R4F"
Rscript cov_analysis.r F.res

bedtools coverage -a genes.bed -b ../aln/J4N9K.realn.bam -d > J.res
echo "J4N9K"
Rscript cov_analysis.r J.res

bedtools coverage -a genes.bed -b ../aln/O7Q0M.realn.bam -d > O.res
echo "O7Q0M"
Rscript cov_analysis.r O.res

bedtools coverage -a genes.bed -b ../aln/P6K2I.realn.bam -d > P.res
echo "P6K2I"
Rscript cov_analysis.r P.res

bedtools coverage -a genes.bed -b ../aln/V8D3K.realn.bam -d > V.res
echo "V8D3K"
Rscript cov_analysis.r V.res

bedtools coverage -a genes.bed -b ../aln/Z3V0Y.realn.bam -d > Z.res
echo "Z3V0Y"
Rscript cov_analysis.r Z.res
