# awk 'sub(/^>/, "")' IDT405amplicons.fa > amplicons.bed

setwd("C:/Users/Delaram/Desktop/cfdna_anal1/")
amplicons = read.table("amplicons.bed")

library(stringr)
df = as.data.frame(str_split_fixed(amplicons$V1, "_", 3))
head(df)
colnames(df) = c("chr", "start", "end")
write.table(df, file = "C:/Users/Delaram/Desktop/cfdna_anal1/amplicons.bed",  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t" )


x <- GRanges("chr1", IRanges(c(2, 2) , c(7, 19)), strand=c("+", "-"))
y <- GRanges("chr1", IRanges(1, 10), strand="-") 
#library(GenomicRanges)
x
y
intersect(x, y, ignore.strand=TRUE)


c = as.factor(0.4)
as.numeric(as.character(c))
as.integer(c)
