setwd("/home/delaram/delaram/bamsplit")
source("functions.R")
source("Visualize.R")


gnomeGenome = fread("/home/delaram/delaram/bamsplit/gnomad.txt", na.strings = ".", header = T, data.table = F)
colnames(gnomeGenome)[1] <- "Chr"
mybed = fread("amplicons.bed", data.table=F)
colnames(mybed) = c("chr", "start", "end")
mybed$range = mybed$end - mybed$start
bed.GRanges = bed_to_granges("amplicons.bed")
param.bed <- ScanBamParam( which=bed.GRanges )
p_param = PileupParam(max_depth = 10000, min_base_quality=13, min_mapq=5,
                      min_nucleotide_depth=1, min_minor_allele_depth=0,
                      distinguish_strands=TRUE, distinguish_nucleotides=TRUE,
                      ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=FALSE,
                      left_bins=NULL, query_bins=NULL, cycle_bins=NULL)

#save(gnomeGenome, mybed, param.bed, bed.GRanges, p_param, file="vars.RData")
gnomeGenome$End <- gnomeGenome$Start
gnomeGenome$Ref <- gnomeGenome$Alt <- ""

overlap = MakeOverlap(mybed ,gnomeGenome)

l <- list.files(".",".bam$")
bf <- mclapply(l, function(i) BamFile(i, index=paste0(i, ".bai")), mc.cores=detectCores())
names(bf) <- sub(".*_","", sub(".bam","",l))
bed.pileup <- mclapply(bf, pileup, pileupParam = p_param, scanBamParam=param.bed, mc.cores = detectCores())
freq <-  mclapply(bed.pileup, pileupFreq, mc.cores = detectCores())
names(freq) <- paste(names(bf), ".freq", sep= "")
#all.pileup <- mclapply(bf, pileupParam = p_param, mc.cores = detectCores())
freq <- mclapply(freq, clean.df, mc.cores = detectCores())
freq <- mclapply(freq, AddAllelFreq, overlap = overlap ,mc.cores = detectCores())
freq <- mclapply(freq, function(i){ i[is.na(i)] = 0 ; return(i)}, mc.cores = detectCores())


low.alf <- mclapply(freq, LowAllele.alt ,  mc.cores = detectCores())
er <- sapply(low.alf, function(i) i[length(i)*0.95]) ##  error cut-off 
low.alf.m <- setNames(melt(low.alf),c("ALT.freq","Sample"))
low.alf.m$Sample <- sub(".freq","",low.alf.m$Sample)


target <- lapply(freq, function(i) i[i$AF.exome.All<0.6 & i$AF.exome.All>0.4 & i$depth>800,] )
target.m <- setNames(melt(lapply(target,function(x)x$F2)),c("AF2","Sample"))
target.m$Sample <- sub(".freq","",target.m$Sample)

pdf("Target_Density.pdf") ## final results
ggplot(target.m, aes(AF2+1e-5, color=Sample))+geom_density()+theme_bw()#+theme(legend.position=c(1,1), legend.justification = c(1,1))
ggplot(target.m, aes(AF2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()#+theme(legend.position=c(1,1), legend.justification = c(1,1))
dev.off()


#### Trying to find distribution of the error ###############
x <- unlist(low.alf)
pdf("Mutation_Distribution_Type_allALT.pdf")
descdist(x)
descdist(log10(x[x>0]))

d <- MASS::fitdistr(log10(x[x>0]), "normal")
zero.rate <- sum(x==0)/length(x)
cat("The share of data with zero substitution error:", zero.rate * 100, "%")

tmp <- 10^rnorm(100000, d$estimate[1], d$estimate[2])
tmp.m <- setNames(melt(list(low.alf, tmp)),c("ALT.freq","Sample"))
tmp.m[is.na(tmp.m)] <- paste0("N(",round(10^d$estimate[1],4),")")
tmp.m$Sample <- sub(".freq","",tmp.m$Sample)

ggplot(tmp.m, aes(ALT.freq+1e-5, color=Sample))+geom_density()+theme_bw()+theme(legend.position=c(1,1), legend.justification = c(1,1)) + ylab("Frequency")
ggplot(tmp.m, aes(ALT.freq+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+theme(legend.position=c(1,1), legend.justification = c(1,1))
ggplot(tmp.m, aes(Sample, ALT.freq+1e-5))+geom_violin(aes(fill=Sample))+geom_boxplot(width=0.2)+theme_bw()+scale_y_log10()+theme_bw()
dev.off()

pdf("Substitution_Error_Density.pdf")
ggplot(low.alf.m, aes(AF2+1e-5, color=Sample))+geom_density()+theme_bw()+theme(legend.position=c(1,1), legend.justification = c(1,1)) + ylab("Frequency")
ggplot(low.alf.m, aes(AF2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+theme(legend.position=c(1,1), legend.justification = c(1,1))
ggplot(low.alf.m, aes(Sample, AF2+1e-5))+geom_violin(aes(fill=Sample))+geom_boxplot(width=0.2)+theme_bw()+scale_y_log10()+theme_bw()
dev.off()

#############################################################

