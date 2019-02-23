setwd("/home/delaram/delaram/bamsplit")

library(MASS)
library(data.table)
library(fitdistrplus)
library(GenomicRanges)
library(IRanges)
library(rlang)
library(ggplot2)
library(Rsamtools)
library(BiocManager)
library(parallel)



############################################## 
pileupFreq <- function(pileupres) {
  nucleotides <- levels(pileupres$nucleotide)
  res <- split(pileupres, pileupres$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tablecounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(pileupres$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  res
}


############################################## 
bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  library("GenomicRanges")
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}


##############################################
clean.df <- function(df){
  df = df[,-(8:10)] # Ali: -(8:10)
  df$depth = df$A + df$C + df$G + df$T + df$N
  colnames(df) = c("chr" ,"position", "A" , "C" , "G" , "T", "N", "depth" )
  tmp <- t(apply(df[,3:7], 1, sort))
  df$F1 <- tmp[,5]/df$depth # Ali, Maximum frequency allele
  df$F2 <- tmp[,4]/df$depth # Ali, Second maximum frequency allele
  df$ALT.freq = (rowSums(tmp[,-5]))/df$depth 
  return(df)
}


###############################################
MakeOverlap <- function(mybed, gnomeExome){
  mybed$chr = gsub("chr", "", mybed$chr)
  mybed.IR = IRanges(start = as.numeric(mybed$start), 
                     end = as.numeric(mybed$end),
                     names = mybed$chr)
  mybed.GR = GRanges(mybed$chr, mybed.IR)
  
  gnomeExome.IR = IRanges(start = gnomeExome$Start, 
                          end = gnomeExome$End, 
                          names = gnomeExome$Chr)
  gnomeExome.GR = GRanges(gnomeExome$Chr, gnomeExome.IR,
                          REF = gnomeExome$Ref,
                          ALT = gnomeExome$Alt,
                          AF.exome.All = gnomeExome$gnomAD_genome_ALL) 
  
  overlap = subsetByOverlaps(gnomeExome.GR, mybed.GR)
  overlap.df = as.data.frame(overlap, row.names = seq(length(overlap)) )
  overlap.df$isSNP = overlap.df$start == overlap.df$end
  overlap.snp = overlap.df[overlap.df$isSNP==T, ]
  return(overlap.snp)
}

######################################################
AddAllelFreq<- function(df, overlap){
  df$chr = gsub("chr", "", df$chr)
  merged.df = merge(df, overlap, by.x = c('chr', 'position'), by.y = c('seqnames', 'start'), all.x = T,sort=T)
  return(merged.df)
}

###################################################### second max alt
LowAllele <- function(df, min.coverage=800, max.AF=0.0001){
  low = df[which(df$AF.exome.All<max.AF & df$depth > min.coverage),]
  # low = low[order(low$ALT.freq, decreasing = F),] #Ali 
  # y <- sapply(low.alf, sort)
  return(sort(low$F2)) 
}
#################################################### all alt
LowAllele.alt <- function(df, min.coverage=800, max.AF=0.0001){
  low = df[which(df$AF.exome.All<max.AF & df$depth > min.coverage),]
  return(sort(low$ALT.freq)) # delaram > use all alterations for error calculation(not just F2)
}

#####################################################
HighAllele <- function(df){
  high = df[which(df$AF.exome.All>0.1 & df$depth > 800),]
  high = high[order(high$ALT.freq, decreasing = F),] 
  return(high$ALT.freq)
}
######################################################

# Sets ggplot theme journal friendly
# base_size: Default size of the fonts
theme_complete_bw <- function(base_size = 25, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      axis.line =         element_blank(),
      axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "black", vjust = 1, margin=unit(c(3,3,3,3),"pt")),
      axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "black", hjust = 1,  margin=unit(c(3,3,3,3),"pt")),
      axis.ticks =        element_line(colour = "black"),
      axis.title.x =      element_text(size = base_size, vjust = 0.5, margin=unit(c(3,3,3,3),"pt")),
      axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5,  margin=unit(c(3,3,3,3),"pt")),
      axis.ticks.length = unit(0.15, "cm"),
      
      legend.background = element_rect(fill="transparent", colour = "NA"), 
      legend.key =        element_rect(fill ="transparent", colour = "NA", size = 0.25),
      legend.key.size =   unit(1.5, "lines"),
      legend.text =       element_text(size = base_size * 0.7),
      legend.title =      element_text(size = base_size * 0.8),
      legend.position =   "right",
      
      panel.background = element_rect(fill = "transparent", colour = NA), 
      panel.border =     element_rect(fill = NA, colour = "black", size=2), 
      panel.grid.major = element_line(colour = NA, size = 0.2), #"grey"
      panel.grid.minor = element_line(colour = NA, size = 0.5), #"grey"
      panel.spacing =     unit(c(0.3,.3,.3,.3), "lines"),
      
      strip.background = element_rect(fill = NA, colour = NA), 
      strip.text.x =     element_text(colour = "black", size = base_size * 0.8),
      strip.text.y =     element_text(colour = "black", size = base_size * 0.8, angle = +90),
      
      plot.background =  element_rect(colour = NA, fill = "transparent"),
      plot.title =       element_text(size = base_size*.8,margin=unit(c(5,5,5,5),"pt")),
      plot.margin =      unit(c(1,1,1,1), "lines"))
}



########################################

library(reshape2)
library(pheatmap)
library(mixtools)
fillKS <- function(ks){
  colnames(ks) = "sample1"
  ks$sample2 = NA
  ks$p_value = -1
  ks$result = NA
  for (i in 1:8){
    for (j in 1:8){
      ks[(i-1)*8 + j,1] = names(bf)[i]
      ks[(i-1)*8 + j,2] = names(bf)[j]
      ks[(i-1)*8 + j,3] = ks.test(target[[i]]$F2, target[[j]]$F2, alternative = "two.sided")[2]
      ks[(i-1)*8 + j,4] = ks.test(target[[i]]$F2, target[[j]]$F2, alternative = "two.sided")[2]<0.05
    }
  }
  return(ks)
}

