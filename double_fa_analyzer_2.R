setwd("/home/delaram/delaram/data/CareDx/CareDx_2019_Feb_23")
library(stringr)
library(data.table)
library(GenomicRanges)
library(IRanges)
library(rlang)
library(ggplot2)
library(Rsamtools)
library(BiocManager)
library(parallel)
library(plyr)
library(UpSetR)

############################################## functions 
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
  res}


clean.df <- function(df){
  df = df[,-(8:10)] 
  df$depth = df$A + df$C + df$G + df$T + df$N
  colnames(df) = c("chr" ,"position", "A" , "C" , "G" , "T", "N", "depth" )
  tmp <- t(apply(df[,3:7], 1, sort))
  df$F1 <- tmp[,5]/df$depth # Maximum frequency allele
  df$F2 <- tmp[,4]/df$depth # Second maximum frequency allele
  df$ALT.freq = (rowSums(tmp[,-5]))/df$depth 
  return(df)}

Clean<- function(df){
  df$depth = df$A + df$C + df$G + df$T
  tmp <- t(apply(df[,3:6], 1, sort))
  df$F1 <- tmp[,4]/df$depth # Maximum frequency allele
  df$F2 <- tmp[,3]/df$depth # Second maximum frequency allele
  df$ALT.freq = (rowSums(tmp[,-4]))/df$depth 
  return(df)}

####################### gnomad
gnomad_genome <-fread("/media/pgdrive/apps/annovar/dls/gnom_genome_sub.txt", na.strings = ".", header = T, data.table = F)
colnames(gnomad_genome)[1] <- "Chr"
gnomad_genome2=subset(gnomad_genome, Ref%in%c('A','T','C','G')&Alt%in%c('A','T','C','G')&gnomAD_genome_ALL>0.35&gnomAD_genome_ALL<0.65)
colnames(gnomad_genome2)[1] <- "Chr"
gnomad_genome2$Chr=paste0('chr',gnomad_genome2$Chr)

gnom.GR = GRanges(seqnames=gnomad_genome2$Chr, IRanges(start=gnomad_genome2$Start,end=gnomad_genome2$End),
                  Ref=gnomad_genome2$Ref,Alt=gnomad_genome2$Alt,AF=gnomad_genome2$gnomAD_genome_ALL)

######################### amplicon info 

amplicons=fread("/home/delaram/delaram/data/CareDx/CareDx_2019_Jan_08/bamsplit/amplicons.bed",na.strings = '.',header=F,data.table = F)
colnames(amplicons)=c('chr','start','end')
amplicons.GR = GRanges(seqnames=amplicons$chr,ranges=IRanges(start=amplicons$start,end=amplicons$end))


amp=as.data.frame(mergeByOverlaps(amplicons.GR,gnom.GR))
amp= amp[,c(1,2,3,7,14,15,16)]
colnames(amp)[1:4]=c('chr','amp.start','amp.end','gnom.pos')
ori.amp=amp

############
length(unique(paste0(ori.amp$gnom.pos,sep='_',ori.amp$chr))) #583 pos 
dim(amp[!duplicated(paste0(amp$pos,sep='_',amp$chr))|amp$AF>0.01,])# 583 again> no duplicates 
###############
amp$id=paste0(amp$chr,sep='_',amp$gnom.pos)
amp$amplicon= paste0(amp$chr,sep='_',amp$amp.start, sep='_', amp$amp.end)
amp$pos= amp$gnom.pos- amp$amp.start +1
dim(amp)
dim(subset(amp,!duplicated(id))) #583 unique positions in amplicons > AF cut-off is not enough
dim(subset(amp,!duplicated(amplicon)))
amp2 = subset(amp, pos>150 & pos< 350)
dim(amp2)
### 388 positions- cut-offs >> AF:(0.35, 0.65) - pos:(150,350) - unique amplicon names(one SNP in each position)
amp3 = subset(subset(amp2,!duplicated(amplicon)))
dim(amp3)
length(unique(amp3$amplicon))

############################# pileup section 

p_param = PileupParam(max_depth = 50000, min_base_quality=10, min_mapq=5,
                      min_nucleotide_depth=0, min_minor_allele_depth=0,
                      distinguish_strands=TRUE, distinguish_nucleotides=TRUE,
                      ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=FALSE,
                      left_bins=NULL, query_bins=NULL, cycle_bins=NULL)

l <- list.files(".",".bam$")
samples <- sub(".*_","", sub(".bam","",l))
#samples <- c('12707','12708','12709','12710','12711','12712')
bf <- mclapply(l, function(i) BamFile(i, index=paste0(i, ".bai")), mc.cores=detectCores()-2)
names(bf)<-samples
bed.pileup <- mclapply(bf, pileup, pileupParam=p_param, mc.core=detectCores()-2)
freq <- mclapply(bed.pileup, pileupFreq, mc.cores=detectCores()-2) 
names(freq) <- paste(names(bf), ".freq", sep= "")
freq <- mclapply(freq, clean.df, mc.cores = detectCores()-2)
sapply(1:length(freq), function(i)colnames(freq[[i]])[1]<<-'Contig')
tmp <- mclapply(freq, function(i) {x=data.frame(str_split_fixed(i$Contig, "_", 3));colnames(x)=c('chr','start','end');return(x)})
freq <- sapply(1:length(freq), function(i)cbind(tmp[[i]],freq[[i]]),simplify =F)
sapply(1:length(freq), function(i)freq[[i]]$T.pos<<-as.numeric(as.character(freq[[i]]$start))+as.numeric(as.character(freq[[i]]$position))-1)
sapply(1:length(freq), function(i)freq[[i]]$chr<<- as.character(freq[[i]]$chr))
m <- sapply(1:length(freq), function(i)merge(freq[[i]],amp3,by.x =c('chr', 'T.pos'),by.y = c('chr', 'gnom.pos'),all.x=F,all.y=T,sort=F),simplify=F)
names(m)=samples

m <- sapply(1:length(m),function(i) m[[i]][,!colnames(m[[i]]) %in% c('amp.start','pos')],simplify = F)
sapply(1:length(m),function(i){x=as.data.frame(str_split_fixed(m[[i]]$Contig, "_", 4));m[[i]]$fasta<<-x$V4})
sapply(1:length(m),function(i){ m[[i]]$fasta <<-substr(as.character(m[[i]]$fasta),1,2)})
## ~15 positions in amp3 file are not included in freq 
lapply(m, function(i) paste(sum(is.na(i$A)),sum(is.na(i$T)),sum(is.na(i$C)),sum(is.na(i$G)))) 
lapply(1:length(m), function(i){m[[i]]$start<<-as.numeric(as.character(m[[i]]$start));
m[[i]]$end<<-as.numeric(as.character(m[[i]]$end))
m[[i]]$Contig<<-as.numeric(as.character(m[[i]]$Contig))
m[[i]]$position<<-as.numeric(as.character(m[[i]]$position))})
lapply(m, function(i) paste(sum(is.na(i$A)),sum(is.na(i$T)),sum(is.na(i$C)),sum(is.na(i$G)))) 
m = lapply(m, function(x){replace(x, is.na(x), 0)} )
lapply(m, function(i) paste(sum(is.na(i$A)),sum(is.na(i$T)),sum(is.na(i$C)),sum(is.na(i$G)))) 
saveRDS(m,'/media/pgdrive/users/delaram/data/DoubleFastaTab_illumina_Feb23.rds')

#m = readRDS('/media/pgdrive/users/delaram/data/DoubleFastaTab2.rds')

############  #pos with more strict cut-offs 
m1 = lapply(m, subset, AF>0.4 & AF<0.6) ## filters NA as well 
lapply(m1, nrow)
m1 = lapply(m1, subset, position> 180 & position< 220)
lapply(m1, nrow)
test = lapply(1:length(m1), function(i) unique.data.frame(m1[[i]][,c('chr','T.pos')]))
lapply(test, nrow) ## ~350 less but better

################################ ddply + add attribute 
m.s = lapply(m, subset, select=c('chr','T.pos','A','T','C','G'))
lapply(1:length(m.s), function(i){m.s[[i]]$T.pos<<-as.character(m.s[[i]]$T.pos); m.s[[i]]$chr<<-as.character(m.s[[i]]$chr)})
lapply(m.s, nrow)
m.f = lapply(m.s,ddply, .(chr, T.pos), numcolwise(sum))
lapply(m.f, nrow) ## same as unique >> correct > 388
lapply(1:length(m.f), function(i) m.f[[i]]$Sample<<-samples[i])


data=do.call(rbind, m.f)
# before merge: total #pos=2328       after merge: each sample #pos=2328 >> no multiple ALt?
data.f=merge(data,amp3,by.x =c('chr', 'T.pos'),by.y = c('chr', 'gnom.pos'),all.x=T,all.y=F,sort=F) 
paste0(nrow(data.f), sep=' ', sum(table(data.f$Alt))) #### NO NA value(all pos in gnomad)
data.f= Clean(data.f)

data.f$m1=NA; data.f$m2=NA
sapply(1:nrow(data.f), function(i){
  tmp=t(apply(data.f[i,3:6],1,sort)) 
  data.f[i,]$m1<<-colnames(tmp)[4]
  ifelse(data.f[i,'ALT.freq']>0, data.f[i,]$m2<<-colnames(tmp)[3] , data.f[i,]$m2<<-NA)})

data.f$alt.C <- sapply(1:nrow(data.f),function(i){i=data.f[i,];base=c('A','C','G','T'); l =!base%in%c(i$Ref); do.call(sum,i[base[l]]) } ) / data.f$depth
data.f$E <- sapply(1:nrow(data.f),function(i){i=data.f[i,];base=c('A','C','G','T'); l =!base%in%c(i$Alt, i$Ref); do.call(sum,i[base[l]]) } ) / data.f$depth
data.f$ref_error <- !(data.f$m1==data.f$Ref & (data.f$m2==data.f$Alt|is.na(data.f$m2)) )& !(data.f$m1==data.f$Alt & (data.f$m2==data.f$Ref|is.na(data.f$m2) ) )
data.f$c_error <- (data.f$alt.C>0.05 & data.f$alt.C<0.35)|(data.f$alt.C>0.65 & data.f$alt.C<0.95)

saveRDS(data.f,'/media/pgdrive/users/delaram/data/dFa_finalTAB_illumina_feb23.rds')
#data.f = readRDS('/media/pgdrive/users/delaram/data/dFasta_finalTAB.rds')

## many are removed since we do not have all the 388 positions as a result of
### hard filtering in base calling step
data.f2 = subset(data.f, depth>400)  

data.f2$label=sapply(1:nrow(data.f2), function(i){
  if(data.f2[i,'ref_error'] & data.f2[i,'c_error'])return('both')
  else if(data.f2[i,'ref_error']) return('ref-alt error')
  else if(data.f2[i,'c_error']) return('CareDx error')
  else return('non')}) 

as.data.frame(table(data.f2$label))
data.f.er <- data.f2[ !(data.f2$m1==data.f2$Ref & (data.f2$m2==data.f2$Alt|is.na(data.f2$m2))) & !(data.f2$m1==data.f2$Alt & (data.f2$m2==data.f2$Ref|is.na(data.f2$m2))) ,]
data.f.cor <- data.f2[ (data.f2$m1==data.f2$Ref & (data.f2$m2==data.f2$Alt| is.na(data.f2$m2))) | (data.f2$m1==data.f2$Alt & (data.f2$m2==data.f2$Ref|is.na(data.f2$m2))) ,]


pdf('DoubleFastaPlot_illumina_feb23.pdf')
ggplot(data.f2, aes(F2+1e-5, color=Sample))+geom_density()+theme_bw()+ggtitle('total data')
ggplot(data.f2, aes(F2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+ggtitle('total data(log)')
ggplot(data.f.cor, aes(F2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+ggtitle('non.ref.error(log) subset')
#ggplot(data.f.er, aes(F2+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+ggtitle('ref.error(log) subset')
ggplot(data.f2, aes(y=F2+1e-5, x=Sample))+scale_y_log10()+geom_violin(aes(fill=Sample))+geom_boxplot(width=0.17)+theme_bw()
ggplot(data.f2, aes(y=F2+1e-5, x=Sample))+scale_y_log10()+geom_boxplot(width=0.7,aes(fill=Sample))+theme_bw()

ggplot(data.f2, aes(E+1e-5, color=Sample))+geom_density()+theme_bw()+ggtitle('total data')
ggplot(data.f2, aes(E+1e-5, color=Sample))+geom_density()+scale_x_log10()+theme_bw()+ggtitle('total data(log)')


df.sub = split(data.f2,data.f2$Sample)
lapply(df.sub, nrow) ## number of targeted positions ~ 360 remained

df.sub = lapply(df.sub,  function(i) i[order(i$alt.C),])
df.sub = lapply(df.sub, function(i){ i$index=1:nrow(i); i})
#pdf('depth.filter_corrected.pdf')
lapply(df.sub, function(t){ggplot(t, aes(x=index,y=alt.C, color=label))+geom_point(size=0.9)+
    theme_bw()+ggtitle(t[1,'Sample'])+
    geom_text(hjust = 0, nudge_x = 0.05,check_overlap = TRUE,size=1.6,colour='black',
              aes(label=ifelse(c_error,paste0(chr,sep='_',amp.start,sep='_',amp.end),'')))
})

#dev.off()


data.f3=do.call(rbind,df.sub)
data.f3$Sample = as.character(data.f3$Sample)
ggplot(data.f3, aes(x=index,y=alt.C, color=Sample))+geom_point(size=0.9)+theme_bw()+scale_color_brewer(palette="Accent")+ggtitle('total data')


x = lapply(df.sub, function(i) i[i$c_error,])
x = lapply(x, function(i) sapply(1:nrow(i), function(j) paste(i[j,'chr'],i[j,'T.pos'],sep='_',collapse = NULL)))
uni = unique(unlist(x))
tab = as.data.frame(uni)
tab[,2:7] =sapply(1:length(x),function(j) {sapply(1:nrow(tab), function(i) ifelse( tab[i,'uni']%in% x[[j]],1, 0))})
colnames(tab)[1]=c('posID')
colnames(tab)[2:7]= names(x)

upset(tab,sets = names(x) , 
      matrix.color = "#990000",mainbar.y.label = 'CareDx error intersection',
      sets.bar.color = c('darkslategrey','cyan4','cyan3','cyan3','cyan2','cyan1'), 
      sets.x.label = c('(alt.C>0.05 & alt.C<0.35)|(alt.C>0.65 & alt.C<0.95)'),
      keep.order = F,nintersects=NA, point.size = 2.6,line.size = 0.7)

#################
data.f3$ref_error=ifelse(data.f3$ref_error,1,0)
data.f3$c_error=ifelse(data.f3$c_error,1,0)
tmp= data.f3[,colnames(data.f3)%in%c('chr','T.pos','alt.C','E','ref_error','Sample','c_error')]
tmp=split.data.frame(tmp,f=tmp$Sample,drop=FALSE)

n=lapply(tmp, function(i){colnames(i)[4:7]=paste(colnames(i)[4:7],sep='_',i[1,'Sample'],collapse=NULL)})
sapply(1:length(tmp), function(i)colnames(tmp[[i]])[4:7]<<-n[[i]] )
sapply(1:length(tmp),function(i) tmp[[i]]$name<<-paste(tmp[[i]]$chr,sep='_',tmp[[i]]$T.pos,collapse=NULL))
sapply(1:length(tmp), function(i) tmp[[i]]<<-tmp[[i]][,-c(1,2,3)])
sapply(1:length(tmp),function(i)rownames(tmp[[i]])<<-tmp[[i]]$name)
name=sapply(1:length(tmp),function(i)tmp[[i]]$name)
ints=Reduce(intersect,name) # 364 positiosn present in all of the samples 
sapply(1:length(tmp), function(i) tmp[[i]]<<-tmp[[i]][,colnames(tmp[[i]])!='name'])

tmp=do.call(cbind,lapply(1:length(tmp), function(i) tmp[[i]][rownames(tmp[[i]])%in%ints,]))
tmp$ref_error.Sum=rowSums(tmp[,colnames(tmp)%in%paste0('ref_error_',samples)])
tmp$c_error.Sum=rowSums(tmp[,colnames(tmp)%in%paste0('c_error_',samples)])
#saveRDS(object = tmp, file = '/home/delaram/delaram/data/CareDx/DoubleFastaErrorTab.rds')
e=subset(tmp,ref_error.Sum>0)
ec=subset(tmp,c_error.Sum>0)
#saveRDS(ec,file= '/home/delaram/delaram/data/CareDx/careDxErrors_SingleFa_feb12.rds')
#ec = readRDS('/home/delaram/delaram/data/CareDx/careDxErrors.rds')
#ggplot(tmp, aes(x=ref_error.Sum, y=c_error.Sum))+geom_point(alpha = 1/5,colour='blue')#geom_bin2d(bins=8)

tmp.p <- ddply(tmp, .(ref_error.Sum, c_error.Sum), summarise, num=length(ref_error.Sum))
ggplot(tmp.p, aes(ref_error.Sum, c_error.Sum))+geom_point(aes(size=num))+scale_size_continuous(trans="log2")

par(mfrow=c(1,2))
barplot(table(tmp$c_error.Sum),main = 'c_error',col='red')
barplot(table(tmp$ref_error.Sum),main='ref_error',col='blue')
par(mfrow=c(1,1))
e.tmp=tmp[,colnames(tmp)%in%paste0('E_',samples)]
thr=1
#e.tmp=subset(e.tmp,E_12707<thr&E_12708<thr&E_12709<thr&E_12710<thr&E_12711<thr&E_12712<thr)
plot(e.tmp)
alt.c.tmp=tmp[,colnames(tmp)%in%paste0('alt.C_',samples)]
plot(alt.c.tmp,col='dark blue')


dev.off()

####### QC
nrow(amplicons)*2  ## total number of initial Contigs (double fasta)
nrow(amplicons)  ## total number of initial Contigs (single fasta)
lapply(1:length(bed.pileup), function(i) length(unique(bed.pileup[[i]]$seqnames))) ## number of unique contigs in bam file
lapply(1:length(bed.pileup), function(i) nrow(amplicons)*2 -length(unique(bed.pileup[[i]]$seqnames))) ##contigs with 0 aligned read(double)
lapply(1:length(bed.pileup), function(i) nrow(amplicons) -length(unique(bed.pileup[[i]]$seqnames))) ##contigs with 0 aligned read (single)
lapply(1:length(freq), function(i) length(unique(freq[[i]]$Contig))) ## same num as above
lapply(1:length(freq), function(i) nrow(freq[[i]])/length(unique(freq[[i]]$Contig))) ## ~90 positions for each contig 
lapply(1:length(freq), function(i) table(freq[[i]]$Contig) )  ## exact 
lapply(freq, function(i) sum(i$depth==0)) ## num positions with zero depth(not removed)
lapply(1:length(freq), function(i) table(round(freq[[i]]$ALT.freq,2)))  ## many ~0 ALT.freq

lapply(1:length(freq), function(i) paste0(nrow(freq[[i]]),sep=' ',nrow(m[[i]]))) # more pos after merge> many ALT positions in gnomad
sapply(1:length(m),function(i) table(as.character(m[[i]]$fasta))/nrow(m[[i]]))#ratio of aligned read to 2 ref files

lapply(1:length(m), function(i) range(as.numeric(m[[i]]$position)) ) #range of snp positions in amplicons 

test = lapply(m, subset, AF>0.4 & AF<0.6)
lapply(test, dim)
test2 = lapply(test, subset, position> 150 & position< 250)
lapply(test2, dim)
test3 = lapply(1:length(test2), function(i) unique.data.frame(test2[[i]][,c('chr','T.pos')]))
lapply(test3, dim)

