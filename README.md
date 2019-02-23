# ddcfDNA
project on analyzing donor-derived cell-free DNA data of transplant patients
## R scripts description 
bam_analyzer: run a pileup on bam file. makes a frequency table of all nucleotides. need to source functions.R
functions: includes different functions needed in the bam_analyzer. 
check_alt: Quality control on bam files+ defines error+ visualization 
modifyfasta: changes fasta file's specific positions. Ref > ALt
DoubleFasta Analyzer: pipleup analysis made for bam files which are aligned to double fasta reference files(included mixture of bam_analyzer + check_alt)
