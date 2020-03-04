
testing <- function() {

library(Rsamtools) # TODO: manage dependencies in the package for this. good times
# needs; S4vectors, IRanges, GenomeInfoDb. BiocGenerics. parallel. Biostrings
# the dependency management for Bioconductor packages is a bit awful.
# DON"T FORGET: setRepositories(ind=1:5)

#library(GenomicAlignments)

library(readr)
library(Rcpp)

Sys.setenv("USE_CXX11" = "yes")
CXX_STD = "CXX11"


Rcpp::sourceCpp('src/Bammit.cpp')


bamFilename <- 'data/1001.sorted.bam'
bedFile <- 'data/Forenseq.txt'




#bamFilename <- 'data/CAU412M_S3_L001.sorted.softClipped.sorted.bam'
#bedFile <- 'data/PrecisionID.bed'

tib <- readr::read_tsv(bedFile, skip=0, col_names=F)
#tib <- tib[,1:3]
colnames(tib)[1:3] <- c('chrom', 'start', 'stop')

if (nrow(tib) == 0) {
  stop("Your bed file is blank!?") 
}

tibby <- tib
tib <- tib[1:3,]
  

# RSAMTOOLS
fields2get <- c('qname', "strand", "pos", "seq", "qual", "cigar", "mapq", "qwidth")

hapSum <- data.frame("Locus"=rep("", nrow(tib)), Counts=rep(-1,nrow(tib)), stringsAsFactors=FALSE)


for (i in 1:nrow(tib)) {
  row <- tib[i,]
  mid <- as.integer( row[[2]]/2 + row[[3]]/2)
  
  granger <- GRanges(row[[1]], ranges=IRanges(mid, mid) )
  param <- ScanBamParam(what=fields2get, 
                        which=granger,
                        mapqFilter=1, 
                        flag=scanBamFlag()
  )
  
  bam <- scanBam(bamFilename, param=param)[[1]] # 1 at a time, yo!
  seqs <- as.character(bam$seq)
  print(length(seqs))
  print(i)
  
  haps <- estimateHaplotypes( seqs, bam$cigar, bam$pos,
                              bam$qwidth, as.integer(bam$strand),
                              row[[2]], row[[3]],
                              row[[2]], row[[3]])
  hapSum$Locus[i] <- row[[4]]
  hapSum$Counts[i] <- sum(haps$PlusCounts) + sum(haps$NegCounts) 
  print( summary(hapSum)  )
}

chroms <- unique(tib[,1])

positionStarts <- c(16)
positionStops <- c(118)

tchroms <- c("chrM", "chrM", "chrM", "chrM", "chrM", "chr1");
tstarts <- c(100L, 121L, 200L, 2L, 599L, 10L);
tstops <- c(150L, 122L, 250L, 10L, 601L, 20L);

qchrom <- c('chrM')
qstarts <- c(122L)
qstops <- c(122L)

p <- partitionBed(qchrom, qstarts, qstops, tchroms, tstarts, tstops)



interval <- 1:200


for (i in 1:1) {
  row <- tib[i,]
  mid <- as.integer( row[[2]]/2 + row[[3]]/2)
  
  granger <- GRanges(row[[1]], ranges=IRanges(mid, mid) )
  param <- ScanBamParam(what=fields2get, 
                        which=granger,
                        mapqFilter=1, 
                        flag=scanBamFlag()
            )
  
  bam <- scanBam(bamFilename, param=param)[[1]] # 1 at a time, yo!
  seqs <- as.character(bam$seq)
  
  haps <- estimateHaplotypes( seqs, bam$cigar, bam$pos,
                      bam$qwidth, as.integer(bam$strand),
                      positionStarts, positionStops,
                      row[[2]], row[[3]])
  
  haps$Counts <- haps$PlusCounts + haps$NegCounts
  haps$LeafStats <- rep(-1L, nrow(haps) )
  haps$Degree <- rep(-1L, nrow(haps) )
  
  approximateNetworkStats(haps$Haplotype, haps$Counts, haps$LeafStats, haps$Degree)
}

}

