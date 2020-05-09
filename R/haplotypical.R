
#' Reads reference fasta
#' 
#' This reads in a the reference sequence (fasta)
#' it does so by wrapping the seqinr::read.fasta function.
#' This has been hard coded for a DNA alphabet and
#' to get character strings
#'
#' @param f (fasta file)
#' @param seqonly (whether to give just strings (TRUE) or seqinr objects (FALSE))
#' @param ... (passed to seqinr::read.fasta)
readReferenceFasta <- function(f, seqonly = FALSE, ...) {
  fa <- seqinr::read.fasta(f, as.string=FALSE, seqtype="DNA", seqonly = seqonly, ...)
}

#' reads BED format
#' 
#' Reads genomic regions in BED format
#' https://genome.ucsc.edu/FAQ/FAQformat.html#format1
#' 
#' Only the first 3 columns are considered
#' Further, the first 3 columns are defined as: chrom, start, stop
#' and start and stop are forced into integers...
#'
#' @param bedfile (name of bed file)
#' @param ... (xtra arguments passed to readr::read_tsv)
readBed <- function(bedfile, ...) {
  bed <- readr::read_tsv(bedfile, col_types=readr::cols(), ...)
  
  if(ncol(bed) < 3) {
    stop("Not bed format. Not enough columns...")
  }
  
  colnames(bed)[1:3] <- c("chrom", "start", "stop")
  bed$start <- as.integer(bed$start)
  bed$stop <- as.integer(bed$stop)
  return(bed)
}

#' gets reference fragments
#' 
#' Takes in an seqinr object and returns a character string for 
#' a given "chromosome" (including 'chrM'), and coordinates
#' (begin=b, end=e) in the reference genome. It returns
#' the sequences associated with these regions (as a vector).
#' 
#' The implementation is vectorized (c, b and e may be of length>1)
#'
#' @param s (seqinr object)
#' @param chrom (chromosomes/scaffolds/fasta identifiers)
#' @param b (begin coordinate, 0-based)
#' @param e (end coordinate, 1-based)
#' @param bed (is bed format? when FALSE 1-based indexing is assumed)
getReferenceSubstring <- function(s, chrom, b, e, bed=TRUE) {
  nseqs <- length(b)
  if (nseqs != length(e)) {
    stop("Coordinates need to come in pairs...")
  }
  
  nchroms <- length(chrom)

  if (nchroms != 1 && nchroms != nseqs) {
    stop("I need either 1 chromosome, or 1 for each coordinate...")
  }
  # bed format uses 0-bases half-open coordinates.
  # seqinr uses 1-based indexing.
  # so...
  if (bed) {
    b <- b + 1
  }
    
  seqs <- rep("", nseqs)

  for (i in seq_along(1:nseqs)) {
   if (nchroms==1) {
     seqs[[i]] <- seqinr::c2s( seqinr::getFrag(s[[ chrom[[1]] ]], b[[i]], e[[i]], as.string=TRUE) )
   } else {
     seqs[[i]] <- seqinr::c2s( seqinr::getFrag(s[[ chrom[[i]] ]], b[[i]], e[[i]], as.string=TRUE) )
   }
  }
  return(seqs)
  
}

#' Generates a single-locus haplotypes
#' 
#' This extracts ALL haplotypes from 1 region specified
#' chrom, begin and end specify the locus coordinate
#' and the bamFilename is the STRING sequence of the bam FILE
#' (not an Rsamtools object)
#' The the Haplotype, the PlusCounts and the NegCounts
#' are returned (as a data frame)
#' 
#' @param chrom (chromosome, really a chromosome/contig ID from the reference mapped to)
#' @param b (begin coordinate (vector))
#' @param e (end coordinate (vector))
#' @param bamFilename (the name of the bam file)
#' @param bed (are the locus in BED format (0-based b, 1-based e). If set to FALSE 1-based indexing is assumed)
#' @param mapqFilter (minimum mapping quality to consider...)
#' @param ... (additional arguments passed to ScanBamParam)
getSingleHaplotypePileup <- function(chrom, b, e, bamFilename, bed=TRUE, mapqFilter=3, ...) {
  b <- as.integer(b)
  e <- as.integer(e)
  if (bed==FALSE) {
    
    b <- b - 1
  }
  fields2get <- c("strand", "pos", "seq",  "cigar", "qwidth")
 # mid <- as.integer(b/2 + e/2)
  
  granger <- GenomicRanges::GRanges(chrom, ranges=IRanges::IRanges(b,e) )
  param <- Rsamtools::ScanBamParam(what=fields2get, 
                        flag=Rsamtools::scanBamFlag(),
                        mapqFilter=mapqFilter,
                        ...
  )
  
  bam <- Rsamtools::scanBam(bamFilename, param=param)[[1]] # 1 at a time, yo!
  seqs <- as.character(bam$seq)
  haps <- estimateHaplotypes( seqs, bam$cigar, bam$pos,
                              bam$qwidth, as.integer(bam$strand),
                              b, e,
                              b, e)
  return(haps)
}

#' Extacts a single locus' haplotypes
#' 
#' *This is a helper function* Do not use unless you know what you're doing!
#' 
#' This extracts the haplotypes associated with a single locus
#' The locus is pre-extracted (bam object) and the begin (b)
#' and end (e) coordinate is specified. 
#' 
#' @param bam (Rsamtools object: from scanBam, singular)
#' @param b (begin coordinate (singular))
#' @param e (end coordinate (singular!))
haplotypeHelper <- function(bam, b, e) {
  seqs <- as.character(bam$seq)
  haps <- estimateHaplotypes( seqs, bam$cigar, bam$pos,
                              bam$qwidth, as.integer(bam$strand),
                              b, e,
                              b, e)
  return(haps)
}

#' Generates a multi-locus haplotypes
#' 
#' This extracts ALL haplotypes from ALL regions specified
#' chrom, begin and end specify the locus coordinates
#' and the bamFilename is the STRING sequence of the bam FILE
#' (not an Rsamtools object)
#' if madedf==TRUE then 1 data frame is returned with:
#' The Locus, the Haplotype, the PlusCounts and the NegCounts
#' else it returns a LIST of dataframes (1 for each Locus)
#' 
#' additional arguments to this function are passed to ScanBamParam
#' (letter you specify what filtering you would prefer)
#' 
#' @param chrom (chromosome, really a chromosome/contig ID from the reference mapped to)
#' @param b (begin coordinate (vector))
#' @param e (end coordinate (vector))
#' @param bamFilename (the name of the bam file)
#' @param bed (are the locus in BED format (0-based b, 1-based e)? If set to FALSE 1-based indexing is assumed)
#' @param makedf (when FALSE, returns a list; each element corresponds to each locus; otherwise it returns a tibble) 
#' @param mapqFilter (minimum mapping quality to consider...)
#' @param ... (additional arguments passed to ScanBamParam)
getMultiHaplotypePileup <- function(chrom, b, e, bamFilename, bed=TRUE, makedf=TRUE, mapqFilter=3, ...) {
  b <- as.integer(b)
  e <- as.integer(e)
  if (bed==FALSE) {
    b <- b - 1
  }
  fields2get <- c("strand", "pos", "seq",  "cigar", "qwidth")
 # mid <- as.integer(b/2 + e/2)
  
  granger <- GenomicRanges::GRanges(chrom, ranges=IRanges::IRanges(b, e) )
  param <- Rsamtools::ScanBamParam(what=fields2get, 
                        which=granger,
                        flag=Rsamtools::scanBamFlag(), 
                        mapqFilter=mapqFilter,
                        ...
  )
  
  bam <- Rsamtools::scanBam(bamFilename, param=param)
  
  mapply(haplotypeHelper, bam=bam, b=b, e=e, SIMPLIFY=FALSE) -> foo
  if (makedf) {
    return(dplyr::bind_rows(foo, .id='Locus'))
  }
  return(foo)
}

#' bed+bam->haplotypes
#' 
#' This takes a bunch of genomic intervals in BED format
#' and a BAM file, and it extracts every haplotype
#' from each genomic region listed...
#' 
#' @param bedFile (A file in bed-format)
#' @param bamFile (an indexed BAM file)
#' @param ... (additional arguments passed to Rsamtools::ScanBamParam)
bedBamAndBeyond <- function(bedFile, bamFile, ...) {
  bedTib <- readBed(bedFile)
  getMultiHaplotypePileup(bedTib$chrom, bedTib$start, bedTib$stop,
                          bamFile,
                          ...
                          ) -> hapPiles
  return(hapPiles)
}


# TODO: Come back to this! Allow for indels in the NN
makeGraphAndAdjustPositions <- function(alleles, diffPos, diffTypes) {
  
  # 2 => insertion; all subsequent positions must be +1
  # 3 => deletion; all subsequent positions must
#  adj <- c(0,1,-1)
 # posAdj <- adj[ diffTypes ]
  

  tib <- tibble::tibble(All=alleles, Pos=diffPos, Types=diffTypes)
  
  dplyr::arrange(tib, Pos, Types) %>%
    dplyr::group_by(Pos) %>%
    dplyr::mutate(Rank=dplyr::row_number()) %>% 
    dplyr::ungroup() -> foo 
    
  forRefMod <- dplyr::filter(foo, Rank==1) # used to construct the modified reference sequence
  forGraph <-  dplyr::filter(foo, Rank>1) # used in the graph to represent differences...
  
  tib %>%
    dplyr::mutate(Adj=
               dplyr::case_when(
                 diffTypes == 1 ~ 0,  # substitution
                 diffTypes == 2 ~ 1,  # insertion
                 diffTypes == 3 ~ -1, # deletion
                 TRUE           ~ -99999
               ),
               CumAdj=base::cumsum(Adj),
               PosAdj=Pos + CumAdj
             ) -> foo
  
  return(foo$PosAdj)
}

testAlignments <- function() {
  
  ref <- "AAAACCCCCTCCCCATG" #16180 in rCRS
  diffTypes <- c(1,1,1,3)
  diffAlleles <- c("C","T","C", "-")
  diffPos <- c(4,9,10,14)
  
  graphy <- makeSequenceGraph(ref, diffAlleles, diffPos, diffTypes)
  
  getSequenceGraphEditDistance(
    alignSequenceGraph(graphy,
                       "AAACCCCCTCCCCATG"
                       )
  )

  
}


