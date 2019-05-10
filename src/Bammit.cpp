#include <Rcpp.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <functional>
#include <unordered_map>
#include <utility>
#include <list>
#include <algorithm>
#include <cstring>
#include <vector>

using namespace Rcpp;
using namespace std;

#define POSITIVE_STRAND 1
#define NEGATIVE_STRAND 2

#define DELIM ";"


// [[Rcpp::plugins(cpp11)]]

class Haplotype
{


public:

  std::string read="";
  int *offsets=NULL;
  int nOffsets; // the number of regions*2; this is an array of pairs,but I keep it flattened
  int nRegions; // this is the number of regions...

  
  Haplotype(int *o, int n, const char* s) {
    
    
    nRegions = n;
    nOffsets = n + n;
    
    offsets = new int[ nOffsets ];
    memcpy(offsets, o, sizeof(int)*nOffsets);

    int buffSize=nRegions-1;
    
    int *tmp = offsets;
    for (int i = 0; i < nRegions; ++i, tmp += 2)
      buffSize += (tmp[1] - tmp[0]);
    
    read = string(buffSize+1, '\0');

    
    tmp = offsets;
    int outIndex =0;
    // the offsets are in half-open coordinates; 0-based start, 1-based stop..
    // thus no +/- 1s are needed.
    for (int i = 0; i < nRegions; ++i, tmp += 2) {
      read.replace(outIndex, tmp[1]-tmp[0], &s[tmp[0]], tmp[1]-tmp[0]);
      outIndex += tmp[1] - tmp[0];
      
       
      if (i < nRegions-1) {
        read.replace(outIndex, 1, DELIM, 1);
        ++outIndex;
      }   
    }


  }
  
  // copy constructor...
  Haplotype(const Haplotype &other) {
 
    read = other.read;
    nRegions = other.nRegions;
    nOffsets = other.nOffsets;
    offsets = new int[ nOffsets ];
    memcpy(offsets, other.offsets, sizeof(int)*nOffsets);

  }
  
  // copy constructor...
  Haplotype(Haplotype &other) {
    read = other.read;
    nRegions = other.nRegions;
    nOffsets = other.nOffsets;
    offsets = new int[ nOffsets ];
  }
  
  
  ~Haplotype() {
    if (offsets != NULL) {
      delete[] offsets;
    }
  }
  
  // TODO: Allow for padded string retrieval
  // ie, pad the string out to include both the reference and nonreference bases
  std::string
  toPaddedString() {
    return ""; 
  }

  /*
   * Two haplotypes are equal IFF
   */
  bool operator==(const Haplotype& other) const {
  
    if (nOffsets != other.nOffsets) // they have the same number of SNP records
      return false;
    
    return(read == other.read);
  }

};


// the hash function to be used with Haplotype objects
struct HaplotypeHasher
{
  std::size_t operator()(const Haplotype& k) const
  {
    return (std::hash<std::string>()(k.read));
  }
};


// global. from the sam file specificatio. bam codes that consume the reference sequence
// taken from: https://samtools.github.io/hts-specs/SAMv1.pdf
bool consumesReference[] = {true, false, true, true, false, false, false, true};
bool consumesQuery[] = {true, true, false, false, true, false, false, true};

/*
 Parses a CIGAR string
 and returns the next cigar OPERATION
 see: https://genome.sph.umich.edu/wiki/SAM
 */
const char*
cigar2op(const char *cig, int &size, int &op) {
  
  int num = 0;
  while (*cig >= '0' && *cig <= '9') { // converts a digit in a char string to its integer representation
    num = (num*10) + static_cast<int>(*cig-'0'); 
    ++cig;
  }
  size=num;
  // taken from: https://samtools.github.io/hts-specs/SAMv1.pdf, page 7
  // maybe should consider using a LUT. in practice will only go a few steps down (S at the most)...
  if (*cig == 'M')
    op = 0;
  else if (*cig == 'I')
    op = 1;
  else if (*cig == 'D')
    op = 2;
  else if (*cig == 'N')
    op = 3;
  else if (*cig == 'S')
    op = 4;
  else if (*cig == 'H')
    op = 5;
  else if (*cig == 'P')
    op = 6;
  else if (*cig == '=')
    op = 7;
  else if (*cig == 'X')
    op = 8;
  else 
    op = 255;
  
  
  return cig+1;
}

//' this parses a cigar string and an interval (positions 5-100 in the *genome* say)
//' and returns whether or not the interval exists in the string 
//' AND it fills out substringStart and substringStop with the substring
//' positions of said intervals in half-open coordinates.
//' half-open allows for 0-length intervals to be assessed 
//' (e.g., insertions that are polymorphic in populations and not found in the queried string)
bool
getIntervals(const char *ciggy, int queryStart, int queryStop, int refOffset, 
             int &substringStart, int &substringStop) {
    
    int op, size;
    const char *origCiggy = ciggy;
    // coordinates, in the reference sequence, wrt to the current cigar operation
    // these are maintained in half-open coordinates
    int refStart = refOffset-1;
    int refStop = refStart;
    
    // these are the equivalent coordinates in the substring.
    // again BED-style half-open coordinates
    int ssStart = 0;
    int ssStop = 0;
    
    // these are the substring start and stop positions (query sequence)
    // for the coordinates in the reference sequence that are sought
    substringStart = -1; // no such coordinate!
    substringStop = -1;
    
  //  Rcout << ciggy << endl;
    
    
    while (*ciggy) {
      ciggy = cigar2op(ciggy, size, op);
   //   Rcout << size << "\t" << op << "\t"  << std::endl;
      
      if (op > 8) {
        Rcerr << "Unexpected cigar operation from " << *origCiggy << std::endl;
        return false;
      }
      
      if (size < 1) // a little unkind, but let's allow it.
        continue;
      
      if (consumesReference[op]) {
        refStart = refStop;
        refStop += size;
      } else {
        refStart = refStop;
      }
      
      if (consumesQuery[op]) {
        ssStart = ssStop;
        ssStop += size;
      } else {
        ssStart = ssStop; 
      }
      
      // Softclipped bases...
      // these should not be used for determining substrings
      if (op == 4)
        continue;
      
     //  Rcout << queryStart << "\t" << queryStop << "\t" << refStart << "\t" << refStop << endl;
      
      // had originally thought that this *should* happen only once.
      // indels are tricky. Asking for the bases directly after a del
      // causes this statement to be true in the del CIGAR operation
      // and in the following CIGAR operation. You want the latter...
      if (
          refStart <= queryStart &&
          refStop >= queryStart) {
         substringStart = (queryStart - refStart) + ssStart;
      }
      
      // this assignment may happen multiple times
      // we want the maximal span
      if (refStart <= queryStop &&
          refStop >= queryStop) {
        substringStop = (queryStop - refStart) + ssStart;
      } 
      
      if (substringStop != -1 && refStop < queryStart)
        break;
      
    }
    
    return substringStop != -1 && substringStart != -1;  
  }

//' Naive and fast Hamming distance (max dist of 2)
//' 
//' This function computes a naive (and fast)
//' Hamming distance between query and target
//' It is only designed to be correct iff that distance < 3
//' and may use short circuit evaluation when the distance is >= 3
//' It optionally ignores homopolymers (e.g., AAAT == AT == distance of 0)
//' @param query (query sequence)
//' @param target (target sequence)
//' @param ignoreHomopolymers (consecutive sequences of the same letter is treated as just 1 instance of that letter)
//' @export
// [[Rcpp::export]]
int
fastCloseHammingPair(std::string &query, std::string &target, bool ignoreHomopolymers=false) {
  
  int diffs=0;
  const char *q, *t;
  
  // if one of the strings is the empty string, that 
  // makes the meat of the code more complicated.
  // thus let's treat it as a corner case.
  if (query.size() == 0) {
    if (target.size() == 0)
      return 0;
    return target.size();
  } else if (target.size() == 0) 
    return query.size();
  
  q = query.c_str();
  t = target.c_str();
  
  for ( ; *q && *t; ++q, ++t) {
    
    if (ignoreHomopolymers) {
      // if the next character is the same as the current.
      while (*q == *(q+1))
        ++q; // go to the next one...
      
      while (*t == *(t+1))
        ++t;
      
    }
    
    if (*q != *t) {
      ++diffs;
      if (diffs > 2)
        break;
    }
    
  }
  
  while (*q && diffs <= 2) {
    if (!ignoreHomopolymers || *q != *(q-1)) {
      ++diffs;
    }
    ++q;
  }
  
  while (*t && diffs <= 2) {
    if (!ignoreHomopolymers || *t != *(t-1)) {
      ++diffs;
    }
    ++t;
  }
  
  return diffs;
}

//' Returns the size of indel between 2 reads. Only single indels are supported
//'
//' This function computes a naive (and fast)
//' "indel" distance. i.e., if the two reads are separated by a single indel whose
//' length is the distance between the two strings.
//' This function is an optimization of the common case of sequence comparisons
//' between reads from massively parallel sequencing that are from the same PCR/sequencing assay
//' in which case the distance between strings is almost always 1.
//' It optionally ignores homopolymers (e.g., AAAT == AT == distance of 0)
//' The size of the indel is returned iff they are separated by an indel.
//' otherwise -1 is returned
//' @param query (query sequence)
//' @param target (target sequence)
//' @param ignoreHomopolymers (consecutive sequences of the same letter is treated as just 1 instance of that letter)
//' @export
// [[Rcpp::export]]
int
fastCloseIndelPair(std::string &query, std::string &target, bool ignoreHomopolymers=false) {
    
    const char *q, *t, *qEnd, *tEnd;
    
    // if one of the strings is the empty string, that 
    // makes the meat of the code more complicated.
    // thus let's treat it as a corner case.
    if (query.size() == 0) {
      if (target.size() == 0)
        return 0;
      return target.size();
    } else if (target.size() == 0) 
      return query.size();
    
    // first char
    q = query.c_str();
    t = target.c_str();
    // last char...
    qEnd = q + query.size() - 1;
    tEnd = t + target.size() -1;
    
    // walk left to right, stop and the first mismatch
    for ( ; *q && *t; ++q, ++t) {
      
      if (ignoreHomopolymers) {
        // if the next character is the same as the current.
        while (*q == *(q+1))
          ++q; // go to the next one...
        
        while (*t == *(t+1))
          ++t;
        
      }
      
      if (*q != *t) {
          break;
      }
      
    }
    // walk from right to left...
    for ( ; qEnd >= q && tEnd >= t; --qEnd, --tEnd) {
      
      if (ignoreHomopolymers) {
        // if the next character is the same as the current.
        while (*qEnd == *(qEnd-1) && qEnd > q)
          --qEnd; // go to the next one...
        
        while (*tEnd == *(tEnd-1) && tEnd > t)
          --tEnd;
       }
     
       if (qEnd >= q && tEnd >= t) {
         if (*qEnd != *tEnd)
           break;
       } else
         break;
      
    }
    
    /*
     * If qEnd < q
     * that means that means the 
     * the first [0..q) bases matched and
     * the last (qEnd.. END] bases matched.
     * if qEnd < q this means that the WHOLE query matched the target
     * the first q letters and the last substring starting at qEnd
     * This further means that they are separated by an deletion in the query of length
     * tEnd - t + 1
     * The same argument is had for tEnd < t
     */
    
    if (qEnd < q) {
      if (tEnd < t)
        return 0;
      return tEnd - t + 1;
    } else if (tEnd < t) {
      if (qEnd < q)
        return 0;
      return qEnd - q + 1;
    }
    // -1 means that the relationship is more complicated...
    return -1;
}

//' Computes fast string distances when the strings are close.
//'
//' This function computes a bounded, naive (and fast)
//' Hamming distance between query and bunch of sequence targets (seqs)
//' It is only designed to be correct iff that distance < 3
//' and may use short circuit evaluation when the distance is >= 3
//' It optionally ignores homopolymers (e.g., AAAT == AT == distance of 0)
//' this version is vectorized (on seqs)
//' @param query (query sequence)
//' @param seqs (target sequences)
//' @param ignoreHomopolymers (consecutive sequences of the same letter is treated as just 1 instance of that letter)
//' @export
// [[Rcpp::export]]
int
fastCloseDistances(
std::string query, 
std::vector<std::string> seqs, 
std::vector<int> result,
bool ignoreHomopolymers=false) {
  
  int diffs;
  int rIndex = 0;
  int rSize = result.size();
  
  for (auto target : seqs) {
    if (rIndex >= rSize)
      break;
    
    diffs = fastCloseHammingPair(query, target, ignoreHomopolymers);
    if (diffs > 2) {
      if (fastCloseIndelPair(query, target, ignoreHomopolymers) > 0) 
        diffs = 1;
      else
        diffs = -1;
    }
    result[rIndex] = diffs;
    ++rIndex;
  }
 
  return rIndex;  
}


//' A heuristic to estimate some network properties of PCR amplified sequences
//'
//' This function uses a heuristic to 
//' infer some very basic network properties on a bunch of sequences
//' We define a leaf to be some sequence with distance < 3 to some other sequence in the set
//' a sequence is a leaf iff:
//' If this is distance 1, then exactly one sequence is distance 1 away
//' If this is distance 2, then only one sequence is distance 2 away
//' if this property holds then the leafStats vector is populated with the index
//' of the allele (some index i in seqs where seqs[i] is the node that this allele is uniquely
//' connected to.
//' The second statistic is the degree. In networks the degree of some node (sequence in our case)
//' is the number of edges connected to that node. For us, it's the number of sequences where that
//' sequence is the argmin of the distance function (aka, a reverse nearest neighbor). 
//' As there may be ties, if some sequence is, say, distance 1 from two sequences, the 
//' node with the higher count (and lowest order amongst equal counts) is selected.
//' all of these network properties are assuming that the distance is small (Hamming < 3)
//' greater distances are evaluated separately.
//' @param seqs (these are unique sequences found at some loci)
//' @param counts (parallel array to seqs; the number of times a particular allele was observed)
//' @param leafStats (leaf indexes-- see description. This parallel to seqs
//' @param degree (same structure as leafStats, but estimates the degree of each allele)
//' @param ignoreHomopolymers (consecutive sequences of the same letter is treated as just 1 instance of that letter)
//' @export
// [[Rcpp::export]]
int
approximateNetworkStats(
  std::vector<std::string> seqs,
  std::vector<int> counts,
  Rcpp::IntegerVector leafStats,
  Rcpp::IntegerVector degree,
  int minCoverage = 10,
  bool ignoreHomopolymers=false) {
  
  int dist;
  
  struct Hit {
    int minDist;
    int secondMinDist;
    int argMin1;
    int argMin2;
  };
  
  Hit *h, hits[ seqs.size() ];
  
  hits[0].minDist= hits[0].secondMinDist=-1;
  
  for (unsigned i=0; i < seqs.size(); ++i)
    hits[i].argMin1 = hits[i].argMin2 = -1;
  
  for (unsigned i=0; i < seqs.size(); ++i)
    degree[i] = 0;
  
  
  int c;
  for (unsigned i = 0; i < seqs.size()-1; ++i) {
    c = counts[i];
    if (c >= minCoverage) { 
      for (unsigned j = i + 1; j < seqs.size(); ++j) {
        dist = fastCloseHammingPair(seqs[i], seqs[j], ignoreHomopolymers);

        if (dist > 2) {
          if (fastCloseIndelPair(seqs[i], seqs[j], ignoreHomopolymers) > 0) 
            dist = 1; // just one indel apart (regardless of length)
          else
            continue; 
        }

        h = &hits[j];
        if (dist < 3) {
          // first entry
          if (h->argMin1 == -1) {
            h->minDist = dist;
            h->argMin1 = i;
          } else if (dist < h->minDist) { // a NEW (strictly) first (best) entry
            h->secondMinDist = h->minDist;
            h->minDist = dist;
            
            h->argMin2 = h->argMin1;
            h->argMin1 = i;
            
          } else if (h->argMin2 == -1 || dist < h->secondMinDist) { // a new SECOND best entry
            h->secondMinDist = dist;
            h->argMin2 = i;
          }
        }
      }    
    } else
      break;
  }
  

  for (unsigned i = 1; i < seqs.size(); ++i) {
    leafStats[i] = -1;
    // a leaf node is defined to be:
    // EITHER:
    // there is only 1 allele that is with 2 steps from it
    // OR
    // there are two alleles that are 1 or 2 steps 
    // and the first allele is 1 step and the second allele is 2 steps
    if (hits[i].argMin1 != -1) { 
      
      
      ++degree[ hits[i].argMin1 ];
      if (hits[i].argMin2 == -1 ||
          hits[i].minDist < hits[i].secondMinDist) {
        leafStats[i] = hits[i].argMin1+1;
      }
    } else
      leafStats[i] = -2;
    
  }

  
  return seqs.size(); 
}

//' Alignment-assisted STRait Razor in R
//'
//' This is an alignment-assisted version of STRait Razor. See:
//' <https://doi.org/10.1016/j.fsigen.2013.04.005>
//' and for shameless self promotion, see:
//' <https://doi.org/10.1016/j.fsigen.2017.05.008>
//' It is used to query a given locus, usually a short tandem repeat, and for it to return the variation in and around the repeat.
//' It is typically applied to DNA that has been PCR amplified.
//' STRait Razor solves this problem by just grabbing the sequence around where varaition is known to exist. 
//' The previous versions used "anchors" (e.g., primers) to grab the sequences associated with a known locus
//' This version uses genomic coordinates. You can optionally query for the entire genomic region, which gives you functionally
//' equivalent allele calls as the prevous versions of STRait Razor. You can optionally also mask out regions (troublesome homopolymers)
//' or you can use it to associate SNP calls in the flanking regions (which are well-aligned) with the STR region itself.
//' Or even to evaluate the co-association of alleles (e.g., physical phase)
//' As input the algorithm expects DNA sequences (seqs) 
//' that have been aligned (cigars)
//' to a reference genome (seqStarts)
//' And a set of queries are given (starts, stops)
//' and all reads that span all starts,stops are assessed
//' and the co-association of alleles in these blocks is assessed.    
//' @param seqs (query sequences)
//' @param cigars (query cigar operations (as per the Bam file format))
//' @param seqStarts (the base position in 1-based indexing of the first base in the alignment)
//' @param qWidths (the size or width of the alignment in the reference genome)
//' @param starts (of the regions you're intersted in, the 0-based start coordinate)
//' @param stops (of the regions you're intersted in, the 1-based stop coordinate)
//' @param regionStart (currently ignored; specifies the amplicon start/stop)
//' @param regionStop (currently ignored; specifies the amplicon start/stop)
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame
estimateHaplotypes(Rcpp::StringVector seqs, Rcpp::StringVector cigars, Rcpp::IntegerVector seqStarts,
                   Rcpp::IntegerVector qWidths, Rcpp::IntegerVector strands,
                   NumericVector starts, NumericVector stops,
                   int regionStart, int regionStop) {

  int j, variantStart, variantStop,  searchLen, strand, qStop, qStart, len = seqs.size();
  
  if (len != cigars.size() ) {
    Rcerr << "Seqs and cigars are different sizes!" << std::endl;
    return NULL;
  }

  unordered_map<Haplotype, pair<int,int>, HaplotypeHasher> hapCounter;


  searchLen = starts.size();

 
 // for every read 
  for (int i = 0; i < len; ++i) {
    qStart = seqStarts[i];
    if (qStart == NA_INTEGER)
      continue;
    
    strand = strands[i];
   
    std::string s = as<std::string>(seqs[i]);   
    std::string c = as<std::string>(cigars[i]);

    qStop = qWidths[i] + qStart;

    int substringPositions[ searchLen + searchLen];
    std::fill_n(substringPositions, searchLen+searchLen, -1);
    
    // for every (relevant) position in the reference sequence
    for (j= 0; j < searchLen; ++j) {
       variantStart = starts[j];
       variantStop = stops[j];
       if ( variantStart > qStop) // the variant we want is AFTER this read ends
         break;
       // the variant we want *may* be in this read...
       else if (variantStart >= qStart && variantStop <= qStop) {

          int qseqStart, qseqStop;
          if ( getIntervals( c.c_str(), variantStart, variantStop, qStart,
                             qseqStart, qseqStop) ) {

            substringPositions[j+j] = qseqStart;
            substringPositions[j+j+1] = qseqStop;
          } else
            break;

       } else
         break;
    }

    /*
     * construct haplotype record here
     */
    if (j == searchLen) {
      Haplotype h(substringPositions, searchLen, s.c_str() );
      
      if (hapCounter.count(h) == 0) {
        hapCounter[h] = std::make_pair(0,0);
      }

      if (strand==POSITIVE_STRAND)
        ++(hapCounter[h].first);
      else
        ++(hapCounter[h].second);

    }
  }

  
  int nrecs = hapCounter.size();

  IntegerVector plusCounts(nrecs);
  IntegerVector negCounts(nrecs);
  
  
  std::vector<std::string> hapVecs;
  hapVecs.resize(nrecs);
  

  // get the map as a vector
  std::vector< std::pair<const Haplotype*, int> > haps ;
  haps.resize(nrecs);
  
  // avoid making deep copy of the haplotypes...
  int i = 0;
  for (auto itr = hapCounter.begin(); itr != hapCounter.end(); ++itr, ++i) {
    // haplotype associated with the sum of the + and - counts...
     haps[i] = make_pair<const Haplotype *, int>( &(itr->first), itr->second.first + itr->second.second);
  }


  // sort the vector in descending order by the sum of the number of supporting reads
  sort(haps.begin(), haps.end(), 
       [=](std::pair<const Haplotype*, int> &a, std::pair<const Haplotype*, int> &b) {
          return a.second > b.second;
       } );

  
  // and construct the parallel data structures...
  for (i = 0; i < nrecs; ++i) {
    const Haplotype *hap = haps[i].first;
    std::pair<int, int> aCounts = hapCounter[*hap];
    hapVecs[i] = hap->read;
    plusCounts[i] = aCounts.first;
    negCounts[i] = aCounts.second;
  }
  
  StringVector haplotypes(nrecs);
  haplotypes=hapVecs;

  return DataFrame::create(
    _["Haplotype"]=haplotypes, 
    _["PlusCounts"]=plusCounts, 
    _["NegCounts"]=negCounts,
    _["stringsAsFactors"] = false); 

}

struct Region {
  int chrom;
  int start; 
  int stop;
};

// used in the lower-bound
// Regions are sorted by chromosome than  by start coordinate then by stop coordinate
// This callback only uses the start coordinate (and conditioned on the same chromosome)
bool RegionCompare(const Region &a, const Region &b) {
  if (a.chrom != b.chrom)
    return a.chrom < b.chrom; 
  
  return a.start < b.start; 
}


/*
 * 
 * Returns true IFF
 * the query is properly enclosed in the target
 */
bool
isSubregion(const Region &target, const Region &query) {
  if (target.chrom != query.chrom)
    return false;
  
  if (query.start >= target.start && query.stop <= target.stop)
    return true;
  
  return false;
}

//' Partitions a bed file. 
//'
//' This function takes two genomic intervals: one for a query (Q)
//' and another for the target (T).
//' The genomic intervals are of the form: chromosome, start, stop
//' which specifies a genomic range. Genomic ranges are in half-open coordinates.
//' See: http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
//' This creates a dataframe of the genomic intersection: every element in Q
//' is paired with every overlapping element in T.
//' This can be useful for associating a set of SNPs, say, with the amplicons in which they
//' arose.
//' @param chromsQ (query chromosome)
//' @param qStart (query start coordinate. 0-based indexing)
//' @param qStop (query stop coordinate. 1-based indexing)
//' @param chromsT (target chromosome)
//' @param tStart (target start coordinate. 0-based indexing)
//' @param tStop (target stop coordinate. 1-based indexing)
//' @export
// [[Rcpp::export]]
DataFrame
partitionBed(Rcpp::StringVector chromsQ, Rcpp::IntegerVector qStarts, Rcpp::IntegerVector qStops,
             Rcpp::StringVector chromsT, Rcpp::IntegerVector tStarts, Rcpp::IntegerVector tStops) {
  

 
  int i =0 ;
  int nrecs = tStarts.size();
  
  Region targets[ nrecs ];
  

  unordered_map<std::string, int> chrom2int;
  
  int chrom2intMax=0;
  
  // create an array of targets.
  // and a simple hash to associate a chromosome string with its index
  for (i = 0; i < nrecs; ++i) {
    
    std::string s = as<std::string>(chromsT[i]);
    if (chrom2int.count(s)==0) {
      chrom2int[s] = chrom2intMax;
      ++chrom2intMax;
    }
    int chromNumber = chrom2int[s];
    
    targets[i].chrom = chromNumber;
    targets[i].start = tStarts[i];
    targets[i].stop = tStops[i];
    
  }



  // sort the bed record by chromosome...
  sort(targets, targets + nrecs, 
       [=](const Region &a, const Region &b) {

         // sort by chromosome
          if (a.chrom != b.chrom)
           return a.chrom < b.chrom; 
          if (a.start != b.start)// then by start position
            return a.start < b.start; 
          return a.stop < b.stop;// then by stop
       } );
  
    
  int nqueries = qStarts.size();
  
  // this is a sort of 2-step.
  // I need vectors to keep which queries are associated with each targets
  // To effeciently use them with R I don't use the R types
  // eg, IntegerVector
  // because their dynamic memory performance (e.g., push_back)
  // is poor.
  // thus I make vectors.
  // then convert the vectors to their Rcpp types
    std::vector<std::string> chromosomes;
  std::vector<int> starts;
  std::vector<int> stops;
  
  std::vector<int> qstarts;
  std::vector<int> qstops;
  
  
  for (i = 0; i < nqueries; ++i) {
    std::string s = as<std::string>(chromsQ[i]);
    if (chrom2int.count(s)!=0) {
      int chromNumber = chrom2int[s];
      Region q;
      q.chrom = chromNumber;
      q.start = qStarts[i];
      q.stop  = qStops[i];
      
     
      Region *tmp = std::lower_bound(targets, targets+nrecs, q, RegionCompare);
      Region *low = tmp-1;

      
      while (low >= targets) {
        if (isSubregion(*low, q)) {
          chromosomes.push_back(s);
          starts.push_back(low->start);
          stops.push_back(low->stop);
          qstarts.push_back(q.start);
          qstops.push_back(q.stop);
          
        } else if (low->chrom != q.chrom)
          break;
        else if (low->start < q.start)
          break;
        
        --low;
      }
      
      low =tmp;     
      // grab all of the records to the  right..
      // that properly enclose the query.
      // even with a sorted vector of Regions
      // the regions that overlap need not be consecutive...
      while (low < targets+nrecs) {
        if (isSubregion(*low, q)) {
          chromosomes.push_back(s);
          starts.push_back(low->start);
          stops.push_back(low->stop);
          
          qstarts.push_back(q.start);
          qstops.push_back(q.stop);
            
          //Rcout << "To the right " << low->start << "\t" << low->stop << endl;  
        } else if (low->chrom != q.chrom)
          break;
        else if (low->start > q.stop)
          break;
        ++low;
        break; 
      }


    }
    
  }
  
  nrecs = chromosomes.size();
  StringVector chromsR(nrecs);
  chromsR=chromosomes;
  
  IntegerVector startsR(nrecs);
  startsR = starts;
  
  IntegerVector stopsR(nrecs);
  stopsR = stops;
  
  IntegerVector qstartsR(nrecs);
  qstartsR = qstarts;
  
  IntegerVector qstopsR(nrecs);
  qstopsR = qstops;
  
  
  
  return DataFrame::create(_["ChromosomeQ"]=chromsR, 
                           _["TStart"]=startsR,
                              _["TStop"]=stopsR,
                              _["QStart"]=qstartsR,
                              _["QStop"]=qstopsR,
                              _["stringsAsFactors"] = false
                              );
}




