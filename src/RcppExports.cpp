// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "Haplotypical_types.h"
#include <Rcpp.h>

using namespace Rcpp;

// fastBoundedHammingRange1NN
Rcpp::IntegerVector fastBoundedHammingRange1NN(Rcpp::StringVector umis, Rcpp::IntegerVector counts, int tolerance);
RcppExport SEXP _Haplotypical_fastBoundedHammingRange1NN(SEXP umisSEXP, SEXP countsSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type umis(umisSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< int >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(fastBoundedHammingRange1NN(umis, counts, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// fastMixtureBoundedHammingGraphDist
Rcpp::List fastMixtureBoundedHammingGraphDist(Rcpp::String hammingGraph, Rcpp::StringVector haplotypes, int nInMix);
RcppExport SEXP _Haplotypical_fastMixtureBoundedHammingGraphDist(SEXP hammingGraphSEXP, SEXP haplotypesSEXP, SEXP nInMixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type hammingGraph(hammingGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type haplotypes(haplotypesSEXP);
    Rcpp::traits::input_parameter< int >::type nInMix(nInMixSEXP);
    rcpp_result_gen = Rcpp::wrap(fastMixtureBoundedHammingGraphDist(hammingGraph, haplotypes, nInMix));
    return rcpp_result_gen;
END_RCPP
}
// fastBoundedHammingGraphDist
Rcpp::IntegerVector fastBoundedHammingGraphDist(Rcpp::String hammingGraph, Rcpp::StringVector toCompare, int maxDist);
RcppExport SEXP _Haplotypical_fastBoundedHammingGraphDist(SEXP hammingGraphSEXP, SEXP toCompareSEXP, SEXP maxDistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type hammingGraph(hammingGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type toCompare(toCompareSEXP);
    Rcpp::traits::input_parameter< int >::type maxDist(maxDistSEXP);
    rcpp_result_gen = Rcpp::wrap(fastBoundedHammingGraphDist(hammingGraph, toCompare, maxDist));
    return rcpp_result_gen;
END_RCPP
}
// makeSequenceHammingGraph
Rcpp::String makeSequenceHammingGraph(Rcpp::String refRS, Rcpp::StringVector diffAllele, Rcpp::IntegerVector positions);
RcppExport SEXP _Haplotypical_makeSequenceHammingGraph(SEXP refRSSEXP, SEXP diffAlleleSEXP, SEXP positionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type refRS(refRSSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type diffAllele(diffAlleleSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type positions(positionsSEXP);
    rcpp_result_gen = Rcpp::wrap(makeSequenceHammingGraph(refRS, diffAllele, positions));
    return rcpp_result_gen;
END_RCPP
}
// fastCloseHammingPair
int fastCloseHammingPair(std::string& query, std::string& target, bool ignoreHomopolymers);
RcppExport SEXP _Haplotypical_fastCloseHammingPair(SEXP querySEXP, SEXP targetSEXP, SEXP ignoreHomopolymersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type query(querySEXP);
    Rcpp::traits::input_parameter< std::string& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< bool >::type ignoreHomopolymers(ignoreHomopolymersSEXP);
    rcpp_result_gen = Rcpp::wrap(fastCloseHammingPair(query, target, ignoreHomopolymers));
    return rcpp_result_gen;
END_RCPP
}
// fastCloseIndelPair
int fastCloseIndelPair(std::string& query, std::string& target, bool ignoreHomopolymers);
RcppExport SEXP _Haplotypical_fastCloseIndelPair(SEXP querySEXP, SEXP targetSEXP, SEXP ignoreHomopolymersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type query(querySEXP);
    Rcpp::traits::input_parameter< std::string& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< bool >::type ignoreHomopolymers(ignoreHomopolymersSEXP);
    rcpp_result_gen = Rcpp::wrap(fastCloseIndelPair(query, target, ignoreHomopolymers));
    return rcpp_result_gen;
END_RCPP
}
// fastCloseDistances
int fastCloseDistances(std::string query, std::vector<std::string> seqs, std::vector<int> result, bool ignoreHomopolymers);
RcppExport SEXP _Haplotypical_fastCloseDistances(SEXP querySEXP, SEXP seqsSEXP, SEXP resultSEXP, SEXP ignoreHomopolymersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type query(querySEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type seqs(seqsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type result(resultSEXP);
    Rcpp::traits::input_parameter< bool >::type ignoreHomopolymers(ignoreHomopolymersSEXP);
    rcpp_result_gen = Rcpp::wrap(fastCloseDistances(query, seqs, result, ignoreHomopolymers));
    return rcpp_result_gen;
END_RCPP
}
// approximateNetworkStats
int approximateNetworkStats(std::vector<std::string> seqs, std::vector<int> counts, Rcpp::IntegerVector leafStats, Rcpp::IntegerVector degree, int minCoverage, bool ignoreHomopolymers);
RcppExport SEXP _Haplotypical_approximateNetworkStats(SEXP seqsSEXP, SEXP countsSEXP, SEXP leafStatsSEXP, SEXP degreeSEXP, SEXP minCoverageSEXP, SEXP ignoreHomopolymersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type seqs(seqsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type leafStats(leafStatsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< int >::type minCoverage(minCoverageSEXP);
    Rcpp::traits::input_parameter< bool >::type ignoreHomopolymers(ignoreHomopolymersSEXP);
    rcpp_result_gen = Rcpp::wrap(approximateNetworkStats(seqs, counts, leafStats, degree, minCoverage, ignoreHomopolymers));
    return rcpp_result_gen;
END_RCPP
}
// estimateHaplotypes
Rcpp::DataFrame estimateHaplotypes(Rcpp::StringVector seqs, Rcpp::StringVector cigars, Rcpp::IntegerVector seqStarts, Rcpp::IntegerVector qWidths, Rcpp::IntegerVector strands, NumericVector starts, NumericVector stops, int regionStart, int regionStop);
RcppExport SEXP _Haplotypical_estimateHaplotypes(SEXP seqsSEXP, SEXP cigarsSEXP, SEXP seqStartsSEXP, SEXP qWidthsSEXP, SEXP strandsSEXP, SEXP startsSEXP, SEXP stopsSEXP, SEXP regionStartSEXP, SEXP regionStopSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type seqs(seqsSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type cigars(cigarsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type seqStarts(seqStartsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type qWidths(qWidthsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type strands(strandsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type starts(startsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stops(stopsSEXP);
    Rcpp::traits::input_parameter< int >::type regionStart(regionStartSEXP);
    Rcpp::traits::input_parameter< int >::type regionStop(regionStopSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateHaplotypes(seqs, cigars, seqStarts, qWidths, strands, starts, stops, regionStart, regionStop));
    return rcpp_result_gen;
END_RCPP
}
// partitionBed
DataFrame partitionBed(Rcpp::StringVector chromsQ, Rcpp::IntegerVector qStarts, Rcpp::IntegerVector qStops, Rcpp::StringVector chromsT, Rcpp::IntegerVector tStarts, Rcpp::IntegerVector tStops);
RcppExport SEXP _Haplotypical_partitionBed(SEXP chromsQSEXP, SEXP qStartsSEXP, SEXP qStopsSEXP, SEXP chromsTSEXP, SEXP tStartsSEXP, SEXP tStopsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type chromsQ(chromsQSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type qStarts(qStartsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type qStops(qStopsSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type chromsT(chromsTSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type tStarts(tStartsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type tStops(tStopsSEXP);
    rcpp_result_gen = Rcpp::wrap(partitionBed(chromsQ, qStarts, qStops, chromsT, tStarts, tStops));
    return rcpp_result_gen;
END_RCPP
}
// makeSequenceGraph
Rcpp::XPtr<GraphlineGraph> makeSequenceGraph(Rcpp::String refRS, Rcpp::StringVector diffAllele, Rcpp::IntegerVector position, Rcpp::IntegerVector etype);
RcppExport SEXP _Haplotypical_makeSequenceGraph(SEXP refRSSEXP, SEXP diffAlleleSEXP, SEXP positionSEXP, SEXP etypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type refRS(refRSSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type diffAllele(diffAlleleSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type position(positionSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type etype(etypeSEXP);
    rcpp_result_gen = Rcpp::wrap(makeSequenceGraph(refRS, diffAllele, position, etype));
    return rcpp_result_gen;
END_RCPP
}
// alignSequenceGraph
Rcpp::XPtr<AlignmentResult> alignSequenceGraph(Rcpp::XPtr<GraphlineGraph> sgrap, Rcpp::String query);
RcppExport SEXP _Haplotypical_alignSequenceGraph(SEXP sgrapSEXP, SEXP querySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<GraphlineGraph> >::type sgrap(sgrapSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type query(querySEXP);
    rcpp_result_gen = Rcpp::wrap(alignSequenceGraph(sgrap, query));
    return rcpp_result_gen;
END_RCPP
}
// getSequenceGraphEditDistance
int getSequenceGraphEditDistance(Rcpp::XPtr<AlignmentResult> alnTable);
RcppExport SEXP _Haplotypical_getSequenceGraphEditDistance(SEXP alnTableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<AlignmentResult> >::type alnTable(alnTableSEXP);
    rcpp_result_gen = Rcpp::wrap(getSequenceGraphEditDistance(alnTable));
    return rcpp_result_gen;
END_RCPP
}
// getAlignmentMissedNodes
Rcpp::IntegerVector getAlignmentMissedNodes(Rcpp::XPtr<GraphlineGraph> sgrap, Rcpp::XPtr<AlignmentResult> alnTable);
RcppExport SEXP _Haplotypical_getAlignmentMissedNodes(SEXP sgrapSEXP, SEXP alnTableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<GraphlineGraph> >::type sgrap(sgrapSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<AlignmentResult> >::type alnTable(alnTableSEXP);
    rcpp_result_gen = Rcpp::wrap(getAlignmentMissedNodes(sgrap, alnTable));
    return rcpp_result_gen;
END_RCPP
}
// graphToString
std::string graphToString(Rcpp::XPtr<GraphlineGraph> sgrap);
RcppExport SEXP _Haplotypical_graphToString(SEXP sgrapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<GraphlineGraph> >::type sgrap(sgrapSEXP);
    rcpp_result_gen = Rcpp::wrap(graphToString(sgrap));
    return rcpp_result_gen;
END_RCPP
}
// ksw2_gg_align
int ksw2_gg_align(std::string Tseq, std::string Qseq, Rcpp::IntegerVector opPos, Rcpp::IntegerVector ops, int sc_mch, int sc_mis, int gapo, int gape, bool extended);
RcppExport SEXP _Haplotypical_ksw2_gg_align(SEXP TseqSEXP, SEXP QseqSEXP, SEXP opPosSEXP, SEXP opsSEXP, SEXP sc_mchSEXP, SEXP sc_misSEXP, SEXP gapoSEXP, SEXP gapeSEXP, SEXP extendedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type Tseq(TseqSEXP);
    Rcpp::traits::input_parameter< std::string >::type Qseq(QseqSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type opPos(opPosSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ops(opsSEXP);
    Rcpp::traits::input_parameter< int >::type sc_mch(sc_mchSEXP);
    Rcpp::traits::input_parameter< int >::type sc_mis(sc_misSEXP);
    Rcpp::traits::input_parameter< int >::type gapo(gapoSEXP);
    Rcpp::traits::input_parameter< int >::type gape(gapeSEXP);
    Rcpp::traits::input_parameter< bool >::type extended(extendedSEXP);
    rcpp_result_gen = Rcpp::wrap(ksw2_gg_align(Tseq, Qseq, opPos, ops, sc_mch, sc_mis, gapo, gape, extended));
    return rcpp_result_gen;
END_RCPP
}
// seqdiffs2seq
std::string seqdiffs2seq(std::string Tseq, Rcpp::IntegerVector positions, Rcpp::IntegerVector types, Rcpp::StringVector events, int initBuff);
RcppExport SEXP _Haplotypical_seqdiffs2seq(SEXP TseqSEXP, SEXP positionsSEXP, SEXP typesSEXP, SEXP eventsSEXP, SEXP initBuffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type Tseq(TseqSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type positions(positionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type types(typesSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< int >::type initBuff(initBuffSEXP);
    rcpp_result_gen = Rcpp::wrap(seqdiffs2seq(Tseq, positions, types, events, initBuff));
    return rcpp_result_gen;
END_RCPP
}
// ksw2_seqs2seqdiffs
Rcpp::DataFrame ksw2_seqs2seqdiffs(std::string Tseq, std::string Qseq, int sc_mch, int sc_mis, int gapo, int gape);
RcppExport SEXP _Haplotypical_ksw2_seqs2seqdiffs(SEXP TseqSEXP, SEXP QseqSEXP, SEXP sc_mchSEXP, SEXP sc_misSEXP, SEXP gapoSEXP, SEXP gapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type Tseq(TseqSEXP);
    Rcpp::traits::input_parameter< std::string >::type Qseq(QseqSEXP);
    Rcpp::traits::input_parameter< int >::type sc_mch(sc_mchSEXP);
    Rcpp::traits::input_parameter< int >::type sc_mis(sc_misSEXP);
    Rcpp::traits::input_parameter< int >::type gapo(gapoSEXP);
    Rcpp::traits::input_parameter< int >::type gape(gapeSEXP);
    rcpp_result_gen = Rcpp::wrap(ksw2_seqs2seqdiffs(Tseq, Qseq, sc_mch, sc_mis, gapo, gape));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Haplotypical_fastBoundedHammingRange1NN", (DL_FUNC) &_Haplotypical_fastBoundedHammingRange1NN, 3},
    {"_Haplotypical_fastMixtureBoundedHammingGraphDist", (DL_FUNC) &_Haplotypical_fastMixtureBoundedHammingGraphDist, 3},
    {"_Haplotypical_fastBoundedHammingGraphDist", (DL_FUNC) &_Haplotypical_fastBoundedHammingGraphDist, 3},
    {"_Haplotypical_makeSequenceHammingGraph", (DL_FUNC) &_Haplotypical_makeSequenceHammingGraph, 3},
    {"_Haplotypical_fastCloseHammingPair", (DL_FUNC) &_Haplotypical_fastCloseHammingPair, 3},
    {"_Haplotypical_fastCloseIndelPair", (DL_FUNC) &_Haplotypical_fastCloseIndelPair, 3},
    {"_Haplotypical_fastCloseDistances", (DL_FUNC) &_Haplotypical_fastCloseDistances, 4},
    {"_Haplotypical_approximateNetworkStats", (DL_FUNC) &_Haplotypical_approximateNetworkStats, 6},
    {"_Haplotypical_estimateHaplotypes", (DL_FUNC) &_Haplotypical_estimateHaplotypes, 9},
    {"_Haplotypical_partitionBed", (DL_FUNC) &_Haplotypical_partitionBed, 6},
    {"_Haplotypical_makeSequenceGraph", (DL_FUNC) &_Haplotypical_makeSequenceGraph, 4},
    {"_Haplotypical_alignSequenceGraph", (DL_FUNC) &_Haplotypical_alignSequenceGraph, 2},
    {"_Haplotypical_getSequenceGraphEditDistance", (DL_FUNC) &_Haplotypical_getSequenceGraphEditDistance, 1},
    {"_Haplotypical_getAlignmentMissedNodes", (DL_FUNC) &_Haplotypical_getAlignmentMissedNodes, 2},
    {"_Haplotypical_graphToString", (DL_FUNC) &_Haplotypical_graphToString, 1},
    {"_Haplotypical_ksw2_gg_align", (DL_FUNC) &_Haplotypical_ksw2_gg_align, 9},
    {"_Haplotypical_seqdiffs2seq", (DL_FUNC) &_Haplotypical_seqdiffs2seq, 5},
    {"_Haplotypical_ksw2_seqs2seqdiffs", (DL_FUNC) &_Haplotypical_ksw2_seqs2seqdiffs, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_Haplotypical(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
