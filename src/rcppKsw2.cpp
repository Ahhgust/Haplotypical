#include <string.h>
#include <string>
#include <stdio.h>
#include <Rcpp.h>
#include <stdint.h>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

/*
 *
 * This is the header file KSW2.h
 *
 */

#define KSW_NEG_INF -0x40000000

#define KSW_EZ_SCORE_ONLY  0x01 // don't record alignment path/cigar
#define KSW_EZ_RIGHT       0x02 // right-align gaps
#define KSW_EZ_GENERIC_SC  0x04 // without this flag: match/mismatch only; last symbol is a wildcard
#define KSW_EZ_APPROX_MAX  0x08 // approximate max; this is faster with sse
#define KSW_EZ_APPROX_DROP 0x10 // approximate Z-drop; faster with sse
#define KSW_EZ_EXTZ_ONLY   0x40 // only perform extension
#define KSW_EZ_REV_CIGAR   0x80 // reverse CIGAR in the output
#define KSW_EZ_SPLICE_FOR  0x100
#define KSW_EZ_SPLICE_REV  0x200
#define KSW_EZ_SPLICE_FLANK 0x400

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    uint32_t max:31, zdropped:1;
    int max_q, max_t;      // max extension coordinate
    int mqe, mqe_t;        // max score when reaching the end of query
    int mte, mte_q;        // max score when reaching the end of target
    int score;             // max score reaching both ends; may be KSW_NEG_INF
    int m_cigar, n_cigar;
    int reach_end;
    uint32_t *cigar;
  } ksw_extz_t;

  /**
  * NW-like extension
  *
  * @param km        memory pool, when used with kalloc
  * @param qlen      query length
  * @param query     query sequence with 0 <= query[i] < m
  * @param tlen      target length
  * @param target    target sequence with 0 <= target[i] < m
  * @param m         number of residue types
  * @param mat       m*m scoring mattrix in one-dimension array
  * @param gapo      gap open penalty; a gap of length l cost "-(gapo+l*gape)"
  * @param gape      gap extension penalty
  * @param w         band width (<0 to disable)
  * @param zdrop     off-diagonal drop-off to stop extension (positive; <0 to disable)
  * @param flag      flag (see KSW_EZ_* macros)
  * @param ez        (out) scores and cigar
  */
  void ksw_extz(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);

  void ksw_extz2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                     int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

  void ksw_extd(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2, int w, int zdrop, int flag, ksw_extz_t *ez);

  void ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                     int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);

  void ksw_exts2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                     int8_t gapo, int8_t gape, int8_t gapo2, int8_t noncan, int zdrop, int flag, ksw_extz_t *ez);

  void ksw_extf2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t mch, int8_t mis, int8_t e, int w, int xdrop, ksw_extz_t *ez);

  /**
  * Global alignment
  *
  * (first 10 parameters identical to ksw_extz_sse())
  * @param m_cigar   (modified) max CIGAR length; feed 0 if cigar==0
  * @param n_cigar   (out) number of CIGAR elements
  * @param cigar     (out) BAM-encoded CIGAR; caller need to deallocate with kfree(km, )
  *
  * @return          score of the alignment
  */
  int ksw_gg(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);
  int ksw_gg2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);
  int ksw_gg2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

  void *ksw_ll_qinit(void *km, int size, int qlen, const uint8_t *query, int m, const int8_t *mat);
  int ksw_ll_i16(void *q, int tlen, const uint8_t *target, int gapo, int gape, int *qe, int *te);

#ifdef __cplusplus
}
#endif

/************************************
*** Private macros and functions ***
************************************/

#ifdef HAVE_KALLOC
#include "kalloc.h"
#else
#include <stdlib.h>
#define kmalloc(km, size) malloc((size))
#define kcalloc(km, count, size) calloc((count), (size))
#define krealloc(km, ptr, size) realloc((ptr), (size))
#define kfree(km, ptr) free((ptr))
#endif

static inline uint32_t *ksw_push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len)
{
  if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
    if (*n_cigar == *m_cigar) {
      *m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
      cigar = (uint32_t*)krealloc(km, cigar, (*m_cigar) << 2);
    }
    cigar[(*n_cigar)++] = len<<4 | op;
  } else cigar[(*n_cigar)-1] += len<<4;
  return cigar;
}

// In the backtrack matrix, value p[] has the following structure:
//   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F, 3 for \tilde{E} and 4 for \tilde{F}
//   bit 3/0x08: 1 if a continuation on the E state (bit 5/0x20 for a continuation on \tilde{E})
//   bit 4/0x10: 1 if a continuation on the F state (bit 6/0x40 for a continuation on \tilde{F})
static inline void ksw_backtrack(void *km, int is_rot, int is_rev, int min_intron_len, const uint8_t *p, const int *off, const int *off_end, int n_col, int i0, int j0,
                                 int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{ // p[] - lower 3 bits: which type gets the max; bit
  int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
  uint32_t *cigar = *cigar_, tmp;
  while (i >= 0 && j >= 0) { // at the beginning of the loop, _state_ tells us which state to check
    int force_state = -1;
    if (is_rot) {
      r = i + j;
      if (i < off[r]) force_state = 2;
      if (off_end && i > off_end[r]) force_state = 1;
      tmp = force_state < 0? p[(size_t)r * n_col + i - off[r]] : 0;
    } else {
      if (j < off[i]) force_state = 2;
      if (off_end && j > off_end[i]) force_state = 1;
      tmp = force_state < 0? p[(size_t)i * n_col + j - off[i]] : 0;
    }
    if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
    else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
    if (state == 0) state = tmp & 7; // TODO: probably this line can be merged into the "else if" line right above; not 100% sure
    if (force_state >= 0) state = force_state;
    if (state == 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
    else if (state == 1 || (state == 3 && min_intron_len <= 0)) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i; // deletion
    else if (state == 3 && min_intron_len > 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 3, 1), --i; // intron
    else cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j; // insertion
  }
  if (i >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, min_intron_len > 0 && i >= min_intron_len? 3 : 2, i + 1); // first deletion
  if (j >= 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
  if (!is_rev)
    for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
      tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
  *m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}

static inline void ksw_reset_extz(ksw_extz_t *ez)
{
  ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
  ez->max = 0, ez->score = ez->mqe = ez->mte = KSW_NEG_INF;
  ez->n_cigar = 0, ez->zdropped = 0, ez->reach_end = 0;
}

static inline int ksw_apply_zdrop(ksw_extz_t *ez, int is_rot, int32_t H, int a, int b, int zdrop, int8_t e)
{
  int r, t;
  if (is_rot) r = a, t = b;
  else r = a + b, t = a;
  if (H > (int32_t)ez->max) {
    ez->max = H, ez->max_t = t, ez->max_q = r - t;
  } else if (t >= ez->max_t && r - t >= ez->max_q) {
    int tl = t - ez->max_t, ql = (r - t) - ez->max_q, l;
    l = tl > ql? tl - ql : ql - tl;
    if (zdrop >= 0 && ez->max - H > zdrop + l * e) {
      ez->zdropped = 1;
      return 1;
    }
  }
  return 0;
}


/*
 * And this is the source file for ksw2_extz
 * This performs global  alignment and extension (e.g., BLAST-like alignment)
 */

typedef struct { int32_t h, e; } eh_t;

void ksw_extz(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int zdrop, int flag, ksw_extz_t *ez)
{
  eh_t *eh;
  int8_t *qp; // query profile
  int32_t i, j, k, max_j = 0, gapoe = gapo + gape, n_col, *off = 0, with_cigar = !(flag&KSW_EZ_SCORE_ONLY);
  uint8_t *z = 0; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be more complex

  ksw_reset_extz(ez);

  // allocate memory
  if (w < 0) w = tlen > qlen? tlen : qlen;
  n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
  qp = (int8_t*)kmalloc(km, qlen * m);
  eh = (eh_t*)kcalloc(km, qlen + 1, 8);
  if (with_cigar) {
    z = (uint8_t*)kmalloc(km, (size_t)n_col * tlen);
    off = (int32_t*)kcalloc(km, tlen, 4);
  }

  // generate the query profile
  for (k = i = 0; k < m; ++k) {
    const int8_t *p = &mat[k * m];
    for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
  }

  // fill the first row
  eh[0].h = 0, eh[0].e = -gapoe - gapoe;
  for (j = 1; j <= qlen && j <= w; ++j)
    eh[j].h = -(gapoe + gape * (j - 1)), eh[j].e = -(gapoe + gapoe + gape * j);
  for (; j <= qlen; ++j) eh[j].h = eh[j].e = KSW_NEG_INF; // everything is -inf outside the band

  // DP loop
  for (i = 0; i < tlen; ++i) { // target sequence is in the outer loop
    int32_t f, h1, st, en, max = KSW_NEG_INF;
    int8_t *q = &qp[target[i] * qlen];
    st = i > w? i - w : 0;
    en = i + w < qlen - 1? i + w : qlen - 1;
    h1 = st > 0? KSW_NEG_INF : -(gapoe + gape * i);
    f  = st > 0? KSW_NEG_INF : -(gapoe + gapoe + gape * i);
    if (!with_cigar) {
      for (j = st; j <= en; ++j) {
        eh_t *p = &eh[j];
        int32_t h = p->h, e = p->e;
        p->h = h1;
        h += q[j];
        h = h >= e? h : e;
        h = h >= f? h : f;
        h1 = h;
        max_j = max > h? max_j : j;
        max   = max > h? max   : h;
        h -= gapoe;
        e -= gape;
        e  = e > h? e : h;
        p->e = e;
        f -= gape;
        f  = f > h? f : h;
      }
    } else if (!(flag&KSW_EZ_RIGHT)) {
      uint8_t *zi = &z[(long)i * n_col];
      off[i] = st;
      for (j = st; j <= en; ++j) {
        eh_t *p = &eh[j];
        int32_t h = p->h, e = p->e;
        uint8_t d; // direction
        p->h = h1;
        h += q[j];
        d = h >= e? 0 : 1;
        h = h >= e? h : e;
        d = h >= f? d : 2;
        h = h >= f? h : f;
        h1 = h;
        max_j = max > h? max_j : j;
        max   = max > h? max   : h;
        h -= gapoe;
        e -= gape;
        d |= e > h? 0x08 : 0;
        e  = e > h? e    : h;
        p->e = e;
        f -= gape;
        d |= f > h? 0x10 : 0; // if we want to halve the memory, use one bit only, instead of two
        f  = f > h? f    : h;
        zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
      }
    } else {
      uint8_t *zi = &z[(long)i * n_col];
      off[i] = st;
      for (j = st; j <= en; ++j) {
        eh_t *p = &eh[j];
        int32_t h = p->h, e = p->e;
        uint8_t d; // direction
        p->h = h1;
        h += q[j];
        d = h > e? 0 : 1;
        h = h > e? h : e;
        d = h > f? d : 2;
        h = h > f? h : f;
        h1 = h;
        max_j = max >= h? max_j : j;
        max   = max >= h? max   : h;
        h -= gapoe;
        e -= gape;
        d |= e >= h? 0x08 : 0;
        e  = e >= h? e    : h;
        p->e = e;
        f -= gape;
        d |= f >= h? 0x10 : 0; // if we want to halve the memory, use one bit only, instead of two
        f  = f >= h? f    : h;
        zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
      }
    }
    eh[j].h = h1, eh[j].e = KSW_NEG_INF;
    // update ez
    if (en == qlen - 1 && eh[qlen].h > ez->mqe)
      ez->mqe = eh[qlen].h, ez->mqe_t = i;
    if (i == tlen - 1)
      ez->mte = max, ez->mte_q = max_j;
    if (ksw_apply_zdrop(ez, 0, max, i, max_j, zdrop, gape)) break;
    if (i == tlen - 1 && en == qlen - 1)
      ez->score = eh[qlen].h;
  }
  kfree(km, qp); kfree(km, eh);
  if (with_cigar) {
    int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
    if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY))
      ksw_backtrack(km, 0, rev_cigar, 0, z, off, 0, n_col, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    else if (ez->max_t >= 0 && ez->max_q >= 0)
      ksw_backtrack(km, 0, rev_cigar, 0, z, off, 0, n_col, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    kfree(km, z); kfree(km, off);
  }
}

/* True global alignment
 * From ksw_gg.c
 *
 */

int ksw_gg(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
  eh_t *eh;
  int8_t *qp; // query profile
  int32_t i, j, k, gapoe = gapo + gape, score, n_col, *off = 0;
  uint8_t *z = 0; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be more complex

  // allocate memory
  if (w < 0) w = tlen > qlen? tlen : qlen;
  n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
  qp = (int8_t*)kmalloc(km, qlen * m);
  eh = (eh_t*)kcalloc(km, qlen + 1, 8);
  if (m_cigar_ && n_cigar_ && cigar_) {
    *n_cigar_ = 0;
    z = (uint8_t*)kmalloc(km, (size_t)n_col * tlen);
    off = (int32_t*)kcalloc(km, tlen, 4);
  }

  // generate the query profile
  for (k = i = 0; k < m; ++k) {
    const int8_t *p = &mat[k * m];
    for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
  }

  // fill the first row
  eh[0].h = 0, eh[0].e = -gapoe - gapoe;
  for (j = 1; j <= qlen && j <= w; ++j)
    eh[j].h = -(gapoe + gape * (j - 1)), eh[j].e = -(gapoe + gapoe + gape * j);
  for (; j <= qlen; ++j) eh[j].h = eh[j].e = KSW_NEG_INF; // everything is -inf outside the band

  // DP loop
  for (i = 0; i < tlen; ++i) { // target sequence is in the outer loop
    int32_t f = KSW_NEG_INF, h1, st, en;
    int8_t *q = &qp[target[i] * qlen];
    st = i > w? i - w : 0;
    en = i + w + 1 < qlen? i + w + 1 : qlen;
    h1 = st > 0? KSW_NEG_INF : -(gapoe + gape * i);
    f  = st > 0? KSW_NEG_INF : -(gapoe + gapoe + gape * i);
    if (m_cigar_ && n_cigar_ && cigar_) {
      uint8_t *zi = &z[(long)i * n_col];
      off[i] = st;
      for (j = st; j < en; ++j) {
        // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
        // Cells are computed in the following order:
        //   H(i,j)   = max{H(i-1,j-1) + S(i,j), E(i,j), F(i,j)}
        //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
        //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
        eh_t *p = &eh[j];
        int32_t h = p->h, e = p->e;
        uint8_t d; // direction
        p->h = h1;
        h += q[j];
        d = h >= e? 0 : 1;
        h = h >= e? h : e;
        d = h >= f? d : 2;
        h = h >= f? h : f;
        h1 = h;
        h -= gapoe;
        e -= gape;
        d |= e > h? 0x08 : 0;
        e  = e > h? e    : h;
        p->e = e;
        f -= gape;
        d |= f > h? 0x10 : 0; // if we want to halve the memory, use one bit only, instead of two
        f  = f > h? f    : h;
        zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
      }
    } else {
      for (j = st; j < en; ++j) {
        eh_t *p = &eh[j];
        int32_t h = p->h, e = p->e;
        p->h = h1;
        h += q[j];
        h = h >= e? h : e;
        h = h >= f? h : f;
        h1 = h;
        h -= gapoe;
        e -= gape;
        e  = e > h? e : h;
        p->e = e;
        f -= gape;
        f  = f > h? f : h;
      }
    }
    eh[en].h = h1, eh[en].e = KSW_NEG_INF;
  }

  // backtrack
  score = eh[qlen].h;
  kfree(km, qp); kfree(km, eh);
  if (m_cigar_ && n_cigar_ && cigar_) {
    ksw_backtrack(km, 0, 0, 0, z, off, 0, n_col, tlen-1, qlen-1, m_cigar_, n_cigar_, cigar_);
    kfree(km, z);
    kfree(km, off);
  }
  return score;
}


void
ksw2_extz_align(std::string Tseq, std::string Qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
    int i;
    int8_t a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    const char *tseq = Tseq.c_str();
    const char *qseq = Qseq.c_str();

    int tl = Tseq.size(), ql = Qseq.size();
    uint8_t *ts, *qs, c[256];

    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
    ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
    for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
      printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
    putchar('\n');
    free(ez.cigar); free(ts); free(qs);

  }


//' Global pairwise alignment from the ksw2 library
//' 
//' This performs global pairwise alignment between a pair of sequences.
//' It returns the number of cigar operations
//' and the operations (M, I, D, optionally =/X if extended CIGARs are used)
//' (ops vector)
//' and the operation positions (opPos vector)
//'
//' Written by Heng Li with small tweaks by August Woerner
//'
//' @param Tseq (a string from the DNA alphabet; this is your target)
//' @param Qseq (like Tseq, but this is your query)
//' @param opPos (integer vector that is modified. It is the position of the cigar operation)
//' @param ops (these are the cigar operations themselves. there encoded as their integer representations, e.g., int('M')
//' @param sc_mch (the score for a match)
//' @param sc_mis (the penalty for a mismatch)
//' @param gapo (gap open penalty)
//' @param gape (gap extend penalty)
//' @param extended (changes M [Match or Mismatch] to =/X [Match or Mismatch, respectively] in the CIGAR)
//' @export
// [[Rcpp::export]]
int
ksw2_gg_align(std::string Tseq,
              std::string Qseq,
              Rcpp::IntegerVector opPos,
              Rcpp::IntegerVector ops,
              int sc_mch=1,
              int sc_mis=-2,
              int gapo=2,
              int gape=1,
              bool extended=true
              )
{
    int i, j, k, outIndex;
    int8_t a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    const char *tseq = Tseq.c_str();
    const char *qseq = Qseq.c_str();

    int tl = Tseq.size(), ql = Qseq.size();
    uint8_t *ts, *qs, c[256];

    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
    ksw_gg(0,
           ql, qs,
           tl, ts,
           5,mat,
           gapo, gape, -1,
           &ez.m_cigar, &ez.n_cigar, &ez.cigar);

    outIndex = 0;
    for (i = 0; i < ez.n_cigar; ++i, ++outIndex) {

      opPos[outIndex] = (int) ez.cigar[i]>>4;
      ops[outIndex] = (int)"MID"[ez.cigar[i]&0xf];

      if (extended) {
        /*
         * This translates the M tag (match or mismatch)
         * into a combination of match tags (=)
         * and mismatch tags ( encoded as 'X')
        */
        if (ops[outIndex] == (int)'M') {
          k = opPos[outIndex]; // the length of the match/mismatch run
          ops[outIndex] = (int)'='; // and use '=' (match) nomenclature

          int previousState = (int) '=';
          if (*tseq != *qseq) { // the first base is a mismatch
            previousState = (int) 'X';
            ops[outIndex] = (int) 'X'; // we know the type of event (mismatch) but not its length
          }

          int thisState;
          int previousRun = 0;
          for (j = 0; j < k; ++j, ++tseq, ++qseq) {
           thisState = (int) '=';
           if (*tseq != *qseq)
             thisState = (int) 'X';

           if (thisState != previousState) {
              opPos[outIndex] = previousRun; // we know how long the previous event was
              ++outIndex;
              ops[outIndex] = thisState; // and we know what type (match XOR mismatch) the current event is
              previousRun = 1; // and we know it's one base long (at least)
            } else {
             ++previousRun; // grow the run.
            }
            previousState = thisState; // updte the previous state
          }
          if (previousRun < k) { // true iff we detected at least one state transition (mismatch)
            opPos[outIndex] = previousRun; // need to update the length of the last event
          }
        } else {
          tseq += opPos[i];
          qseq += opPos[i];
        }
      } else if (ops[i] == (int)'D') { // opPos bases were deleted in the target
         tseq += opPos[i];
      } else { // or inserted in the query.
         qseq += opPos[i];
      }
    }

    free(ez.cigar); free(ts); free(qs);
    return outIndex;
  }


//' String reconstruction from a difference encoding...
//'
//' This takes the output from: ksw2_gg_align_df(Tseq, Qseq, ... )
//' and the reference (Tseq)
//' and uses it to re-create Qseq.
//' The main advantage of this approach is that particular types
//' of mutations can be filtered, as well as variants that fall in
//' particular regions.
//'
//' @param Tseq (the target sequence; e.g., the reference genome)
//' @param positions (1-based positions where the sequence differences are)
//' @param types (the types of events (0,1,2 for mismatch, deletions and insertions, respectively)
//' @param events (the nucleotides involved with the event; ignored if deletion)
//' @param initBuff (guess as to the final size of the query sequence. Overestimating is better than under)
//' @export
// [[Rcpp::export]]
std::string
seqdiffs2seq(std::string Tseq,
             Rcpp::IntegerVector positions,
             Rcpp::IntegerVector types,
             Rcpp::StringVector events,
             int initBuff=-1) {


  
  int nEvents = positions.size();
  if (nEvents==0) // optimize for exact matching...
    return Tseq;
  
  // guess how big the output string will be...
  int i;
  if (initBuff < 1)
    initBuff = Tseq.size() + (int)(Tseq.size()*0.3) + 100;

  // and make a nul string of that size
  // the query is a reconstruction of the query sequence (that generated the alignment ops)
  std::string query(initBuff, '\0');

  char eventCode;
  int eventPos;
  std::string event;

  int previousPos = -1;
  int queryIndex = 0;

  const char* tseq = Tseq.c_str();
  const char* tseqHead = tseq;

  // handle the case of insertions PRIOR to the first base separately!
  i = 0;
  if (positions[0]==0) {
    if (types[0] == 3) {
      event = events[0];
      query.replace(queryIndex, event.size(), event);
      queryIndex += event.size();
    } else {
     Rcerr << "Non-insertion event occurring before the start of the string... cannot be!\n"; 
    }
    ++i;
  }
  
  for(; i < nEvents; ++i) {
    eventPos = positions[i]-1; // convert back to 0-based indexing...
    
    eventCode = 'X'; // mismatch
    if (types[i]==2)
      eventCode= 'D';
    else if (types[i] == 3) {
      eventCode='I';
      ++eventPos; // the index of insertions is tricky.
      // the position refers to the base *after* which the insertion occurs
      // position[0] == 0 -> eventPos==-1; gets changed back to a
    }
    
    event = events[i];

    // fill in the string for all bases prior to the variation
    if (previousPos < eventPos-1) {
      int prevMatches = eventPos - previousPos - 1;
     query.replace(queryIndex, prevMatches, tseq);
      tseq += prevMatches;
      queryIndex += prevMatches;
    }

    // and then fill in the string for all bases of the event
    if (eventCode == 'X') {
      query.replace(queryIndex, 1, event);
      ++tseq;
      ++queryIndex;
    } else if (eventCode == 'D') { // deletion in the query...
      tseq += event.size(); // by construction these are 1 unit long... but just to be safe.
    } else { // bases inserted
      query.replace(queryIndex, event.size(), event);
      queryIndex += event.size();
    }

    previousPos = positions[i]-1;
  }
 // record the bases *after* the last mutational event
  int basesRemaining = Tseq.size() - (tseq-tseqHead);
  if (basesRemaining > 0) {
    query.replace(queryIndex, basesRemaining, tseq);
    queryIndex += basesRemaining;
  }
// and adjust the memory buffer...
  query.resize(queryIndex);
  return query;
}



//' Generates a diff between two sequences.
//' 
//' This is written by Heng Li with small tweaks by August Woerner
//' It performs global alignment (parameterized below)
//' and returns a data frame  of the sequence differences
//'
//' @param Tseq (the target sequence; e.g., the reference genome)
//' @param Qseq (the query sequence)
//' @param sc_mch (the score for a match)
//' @param sc_mis (the penalty for a mismatch)
//' @param gapo (gap open penalty)
//' @param gape (gap extend penalty)
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame
ksw2_seqs2seqdiffs(std::string Tseq,
                std::string Qseq,
                int sc_mch=1,
                int sc_mis=-2,
                int gapo=2,
                int gape=1)
  {
    int i, j, k, outIndex;
    int8_t a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    const char *tseq = Tseq.c_str();
    const char *qseq = Qseq.c_str();

    const char *tseqHead = tseq;

    int tl = Tseq.size(), ql = Qseq.size();

    std::vector<int> mismatchTypes; // I (insertion) D (deletion) X (Mismatch).
    std::vector<int> mismatchPositions; // the where, in the reference sequence, for the event
    std::vector<std::string> mismatchEvents; // the what. the bases inserted, deleted,or mismatched

    uint8_t *ts, *qs, c[256];

    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
    ksw_gg(0,
           ql, qs,
           tl, ts,
           5,mat,
           gapo, gape, -1,
           &ez.m_cigar, &ez.n_cigar, &ez.cigar);

    outIndex = 0;

    int op, opPos;



    for (i = 0; i < ez.n_cigar; ++i, ++outIndex) {

      opPos = (int) ez.cigar[i]>>4;
      op= (int)"MID"[ez.cigar[i]&0xf];

        /*
        * This evaluates the M tag (match or mismatch)
        * encoding the mismatches as 'X' events
        */
        if (op == (int)'M') {
          k = opPos; // the length of the match/mismatch run
          for (j = 0; j < k; ++j, ++tseq, ++qseq) {
            if (*tseq != *qseq) {
              mismatchTypes.push_back( 1 );
              mismatchPositions.push_back( static_cast<int>(tseq - tseqHead) + 1); // 1-based indexing...
              mismatchEvents.push_back( std::string(qseq, 1) );
            }
          }
      } else if (op == (int)'D') { // opPos bases were deleted in the query
// each base of the deletion is encoded...
// this makes for easy subsetting
        for (j=0; j < opPos; ++j, ++tseq) {
          mismatchTypes.push_back( 2 );
          mismatchPositions.push_back( static_cast<int>(tseq - tseqHead) + 1);
          mismatchEvents.push_back( std::string(tseq, 1) );
        }
      } else { // or inserted in the query.
  // insertions get encoded as a single event.
        mismatchTypes.push_back( 3 );
        mismatchPositions.push_back( static_cast<int>(tseq - tseqHead) + 1);
        mismatchEvents.push_back( std::string(qseq, opPos) );
        qseq += opPos;

      }
    }

    free(ez.cigar); free(ts); free(qs);
    int nOps = mismatchTypes.size();

    /*
     * Convert from c++ types to R-friendly types
     * First, the difference types (encoded as a factor)
     */
    IntegerVector typesR(nOps);
    typesR = mismatchTypes;
    typesR.attr("class") = "factor";
    CharacterVector diffTypes = CharacterVector::create("X", "D", "I");
    typesR.attr("levels") = diffTypes;

    // and the difference locations
    IntegerVector posR(nOps);
    posR = mismatchPositions;

    // and the bases involved in the
    StringVector eventsR(nOps);
    eventsR = mismatchEvents;

// and make a nice data frame...
    return DataFrame::create(_["Positions"]=posR,
                             _["Types"]=typesR,
                             _["Events"]=eventsR,
                             _["stringsAsFactors"] = false);
}



