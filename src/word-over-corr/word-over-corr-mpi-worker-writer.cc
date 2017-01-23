//% generate all word overlap correlations
// This file is part of the Word-overlap collection.
// Copyright (C) 2013-2017 Paul Leopardi
// Parts of this code are based on code by Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the src directory.

#include "mpi.h"
#include <stdlib.h>

#include "fxt/bits-all.h"
//#include "fxt/comb/setpart-p-rgs-lex.h"
#include "setpart-p-rgs-lex-exp.h"

#include "fxt/fxtio.h"

#include "fxt/fxttypes.h"  // ulong

#include "fxt/jjassert.h"

#include <string>
#include <iosfwd>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <cmath>

#include <sys/gmon.h>

// Uncomment to disable detailed printing
// #define WOC_TOTALS_ONLY

typedef std::string string;

const int END_TAG  = 0;
const int MORE_TAG = 1;

#ifndef WORD_LEN
#define WORD_LEN 2
#endif

// Number of characters in a word
ulong T = WORD_LEN;

// Number of different characters in the alphabet
ulong alpha;

// Number of worker processes
ulong nbr_workers;

#ifndef WOC_DELTA_OFFSET
#define WOC_DELTA_OFFSET 1
#endif

// Maximum number of bits needed to represent a character
ulong log2_alpha;

#ifdef WOC_DEBUG
// Maximum number of bits needed to represent a word
ulong nbw;
#endif

// Number of bits in mask
ulong nbm;

// Initial mask
ulong m0;

// Maximum number of different correlation values
ulong N;

struct hist_entry
{
  ulong m_beta;
  ulong m_a;
  ulong m_b;
  ulong m_count;
};

typedef hist_entry*     pu;  // Pointer to hist_entry
typedef hist_entry**    ppu; // Pointer to pointer to hist_entry
typedef hist_entry***   p3u; // Pointer to pointer to pointer to hist_entry
typedef hist_entry****  p4u; // Pointer to pointer to pointer to pointer to hist_entry

// Histogram
p4u hist;

// MPI variables
int numprocs;
int myid;
char processor_name[MPI_MAX_PROCESSOR_NAME];
int namelen;

struct word_corr_pair
{
  ulong m_beta;
  ulong m_a;
  ulong m_b;
  ulong m_aa;
  ulong m_bb;
  ulong m_ab;
  ulong m_ba;
  ulong m_count;
};

MPI_Datatype pair_mpi_struct_type;

void create_pair_mpi_struct_type()
{
  // MPI variables
  int pair_mpi_blocklengths[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
  MPI_Aint  pair_mpi_displacements[8];
  MPI_Datatype pair_mpi_types[8] = { MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,
                                     MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,
                                     MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,
                                     MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG };
  word_corr_pair dummy_pair;
  word_corr_pair* buff = &dummy_pair;
  MPI_Aint buff_address;
  MPI_Get_address(buff, &buff_address);
  MPI_Get_address(&(buff->m_beta),  &(pair_mpi_displacements[0]));
  MPI_Get_address(&(buff->m_a),     &(pair_mpi_displacements[1]));
  MPI_Get_address(&(buff->m_b),     &(pair_mpi_displacements[2]));
  MPI_Get_address(&(buff->m_aa),    &(pair_mpi_displacements[3]));
  MPI_Get_address(&(buff->m_bb),    &(pair_mpi_displacements[4]));
  MPI_Get_address(&(buff->m_ab),    &(pair_mpi_displacements[5]));
  MPI_Get_address(&(buff->m_ba),    &(pair_mpi_displacements[6]));
  MPI_Get_address(&(buff->m_count), &(pair_mpi_displacements[7]));

  for (int k=0; k!=8; k++)
    pair_mpi_displacements[k] -= buff_address;

  MPI_Type_create_struct(8, pair_mpi_blocklengths, pair_mpi_displacements, pair_mpi_types, &pair_mpi_struct_type);
  MPI_Type_commit(&pair_mpi_struct_type);
}

#ifndef WOC_PAIR_BUFFER_SIZE
#define WOC_PAIR_BUFFER_SIZE 10000
#endif

const int pair_buffer_size = WOC_PAIR_BUFFER_SIZE;
word_corr_pair* pair_buffer;

//: Default is to print zeros as dots.
static const char n01[] = {'.', '1'};

void fprint_bin(const char *bla, unsigned long long x, int pd/*=0*/, const char *c01/*=0*/)
// Print x to radix-2.
// pd: how many bits to print.
{
    cout << bla;

    if ( 0==pd )  pd = BITS_PER_LONG;
    unsigned long long m = 1ULL << (BITS_PER_LONG_LONG-1);
    m >>= (BITS_PER_LONG_LONG-pd);

    const char *d = ( 0==c01 ?  n01 : c01 );

    for(  ; 0!=m; m>>=1)  cout << d[ 0!=(x & m) ];
}

void fprint_char(const char *bla, unsigned long long x, const ulong T, const ulong log2_alpha, const char c='A')
// Print x to radix-2.
// pd: how many bits to print.
{
    cout << bla;

    unsigned long long m0 = (1ULL<<log2_alpha)-1;
    //cout << endl << "m0==" << m0 << endl;
    unsigned long long m = m0 << log2_alpha*(T-1);
    //cout << endl << "m==" << m << endl;
    for(int k=T-1; 0!=m; --k, m>>=log2_alpha)
    {
      //cout << "k==" << k << endl;
      //cout << "(x & m)==" << (x & m) << endl;
      //cout << "(x & m)>>(k*log2_alpha)==" << ((x & m)>>(k*log2_alpha)) << endl;
      cout << (char)(c+((x & m)>>(k*log2_alpha))); // << endl;
    }
}

string poly(const ulong T, ulong x)
{
  // Translate a ulong representing a correlation into a string repesenting a correlation polynomial
  std::ostringstream op;
  op << "z->";
  if (x == 0)
    op << "0";
  else
  {
    int i;
    for (i=(T-1); i >= 0; --i, x >>= 1)
      if (x % 2)
      {
        if (i == 0)
          op << "1";
        else if (i == 1)
          op << "z";
        else
          op << "z^" << i;
        break;
      }
    --i;
    x >>= 1;
    for (; i >= 0; --i, x >>= 1)
      if (x % 2)
      {
        op << "+";
        if (i == 0)
          op << "1";
        else if (i == 1)
          op << "z";
        else
          op << "z^" << i;
      }
  }
  return op.str();
}

void send_hist()
{
  int buff_index = 0;
  for (ulong ab=0; ab<N; ++ab)
  {
    const p3u hab = hist[ab];
    if (hab != NULL)
      for (ulong ba=0; ba<N; ++ba)
      {
        const ppu habba = hab[ba];
        if (habba != NULL)
          for (ulong aa=0; aa<N; ++aa)
          {
            const pu habbaaa = habba[aa];
            if (habbaaa != NULL)
              for (ulong bb=0; bb<N; ++bb)
              {
                const hist_entry& habbaaabb = habbaaa[bb];
                if (habbaaabb.m_count)
                {
                  word_corr_pair* buff = &(pair_buffer[buff_index]);
                  buff->m_beta = habbaaabb.m_beta;
                  buff->m_a    = habbaaabb.m_a;
                  buff->m_b    = habbaaabb.m_b;
                  buff->m_aa   = aa;  buff->m_bb = bb;
                  buff->m_ab   = ab;  buff->m_ba = ba;
                  buff->m_count = habbaaabb.m_count;
                  buff_index++;
                  if (buff_index == pair_buffer_size)
                  {
                    MPI_Send((void*) &(pair_buffer[0]), buff_index, pair_mpi_struct_type, 0, MORE_TAG, MPI_COMM_WORLD);
#ifdef WOC_MPI_DEBUG
                    cout << "Process " << myid << " sent  : " << buff_index << endl;
#endif
                    buff_index = 0;
                  }
                }
              }
          }
      }
  }
  MPI_Send((void*) &(pair_buffer[0]), buff_index, pair_mpi_struct_type, 0, END_TAG, MPI_COMM_WORLD);
#ifdef WOC_MPI_DEBUG
  cout << "Process " << myid << " sent  : " << buff_index << " after loop " << endl;
#endif
}

void output_hist(const ulong beta)
{
  ulong pairs = 0;    // count hist
  ulong classes = 0;  // count nonzero values
  for (ulong ab=0; ab<N; ++ab)
  {
    const p3u hab = hist[ab];
    if (hab != NULL)
      for (ulong ba=0; ba<N; ++ba)
      {
        const ppu habba = hab[ba];
        if (habba != NULL)
          for (ulong aa=0; aa<N; ++aa)
          {
            const pu habbaaa = habba[aa];
            if (habbaaa != NULL)
              for (ulong bb=0; bb<N; ++bb)
              {
                const hist_entry& habbaaabb = habbaaa[bb];
                if (habbaaabb.m_count)
                {
                  pairs += habbaaabb.m_count;
                  classes++;
#ifdef WOC_PRINT_NEW_CLASSES
                  if (beta > 3 && beta == habbaaabb.m_beta)
                  {
                    fprint_char(" X==", habbaaabb.m_a, T, log2_alpha);
                    cout << "; ";
                    fprint_char(" Y==", habbaaabb.m_b, T, log2_alpha);
                    cout << "; ";
                    cout << endl;
                    fprint_bin("XX==", N+aa, T, "01");
                    cout << "; ";
                    fprint_bin("YY==", N+bb, T, "01");
                    cout << "; ";
                    cout << endl;
                    fprint_bin("XY==", ab, T, "01");
                    cout << "; ";
                    fprint_bin("YX==", ba, T, "01");
                    cout << "; ";
                    cout << endl;
                    cout << endl;
                  }
#endif
#ifndef WOC_TOTALS_ONLY
                  if (beta == alpha)
                  {
                    //print_bin("", ab, T);
                    cout <<   "pXY[" << myid << "," << classes << "]:=" << poly(T, ab) << ":";
                    //print_bin("", ba, T);
                    cout <<   "pYX[" << myid << "," << classes << "]:=" << poly(T, ba) << ":";
                    //print_bin("", N+aa, T);
                    cout <<   "pXX[" << myid << "," << classes << "]:=" << poly(T, N+aa) << ":";
                    //print_bin("", N+bb, T);
                    cout <<   "pYY[" << myid << "," << classes << "]:=" << poly(T, N+bb) << ":";
                    cout << "count[" << myid << "," << classes << "]:=" << habbaaabb.m_count << ":" << endl;
                  }
#endif
                }
              }
          }
      }
  }
  cout << "classes    ==" << setw(20) << classes
       << " (number of correlation classes so far)" << endl;
  cout << "pairs      ==" << setw(20) << pairs
       << " (number of different pairs of unequal words so far)" << endl;
  cout << "-------------" << endl;
}

void free_hist()
{
#ifdef WOC_MALLOC_DEBUG
  cout << "Process " << myid << " free_hist() called." << endl;
#endif
  for (ulong ab=0; ab<N; ++ab)
  {
    const p3u hab = hist[ab];
    if (hab != NULL)
    {
      for (ulong ba=0; ba<N; ++ba)
      {
        const ppu habba = hab[ba];
        if (habba != NULL)
        {
          for (ulong aa=0; aa<N; ++aa)
          {
            const pu habbaaa = habba[aa];
            if (habbaaa != NULL)
#ifdef WOC_MALLOC_DEBUG
              cout << "Process " << myid << " calling free(habba[aa]); aa==" << aa << endl;
#endif
              free(habba[aa]);
          }
#ifdef WOC_MALLOC_DEBUG
          cout << "Process " << myid << " calling free(hab[ba]); ba==" << ba << endl;
#endif
          free(hab[ba]);
        }
      }
#ifdef WOC_MALLOC_DEBUG
      cout << "Process " << myid << " calling free(hist[ab]); ab==" << ab << endl;
#endif
      free(hist[ab]);
      hist[ab] = NULL;
    }
  }
}

void swap(ulong& a, ulong& b)
{
  const ulong tmp = a;
  a = b;
  b = tmp;
}

inline ulong autocorr(ulong m, const ulong x)
{
  ulong xx = 0;  // Correlation bit vector as ulong
  ulong xt = x;  // Copy of word, to be shifted
  do
  {
    xt >>= log2_alpha;
    xx <<= 1;
    xx |= ( (m &( x^xt )) == 0 );
    m >>= log2_alpha;
  }
  while ( m );
  return xx;
}

void calc_ab_ba(ulong *ab, ulong *ba, ulong m, ulong a, ulong b)
{
  ulong bt = b;  // aux, to be shifted
  ulong at = a;  // aux, to be shifted
  do
  {
    bt >>= log2_alpha;
    *ab <<= 1;
    *ab |= ( (m &( a^bt )) == 0 );
    at >>= log2_alpha;
    *ba <<= 1;
    *ba |= ( (m &( b^at )) == 0 );
    m >>= log2_alpha;
  }
  while ( m );
}

void insert_hist(const ulong beta,
                 const ulong a, const ulong b, ulong aa, ulong bb, ulong ab, ulong ba, const ulong increment)
{
  bool swapped_ab_ba = false;
  if (ab < ba)
  {
    swap(ab, ba);
    swapped_ab_ba = true;
  }
  p3u hab = hist[ab];
  if (hab == NULL)
  {
#ifdef WOC_MALLOC_DEBUG
    cout << "Process " << myid << " calling malloc(N*sizeof(ppu)); ab==" << ab << endl;
#endif
    hist[ab] = (p3u) malloc(N*sizeof(ppu));
    hab = hist[ab];
    for (ulong k=0; k<N; ++k)
      hab[k] = NULL;
  }
  ppu habba = hab[ba];
  if (habba == NULL)
  {
#ifdef WOC_MALLOC_DEBUG
    cout << "Process " << myid << " calling malloc(N*sizeof(pu)); ba==" << ba << endl;
#endif
    hab[ba] = (ppu) malloc(N*sizeof(pu));
    habba = hab[ba];
    for (ulong k=0; k<N; ++k)
      habba[k] = NULL;
  }
  bool swapped_aa_bb = false;
  if (aa < bb)
  {
    swap(aa, bb);
    swapped_aa_bb = true;
  }
  pu habbaaa = habba[aa];
  if (habbaaa == NULL)
  {
#ifdef WOC_MALLOC_DEBUG
    cout << "Process " << myid << " calling malloc(N*sizeof(hist_entry)); aa==" << aa << endl;
#endif
    habba[aa] = (pu) malloc(N*sizeof(hist_entry));
    habbaaa = habba[aa];
    for (ulong k=0; k<N; ++k)
      habbaaa[k].m_count = 0;
  }
  if (habbaaa[bb].m_count == 0 || habbaaa[bb].m_beta > beta)
  {
    habbaaa[bb].m_beta = beta;
    habbaaa[bb].m_a = a;
    habbaaa[bb].m_b = b;
  }
  habbaaa[bb].m_count += increment;
#ifdef WOC_PER_PAIR_DEBUG
  cout << "Process " << myid
                     << " ab==" << ab
                     << " ba==" << ba
                     << " aa==" << aa
                     << " bb==" << bb
                     << " increment==" << increment
                     << "; habbaaa[" << bb << "].m_count==" << habbaaa[bb].m_count << endl;
#endif
}

bool arrays_are_equal(ulong n, const ulong* s, const ulong *t)
{
  for (ulong k=0; k!=n; k++)
    if (s[k] != t[k])
      return false;
  return true;
}

ulong* new_setpart_p_rgs_array(const ulong beta, const ulong gamma)
{
    ulong* s = new ulong[2*T];
    if (s == NULL)
      return s;
    for (ulong k = 0; k != 2*T; k++)
      s[k] = 0UL;
    if (gamma > beta)
      s[0] = ~0UL;
    else
    {
      for (ulong k=0; k != gamma; k++)
        s[k] = k;
      for (ulong k=0; gamma + k < beta; k++)
        s[2*T-beta+gamma+k] = gamma+k;
      if (!is_setpart_p_rgs_array(2*T, beta, s))
      {
        cout << "Failure in new_setpart_p_rgs_array for ";
        for (ulong k=0; k < 2*T; k++)
          cout << s[k] << " ";
        cout << endl;
        return (ulong *) NULL;
      }
    }
    return s;
}

ulong check_setpart_p_rgs_range(const ulong beta, const ulong orbit_size,
                                const ulong* rgs_begin,const ulong* rgs_end)
{
#ifdef WOC_MPI_DEBUG
  cout << "Process " << myid << " check_setpart_p_rgs_range" << endl;
  cout << "Process " << myid << " rgs_begin ==";
  for (ulong k=0; k < 2*T; k++)
    cout << rgs_begin[k] << " ";
  cout << endl;
  cout << "Process " << myid << " rgs_end   ==";
  for (ulong k=0; k < 2*T; k++)
    cout << rgs_end[k] << " ";
  cout << endl;
#endif
  if (!is_setpart_p_rgs_array(2*T, beta, rgs_begin))
  {
    cout << "Process " << myid << " failure in rgs_begin" << endl;
    cout << "Process " << myid << " rgs_begin ==";
    for (ulong k=0; k < 2*T; k++)
      cout << rgs_begin[k] << " ";
    cout << endl;
    exit(-1);
  }
  setpart_p_rgs_lex p = setpart_p_rgs_lex(2*T,beta,rgs_begin);
  ulong olda = 0;
  ulong aa = autocorr(m0, olda);
  ulong steps = 0;
  bool more;
  for (more=true; more; more=p.next())
  {
    const ulong* rgs_data = p.data();

    if (arrays_are_equal(2*T, rgs_data, rgs_end))
      break;

    steps++;

    ulong a = rgs_data[0];
    ulong b = rgs_data[T];
    for (ulong k=1; k<T; ++k)
    {
      a <<= log2_alpha;
      b <<= log2_alpha;
      a |= rgs_data[k];
      b |= rgs_data[T+k];
    }
    if (a != olda)
    {
      aa = autocorr(m0, a);
      olda = a;
    }
    if (a != b)
    {
      ulong bb = autocorr(m0, b);
      ulong ab = 0;  // corr
      ulong ba = 0;  // corr
      calc_ab_ba(&ab, &ba, m0, a, b);
      insert_hist(beta, a, b, aa, bb, ab, ba, orbit_size);
    }
  }
  return steps;
}

int
main(int argc, char **argv)
{
  extern void* _start;
  extern void* etext;
  monstartup ((u_long) &_start, (u_long) &etext);

  MPI_Init(&argc,&argv);
  create_pair_mpi_struct_type();
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);

  cerr << "Process " << myid << " of " << numprocs << " on " << processor_name << endl;

  // Number of characters in a word
  if (argc > 1)
    sscanf(argv[1],"%lu",&T);

  // Number of different characters in the alphabet
  alpha = T*2;
  if (argc > 2)
    sscanf(argv[2],"%lu",&alpha);

  // Number of worker  processes
  nbr_workers = numprocs - 1;
  if (nbr_workers < 0)
    return -1;

  // Maximum number of bits needed to represent a character
  log2_alpha = (ulong) std::ceil(std::log(alpha)/std::log(2.0));
#ifdef WOC_DEBUG
  // Maximum number of bits needed to represent a word
  nbw = T*log2_alpha;
#endif
  // Number of bits in mask
  nbm = (T-1)*log2_alpha;
  // Initial mask
  m0 = (1UL << nbm) - 1UL;
  // Maximum number of different correlation values
  N = 1UL<<(T-1);

  //jjassert(nbw <= BITS_PER_LONG);
  //jjassert(T > 1);
  //jjassert(alpha <= 2*T);
  if (myid == 0)
  {
    cout << "T          ==" << setw(20) << T
         << " (number of characters in a word)" << endl;
    cout << "alpha      ==" << setw(20) << alpha
         << " (number different characters in the alphabet)" << endl;
    cout << "nbr_workers==" << setw(20) << nbr_workers
         << " (number of worker threads to use)" << endl;
  #ifdef WOC_DEBUG
    cout << "log2_alpha ==" << setw(20) << log2_alpha
         << " (number of bits needed to represent a character)" << endl;
    cout << "nbw        ==" << setw(20) << nbw
         << " (number of bits needed to represent a word)" << endl;
    cout << "nbm        ==" << setw(20) << nbm
         << " (number of bits in mask)" << endl;
    cout << "N          ==" << setw(20) << N
         << " (maximum number of different possible correlation values)" << endl;
  #endif
    cout << "-------------" << endl;
  }
  typedef word_corr_pair* pword_pair;
#ifdef WOC_MALLOC_DEBUG
    cout << "Process " << myid << " calling malloc(N*sizeof(p3u));" << endl;
#endif
  hist = (p4u) malloc(N*sizeof(p3u));
  for (ulong k=0; k < N; ++k)
    hist[k] = NULL;
  pair_buffer = new word_corr_pair[pair_buffer_size];
  for (ulong beta=2; beta<=alpha; ++beta)
  {
    ulong orbit_size=1;
    for (ulong k=alpha; k>alpha-beta; --k)
      orbit_size *=k;
    if (myid == 0)
    {
      cout << "beta       ==" << setw(20) << beta
           << " (number of different characters in the word pair)" << endl;
      cout << "orbit_size ==" << setw(20) << orbit_size
           << " (size of an orbit under permutation of the alphabet)" << endl;
    }
    ulong partitions = 0;
    if (myid == 0 && nbr_workers == 0)
    {
      ulong* rgs_begin = new_setpart_p_rgs_array(beta, 1);
      if (rgs_begin == (ulong *) NULL)
        return -1;

      ulong* rgs_end = new_setpart_p_rgs_array(beta, beta+1);
      if (rgs_end == (ulong *) NULL)
        return -1;

      partitions = check_setpart_p_rgs_range(beta, orbit_size, rgs_begin, rgs_end);
      delete [] rgs_begin;
      delete [] rgs_end;
    }
    else
    {
      ulong steps = 0;
      if (myid > 0)
      {
        ulong delta = (beta+WOC_DELTA_OFFSET > 2*T) ? beta : beta+WOC_DELTA_OFFSET;

        setpart_p_rgs_lex r = setpart_p_rgs_lex(2*T,beta);
        const ulong* r_data = r.data();
        ulong* rgs_begin = new ulong[2*T];
        if (rgs_begin == (ulong *) NULL)
          return -1;
        for (int i = 0; i != 2*T; i++)
          rgs_begin[i] = r_data[i];

        ulong* rgs_end = new ulong[2*T];
        if (rgs_end == (ulong *) NULL)
          return -1;
        for (int i = delta; i < 2*T; i++)
          rgs_end[i] = 0;

        setpart_p_rgs_lex q = setpart_p_rgs_lex(delta,beta);
        ulong delta_index = 0;
        for (bool more = true; more; more = q.next(), delta_index++)
        {
          const ulong* q_data = q.data();
          for (int i = 0; i != delta; i++)
            rgs_end[i] = q_data[i];
          if (delta_index % nbr_workers + 1 == (ulong) myid)
          {
#ifdef WOC_MPI_DEBUG
            cout << "Process " << myid
                 << " delta_index % nbr_workers + 1 ==" << delta_index % nbr_workers + 1
                 << " delta_index ==" << delta_index
                 << endl;
#endif
            steps += check_setpart_p_rgs_range(beta, orbit_size, rgs_begin, rgs_end);
          }
          for (int i = 0; i != 2*T; i++)
            rgs_begin[i] = rgs_end[i];
        }
        rgs_end[0] = ~0UL;
        if (delta_index % nbr_workers + 1 == (ulong) myid)
        {
#ifdef WOC_MPI_DEBUG
          cout << "Process " << myid
               << " delta_index % nbr_workers + 1 ==" << delta_index % nbr_workers + 1
               << " delta_index ==" << delta_index
               << endl;
#endif
          steps += check_setpart_p_rgs_range(beta, orbit_size, rgs_begin, rgs_end);
        }
#ifdef WOC_MPI_DEBUG
        cout << "Process " << myid << " calling send_hist() " << endl;
#endif
        send_hist();
        delete [] rgs_begin;
        delete [] rgs_end;
      }
      else
      {
        MPI_Status status;
        for (ulong nbr_end_tags = 0; nbr_end_tags != nbr_workers;)
        {
          int nbr_word_pairs;
          MPI_Recv((void*) &(pair_buffer[0]), pair_buffer_size, pair_mpi_struct_type, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          MPI_Get_count(&status, pair_mpi_struct_type, &nbr_word_pairs);
#ifdef WOC_MPI_DEBUG
          cout << "Process " << myid << " received " << nbr_word_pairs << endl;
#endif
          for (int j = 0; j != nbr_word_pairs; j++)
          {
            word_corr_pair* buff = &(pair_buffer[j]);
#ifdef WOC_MPI_PER_PAIR_DEBUG
            cout << "Process " << myid << " calling insert_hist(); j==" << j << endl;
#endif
            insert_hist(buff->m_beta,
                        buff->m_a,  buff->m_b,
                        buff->m_aa, buff->m_bb,
                        buff->m_ab, buff->m_ba,
                        buff->m_count);
          }
          if (status.MPI_TAG == END_TAG)
          {
#ifdef WOC_MPI_DEBUG
            cout << "Process " << myid << " received END_TAG" << endl;
#endif
            nbr_end_tags++;
          }
        }
      }
#ifdef WOC_MPI_DEBUG
      cout << "Process " << myid << " approached Barrier 1" << endl;
#endif
      MPI_Barrier(MPI_COMM_WORLD);
#ifdef WOC_MPI_DEBUG
      cout << "Process " << myid << " passed Barrier 1" << endl;
#endif
      MPI_Reduce(&steps, &partitions, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef WOC_MPI_VERBOSE
      if (myid > 0)
      {
        cout << "P" << setw(2) << myid;
        cout << " steps  ==" << setw(20) << steps
             << " (number of restricted growth strings generated for this value of beta, myid)"
             << endl;
      }
#endif
    }
    if (myid == 0)
    {
      cout << "partitions ==" << setw(20) << partitions
           << " (number of restricted growth strings generated for this value of beta)"
           << endl;
    }
#ifdef WOC_MPI_DEBUG
    cout << "Process " << myid << " approached Barrier 2" << endl;
#endif
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef WOC_MPI_DEBUG
    cout << "Process " << myid << " passed Barrier 2" << endl;
#endif
    if (myid == 0)
    {
      output_hist(beta);
      if (nbr_workers > 0)
        free_hist();
    }
#ifdef WOC_MPI_DEBUG
    cout << "Process " << myid << " approached Barrier 3" << endl;
#endif
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef WOC_MPI_DEBUG
    cout << "Process " << myid << " passed Barrier 3" << endl;
#endif
  }
  MPI_Finalize();
  return 0;
}
