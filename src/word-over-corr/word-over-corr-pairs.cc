//% generate all word overlap correlations
// This file is part of the Word-overlap collection.
// Copyright (C) 2010-2017 Paul Leopardi
// Parts of this code are based on code by Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the src directory.

#include <stdlib.h>

#include "fxt/bits-all.h"
#include "fxt/comb/setpart-p-rgs-lex.h"

#include "fxt/fxtiomanip.h"

#include "fxt/fxttypes.h"  // ulong

#include "fxt/jjassert.h"

#include <string>
#include <iosfwd>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <cmath>
// Uncomment to disable detailed printing
// #define TOTALS_ONLY

typedef struct
{
  ulong count;
  ulong beta;
  ulong X;
  ulong Y;
} count_pair;

typedef count_pair*     pu;  // Pointer to count_pair
typedef count_pair**    ppu; // Pointer to pointer to count_pair
typedef count_pair***   p3u; // Pointer to pointer to pointer to count_pair
typedef count_pair****  p4u; // Pointer to pointer to pointer to pointer to count_pair

typedef std::string string;

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

void output_hist(const ulong T, const ulong N, const p4u hist)
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
                const count_pair habbaaabb = habbaaa[bb];
                if ( habbaaabb.count )
                {
                    pairs += habbaaabb.count;
                    classes++;
#ifndef TOTALS_ONLY
                    //print_bin("", ab, T);
                    cout << "pXY[" << classes << "]:=" << poly(T, ab) << ":";
                    //print_bin("", ba, T);
                    cout << "pYX[" << classes << "]:=" << poly(T, ba) << ":";
                    //print_bin("", N+aa, T);
                    cout << "pXX[" << classes << "]:=" << poly(T, N+aa) << ":";
                    //print_bin("", N+bb, T);
                    cout << "pYY[" << classes << "]:=" << poly(T, N+bb) << ":";
                    cout << "count[" << classes << "]:=" << habbaaabb.count << ":" << endl;
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

void swap(ulong& a, ulong& b)
{
  const ulong tmp = a;
  a = b;
  b = tmp;
}

inline ulong autocorr(const ulong log2_alpha, const ulong m0, const ulong x)
{
  ulong xx = 0;  // Correlation bit vector as ulong
  ulong xt = x;  // Copy of word, to be shifted
  ulong m = m0;  // Copy of mask, to be shifted
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

#ifndef WORD_LEN
#define WORD_LEN 2
#endif

int
main(int argc, char **argv)
{
  // Number of characters in a word
  ulong T = WORD_LEN;
  if (argc > 1)
    sscanf(argv[1],"%uld",&T);

  // Number of different characters in the alphabet
  ulong alpha = T*2;
  if (argc > 2)
    sscanf(argv[2],"%uld",&alpha);

  // Maximum number of bits needed to represent a character
  const ulong log2_alpha = (ulong) std::ceil(std::log(alpha)/std::log(2.0));
  // Maximum number of bits needed to represent a word
  const ulong nbw = T*log2_alpha;
  // Number of bits in mask
  const ulong nbm = (T-1)*log2_alpha;
  // Initial mask
  const ulong m0 = (1UL << nbm) - 1UL;
  // Maximum number of different correlation values
  const ulong N = 1UL<<(T-1);

  cout << "T          ==" << setw(20) << T
            << " (number of characters in a word)" << endl;
  cout << "alpha      ==" << setw(20) << alpha
            << " (number different characters in the alphabet)" << endl;
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

  if (!(nbw <= BITS_PER_LONG))
  {
    cout << "Must have (nbw <= BITS_PER_LONG) but nbw==" << nbw
         << "; BITS_PER_LONG==" << BITS_PER_LONG << endl;
    exit(1);
  }
  if (!(T > 1))
  {
    cout << "Must have (T > 1) but T==" << T << endl;
    exit(1);
  }
  if (!(alpha <= 2*T))
  {
    cout << "Must have (alpha <= 2*T) but alpha==" << alpha
         << "; 2*T==" << (2*T) << endl;
    exit(1);
  }

  const p4u hist = (p4u) malloc(N*sizeof(p3u));
  for (ulong k=0; k<N; ++k)
    hist[k] = NULL;

  for (ulong beta=2; beta<=alpha; ++beta)
  {
    ulong orbit_size=1;
    for (ulong k=alpha; k>alpha-beta; --k)
      orbit_size *=k;

    cout << "beta       ==" << setw(20) << beta
              << " (number of different characters in the word pair)" << endl;
    cout << "orbit_size ==" << setw(20) << orbit_size
              << " (size of an orbit under permutation of the alphabet)" << endl;

    setpart_p_rgs_lex p = setpart_p_rgs_lex(2*T,beta);
    ulong olda = 0;
    ulong oldaa0 = autocorr(log2_alpha, m0, olda);
    long partitions = 0;
    for (bool more=true; more; more=p.next())
    {
      const ulong *rgs_data = p.data();
      partitions++;

      ulong a = rgs_data[0];
      ulong b = rgs_data[T];
      for (ulong k=1; k<T; ++k)
      {
        a <<= log2_alpha;
        b <<= log2_alpha;
        a |= rgs_data[k];
        b |= rgs_data[T+k];
      }
      if (a!=b)
      {
        ulong aa0;
        if (a == olda)
        {
          aa0 = oldaa0;
        }
        else
        {
          aa0 = autocorr(log2_alpha, m0, a);
          olda = a;
          oldaa0 = aa0;
        }
        ulong ab = 0;  // corr
        ulong ba = 0;  // corr
        ulong bt = b;  // aux, to be shifted
        ulong at = a;  // aux, to be shifted
        ulong m = m0;  // aux, to be shifted
        do
        {
          bt >>= log2_alpha;
          ab <<= 1;
          ab |= ( (m &( a^bt )) == 0 );
          at >>= log2_alpha;
          ba <<= 1;
          ba |= ( (m &( b^at )) == 0 );
          m >>= log2_alpha;
        }
        while ( m );

        bool swapped_ab_ba = false;
        if (ab < ba)
        {
          swap(ab, ba);
          swapped_ab_ba = true;
        }
        p3u hab = hist[ab];
        if (hab == NULL)
        {
          hist[ab] = (p3u) malloc(N*sizeof(ppu));
          hab = hist[ab];
          for (ulong k=0; k<N; ++k)
            hab[k] = NULL;
        }
        ppu habba = hab[ba];
        if (habba == NULL)
        {
          hab[ba] = (ppu) malloc(N*sizeof(pu));
          habba = hab[ba];
          for (ulong k=0; k<N; ++k)
            habba[k] = NULL;
        }
        ulong bb = autocorr(log2_alpha, m0, b);
        ulong aa = aa0;
        bool swapped_aa_bb = false;
        if (aa < bb)
        {
          swap(aa, bb);
          swapped_aa_bb = true;
        }
        pu habbaaa = habba[aa];
        if (habbaaa == NULL)
        {
          habba[aa] = (pu) malloc(N*sizeof(count_pair));
          habbaaa = habba[aa];
          for (ulong k=0; k<N; ++k)
            habbaaa[k].count = 0;
        }
        if (beta > 3 && (habbaaa[bb].count == 0 || (beta > habbaaa[bb].beta && habbaaa[bb].beta > 2)))
        {
          fprint_char(" X==", a, T, log2_alpha);
          cout << "; ";
          fprint_char(" Y==", b, T, log2_alpha);
          cout << "; ";
          cout << endl;
          fprint_bin("XX==", N+(swapped_aa_bb ? bb : aa), T, "01");
          cout << "; ";
          fprint_bin("YY==", N+(swapped_aa_bb ? aa : bb), T, "01");
          cout << "; ";
          //if (swapped_aa_bb) cout << " swapped_aa_bb";
          cout << endl;
          fprint_bin("XY==", (swapped_ab_ba ? ba : ab), T, "01");
          cout << "; ";
          fprint_bin("YX==", (swapped_ab_ba ? ab : ba), T, "01");
          cout << "; ";
          //if (swapped_ab_ba) cout << " swapped_ab_ba";
          cout << endl;
          if (habbaaa[bb].count != 0 && (beta > habbaaa[bb].beta && habbaaa[bb].beta > 2))
          {
            fprint_char("pX==", habbaaa[bb].X, T, log2_alpha);
            cout << "; ";
            fprint_char("pY==", habbaaa[bb].Y, T, log2_alpha);
            cout << "; ";
            cout << endl;
          }
          else if (habbaaa[bb].count == 0)
          {
            cout << "*** NEW CORRELATION CLASS ***" << endl;
          }
          cout << endl;
          habbaaa[bb].beta = beta;
        }
        if (habbaaa[bb].count == 0)
        {
          habbaaa[bb].beta = beta;
          habbaaa[bb].X = a;
          habbaaa[bb].Y = b;
        }
        habbaaa[bb].count += orbit_size;
      }
    }
    cout << "partitions ==" << setw(20) << partitions
              << " (number of restricted growth strings generated for this value of beta)" << endl;
    output_hist(T, N, hist);
  }

  return 0;
}
