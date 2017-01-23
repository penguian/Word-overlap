
//% generate all overlap correlations
// This file is part of the Word-overlap collection.
// Copyright (C) 2008-2016 Paul Leopardi
// Parts of this code are based on code by Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the src directory.

#include "fxt/bits-all.h"

#include "fxt/fxtiomanip.h"

#include "fxt/fxttypes.h"  // ulong

#include "fxt/jjassert.h"
#include <iomanip>

//#define TIMING // uncomment to disable printing
#ifndef WORD_LEN
#define WORD_LEN 4
#endif
// Number of symbols
#define n  (WORD_LEN)

#define nm1 (n-1)

#ifndef LOG_CHARS
#define LOG_CHARS 2
#endif
// bits per symbol (2**lk)-ary symbols
#define lk (LOG_CHARS)
// max number of coor values
#define N  (1UL<<nm1)

const ulong nb = n*lk;   // number of bits in words
const ulong nbm = nm1*lk;   // number of bits in mask
const ulong Ns = 1UL << nb;  // number of symbols

typedef ulong      anu[N];
typedef ulong*     pu;
typedef pu         anpu[N];
/*
typedef anu*       panu;
typedef panu       anpanu[N];
typedef anpanu*    panpanu;
typedef panpanu    anpanpanu[N];
typedef anpanpanu* panpanpanu;
typedef panpanpanu histtype[N];
*/

void output_hist(const anpu hist)
{
    ulong ct = 0;  // count hist
    ulong ctnz = 0;  // count nonzero values
    for (ulong j=0; j<N; ++j)
        if (hist[j] != NULL)
        {
            const pu hj = hist[j];
            for (ulong k=0; k<N; ++k)
            {
                const ulong hjk = hj[k];
                ct += hjk;
#ifndef TIMING
                if ( hjk )
                {
                    ctnz+=1;
                    print_bin("", j, n);
                    cout << ":";
                    print_bin("", k, n);
                    cout << " * " << setw(10) << hjk*2 << endl;
                }
#endif
            }
        }
    cout << "ctnz=" << ctnz << endl;

    cout << "ct=" << ct << endl;

    const ulong nc = Ns*Ns;
}

int
main(int argc, char **argv)
{
    jjassert( nb <= BITS_PER_LONG );

    cout << " Ns " << Ns << endl;

    anpu hist;
    for (ulong j=0; j<N; ++j)
        hist[j] = NULL;

    const ulong m0 = (1 << nbm) - 1;  // aux mask
    print_bin(" m0=", m0, 0);  cout << endl;

    cout << setw(8) << Ns << endl;
    for (ulong a=0; a<Ns; ++a)
        for (ulong b=0; b<a; ++b)
        {
            ulong ab = 0;  // corr
            ulong ba = 0;  // corr
            ulong bt = b;  // aux, to be shifted
            ulong at = a;  // aux, to be shifted
            ulong m = m0;  // aux, to be shifted
            do
            {
                bt >>= lk;
                ab <<= 1;
                ab |= ( (m &( a^bt )) == 0 );

                at >>= lk;
                ba <<= 1;
                ba |= ( (m &( b^at )) == 0 );

                m >>= lk;
            }
            while ( m );

            pu hab = hist[ab];
            if (hab == NULL)
            {
              hist[ab] = (pu) malloc(sizeof(anu));
              hab = hist[ab];
              for (ulong k=0; k<N; ++k)
                  hab[k] = 0;
            }
            hab[ba] += 1;

//            jjassert(cv < N);
        }

    output_hist(hist);

    return 0;
}
// -------------------------

/*

HTHTTH
HTTHT  0
 HTTHT 0
  HTTH 1
   HTT 0
    HT 0
     H 1


% for n in $(seq 1 12); do ./bin $n 1 | g ctnz; done
ctnz=2
ctnz=4
ctnz=7
ctnz=11
ctnz=17
ctnz=25
ctnz=35
ctnz=48
ctnz=65
ctnz=86
ctnz=113
ctnz=143
*/

/*
Timing:

*/
