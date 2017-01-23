
//% generate all overlap correlations
// This file is part of the Word-overlap collection.
// Copyright (C) 2008-2016 Paul Leopardi
// Parts of this code are based on code by Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the src directory.

#include "fxt/bits-all.h"

#include "fxt/fxtiomanip.h"

//#include "demo/nextarg.h"
#include "fxt/fxttypes.h"  // ulong

#include "fxt/jjassert.h"
#include <iomanip>

//#define TIMING // uncomment to disable printing
#ifndef WORD_LEN
#define WORD_LEN 4
#endif
#define n  (WORD_LEN)

#ifndef LOG_CHARS
#define LOG_CHARS 2
#endif
#define lk (LOG_CHARS)

#define N  (1UL<<n)

int
main(int argc, char **argv)
{
    //ulong n = 2;
    //NXARG(n, "Number of symbols");
    //const ulong N = (1UL<<n);   // max number of coor values

//    ulong lk = 1;
//    NXARG(lk, "bits per symbol (2**lk)-ary symbols");
    const ulong nb = n*lk;   // number of bits in words
    jjassert( nb <= BITS_PER_LONG );
    const ulong Ns = 1UL << nb;  // number of symbols
    cout << " Ns " << Ns << endl;

    ulong (*(hist[N]))[N];
    for (ulong j=0; j<N; ++j)
        hist[j] = NULL;

    const ulong m0 = (1 << nb) - 1;  // aux mask
    print_bin(" m0=", m0, 0);  cout << endl;

    cout << setw(8) << Ns << endl;
    for (ulong a=0; a<Ns; ++a)
    {
        cout << setw(8) << a << endl;
        for (ulong b=0; b<a; ++b)
        {
            ulong ab = 0;  // corr
            ulong ba = 0;  // corr
            ulong bt = b;  // aux, to be shifted
            ulong at = a;  // aux, to be shifted
            ulong m = m0;  // aux, to be shifted
            do
            {
                ab <<= 1;
                ab |= ( (m &( a^bt )) == 0 );
                bt >>= lk;

                ba <<= 1;
                ba |= ( (m &( b^at )) == 0 );
                at >>= lk;

                m >>= lk;
            }
            while ( m );

            if (hist[ab] == NULL)
            {
              hist[ab] = (ulong(*)[N]) malloc(N*sizeof(ulong));
              for (ulong k=0; k<N; ++k)
                  (*(hist[ab]))[k] = 0;
            }
            (*(hist[ab]))[ba] += 1;

//            jjassert(cv < N);
        }
    }
    ulong ct = 0;  // count hist
    ulong ctnz = 0;  // count nonzero values
    for (ulong j=0; j<N; ++j)
        if (hist[j] != NULL)
            for (ulong k=0; k<N; ++k)
            {
                const ulong hjk = (*(hist[j]))[k];
                ct += hjk;
#ifndef TIMING
                if ( hjk )
                {
                    ctnz+=1;
                    print_bin("", j, n);
                    cout << ":";
                    print_bin("", k, n);
        //            cout << j << ", ";
                    cout << " * " << setw(10) << hjk*2 << endl;
                }
#endif
            }
    cout << "ctnz=" << ctnz << endl;

    cout << "ct=" << ct << endl;

    const ulong nc = Ns*Ns;
    //jjassert( ct == nc );

    //delete [] hist;

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
