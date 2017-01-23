
//% generate all overlap correlations
// This file is part of the Word-overlap collection.
// Copyright (C) 2008-2016 Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the src directory.


#include "fxt/bits-all.h"

#include "fxt/fxtiomanip.h"

//#include "demo/nextarg.h"
#include "fxt/fxttypes.h"  // ulong

#include "fxt/jjassert.h"


//#define TIMING // uncomment to disable printing

#define lk 2

int
main(int argc, char **argv)
{
    ulong n = 4;
    //NXARG(n, "Number of symbols");
    const ulong N = (1UL<<n);   // max number of coor values

//    ulong lk = 1;
//    NXARG(lk, "bits per symbol (2**lk)-ary symbols");
    const ulong nb = n*lk;   // number of bits in words
    jjassert( nb <= BITS_PER_LONG );
    const ulong Ns = 1UL << nb;  // number of symbols
    cout << " Ns " << Ns << endl;

    ulong *hist = new ulong[N];
    for (ulong j=0; j<N; ++j)  hist[j] = 0;

    const ulong m0 = first_comb(nb);  // aux mask
    print_bin(" m0=", m0, 0);  cout << endl;

    for (ulong a=0; a<Ns; ++a)
    {
        for (ulong b=0; b<Ns; ++b)
        {
            ulong cv = 0;  // corr
            ulong bt = b;  // aux, to be shifted
            ulong m = m0;  // aux, to be shifted
            do
            {
                cv <<= 1;
                cv |= ( (m&( a^bt )) == 0 );
                bt >>= lk;
                m >>= lk;
            }
            while ( m );

            hist[cv] += 1;

//            jjassert(cv < N);
        }
    }

    ulong ct = 0;  // count hist
    ulong ctnz = 0;  // count nonzero values
    for (ulong j=0; j<N; ++j)
    {
        const ulong hj = hist[j];
        ct += hj;
#ifndef TIMING
        if ( hj )
        {
            ctnz+=1;
            print_bin("", j, n);
//            cout << j << ", ";
            cout << " * " << hj << endl;
        }
#endif
    }
    cout << "ctnz=" << ctnz << endl;

    cout << "ct=" << ct << endl;

    const ulong nc = Ns*Ns;
    jjassert( ct == nc );

    delete [] hist;

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
