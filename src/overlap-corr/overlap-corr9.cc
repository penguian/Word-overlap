//% generate all overlap correlations

#include "fxt/bits-all.h"

#include "fxt/fxtiomanip.h"

#include "fxt/fxttypes.h"  // ulong

#include "fxt/jjassert.h"

#include <string>
#include <iosfwd>
#include <iostream>
#include <iomanip>
#include <sstream>

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
typedef ulong**    ppu;
typedef ppu        anppu[N];
typedef ulong***   pppu;
typedef pppu       anpppu[N];

typedef std::string string;

string poly(ulong x)
{
  std::ostringstream op;
  op << "z->";
  if (x == 0)
    op << "0";
  else
  {
    int i;
    for (i=nm1; i >= 0; --i, x >>= 1)
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

void output_hist(const anpppu hist)
{
  ulong ct = 0;  // count hist
  ulong ctnz = 0;  // count nonzero values

  for (ulong ab=0; ab<N; ++ab)
  {
    const pppu hab = hist[ab];
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
                const ulong habbaaabb = habbaaa[bb];
                if ( habbaaabb )
                {
                    ct += habbaaabb;
                    ctnz+=1;
#ifndef TIMING
                    //print_bin("", ab, n);
                    cout << "pXY[" << ctnz << "]:=" << poly(ab) << ":";
                    //print_bin("", ba, n);
                    cout << "pYX[" << ctnz << "]:=" << poly(ba) << ":";
                    //print_bin("", N+aa, n);
                    cout << "pXX[" << ctnz << "]:=" << poly(N+aa) << ":";
                    //print_bin("", N+bb, n);
                    cout << "pYY[" << ctnz << "]:=" << poly(N+bb) << ":";
                    cout << "count[" << ctnz << "]:=" << habbaaabb*2 << ":" << endl;
#endif
                }
              }
          }
        }
    }  
    //cout << "ctnz=" << ctnz << endl;
    cout << "ct=" << ct << endl;
}

int
main(int argc, char **argv)
{
  jjassert( nb <= BITS_PER_LONG );
  //cout << " Ns " << setw(10) << Ns << endl;

  const ulong m0 = (1 << nbm) - 1;  // aux mask
  //print_bin(" m0=", m0, 0);  cout << endl;

  ulong autohist[Ns];
  for (ulong x=0; x<Ns; ++x)
  {
      ulong xx = 0;  // corr
      ulong xt = x;  // aux, to be shifted
      ulong m = m0;  // aux, to be shifted
      do
      {
          xt >>= lk;
          xx <<= 1;
          xx |= ( (m &( x^xt )) == 0 );
          m >>= lk;
      }
      while ( m );
      autohist[x] = xx;
  }

  anpppu hist;
  for (ulong ab=0; ab<N; ++ab)
      hist[ab] = NULL;  

  for (ulong a=0; a<Ns; ++a)
  {
    const ulong aa = autohist[a];
    for (ulong b=0; b<a; ++b)
    {
      const ulong bb = autohist[b];
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

      pppu hab = hist[ab];
      if (hab == NULL)
      {
        hist[ab] = (pppu) malloc(sizeof(anppu));
        hab = hist[ab];
        for (ulong k=0; k<N; ++k)
            hab[k] = NULL;
      }
      ppu habba = hab[ba];
      if (habba == NULL)
      {
        hab[ba] = (ppu) malloc(sizeof(anpu));
        habba = hab[ba];
        for (ulong k=0; k<N; ++k)
            habba[k] = NULL;
      }
      pu habbaaa = habba[aa];
      if (habbaaa == NULL)
      {
        habba[aa] = (pu) malloc(sizeof(anu));
        habbaaa = habba[aa];
        for (ulong k=0; k<N; ++k)
            habbaaa[k] = 0;
      }
      habbaaa[bb] += 1;
    }
  }
  output_hist(hist);

  return 0;
}
