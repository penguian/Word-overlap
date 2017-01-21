//% generate all word overlap correlations

#include "mpi.h"
#include <stdlib.h>

#include "fxt/bits-all.h"
//#include "fxt/comb/setpart-p-rgs-lex.h"
#include "setpart-p-rgs-lex-exp.h"

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

typedef ulong*     pu;  // Pointer to ulong 
typedef ulong**    ppu; // Pointer to pointer to ulong
typedef ulong***   p3u; // Pointer to pointer to pointer to ulong
typedef ulong****  p4u; // Pointer to pointer to pointer to pointer to ulong

typedef std::string string;

const int END_TAG  = 0;
const int MORE_TAG = 1;

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

void output_hist(const ulong N, const ulong T, const p4u hist)
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
                const ulong habbaaabb = habbaaa[bb];
                if ( habbaaabb )
                {
                    pairs += habbaaabb;
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
                    cout << "count[" << classes << "]:=" << habbaaabb << ":" << endl;
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

inline ulong autocorr(const ulong log2_alpha, ulong m, const ulong x)
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

#ifndef WORD_LEN
#define WORD_LEN 2
#endif

void calc_ab_ba(ulong *ab, ulong *ba, const ulong log2_alpha, ulong m, ulong a, ulong b)
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

void insert_hist(const p4u hist, const ulong N, const ulong T, const ulong log2_alpha, const ulong beta, const ulong orbit_size, const ulong m0, 
                 const ulong a, const ulong b, ulong aa, ulong bb, ulong ab, ulong ba)
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
  bool swapped_aa_bb = false;
  if (aa < bb)
  {
    swap(aa, bb);
    swapped_aa_bb = true;
  }
  pu habbaaa = habba[aa];
  if (habbaaa == NULL)
  {
    habba[aa] = (pu) malloc(N*sizeof(ulong));
    habbaaa = habba[aa];
    for (ulong k=0; k<N; ++k)
      habbaaa[k] = 0;
  }
  if (beta>3 && habbaaa[bb] == 0)
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
    cout << endl;
  }
  habbaaa[bb] += orbit_size;
}

bool arrays_are_equal(ulong n, const ulong* s, const ulong *t)
{
  for (ulong k=0; k!=n; k++)
    if (s[k] != t[k])
      return false;
  return true;
}

ulong* new_setpart_p_rgs_array(const ulong T, const ulong beta, const ulong gamma)
{
    ulong* s = new ulong[2*T];
    for (ulong k=0; k < 2*T; k++)
      s[k] = 0UL;
    if (gamma > beta)
      s[0] = ~0UL;
    else
    {
      for (ulong k=0; k < gamma; k++)
        s[k] = k;
      for (ulong k=0; gamma+k < beta; k++)
        s[2*T-beta+gamma+k] = gamma+k;
      //for (ulong k=0; k < 2*T; k++)
        //cout << s[k] << " ";
      //cout << endl;
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

ulong check_setpart_p_rgs_range(const ulong T, const ulong log2_alpha, const ulong m0, const ulong N,
                                const ulong beta, const ulong orbit_size, 
                                const ulong* rgs_begin,const ulong* rgs_end, const p4u hist, const int myid)
{
  setpart_p_rgs_lex p = setpart_p_rgs_lex(2*T,beta,rgs_begin);
  ulong olda = 0;
  ulong aa = autocorr(log2_alpha, m0, olda);
  ulong partitions = 0;
  bool more;
  ulong buffer[6];
  
  for (more=true; more; more=p.next())
  {
    const ulong* rgs_data = p.data();

    if (arrays_are_equal(2*T, rgs_data, rgs_end))
    {
      if (myid > 0)
      {
        buffer[0] = ~0UL;
        MPI_Send(buffer, 6, MPI_UNSIGNED_LONG, 0, END_TAG, MPI_COMM_WORLD);
      }
      break;
    }
    //for (ulong k=0;k!=2*T;k++)
     //cout << rgs_data[k];
    //cout << endl;
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
    if (a != olda)
    {
      aa = autocorr(log2_alpha, m0, a);
      olda = a;
    }
    if (a!=b)
    {
      ulong bb = autocorr(log2_alpha, m0, b);
      ulong ab = 0;  // corr
      ulong ba = 0;  // corr
      calc_ab_ba(&ab, &ba, log2_alpha, m0, a, b);
      if (myid == 0)
        insert_hist(hist, N, T, log2_alpha, beta, orbit_size, m0, a, b, aa, bb, ab, ba); 
      else
      {
        buffer[0] = a;   buffer[1] = b;
        buffer[2] = aa;  buffer[3] = bb;
        buffer[4] = ab;  buffer[5] = ba;
        MPI_Send(buffer, 6, MPI_UNSIGNED_LONG, 0, MORE_TAG, MPI_COMM_WORLD);
      }
    }
  }
  if (!more && (myid > 0))
  {
    buffer[0] = ~0UL;
    MPI_Send(buffer, 6, MPI_UNSIGNED_LONG, 0, END_TAG, MPI_COMM_WORLD);
  }
  return partitions;
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  
  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  MPI_Get_processor_name(processor_name,&namelen);
  fprintf(stderr,"Process %d of %d on %s\n", myid, numprocs, processor_name);

  // Number of characters in a word
  ulong T = WORD_LEN;
  if (argc > 1)
    sscanf(argv[1],"%lu",&T);

  // Number of different characters in the alphabet
  ulong alpha = T*2;
  if (argc > 2)
    sscanf(argv[2],"%lu",&alpha);

  // Maximum number of ranges
  ulong max_ranges = T*3/2;
  if (argc > 3)
    sscanf(argv[3],"%lu",&max_ranges);
  max_ranges = (max_ranges>2*T) ? 2*T : ( (max_ranges > (ulong) numprocs-1) ? (numprocs-1) : max_ranges );
  
  // Maximum number of bits needed to represent a character
  const ulong log2_alpha = (ulong) std::ceil(std::log(alpha)/std::log(2.0));
#ifdef WOC_DEBUG
  // Maximum number of bits needed to represent a word
  const ulong nbw = T*log2_alpha;
#endif
  // Number of bits in mask
  const ulong nbm = (T-1)*log2_alpha;
  // Initial mask
  const ulong m0 = (1UL << nbm) - 1UL;
  // Maximum number of different correlation values
  const ulong N = 1UL<<(T-1);

  p4u hist;
  
  //jjassert(nbw <= BITS_PER_LONG);
  //jjassert(T > 1);
  //jjassert(alpha <= 2*T);
  if (myid == 0)
  {
    cout << "T          ==" << setw(20) << T
              << " (number of characters in a word)" << endl;
    cout << "alpha      ==" << setw(20) << alpha
              << " (number different characters in the alphabet)" << endl;
    cout << "max_ranges ==" << setw(20) << max_ranges
              << " (maximum number of restricted growth string ranges to use)" << endl;
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

    hist = (p4u) malloc(N*sizeof(p3u));
    for (ulong k=0; k<N; ++k)
      hist[k] = NULL;  
  }

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
    if (beta > max_ranges)
      if (myid == 0)
      {
        ulong* rgs_begin = new_setpart_p_rgs_array(T, beta, 1);
        if (rgs_begin == (ulong *) NULL)
          return -1;
        
        ulong* rgs_end = new_setpart_p_rgs_array(T, beta, beta+1);
        if (rgs_end == (ulong *) NULL)
          return -1;

        partitions = check_setpart_p_rgs_range(T, log2_alpha, m0, N, beta, orbit_size, 
                                              rgs_begin, rgs_end, hist, myid);
        delete [] rgs_begin;
        delete [] rgs_end;
      }
    if (beta <= max_ranges)
    {
      ulong steps = 0;
      if (myid == 0)
      {
        ulong nranges = beta;
        ulong buffer[6];
        MPI_Status status;

        while (nranges!=0)
        {
          MPI_Recv(buffer, 6, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          const ulong a  = buffer[0];  const ulong b  = buffer[1];
          const ulong aa = buffer[2];  const ulong bb = buffer[3];
          const ulong ab = buffer[4];  const ulong ba = buffer[5];
          if (a == ~0UL)
            nranges--;

          else
            insert_hist(hist, N, T, log2_alpha, beta, orbit_size, m0, a, b, aa, bb, ab, ba); 
        }
      }
      else
      {
        ulong gamma = myid;
        if (gamma <= beta)
        {
          ulong* rgs_begin = new_setpart_p_rgs_array(T, beta, gamma);
          if (rgs_begin == (ulong *) NULL)
            return -1;
          
          ulong* rgs_end = new_setpart_p_rgs_array(T, beta, gamma+1);
          if (rgs_end == (ulong *) NULL)
            return -1;

          steps = check_setpart_p_rgs_range(T, log2_alpha, m0, N, beta, orbit_size, 
                                           rgs_begin, rgs_end, hist, myid);
          delete [] rgs_begin;
          delete [] rgs_end;
          /*
          if (steps > 0)
          {
            cout << "gamma      ==" << setw(20) << gamma << endl;
            cout << "steps      ==" << setw(20) << steps
                 << " (number of restricted growth strings generated for this value of beta, gamma)" 
                 << endl;
          }
          */
        }
      }
      MPI_Reduce(&steps, &partitions, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if (myid == 0)
    {
      cout << "partitions ==" << setw(20) << partitions
           << " (number of restricted growth strings generated for this value of beta)" 
           << endl;
      output_hist(N, T, hist);
    }
  }
  MPI_Finalize(); 
  return 0;
}
