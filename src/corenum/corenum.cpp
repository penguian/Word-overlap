#include <string>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>

#ifndef NUM_CHARS
#define NUM_CHARS      4
#endif
#ifndef WORD_LEN
#define WORD_LEN       2
#endif

//#define LOG2_NUM_WORDS (2*WORD_LEN)
//#define NUM_WORDS      (1 << LOG2_NUM_WORDS)

typedef unsigned long ulong;
typedef std::string string;
typedef std::bitset<WORD_LEN> bitstring;

const ulong NUM_WORDS = ulong(pow(NUM_CHARS,WORD_LEN));

const string to_str(ulong val)
{
  string word(WORD_LEN,'0');
  for (int i=0; i != WORD_LEN; i++)
  {
    word[i] = '0' + val % NUM_CHARS;
    val /= NUM_CHARS;
  }
  return word;
}

const bitstring corr(const string x, const string y)
{
  bitstring result(false);
  for (int i=1; i != WORD_LEN; i++)
    result[i] = x.substr(i) == y.substr(0,WORD_LEN-i);
  return result;
}

const bitstring autocorr(const string x)
{
  bitstring result(true);
  for (int i=1; i != WORD_LEN; i++)
    result[i] = x.substr(i) == x.substr(0,WORD_LEN-i);
  return result;
}

ulong quad_to_ulong(const bitstring xy, const bitstring yx, const bitstring xx, const bitstring yy)
{
  const ulong wl = WORD_LEN;
  return (xy.to_ulong() << (wl*3)) + (yx.to_ulong() << (wl*2)) + (xx.to_ulong() << wl) + yy.to_ulong();
}

struct corr_pattern
{
  corr_pattern(const bitstring& xyr, const bitstring& yxr, 
               const bitstring& xxr, const bitstring& yyr) :
  xy (xyr), yx (yxr), xx (xxr),  yy (yyr),
  count( 1 )
  {}
  corr_pattern(string xstr, string ystr) :
  xy (corr(xstr,ystr)), yx (corr(ystr,xstr)),
  xx (autocorr(xstr)),  yy (autocorr(ystr)),
  count( 1 )
  {}
  bitstring xy;
  bitstring yx;
  bitstring xx;
  bitstring yy;
  ulong count;
};

int main(int argc, char ** argv)
{

  typedef std::map< ulong, corr_pattern > corr_pattern_map;
  typedef corr_pattern_map::iterator iterator;
  typedef corr_pattern_map::const_iterator const_iterator;

  corr_pattern_map corr_map;

  const std::streamsize width = 10;
  for (ulong x=0; x != NUM_WORDS; ++x)
    for (ulong y=0; y != x; ++y)
    {
      const string xstr = to_str(x);
      const string ystr = to_str(y);
      const bitstring xy = corr(xstr,ystr);
      const bitstring yx = corr(ystr,xstr);
      const bitstring xx = autocorr(xstr);
      const bitstring yy = autocorr(ystr);

      const ulong key = quad_to_ulong(xy, yx, xx, yy);
      const iterator corr_it = corr_map.find(key);
      if (corr_it == corr_map.end())
        corr_map.insert(std::make_pair(key, corr_pattern(xy, yx, yy, xx)));
      else
        corr_it->second.count++;
      /*
      std::cout << std::setw(width)
                << quad_to_ulong(xy, yx, xx, yy)
         << ":" << corr(xstr,ystr)
         << ":" << corr(ystr,xstr)
         << ":" << autocorr(xstr)
         << ":" << autocorr(ystr)
         << ":" << xstr 
         << ":" << ystr 
         << std::endl;
      */
    }
    for (const_iterator 
      corr_it=corr_map.begin(); corr_it != corr_map.end(); ++corr_it)
    {
      const ulong key = corr_it->first;
      const corr_pattern& pattern = corr_it->second;
      std::cout << std::setw(width)
                << key
         << ":" << pattern.xy
         << ":" << pattern.yx
         << ":" << pattern.xx
         << ":" << pattern.yy
         << ":" << std::setw(width) 
                << pattern.count
         << std::endl;
      if ((pattern.xx & pattern.xy).to_ulong())
        std::cout << ": xx&xy == " << pattern.xx & pattern.xy
        << std::endl;
      if ((pattern.yx & pattern.yy).to_ulong())
        std::cout << ": yx&yy == " << pattern.yx & pattern.yy
        << std::endl;
    }
  return 0; 
}
