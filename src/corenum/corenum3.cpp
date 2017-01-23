#include <string>
#include <iosfwd>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <bitset>
#include <map>
// This file is part of the Word-overlap collection.
// Copyright (C) 2008-2016 Paul Leopardi
// Parts of this code are based on code by Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the src directory.

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

string poly(const bitstring x)
{
  std::ostringstream op;
  op << "z->";
  if (x.count() == 0)
    op << "0";
  else
  {
    int i;
    for (i=0; i != WORD_LEN; i++)
      if (x[i])
      {
        if (i == 0)
          op << "1";
        else
          op << "z^" << i;
        break;
      }
    for (i++; i <= WORD_LEN; i++)
      if (x[i])
        op << "+z^" << i;
  }
  return op.str();
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
  bitstring xy;
  bitstring yx;
  bitstring xx;
  bitstring yy;
  ulong count;
};

void swap(bitstring& a, bitstring& b)
{
  const bitstring tmp = a;
  a = b;
  b = tmp;
}

int main(int argc, char ** argv)
{

  typedef std::map< ulong, corr_pattern > corr_pattern_map;
  typedef corr_pattern_map::iterator iterator;
  typedef corr_pattern_map::const_iterator const_iterator;

  corr_pattern_map corr_map;

  const std::streamsize width = 10;
  for (ulong x=0; x != NUM_WORDS; ++x)
  {
    const string xstr = to_str(x);
    bitstring xxx = autocorr(xstr);
    for (ulong y=0; y != x; ++y)
    {
      const string ystr = to_str(y);
      bitstring xy = corr(xstr,ystr);
      bitstring yx = corr(ystr,xstr);
      bitstring xx = xxx;
      bitstring yy = autocorr(ystr);
      if (yx.to_ulong() < xy.to_ulong())
        swap(yx, xy);
      if (yy.to_ulong() < xx.to_ulong())
        swap(yy, xx);
      const ulong key = quad_to_ulong(xy, yx, xx, yy);
      const iterator corr_it = corr_map.find(key);
      if (corr_it == corr_map.end())
        corr_map.insert(std::make_pair(key, corr_pattern(xy, yx, yy, xx)));
      else
        corr_it->second.count++;
    }
  }
  ulong i=0;
  for (const_iterator
    corr_it=corr_map.begin(); corr_it != corr_map.end(); ++corr_it, ++i)
  {
    const ulong key = corr_it->first;
    const corr_pattern& pattern = corr_it->second;
    std::cout
        <<   "pXY[" << i << "]:=" << poly(pattern.xy)
        << ": pYX[" << i << "]:=" << poly(pattern.yx)
        << ": pXX[" << i << "]:=" << poly(pattern.xx)
        << ": pYY[" << i << "]:=" << poly(pattern.yy)
        << ": count[" << i << "]:= " << pattern.count * 2
        << " :" << std::endl;
  }
  return 0;
}
