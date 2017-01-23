#include <string>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <cmath>
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
    word[WORD_LEN-1-i] = '0' + val % NUM_CHARS;
    val /= NUM_CHARS;
  }
  return word;
}

const bitstring corr(const string x, const string y)
{
  bitstring result(false);
  for (int i=1; i != WORD_LEN; i++)
    result[WORD_LEN-1-i] = x.substr(i) == y.substr(0,WORD_LEN-i);
  return result;
}

const bitstring autocorr(const string x)
{
  bitstring result(false);
  result[WORD_LEN-1] = true;
  for (int i=1; i != WORD_LEN; i++)
    result[WORD_LEN-1-i] = x.substr(i) == x.substr(0,WORD_LEN-i);
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
               const bitstring& xxr, const bitstring& yyr,
               const string& xstrr, const string& ystrr) :
  xy (xyr), yx (yxr), xx (xxr),  yy (yyr),
  xstr (xstrr), ystr (ystrr),
  count( 1 )
  {}
  corr_pattern(string xstrr, string ystrr) :
  xy (corr(xstrr,ystrr)), yx (corr(ystrr,xstrr)),
  xx (autocorr(xstrr)),  yy (autocorr(ystrr)),
  xstr (xstrr), ystr (ystrr),
  count( 1 )
  {}
  bitstring xy;
  bitstring yx;
  bitstring xx;
  bitstring yy;
  string xstr;
  string ystr;
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
  for (ulong x=0; x != NUM_WORDS/NUM_CHARS; ++x)
  {
    const string xstr = to_str(x);
    bitstring xx = autocorr(xstr);
    for (ulong y=0; y != NUM_WORDS; ++y)
    if (y != x)
    {
      const string ystr = to_str(y);
      bitstring xy = corr(xstr,ystr);
      bitstring yx = corr(ystr,xstr);
      bitstring xxx = xx;
      bitstring yy = autocorr(ystr);
      if (yx.to_ulong() < xy.to_ulong())
        swap(yx, xy);
      if (yy.to_ulong() < xxx.to_ulong())
        swap(yy, xxx);
      const ulong key = quad_to_ulong(xy, yx, xxx, yy);
      const iterator corr_it = corr_map.find(key);
      if (corr_it == corr_map.end())
      {
        const corr_pattern& pattern = corr_pattern(xy, yx, xxx, yy, xstr, ystr);
        corr_map.insert(std::make_pair(key, pattern));
      }
      else
        corr_it->second.count++;
    }
  }
  for (const_iterator
    corr_it=corr_map.begin(); corr_it != corr_map.end(); ++corr_it)
  {
    const ulong key = corr_it->first;
    const corr_pattern& pattern = corr_it->second;
    std::cout
        << ":" << pattern.xy
        << ":" << pattern.yx
        << ":" << pattern.xx
        << ":" << pattern.yy
        << ":" << pattern.xstr
        << ":" << pattern.ystr
        << ":" << std::setw(width)
              << pattern.count * 2
        << std::endl;
  }
  return 0;
}
