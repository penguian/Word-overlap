#if !defined  HAVE_SETPART_P_RGS_LEX_H__
#define       HAVE_SETPART_P_RGS_LEX_H__

//#include "sort/minmaxmed23.h"  // max2()

#include "fxttypes.h"

class setpart_p_rgs_lex
// Set partitions of the n-set into p parts (where 2<=p<=n)
// as restricted growth strings (RGS).
// Lexicograpic order.
{
public:
    ulong n_;    // Number of elements of set (set = {1,2,3,...,n})
    ulong p_;    // Exactly p subsets
    ulong *m_;   // m[k+1] = max(s[0], s[1],..., s[k]) + 1
    ulong *s_;   // RGS

public:
    setpart_p_rgs_lex(ulong n, ulong p)
    {
        n_ = n;
        s_ = new ulong[n_];
        m_ = new ulong[n_+1];
        m_[0] = ~0UL;    // sentinel m_[0] = infinity
        first(p);
    }
    
    setpart_p_rgs_lex(ulong n, ulong p, const ulong *s)
    {
        n_ = n;
        s_ = new ulong[n_];
        for (ulong k=0; k<n_; ++k) s_[k] = s[k];
        m_ = new ulong[n_+1];
        m_[0] = ~0UL;    // sentinel m_[0] = infinity
        m_[1] = s_[0]+1;
        for (ulong k=1; k<n_; ++k)  
        {
            ulong mm = m_[k];
            mm += (s_[k]>=mm);
            m_[k+1] = mm;  // == max2(m_[k], s_[k]+1);
        }
        p_ = p;
    }

    ~setpart_p_rgs_lex()
    {
        delete [] m_;
        delete [] s_;
    }

    void first(ulong p)
    // Must have  2<=p<=n
    {
        for (ulong k=0; k<n_; ++k)  s_[k] = 0;
        for (ulong k=n_-p+1, j=1; k<n_; ++k, ++j)  s_[k] = j;

        for (ulong k=1; k<=n_; ++k)  m_[k] = s_[k-1]+1;
        p_ = p;
    }

    bool next()
    {
        // if ( 1==p_ )  return false;  // make things work with p==1

        ulong k = n_;
        bool q;
        do
        {
            --k;
            const ulong sk1 = s_[k] + 1;
            q = (sk1 > m_[k]);   // greater max
            q |= (sk1 >= p_);    // more than p parts
        }
        while ( q && (k !=0) );

        if ( k == 0 )  return false;

        s_[k] += 1UL;
        ulong mm = m_[k];
        mm += (s_[k]>=mm);
        m_[k+1] = mm;  // == max2(m_[k], s_[k]+1);

        while ( ++k<n_ )
        {
            s_[k] = 0;
            m_[k+1] = mm;
        }

        ulong p = p_;
        if ( mm<p )  // repair tail
        {
            do  { m_[k] = p; --k; --p; s_[k] = p; }
            while ( m_[k] < p );
        }

        return true;
    }

    const ulong* data()  const  { return s_; }

    void print()  const;
    void print_set(ulong off=1)  const;
};

bool is_setpart_p_rgs_array(ulong n, ulong p, const ulong *s)
{
    ulong* m = new ulong[n];
    m[0] = s[0]+1;
    for (ulong k=1; k<n; ++k)  
    {
        ulong mm = m[k-1];
        mm += (s[k]>=mm);
        m[k] = mm;  // == max2(m[k-1], s[k]+1);
    }
    delete [] m;
    return (m[n-1] == p);
}
// -------------------------


#endif  // !defined HAVE_SETPART_P_RGS_LEX_H__
