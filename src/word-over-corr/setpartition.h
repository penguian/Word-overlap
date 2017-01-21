#if !defined  HAVE_SETPARTITION_H__
#define       HAVE_SETPARTITION_H__

#include "fxttypes.h"


class set_partition
// Set partitions of the set {1,2,3,...,n}
// By default in minimal-change order
{
public:
    ulong n_;   // Number of elements of set (set = {1,2,3,...,n})
    int *p_;    // p[] contains set partitions of length 1,2,3,...,n
    int **pp_;  // pp[k] points to start of set partition k
    int *ns_;   // ns[k] Number of Sets in set partition k
    int *as_;   // element k attached At Set (0<=as[k]<=k) of set(k-1)
    int *d_;    // direction with recursion (+1 or -1)
//    int *h_;    // auxiliary array
    int *x_;    // current set partition (==pp[n])
    bool xdr_;  // whether to change direction in recursion (==> minimal-change order)
    int dr0_;   // dr0: starting direction in each recursive step:
    //   dr0=+1  ==> start with partition  {{1,2,3,...,n}}
    //   dr0=-1  ==> start with partition  {{1},{2},{3},...,{n}}}

public:
    set_partition(ulong n, bool xdr=true, int dr0=+1)
        : n_(n)
    {
        ulong np = (n_*(n_+1))/2;  // == \sum_{k=1}^{n}{k}
        p_ = new int[np];

        pp_ = new int *[n_+1];
        pp_[0] = 0;  // unused
        pp_[1] = p_;
        for (ulong k=2; k<=n_; ++k)  pp_[k] = pp_[k-1] + (k-1);

        ns_ = new int[n_+1];
        as_ = new int[n_+1];
        d_ = new int[n_+1];
//        h_ = new int[n_+1];
        x_ = pp_[n_];

        init(xdr, dr0);
    }

    ~set_partition()
    {
        delete [] p_;
        delete [] pp_;
        delete [] ns_;
        delete [] as_;
        delete [] d_;
//        delete [] h_;
    }

    void init(bool xdr, int dr0);

    bool next()  { return next_rec(n_); }

    const int* data()  const  { return x_; }

    ulong print()  const
    // Print current set partition
    // Return number of chars printed
    { return print_p(n_); }

    ulong print_p(ulong k)  const;
    void print_internal()  const;  // print internal state

protected:
    int cp_append(const int *src, int *dst, ulong k, ulong a);
    int next_rec(ulong k);
};
// -------------------------


#endif  // !defined HAVE_SETPARTITION_H__
