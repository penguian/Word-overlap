export CFLAGS="-O3 -pipe -march=native -floop-interchange -floop-strip-mine -floop-block -ftree-vectorize -flto"
export LDFLAGS="-Wl,-O1 -Wl,--as-needed -Wl,--hash-style=gnu -Wl,--sort-common -Wl,-flto"
mpic++ $CFLAGS -DWOC_PAIR_BUFFER_SIZE=100000 -DWOC_DELTA_OFFSET=2 -DTOTALS_ONLY -I. -I$LOCAL/include -I$LOCAL/include/fxt $LDFLAGS -L$LOCAL/lib -lfxt -lmpi word-over-corr-mpi-worker-hist.cc -o word-over-corr-mpi-worker-hist-opt.exe
