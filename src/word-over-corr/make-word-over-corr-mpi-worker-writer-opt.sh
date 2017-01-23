export CFLAGS="-O3 -pipe -march=native -floop-interchange -floop-strip-mine -floop-block -ftree-vectorize -flto"
export DEFINES="-DWOC_PAIR_BUFFER_SIZE=100000 -DWOC_DELTA_OFFSET=2 -DWOC_PRINT_NEW_CLASSES -DWOC_TOTALS_ONLY"
export INCLUDES="-I. -I$LOCAL/include -I$LOCAL/include/fxt"
export LDFLAGS="-Wl,-O1 -Wl,--as-needed -Wl,--hash-style=gnu -Wl,--sort-common -Wl,-flto"
export LIBS="-L$LOCAL/lib -lfxt -lmpi"
mpic++ $CFLAGS $DEFINES $INCLUDES $LDFLAGS $LIBS word-over-corr-mpi-worker-writer.cc -o word-over-corr-mpi-worker-writer-opt.exe
