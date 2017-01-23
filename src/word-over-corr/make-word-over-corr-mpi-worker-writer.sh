export CFLAGS="-O1 -g"
export DEFINES="-DWOC_PRINT_NEW_CLASSES" # "-DWOC_TOTALS_ONLY"
export INCLUDES="-I. -I$LOCAL/include -I$LOCAL/include/fxt"
export LIBS="-L$LOCAL/lib -lfxt -lmpi"
mpic++ $CFLAGS $DEFINES $INCLUDES $LIBS word-over-corr-mpi-worker-writer.cc -o word-over-corr-mpi-worker-writer.exe
