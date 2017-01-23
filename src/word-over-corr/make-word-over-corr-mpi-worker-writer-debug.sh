export CFLAGS="-O0 -g3"
export DEFINES="-DWOC_DEBUG -DWOC_MALLOC_DEBUG -DWOC_MPI_DEBUG -DWOC_MPI_VERBOSE -DWOC_PRINT_NEW_CLASSES"
export INCLUDES="-I. -I$LOCAL/include -I$LOCAL/include/fxt"
export LIBS="-L$LOCAL/lib -lfxt -lmpi"
mpic++ $CFLAGS $DEFINES $INCLUDES $LIBS word-over-corr-mpi-worker-writer.cc -o word-over-corr-mpi-worker-writer-debug.exe
