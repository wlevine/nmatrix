//definitions copied from LAPACKE
//need to add license info

#ifndef LAPACKE_H
#define LAPACKE_H


#ifndef lapack_int
#define lapack_int     int
#endif

#define LAPACK_ROW_MAJOR               101
#define LAPACK_COL_MAJOR               102

#define LAPACK_WORK_MEMORY_ERROR       -1010
#define LAPACK_TRANSPOSE_MEMORY_ERROR  -1011

#ifndef LAPACKE_malloc
#define LAPACKE_malloc( size ) malloc( size )
#endif
#ifndef LAPACKE_free
#define LAPACKE_free( p )      free( p )
#endif

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

//these are the functions provided by liblapack.so
extern "C" {
  void sgetri_(int* n, float* a, int* lda, const int* ipiv, float* work, int* lwork, int* info);
}

//functions provided by LAPACKE

void LAPACKE_xerbla( const char *name, lapack_int info );

void LAPACKE_sge_trans( int matrix_order, lapack_int m, lapack_int n,
                        const float* in, lapack_int ldin,
                        float* out, lapack_int ldout );

lapack_int LAPACKE_sgetri_work( int matrix_order, lapack_int n, float* a,
                                lapack_int lda, const lapack_int* ipiv,
                                float* work, lapack_int lwork );

lapack_int LAPACKE_sgetri( int matrix_order, lapack_int n, float* a, 
                           lapack_int lda, const lapack_int* ipiv );

#endif // LAPACKE_H
