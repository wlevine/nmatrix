//function implementations copied from LAPACKE
//put license info here

#include "lapacke.h"

void LAPACKE_xerbla( const char *name, lapack_int info )
{
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        printf( "Not enough memory to allocate work array in %s\n", name );
    } else if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
        printf( "Not enough memory to transpose matrix in %s\n", name );
    } else if( info < 0 ) {
        printf( "Wrong parameter %d in %s\n", -(int) info, name );
    }
}

void LAPACKE_sge_trans( int matrix_order, lapack_int m, lapack_int n,
                        const float* in, lapack_int ldin,
                        float* out, lapack_int ldout )
{
    lapack_int i, j, x, y;

    if( in == NULL || out == NULL ) return;

    if( matrix_order == LAPACK_COL_MAJOR ) {
        x = n;
        y = m;
    } else if ( matrix_order == LAPACK_ROW_MAJOR ) {
        x = m;
        y = n;
    } else {
        /* Unknown input layout */
        return;
    }

    /* In case of incorrect m, n, ldin or ldout the function does nothing */
    for( i = 0; i < MIN( y, ldin ); i++ ) {
        for( j = 0; j < MIN( x, ldout ); j++ ) {
            out[ (size_t)i*ldout + j ] = in[ (size_t)j*ldin + i ];
        }
    }
}

lapack_int LAPACKE_sgetri_work( int matrix_order, lapack_int n, float* a,
                                lapack_int lda, const lapack_int* ipiv,
                                float* work, lapack_int lwork )
{

  std::cout << matrix_order << std::endl;
  std::cout << n << std::endl;
  std::cout << a[0] << " " << a[1] << std::endl;
  std::cout << lda << std::endl;
  std::cout << ipiv[0] << " " << ipiv[1] << " " << ipiv[2] << std::endl;
    lapack_int info = 0;
    if( matrix_order == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        sgetri_( &n, a, &lda, ipiv, work, &lwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_order == LAPACK_ROW_MAJOR ) {
        lapack_int lda_t = MAX(1,n);
        float* a_t = NULL;
        /* Check leading dimension(s) */
        if( lda < n ) {
            info = -4;
            LAPACKE_xerbla( "LAPACKE_sgetri_work", info );
            return info;
        }
        /* Query optimal working array(s) size if requested */
        if( lwork == -1 ) {
            sgetri_( &n, a, &lda_t, ipiv, work, &lwork, &info );
            return (info < 0) ? (info - 1) : info;
        }
        /* Allocate memory for temporary array(s) */
        a_t = (float*)LAPACKE_malloc( sizeof(float) * lda_t * MAX(1,n) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        /* Transpose input matrices */
        LAPACKE_sge_trans( matrix_order, n, n, a, lda, a_t, lda_t );
        /* Call LAPACK function and adjust info */
        sgetri_( &n, a_t, &lda_t, ipiv, work, &lwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, a_t, lda_t, a, lda );
        /* Release memory and exit */
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_sgetri_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_sgetri_work", info );
    }
    return info;
}

lapack_int LAPACKE_sgetri( int matrix_order, lapack_int n, float* a, 
                           lapack_int lda, const lapack_int* ipiv ) 
{ 
    lapack_int info = 0; 
    lapack_int lwork = -1; 
    float* work = NULL; 
    float work_query; 
    if( matrix_order != LAPACK_COL_MAJOR && matrix_order != LAPACK_ROW_MAJOR ) { 
        LAPACKE_xerbla( "LAPACKE_sgetri", -1 ); 
        return -1; 
    } 
#ifndef LAPACK_DISABLE_NAN_CHECK 
    /* Optionally check input matrices for NaNs */ 
    /*if( LAPACKE_sge_nancheck( matrix_order, n, n, a, lda ) ) { 
        return -3; 
    }*/
#endif 
    /* Query optimal working array(s) size */ 
    info = LAPACKE_sgetri_work( matrix_order, n, a, lda, ipiv, &work_query, 
                                lwork ); 
    if( info != 0 ) { 
        goto exit_level_0; 
    } 
    lwork = (lapack_int)work_query; 
    /* Allocate memory for work arrays */ 
    work = (float*)LAPACKE_malloc( sizeof(float) * lwork ); 
    if( work == NULL ) { 
        info = LAPACK_WORK_MEMORY_ERROR; 
        goto exit_level_0; 
    } 
    /* Call middle-level interface */ 
    info = LAPACKE_sgetri_work( matrix_order, n, a, lda, ipiv, work, lwork ); 
    /* Release memory and exit */ 
    LAPACKE_free( work ); 
exit_level_0: 
    if( info == LAPACK_WORK_MEMORY_ERROR ) { 
        LAPACKE_xerbla( "LAPACKE_sgetri", info ); 
    } 
    return info; 
}
