/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2014, Ruby Science Foundation
// NMatrix is Copyright (c) 2012 - 2014, John Woods and the Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == getri.h
//
// getri function in native C++.
//

/*
 *             Automatically Tuned Linear Algebra Software v3.8.4
 *                    (C) Copyright 1999 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef GETRI_H
#define GETRI_H

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

extern "C" {
  void sgetri_(int* n, float* a, int* lda, const int* ipiv, float* work, int* lwork, int* info);
}

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

namespace nm { namespace math {

template <typename DType>
inline int getri(const enum CBLAS_ORDER order, const int n, DType* a, const int lda, const int* ipiv) {
  rb_raise(rb_eNotImpError, "getri not yet implemented for non-BLAS dtypes");
  return 0;
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, float* a, const int lda, const int* ipiv) {
  return LAPACKE_sgetri(order, n, a, lda, ipiv);
}

#if defined (HAVE_CLAPACK_H) || defined (HAVE_ATLAS_CLAPACK_H)
template <>
inline int getri(const enum CBLAS_ORDER order, const int n, double* a, const int lda, const int* ipiv) {
  return clapack_dgetri(order, n, a, lda, ipiv);
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, Complex64* a, const int lda, const int* ipiv) {
  return clapack_cgetri(order, n, reinterpret_cast<void*>(a), lda, ipiv);
}

template <>
inline int getri(const enum CBLAS_ORDER order, const int n, Complex128* a, const int lda, const int* ipiv) {
  return clapack_zgetri(order, n, reinterpret_cast<void*>(a), lda, ipiv);
}
#endif

/*
 * Function signature conversion for calling LAPACK's getri functions as directly as possible.
 *
 * For documentation: http://www.netlib.org/lapack/double/dgetri.f
 *
 * This function should normally go in math.cpp, but we need it to be available to nmatrix.cpp.
 */
template <typename DType>
inline int clapack_getri(const enum CBLAS_ORDER order, const int n, void* a, const int lda, const int* ipiv) {
  return getri<DType>(order, n, reinterpret_cast<DType*>(a), lda, ipiv);
}


} } // end nm::math

#endif // GETRI_H
