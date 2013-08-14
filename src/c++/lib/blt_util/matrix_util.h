// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// CodeAxe : phylogenetic analysis and simulation tools
//
//   http://www.phrap.org
//
//
// Copyright 2007 Christopher T Saunders (ctsa@u.washington.edu)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
//
//

// $Id: matrix_util.h 1055 2007-12-08 04:27:31Z ctsa $

/// \file
///
/// \brief small matrix ops
///

#ifndef __MATRIX_UTIL_H
#define __MATRIX_UTIL_H

#include "array_util.h"


#include <cmath>

#include <complex>
#include <limits>


/// \brief test square matrix for symmetry
///
template <typename FloatType>
inline
bool
matrix_is_symm(FloatType m[],
               const unsigned N)
{
    for (unsigned i(0); i<(N-1); ++i)
        for (unsigned j(i+1); j<N; ++j)
            if ( std::fabs(m[j+i*N]-m[i+j*N]) > std::numeric_limits<FloatType>::epsilon() ) return false;
    return true;
}


/// \brief get NxN zero matrix
///
template <typename FloatType>
inline
void
matrix_zero(FloatType m[],
            const unsigned N)
{
    array_zero(m,N*N);
}


/// \brief get NxN identity matrix
///
template <typename FloatType>
inline
void
matrix_identity(FloatType m[],
                const unsigned N)
{
    static const FloatType one(static_cast<FloatType>(1));
    matrix_zero<FloatType>(m,N);
    for (unsigned i(0); i<N; ++i) m[i*(N+1)] = one;
}


/// \brief set matrix diagonal so that "from" (index1) sums to 0
///
template <typename FloatType>
inline
void
matrix_balance_diagonal(FloatType m[],
                        const unsigned N)
{
    static const FloatType zero(static_cast<FloatType>(0));
    for (unsigned i(0); i<N; ++i)
    {
        FloatType sum(zero);
        for (unsigned j(0); j<N; ++j)
        {
            if ( i != j ) sum -= m[j+i*N];
        }
        m[i*(N+1)] = sum;
    }
}


/// \brief copy one NxN matrix to another
///
template <typename FloatType>
inline
void
matrix_copy(const FloatType from[],
            FloatType to[],
            const unsigned N)
{
    std::copy(from,from+N*N,to);
}


template <typename FloatType>
inline
void
matrix_transpose_inplace(FloatType m[],
                         const unsigned N)
{
    for (unsigned i(0); i<(N-1); ++i)
        for (unsigned j(i+1); j<N; ++j)
            std::swap(m[j+i*N],m[i+j*N]);
}

/// \brief copy one NxN matrix to its transpose
///
template <typename FloatType>
inline
void
matrix_transpose_copy(const FloatType from[],
                      FloatType to[],
                      const unsigned N)
{
    for (unsigned i(0); i<N; ++i)
        for (unsigned j(0); j<N; ++j)
            to[j+i*N] = from[i+j*N];
}


template <typename FloatType1,
         typename FloatType2>
inline
void
matrix_scale(FloatType1 m[],
             const unsigned N,
             const FloatType2 s)
{
    array_scale(m,N*N,s);
}


/// \brief return to = to + from
///
template <typename FloatType>
inline
void
matrix_sum(const FloatType from[],
           FloatType to[],
           const unsigned N)
{
    const unsigned N2(N*N);
    for (unsigned i(0); i<N2; ++i) to[i] += from[i];
}


template <typename FloatType>
inline
void
matrix_scale_sum(const FloatType from[],
                 const FloatType s,
                 FloatType to[],
                 const unsigned N)
{
    array_scale_sum(from,s,to,N*N);
}


template <typename FloatType>
inline
void
matrix_mult(const FloatType from1[],
            const FloatType from2[],
            FloatType to[],
            const unsigned N,
            const bool is_symm=false)
{
    scale_matrix_mult(FloatType(1),from1,from2,to,N,is_symm);
}

// scale*from1*from2
//
template <typename FloatType>
inline
void
scale_matrix_mult(const FloatType scale,
                  const FloatType from1[],
                  const FloatType from2[],
                  FloatType to[],
                  const unsigned N,
                  const bool is_symm=false)
{
    matrix_zero<FloatType>(to,N);
    for (unsigned i(0); i<N; ++i)
    {
        for (unsigned j(0); j<N; ++j)
        {
            for (unsigned k(0); k<N; ++k)
            {
                to[j+i*N] += from1[k+i*N]*from2[j+k*N];
            }
            to[j+i*N] *= scale;
        }
    }
}

#ifdef USE_BLAS
template <>
inline
void
scale_matrix_mult(const float scale,
                  const float from1[],
                  const float from2[],
                  float to[],
                  const unsigned N,
                  const bool is_symm)
{
#ifdef FORTRAN_BLAS
    static float fzero(0.);
    static char nchar('N');
    static char lchar('L');
    static char uchar('U');
    int ncopy(N);
    if (is_symm)
    {
        FBLAS::ssymm(&lchar,&uchar,ncopy,ncopy,const_cast<float>(scale),const_cast<float*>(from1),
                     ncopy,const_cast<float*>(from2),ncopy,fzero,to,ncopy);
    }
    else
    {
        FBLAS::sgemm(&nchar,&nchar,ncopy,ncopy,ncopy,const_cast<float>(scale),const_cast<float*>(from1),
                     ncopy,const_cast<float*>(from2),ncopy,fzero,to,ncopy);
    }
#else
    if (is_symm)
    {
        cblas_ssymm(CblasRowMajor,CblasLeft,CblasUpper,N,N,scale,from1,N,from2,N,0.,to,N);
    }
    else
    {
        cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,N,scale,from1,N,from2,N,0.,to,N);
    }
#endif
}


template <>
inline
void
scale_matrix_mult(const double scale,
                  const double from1[],
                  const double from2[],
                  double to[],
                  const unsigned N,
                  const bool is_symm)
{
#ifdef FORTRAN_BLAS
    static double dzero(0.);
    static char nchar('N');
    static char lchar('L');
    static char uchar('U');
    int ncopy(N);
    if (is_symm)
    {
        FBLAS::dsymm(&lchar,&uchar,ncopy,ncopy,const_cast<double>(scale),const_cast<double*>(from1),
                     ncopy,const_cast<double*>(from2),ncopy,dzero,to,ncopy);
    }
    else
    {
        FBLAS::dgemm(&nchar,&nchar,ncopy,ncopy,ncopy,const_cast<double>(scale),const_cast<double*>(from1),
                     ncopy,const_cast<double*>(from2),ncopy,dzero,to,ncopy);
    }
#else
    if (is_symm)
    {
        cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,N,N,scale,from1,N,from2,N,0.,to,N);
    }
    else
    {
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,N,scale,from1,N,from2,N,0.,to,N);
    }
#endif
}


template <>
inline
void
scale_matrix_mult(const std::complex<double> scale,
                  const std::complex<double> from1[],
                  const std::complex<double> from2[],
                  std::complex<double> to[],
                  const unsigned N,
                  const bool is_symm)
{
#ifdef FORTRAN_BLAS
    static char nchar('N');
    static char lchar('L');
    static char uchar('U');
    int ncopy(N);
    if (is_symm)
    {
        FBLAS::zsymm(&lchar,&uchar,&ncopy,&ncopy,const_cast<void*>(sc),const_cast<void*>(f1),
                     &ncopy,const_cast<void*>(f2),&ncopy,const_cast<void*>(z),t,&ncopy);
    }
    else
    {
        FBLAS::zgemm(&nchar,&nchar,&ncopy,&ncopy,&ncopy,const_cast<void*>(sc),const_cast<void*>(f1),
                     &ncopy,const_cast<void*>(f2),&ncopy,const_cast<void*>(z),t,&ncopy);
    }
#else
    static const std::complex<double> zero(0.);
    const void* sc = reinterpret_cast<const void*>(&scale);
    const void* f1 = reinterpret_cast<const void*>(from1);
    const void* f2 = reinterpret_cast<const void*>(from2);
    const void* z = reinterpret_cast<const void*>(&zero);
    void* t = reinterpret_cast<void*>(to);

    if (is_symm)
    {
        cblas_zsymm(CblasRowMajor,CblasLeft,CblasUpper,N,N,sc,f1,N,f2,N,z,t,N);
    }
    else
    {
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,N,sc,f1,N,f2,N,z,t,N);
    }
#endif
}
#endif


// perform Blas 2 _gemv with alpha,beta=1.
//
template <typename FloatType>
inline
void
matrix_vector_mult_sum(const FloatType A[], // MxN
                       const FloatType x[],  // N
                       FloatType y[],        // M
                       const unsigned M,
                       const unsigned N)
{
    for (unsigned i(0); i<N; ++i)
    {
        array_scale_sum(A+i*M,x[i],y,M);
    }
}

#ifdef USE_BLAS
template <>
inline
void
matrix_vector_mult_sum(const float A[],
                       const float x[],
                       float y[],
                       const unsigned M,
                       const unsigned N)
{
#ifdef FORTRAN_BLAS
    static char nchar('N');
    static float fone(1.);
    static int one(1);
    int mcopy(M);
    int ncopy(N);
    FBLAS::sgemv(&nchar,mcopy,ncopy,fone,const_cast<float*>(A),mcopy,const_cast<float*>(x),one,fone,y,one);
#else
    cblas_sgemv(CblasColMajor,CblasNoTrans,M,N,1.,A,M,x,1,1.,y,1);
#endif
}

template <>
inline
void
matrix_vector_mult_sum(const double A[],
                       const double x[],
                       double y[],
                       const unsigned M,
                       const unsigned N)
{
#ifdef FORTRAN_BLAS
    static char nchar('N');
    static double done(1.);
    static int one(1);
    int mcopy(M);
    int ncopy(N);
    FBLAS::dgemv(&nchar,mcopy,ncopy,done,const_cast<double*>(A),mcopy,const_cast<double*>(x),one,done,y,one);
#else
    cblas_dgemv(CblasColMajor,CblasNoTrans,M,N,1.,A,M,x,1,1.,y,1);
#endif
}
#endif


// perform Blas 2 _gemv with alpha,beta=1.
//
template <typename FloatType>
inline
void
transposed_matrix_vector_mult_sum(const FloatType A[], // MxN
                                  const FloatType x[],  // N
                                  FloatType y[],        // M
                                  const unsigned M,
                                  const unsigned N)
{
    for (unsigned i(0); i<M; ++i)
    {
        y[i] += array_dot(A+i*N,x,N);
    }
}

#ifdef USE_BLAS
template <>
inline
void
transposed_matrix_vector_mult_sum(const float A[],
                                  const float x[],
                                  float y[],
                                  const unsigned M,
                                  const unsigned N)
{
#ifdef FORTRAN_BLAS
    die("undefined");
#else
    cblas_sgemv(CblasColMajor,CblasTrans,N,M,1.,A,N,x,1,1.,y,1);
#endif
}

template <>
inline
void
transposed_matrix_vector_mult_sum(const double A[],
                                  const double x[],
                                  double y[],
                                  const unsigned M,
                                  const unsigned N)
{
#ifdef FORTRAN_BLAS
    die("undefined");
#else
    cblas_dgemv(CblasColMajor,CblasTrans,N,M,1.,A,N,x,1,1.,y,1);
#endif
}
#endif


template <typename FloatType>
inline
void
matrix_square(const FloatType from[],
              FloatType to[],
              const unsigned N)
{
    matrix_mult(from,from,to,N);
}


template <unsigned N,typename FloatType>
inline
void
matrix_mult_inplace(const FloatType from[],
                    FloatType to[])
{

    FloatType to_copy[N*N];

    matrix_copy<N,FloatType>(to,to_copy,N);
    matrix_mult<N,FloatType>(from,to_copy,to,N);
}


template <typename FloatType>
inline
void
matrix_net_asym(FloatType asym[],
                const FloatType mat[],
                const unsigned N)
{
    static const FloatType zero(static_cast<FloatType>(0));
    for (unsigned i(0); i<N; ++i)
    {
        for (unsigned j(0); j<N; ++j)
        {
            asym[j+i*N] = mat[j+i*N]-mat[i+j*N];
            if (asym[j+i*N]<zero) asym[j+i*N] = zero;
        }
    }
}


template <typename FloatType>
inline
void
matrix_relative_asym(FloatType asym[],
                     const FloatType mat[],
                     const unsigned N)
{
    static const FloatType zero(static_cast<FloatType>(0));
    for (unsigned i(0); i<N; ++i)
    {
        for (unsigned j(0); j<N; ++j)
        {
            FloatType mdiff = mat[j+i*N]-mat[i+j*N];
            FloatType msum = mat[j+i*N]+mat[i+j*N];
            if ( mdiff < zero || msum == zero)
            {
                asym[j+i*N] = zero;
            }
            else
            {
                asym[j+i*N] = mdiff/msum;
            }
        }
    }
}


template <typename FloatType>
void
matrix_l2lr_asym(FloatType asym[],
                 const FloatType mat[],
                 const unsigned N)
{
    static const FloatType zero(static_cast<FloatType>(0));
    for (unsigned i(0); i<N; ++i)
    {
        for (unsigned j(0); j<N; ++j)
        {
            if ( i == j ||
                 std::fabs(mat[i+j*N]) < std::numeric_limits<FloatType>::epsilon() ||
                 mat[j+i*N]<mat[i+j*N] )
            {
                asym[j+i*N] = zero;
            }
            else
            {
                asym[j+i*N] = std::log(mat[j+i*N]/mat[i+j*N])/M_LN2;
            }
        }
    }
}




// x(i,j) i,j in N in the full matrix are averaged/summed to
// y(k,l) k,l in M given a mapping from N to M
//
template <typename FloatType>
void
matrix_state_reduction(FloatType* reduced_mat,
                       const unsigned reduced_size,
                       const FloatType* full_mat,
                       const unsigned full_size,
                       const unsigned* reduction_map,
                       const bool is_average = false,
                       const bool is_exclude_diag = true); // default is sum


#include "matrix_util.hh"

#endif
