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

// $Id: array_util.h 906 2007-10-10 23:30:33Z ctsa $

/// \file
///
/// \brief simple vector/array ops
///


#ifndef __ARRAY_UTIL_H
#define __ARRAY_UTIL_H

#ifdef USE_BLAS
extern "C" {
#include <cblas.h>
}
#ifdef XLC_HACK
namespace ESSL
{
extern "C" {
#define sdot      esvsdot
#define ddot      esvddot
    float  esvsdot(int,  float*, int,  float*, int);
    double esvddot(int, double*, int, double*, int);
}
}
#endif
#endif

#include <cmath>

#include <algorithm>
#include <limits>


/// adjust expectation of array a to target_expect, where each element
/// in a has probability given by a_prob, returns final scale factor
///
template <typename ForwardIterator,
         typename InputIterator,
         typename T>
T
array_scale_expectation(ForwardIterator a,
                        const ForwardIterator a_end,
                        InputIterator a_prob,
                        const T target_expect);


/// adjust expectation of array a to target_expect, where each element
/// in a has probability given by a_prob. adjust expectation by only
/// changing elements corresponding to true in a_selected. if no
/// elements selected, adjust all elements, returns final scale factor
///
template <typename ForwardIterator1,
         typename InputIterator,
         typename ForwardIterator2,
         typename T>
T
array_scale_expectation_selected_only(ForwardIterator1 a,
                                      const ForwardIterator1 a_end,
                                      InputIterator a_prob,
                                      ForwardIterator2 a_selected,
                                      const T target_expect);


template <typename FloatType>
inline
FloatType
array_abs_max(const FloatType v[],
              const unsigned N)
{
    if (N==0) return std::numeric_limits<FloatType>::min();
    FloatType m = std::abs(*v);
    for (unsigned i(1); i<N; ++i) if (m<std::abs(v[i])) m=std::abs(v[i]);
    return m;
}


template <typename FloatType>
inline
FloatType
array_dot(const FloatType v1[],
          const FloatType v2[],
          const unsigned N)
{
    FloatType d = 0;
    for (unsigned i(0); i<N; ++i)
    {
        d += v1[i]*v2[i];
    }
    return d;
}


#if USE_BLAS
template <>
inline
float
array_dot(const float v1[],
          const float v2[],
          const unsigned N)
{
#ifndef XLC_HACK
    return cblas_sdot(N,v1,1,v2,1);
#else
    int ncopy(N);
    return ESSL::sdot(ncopy,const_cast<float*>(v1),1,const_cast<float*>(v2),1);
#endif
}

template <>
inline
double
array_dot(const double v1[],
          const double v2[],
          const unsigned N)
{
#ifndef XLC_HACK
    return cblas_ddot(N,v1,1,v2,1);
#else
    int ncopy(N);
    return ESSL::ddot(ncopy,const_cast<double*>(v1),1,const_cast<double*>(v2),1);
#endif

}
#endif


// get zero N vector
template <typename FloatType>
inline
void
array_zero(FloatType v[],
           const unsigned N)
{
    for (unsigned i(0); i<N; ++i) v[i] = FloatType(0);
}

template <typename FloatType>
inline
FloatType
array_sum(const FloatType* v,
          const unsigned N)
{

    FloatType s(0);
    for (unsigned i(0); i<N; ++i) s += v[i];
    return s;
}

template <typename FloatType1,
         typename FloatType2>
inline
void
array_scale(FloatType1* v,
            const unsigned N,
            const FloatType2 s)
{
    for (unsigned i(0); i<N; ++i) v[i] = v[i] * s;
}

template <typename FloatType>
inline
void
array_scale_plus_vector(const FloatType* v1,
                        const FloatType s,
                        FloatType* v2,
                        const unsigned N)
{
    for (unsigned i(0); i<N; ++i) v2[i] += s*v1[i];
}
#if USE_BLAS
template <>
inline
void
array_scale_plus_vector(const float* v1,
                        const float s,
                        float* v2,
                        const unsigned N)
{
#ifdef FORTRAN_BLAS
    static int one(1);
    int ncopy(N);
    return FBLAS::saxpy(ncopy,const_cast<float>(s),const_cast<float*>(v1),one,const_cast<float*>(v2),one);
#else
    return cblas_saxpy(N,s,v1,1,v2,1);
#endif
}

template <>
inline
void
array_scale_plus_vector(const double* v1,
                        const double s,
                        double* v2,
                        const unsigned N)
{
#ifdef FORTRAN_BLAS
    static int one(1);
    int ncopy(N);
    return FBLAS::daxpy(ncopy,const_cast<double>(s),const_cast<double*>(v1),one,const_cast<double*>(v2),one);
#else
    return cblas_daxpy(N,s,v1,1,v2,1);
#endif
}
#endif


template <typename FloatType>
inline
FloatType
array_norm2(const FloatType* v,
            const unsigned N)
{
    FloatType sum(0);
    for (unsigned i(0); i<N; ++i) sum += v[i]*v[i];
    return std::sqrt(sum);
}

#if USE_BLAS
template <>
inline
float
array_norm2(const float* v,
            const unsigned N)
{
#ifdef FORTRAN_BLAS
    static int one(1);
    int ncopy(N);
    return FBLAS::snrm2(ncopy,const_cast<float*>(v),one);
#else
    return cblas_snrm2(N,v,1);
#endif
}

template <>
inline
double
array_norm2(const double* v,
            const unsigned N)
{
#ifdef FORTRAN_BLAS
    static int one(1);
    int ncopy(N);
    return FBLAS::dnrm2(ncopy,const_cast<double*>(v),one);
#else
    return cblas_dnrm2(N,v,1);
#endif
}
#endif


template <typename FloatType>
inline
FloatType
normalized_dot(const FloatType v1[],
               const FloatType v2[],
               const unsigned N)
{
    return array_dot(v1,v2,N)/(array_norm2(v1,N)*array_norm2(v2,N));
}


template <typename FloatType>
inline
void
array_log(FloatType* v,
          const unsigned N)
{
    for (unsigned i(0); i<N; ++i) v[i] = std::log(v[i]);
}

template <typename FloatType>
inline
void
array_log_lowres(FloatType* v,
                 const unsigned N)
{
    for (unsigned i(0); i<N; ++i) v[i] = std::log(static_cast<float>(v[i]));
}


#ifdef USE_VML
extern "C" {
    void vsLn(int,const float*,float*);
    void vdLn(int,const double*,double*);
}

template <>
inline
void
array_log(float* v,
          const unsigned N)
{
    vsLn(N,v,v);
}

template <>
inline
void
array_log(double* v,
          const unsigned N)
{
    vdLn(N,v,v);
}
#endif


template <typename FloatType>
inline
void
array_sqrt(FloatType* v,
           const unsigned N)
{
    for (unsigned i(0); i<N; ++i) v[i] = std::sqrt(v[i]);
}

#ifdef USE_VML
extern "C" {
    void vsSqrt(int,const float*,float*);
    void vdSqrt(int,const double*,double*);
}

template <>
inline
void
array_sqrt(float* v,
           const unsigned N)
{
    vsSqrt(N,v,v);
}

template <>
inline
void
array_sqrt(double* v,
           const unsigned N)
{
    vdSqrt(N,v,v);
}
#endif


// lifts all values to at least minval
template <typename FloatType>
inline
void
array_set_min(FloatType* v,
              const FloatType minval,
              const unsigned N)
{
    for (unsigned i(0); i<N; ++i) v[i] = std::max(v[i],minval);
}

template <typename FloatType>
inline
void
array_scale_sum(const FloatType from[],
                const FloatType s,
                FloatType to[],
                const unsigned N)
{
    for (unsigned i(0); i<N; ++i) to[i] += from[i]*s;
}

#ifdef USE_BLAS
template <>
inline
void
array_scale_sum(const float from[],
                const float s,
                float to[],
                const unsigned N)
{
#ifdef FORTRAN_BLAS
    static int one(1);
    int ncopy(N);
    FBLAS::saxpy(ncopy,const_cast<float>(s),const_cast<float*>(from),one,to,one);
#else
    cblas_saxpy(N,s,from,1,to,1);
#endif
}


template <>
inline
void
array_scale_sum(const double from[],
                const double s,
                double to[],
                const unsigned N)
{
#ifdef FORTRAN_BLAS
    static int one(1);
    int ncopy(N);
    FBLAS::daxpy(ncopy,const_cast<double>(s),const_cast<double*>(from),one,to,one);
#else
    cblas_daxpy(N,s,from,1,to,1);
#endif
}
#endif

#include "array_util.hh"

#endif
