#ifndef _BASICMATH_H
#define _BASICMATH_H

#define __SSE_4_2__
#include <immintrin.h>

#include <assert.h>

#include "dsReal.h"
#include "dsComplex.h"
//#include "spinor.h"
#include "su3v.h"
#include "su3m.h"

namespace delsquared {

template <typename Tout, typename Tin1, typename Tin2>
inline void cXc(dsComplex<Tout> &cout, const dsComplex<Tin1> & cin1, const dsComplex<Tin2> & cin2, const unsigned int vlength){
  for(unsigned int i = 0; i < vlength; ++i){
    cout.r[i] = cin1.r(i) * cin2.r(i) - cin1.i(i) * cin2.i(i);
    cout.i[i] = cin1.r(i) * cin2.i(i) + cin1.i(i) * cin2.r(i);
  }
}

template <typename Tout, typename Tin1, typename Tin2>
inline void rXc(dsComplex<Tout> &cout, dsReal<Tin1> & rin1, dsComplex<Tin2> & cin2, const unsigned int vlength){
  for(unsigned int i = 0; i < vlength; ++i){
    cout.r[i] = rin1.r(i) * cin2.r(i);
    cout.i[i] = rin1.r(i) * cin2.i(i);
  }
}

#define i_cXc( out, in1, in2 )\
  (out).r[i] = (in1).r(i) * (in2).r(i) - (in1).i(i) * (in2).i(i);\
  (out).i[i] = (in1).r(i) * (in2).i(i) + (in1).i(i) * (in2).r(i);

#define cassign( out, in )\
  (out).r[i] = (in).r(i);\
  (out).i[i] = (in).i(i);

#define gather_3v( out, in0, in1, in2 )\
  (out).r[i] = (in0).r(i) + (in1).r(i) + (in2).r(i);\
  (out).i[i] = (in0).i(i) + (in1).i(i) + (in2).i(i);


template <typename Tout, typename Tin1, typename Tin2>
inline void su3mXsu3v(su3v<Tout>& vout, const su3m<Tin1>& min, const su3v<Tin2>& vin, const unsigned int vlength, su3v<Tout>* const t ){
  for(unsigned int i = 0; i < vlength; ++i){
    i_cXc( t[0].c0 , min.c00 , vin.c0 );
    i_cXc( t[1].c0 , min.c10 , vin.c0 );
    i_cXc( t[2].c0 , min.c20 , vin.c0 );
    
    i_cXc( t[0].c1 , min.c01 , vin.c1 );
    i_cXc( t[1].c1 , min.c11 , vin.c1 );
    i_cXc( t[2].c1 , min.c21 , vin.c1 );
    
    i_cXc( t[0].c2 , min.c02 , vin.c2 );
    i_cXc( t[1].c2 , min.c12 , vin.c2 );
    i_cXc( t[2].c2 , min.c22 , vin.c2 );
  }
  for( unsigned int i = 0; i < vlength; ++i){
    gather_3v( vout.c0, t[0].c0, t[0].c1, t[0].c2 );
    gather_3v( vout.c1, t[1].c0, t[1].c1, t[1].c2 );
    gather_3v( vout.c2, t[2].c0, t[2].c1, t[2].c2 );    
  }
}

//template <typename Tout, typename Tin1, typename Tin2>
//inline void su3mXsu3v(su3vector<Tout> &vout, su3matrix &min<Tin1>, su3vector &vin<Tin2>){
//  assert( min.m[c00].r[0].l == vin.m[c0].r[0].l );


//template void ctc(dsComplex<__m128> &cout, dsComplex<__m128> &cin1, dsComplex<__m128> &cin2);
//{
//  assert( cin2.l == cin2 );
//  for(int i = 0; i < cin1.l; ++i){
//    cout.r[i] = __mm_sub_ps( _mm_mul_ps( cin1.r[i], cin2.r[i] ), _mm_mul_ps( cin1.i[i], cin2.i[i] ) );
//    cout.i[i] = __mm_add_ps( _mm_mul_ps( cin1.r[i], cin2.i[i] ), _mm_mul_ps( cin1.i[i], cin2.r[i] ) );
//  }
}

#endif // _BASICMATH_H



