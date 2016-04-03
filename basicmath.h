#ifndef _BASICMATH_H
#define _BASICMATH_H

#include <assert.h>

#include "dsReal.h"
#include "dsComplex.h"
//#include "spinor.h"
#include "su3v.h"
#include "su3m.h"

#include "immintrin.h"

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

#define i_cXc_direct( out, in1, in2 )\
  (out).r.v[i]  = (in1).r.v[i] * (in2).r.v[i];\
  (out).r.v[i] -= (in1).i.v[i] * (in2).i.v[i];\
  (out).i.v[i]  = (in1).r.v[i] * (in2).i.v[i];\
  (out).i.v[i] += (in1).i.v[i] * (in2).r.v[i];

#define gather_3v_direct( out, in0, in1, in2 )\
  (out).r.v[i]  = (in0).r.v[i] + (in1).r.v[i];\
  (out).r.v[i] += (in2).r.v[i];\
  (out).i.v[i]  = (in0).i.v[i] + (in1).i.v[i];\
  (out).i.v[i] += (in2).i.v[i];

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
//  }
//  for( unsigned int i = 0; i < vlength; ++i){
    gather_3v( vout.c0, t[0].c0, t[0].c1, t[0].c2 );
    gather_3v( vout.c1, t[1].c0, t[1].c1, t[1].c2 );
    gather_3v( vout.c2, t[2].c0, t[2].c1, t[2].c2 );    
  }
}

//template <typename Tout, typename Tin1, typename Tin2>
inline void su3mXsu3v_intrin_float(su3v<float>& vout, const su3m<float>& min, const su3v<float>& vin, const unsigned int vlength){
  __m128 mc0r, mc0i, mc1r, mc1i, mc2r, mc2i;
  __m128 vc0r, vc0i, vc1r, vc1i, vc2r, vc2i;
  __m128 tc0r, tc0i, tc1r, tc1i, tc2r, tc2i;
  __m128 ttc0r, ttc0i, ttc1r, ttc1i, ttc2r, ttc2i;
  __m128 tttc0r, tttc0i, tttc1r, tttc1i, tttc2r, tttc2i;
  for(unsigned int i = 0; (i+3) < vlength; i+=4){
    mc0r = _mm_load_ps( &(min.c00.r(i)) );
    vc0r = _mm_load_ps( &(vin.c0.r(i)) );
    tc0r = _mm_mul_ps( mc0r, vc0r );
    
    mc1r = _mm_load_ps( &(min.c01.r(i)) );
    vc1r = _mm_load_ps( &(vin.c1.r(i)) );
    tc1r = _mm_mul_ps( mc1r, vc1r );
    
    mc2r = _mm_load_ps( &(min.c02.r(i)) );
    vc2r = _mm_load_ps( &(vin.c2.r(i)) );
    tc2r = _mm_mul_ps( mc2r, vc2r );
    
    mc0i = _mm_load_ps( min.c00.i.v+i );
    tc0i = _mm_mul_ps( mc0i, vc0r );

    mc1i = _mm_load_ps( &(min.c01.i(i)) );
    tc1i = _mm_mul_ps( mc1i, vc1r );

    mc2i = _mm_load_ps( &(min.c02.i(i)) );
    tc2i = _mm_mul_ps( mc2i, vc2r );

    vc0i = _mm_load_ps( &(vin.c0.i(i)) );
    ttc0r = _mm_mul_ps( mc0i, vc0i );
    vc1i = _mm_load_ps( &(vin.c1.i(i)) );
    tttc0r = _mm_sub_ps( tc0r, ttc0r );
    ttc1r = _mm_mul_ps( mc1i, vc1i );
    vc2i = _mm_load_ps( &(vin.c2.i(i)) );
    tttc1r = _mm_sub_ps( tc1r, ttc1r );
    ttc2r = _mm_mul_ps( mc2i, vc2i );
    tttc2r = _mm_sub_ps( tc2r, ttc2r );

    tc1r = _mm_add_ps( tttc0r, tttc1r );
    tc0r = _mm_add_ps( tc1r, tttc2r );
    _mm_store_ps( &(vout.c0.r[i]), tc0r );

    ttc0i = _mm_mul_ps( mc0r, vc0i );
    tttc0i = _mm_add_ps( tc0i, ttc0i );

    ttc1i = _mm_mul_ps( mc1r, vc1i );
    tttc1i = _mm_add_ps( tc1i, ttc1i );
    tc1i = _mm_add_ps( tttc1i, tttc0i );

    ttc2i = _mm_mul_ps( mc2r, vc2i );
    tttc2i = _mm_add_ps( tc2i, ttc2i );
    tc0i = _mm_add_ps( tc1i, tttc2i );
    _mm_store_ps( &(vout.c0.i[i]), tc0i );

    /////////////////////////////////////

    mc0r = _mm_load_ps( &(min.c10.r(i)) );
    tc0r = _mm_mul_ps( mc0r, vc0r );
    
    mc1r = _mm_load_ps( &(min.c11.r(i)) );
    tc1r = _mm_mul_ps( mc1r, vc1r );
    
    mc2r = _mm_load_ps( &(min.c12.r(i)) );
    tc2r = _mm_mul_ps( mc2r, vc2r );
    
    mc0i = _mm_load_ps( &(min.c10.i(i)) );
    tc0i = _mm_mul_ps( mc0i, vc0r );

    mc1i = _mm_load_ps( &(min.c11.i(i)) );
    tc1i = _mm_mul_ps( mc1i, vc1r );

    mc2i = _mm_load_ps( &(min.c12.i(i)) );
    tc2i = _mm_mul_ps( mc2i, vc2r );

    ttc0r = _mm_mul_ps( mc0i, vc0i );
    tttc0r = _mm_sub_ps( tc0r, ttc0r );
    ttc1r = _mm_mul_ps( mc1i, vc1i );
    tttc1r = _mm_sub_ps( tc1r, ttc1r );
    ttc2r = _mm_mul_ps( mc2i, vc2i );
    tttc2r = _mm_sub_ps( tc2r, ttc2r );

    tc1r = _mm_add_ps( tttc0r, tttc1r );
    tc0r = _mm_add_ps( tc1r, tttc2r );
    _mm_store_ps( &(vout.c1.r[i]), tc0r );

    ttc0i = _mm_mul_ps( mc0r, vc0i );
    tttc0i = _mm_add_ps( tc0i, ttc0i );

    ttc1i = _mm_mul_ps( mc1r, vc1i );
    tttc1i = _mm_add_ps( tc1i, ttc1i );
    tc1i = _mm_add_ps( tttc1i, tttc0i );

    ttc2i = _mm_mul_ps( mc2r, vc2i );
    tttc2i = _mm_add_ps( tc2i, ttc2i );
    tc0i = _mm_add_ps( tc1i, tttc2i );
    _mm_store_ps( &(vout.c1.i[i]), tc0i );

    //////////////////////////////////////////////////////////////////
    //TODO: intersperse with load instructions for next iteration...//
    //////////////////////////////////////////////////////////////////

    mc0r = _mm_load_ps( &(min.c20.r(i)) );
    tc0r = _mm_mul_ps( mc0r, vc0r );
    
    mc1r = _mm_load_ps( &(min.c21.r(i)) );
    tc1r = _mm_mul_ps( mc1r, vc1r );
    
    mc2r = _mm_load_ps( &(min.c22.r(i)) );
    tc2r = _mm_mul_ps( mc2r, vc2r );
    
    mc0i = _mm_load_ps( &(min.c20.i(i)) );
    tc0i = _mm_mul_ps( mc0i, vc0r );

    mc1i = _mm_load_ps( &(min.c21.i(i)) );
    tc1i = _mm_mul_ps( mc1i, vc1r );

    mc2i = _mm_load_ps( &(min.c22.i(i)) );
    tc2i = _mm_mul_ps( mc2i, vc2r );

    ttc0r = _mm_mul_ps( mc0i, vc0i );
    tttc0r = _mm_sub_ps( tc0r, ttc0r );
    ttc1r = _mm_mul_ps( mc1i, vc1i );
    tttc1r = _mm_sub_ps( tc1r, ttc1r );
    ttc2r = _mm_mul_ps( mc2i, vc2i );
    tttc2r = _mm_sub_ps( tc2r, ttc2r );

    tc1r = _mm_add_ps( tttc0r, tttc1r );
    tc0r = _mm_add_ps( tc1r, tttc2r );
    _mm_store_ps( &(vout.c2.r[i]), tc0r );

    ttc0i = _mm_mul_ps( mc0r, vc0i );
    tttc0i = _mm_add_ps( tc0i, ttc0i );

    ttc1i = _mm_mul_ps( mc1r, vc1i );
    tttc1i = _mm_add_ps( tc1i, ttc1i );
    tc1i = _mm_add_ps( tttc1i, tttc0i );

    ttc2i = _mm_mul_ps( mc2r, vc2i );
    tttc2i = _mm_add_ps( tc2i, ttc2i );
    tc0i = _mm_add_ps( tc1i, tttc2i );
    _mm_store_ps( &(vout.c2.i[i]), tc0i );
  }
}

template <typename Tout, typename Tin1, typename Tin2>
inline void su3mXsu3v_direct(su3v<Tout>& vout, const su3m<Tin1>& min, const su3v<Tin2>& vin, const unsigned int vlength, su3v<Tout>* const t ){
  for(unsigned int i = 0; i < vlength; ++i){
    i_cXc_direct( t[0].c0 , min.c00 , vin.c0 );
    i_cXc_direct( t[1].c0 , min.c10 , vin.c0 );
    i_cXc_direct( t[2].c0 , min.c20 , vin.c0 );
    
    i_cXc_direct( t[0].c1 , min.c01 , vin.c1 );
    i_cXc_direct( t[1].c1 , min.c11 , vin.c1 );
    i_cXc_direct( t[2].c1 , min.c21 , vin.c1 );
    
    i_cXc_direct( t[0].c2 , min.c02 , vin.c2 );
    i_cXc_direct( t[1].c2 , min.c12 , vin.c2 );
    i_cXc_direct( t[2].c2 , min.c22 , vin.c2 );
//  }
//  for( unsigned int i = 0; i < vlength; ++i){
    gather_3v_direct( vout.c0, t[0].c0, t[0].c1, t[0].c2 );
    gather_3v_direct( vout.c1, t[1].c0, t[1].c1, t[1].c2 );
    gather_3v_direct( vout.c2, t[2].c0, t[2].c1, t[2].c2 );    
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



