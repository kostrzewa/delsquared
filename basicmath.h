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

/* #define i_cXc( out, in1, in2 )\
  (out).r.v[i  ]  = (in1).r.v[i  ] * (in2).r.v[i  ];\
  (out).r.v[i+1]  = (in1).r.v[i+1] * (in2).r.v[i+1];\
  (out).r.v[i+2]  = (in1).r.v[i+2] * (in2).r.v[i+2];\
  (out).r.v[i+3]  = (in1).r.v[i+3] * (in2).r.v[i+3];\
  (out).r.v[i  ] -= (in1).i.v[i  ] * (in2).i.v[i  ];\
  (out).r.v[i+1] -= (in1).i.v[i+1] * (in2).i.v[i+1];\
  (out).r.v[i+2] -= (in1).i.v[i+2] * (in2).i.v[i+2];\
  (out).r.v[i+3] -= (in1).i.v[i+3] * (in2).i.v[i+2];\
  (out).i.v[i  ]  = (in1).r.v[i  ] * (in2).i.v[i  ];\
  (out).i.v[i+1]  = (in1).r.v[i+1] * (in2).i.v[i+1];\
  (out).i.v[i+2]  = (in1).r.v[i+2] * (in2).i.v[i+2];\
  (out).i.v[i+3]  = (in1).r.v[i+3] * (in2).i.v[i+3];\
  (out).i.v[i]   += (in1).i.v[i  ] * (in2).r.v[i  ];\
  (out).i.v[i+1] += (in1).i.v[i+1] * (in2).r.v[i+1];\
  (out).i.v[i+2] += (in1).i.v[i+2] * (in2).r.v[i+2];\
  (out).i.v[i+3] += (in1).i.v[i+3] * (in2).r.v[i+3]; */

//#define i_cXc( out, in1, in2 )\
//  (out).r.v[i  ]  = (in1).r.v[i  ] * (in2).r.v[i  ];\
//  (out).r.v[i  ] -= (in1).i.v[i  ] * (in2).i.v[i  ];\
//  (out).i.v[i  ]  = (in1).r.v[i  ] * (in2).i.v[i  ];\
//  (out).i.v[i]   += (in1).i.v[i  ] * (in2).r.v[i  ];

#define i_cXc( out, in1, in2 )\
  (out).r[i  ]  = (in1).r(i  ) * (in2).r(i  );\
  (out).r[i  ] -= (in1).i(i  ) * (in2).i(i  );\
  (out).i[i  ]  = (in1).r(i  ) * (in2).i(i  );\
  (out).i[i]   += (in1).i(i  ) * (in2).r(i  );

//#define i_cXc( out, in1, in2 )\
//  (out).r.v[i  ]  = (in1).r.v[i  ] * (in2).r.v[i  ] \
//  - (in1).i.v[i  ] * (in2).i.v[i  ];\
//  (out).i.v[i  ]  = (in1).r.v[i  ] * (in2).i.v[i  ] \
//  + (in1).i.v[i  ] * (in2).r.v[i  ];

#define cassign( out, in )\
  (out).r[i] = (in).r(i);\
  (out).i[i] = (in).i(i);

/* #define gather_3v( out, in0, in1, in2 )\
  (out).r.v[i  ]  = (in0).r.v[i  ] + (in1).r.v[i  ];\
  (out).r.v[i+1]  = (in0).r.v[i+1] + (in1).r.v[i+1];\
  (out).r.v[i+2]  = (in0).r.v[i+2] + (in1).r.v[i+2];\
  (out).r.v[i+3]  = (in0).r.v[i+3] + (in1).r.v[i+3];\
  (out).r.v[i  ] += (in2).r.v[i  ];\
  (out).r.v[i+1] += (in2).r.v[i+1];\
  (out).r.v[i+2] += (in2).r.v[i+2];\
  (out).r.v[i+3] += (in2).r.v[i+3];\
  (out).i.v[i  ]  = (in0).i.v[i  ] + (in1).i.v[i  ];\
  (out).i.v[i+1]  = (in0).i.v[i+1] + (in1).i.v[i+1];\
  (out).i.v[i+2]  = (in0).i.v[i+2] + (in1).i.v[i+2];\
  (out).i.v[i+3]  = (in0).i.v[i+3] + (in1).i.v[i+3];\
  (out).i.v[i  ] += (in2).i.v[i  ];\
  (out).i.v[i+1] += (in2).i.v[i+1];\
  (out).i.v[i+2] += (in2).i.v[i+2];\
  (out).i.v[i+3] += (in2).i.v[i+3];\ */

//#define gather_3v( out, in0, in1, in2 )\
//  (out).r.v[i  ]  = (in0).r.v[i  ] + (in1).r.v[i  ];\
//  (out).r.v[i  ] += (in2).r.v[i  ];\
//  (out).i.v[i  ]  = (in0).i.v[i  ] + (in1).i.v[i  ];\
//  (out).i.v[i  ] += (in2).i.v[i  ];

#define gather_3v( out, in0, in1, in2 )\
  (out).r.v[i  ]  = (in0).r.v[i  ] + (in1).r.v[i  ];\
  (out).r.v[i  ] += (in2).r.v[i  ];\
  (out).i.v[i  ]  = (in0).i.v[i  ] + (in1).i.v[i  ];\
  (out).i.v[i  ] += (in2).i.v[i  ];

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

#define AL __attribute__ ((aligned (32)))

#define vmul( out, in1, in2 )\
  *((out)+0) = *((in1)+0) * *((in2)+0);\
  *((out)+1) = *((in1)+1) * *((in2)+1);\
  *((out)+2) = *((in1)+2) * *((in2)+2);\
  *((out)+3) = *((in1)+3) * *((in2)+3);

template <typename Tout, typename Tin1, typename Tin2>
inline void su3mXsu3v(su3v<Tout>& vout, const su3m<Tin1>& min, const su3v<Tin2>& vin, const unsigned int vlength, su3v<Tout>* const t ){
//  const unsigned int vincr = 4;
//  Tout t0r[vlength] AL, t0i[vlength] AL, t1r[vlength] AL, t1i[vlength] AL, t2r[vlength] AL, t2i[vlength] AL;
//
//  const Tin1* __restrict__ const mc00r AL = &(min.c00.r(0)); const Tin1* __restrict__ const mc00i AL = &(min.c00.i(0));    
//  const Tin1* __restrict__ const mc01r AL = &(min.c01.r(0)); const Tin1* __restrict__ const mc01i AL = &(min.c01.i(0));
//  const Tin1* __restrict__ const mc02r AL = &(min.c02.r(0)); const Tin1* __restrict__ const mc02i AL = &(min.c02.i(0));
//  const Tin1* __restrict__ const mc10r AL = &(min.c10.r(0)); const Tin1* __restrict__ const mc10i AL = &(min.c10.i(0));
//  const Tin1* __restrict__ const mc11r AL = &(min.c11.r(0)); const Tin1* __restrict__ const mc11i AL = &(min.c11.i(0));
//  const Tin1* __restrict__ const mc12r AL = &(min.c12.r(0)); const Tin1* __restrict__ const mc12i AL = &(min.c12.i(0));
//  const Tin1* __restrict__ const mc20r AL = &(min.c20.r(0)); const Tin1* __restrict__ const mc20i AL = &(min.c20.i(0));
//  const Tin1* __restrict__ const mc21r AL = &(min.c21.r(0)); const Tin1* __restrict__ const mc21i AL = &(min.c21.i(0));
//  const Tin1* __restrict__ const mc22r AL = &(min.c22.r(0)); const Tin1* __restrict__ const mc22i AL = &(min.c22.i(0));
//
//  const Tin2* __restrict__ const vc0r AL = &(vin.c0.r(0));   const Tin2* __restrict__ const vc0i AL = &(vin.c0.r(0));
//  const Tin2* __restrict__ const vc1r AL = &(vin.c1.r(0));   const Tin2* __restrict__ const vc1i AL = &(vin.c1.r(0));
//  const Tin2* __restrict__ const vc2r AL = &(vin.c2.r(0));   const Tin2* __restrict__ const vc2i AL = &(vin.c2.r(0));
//
//  Tout* __restrict__ const oc0r AL = &(vout.c0.r[0]); Tout* __restrict__ const oc0i AL = &(vout.c0.r[0]);
//  Tout* __restrict__ const oc1r AL = &(vout.c1.r[0]); Tout* __restrict__ const oc1i AL = &(vout.c1.r[0]);
//  Tout* __restrict__ const oc2r AL = &(vout.c2.r[0]); Tout* __restrict__ const oc2i AL = &(vout.c2.r[0]);
//
//  for(unsigned int i = 0; i < vlength; ++i){
//    *(oc0r+i) = *(mc00r+i) * *(vc0r+i);
//  } for(unsigned int i = 0; i < vlength; ++i){
//    *(oc1r+i) = *(mc01r+i) * *(vc1r+i);
//    *(oc2r+i) = *(mc02r+i) * *(vc2r+i);
//    *(oc0i+i) = *(mc00i+i) * *(vc0i+i);
//    *(oc1i+i) = *(mc01i+i) * *(vc1i+i);
//    *(oc2i+i) = *(mc02i+i) * *(vc2i+i);

//    vmul( oc0r+i, mc00r+i, vc0r+i );
//    vmul( oc1r+i, mc01r+i, vc1r+i );
//    vmul( oc2r+i, mc02r+i, vc2r+i );
//    vmul( oc0i+i, mc01i+i, vc0i+i );
//    vmul( oc1i+i, mc02i+i, vc1i+i );
//    vmul( oc2i+i, mc02i+i, vc2i+i );

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c0.r[i] = min.c00.r(i) * vin.c0.r(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c0.r[i] = min.c10.r(i) * vin.c0.r(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c0.r[i] = min.c20.r(i) * vin.c0.r(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c1.r[i] = min.c01.r(i) * vin.c1.r(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c1.r[i] = min.c11.r(i) * vin.c1.r(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c1.r[i] = min.c21.r(i) * vin.c1.r(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c2.r[i] = min.c02.r(i) * vin.c2.r(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c2.r[i] = min.c12.r(i) * vin.c2.r(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c2.r[i] = min.c22.r(i) * vin.c2.r(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c0.r[i] -= min.c00.i(i) * vin.c0.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c0.r[i] -= min.c10.i(i) * vin.c0.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c0.r[i] -= min.c20.i(i) * vin.c0.i(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c1.r[i] -= min.c01.i(i) * vin.c1.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c1.r[i] -= min.c11.i(i) * vin.c1.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c2.r[i] -= min.c21.i(i) * vin.c1.i(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c2.r[i] -= min.c02.i(i) * vin.c2.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c2.r[i] -= min.c12.i(i) * vin.c2.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c2.r[i] -= min.c22.i(i) * vin.c2.i(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c0.i[i] = min.c00.r(i) * vin.c0.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c0.i[i] = min.c10.r(i) * vin.c0.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c0.i[i] = min.c20.r(i) * vin.c0.i(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c1.i[i] = min.c01.r(i) * vin.c1.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c1.i[i] = min.c11.r(i) * vin.c1.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c1.i[i] = min.c21.r(i) * vin.c1.i(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c2.i[i] = min.c02.r(i) * vin.c2.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c2.i[i] = min.c12.r(i) * vin.c2.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c2.i[i] = min.c22.r(i) * vin.c2.i(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c0.i[i] += min.c00.r(i) * vin.c0.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c0.i[i] += min.c10.r(i) * vin.c0.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c0.i[i] += min.c20.r(i) * vin.c0.i(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c1.i[i] += min.c01.r(i) * vin.c1.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c1.i[i] += min.c11.r(i) * vin.c1.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c2.i[i] += min.c21.r(i) * vin.c1.i(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    t[0].c2.i[i] += min.c02.r(i) * vin.c2.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[1].c2.i[i] += min.c12.r(i) * vin.c2.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    t[2].c2.i[i] += min.c22.r(i) * vin.c2.i(i);
  }

  for(unsigned int i = 0; i < vlength; ++i){
    vout.c0.r[i] = t[0].c0.r(i) + t[0].c1.r(i) + t[0].c2.r(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    vout.c0.i[i] = t[0].c0.i(i) + t[0].c1.i(i) + t[0].c2.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    vout.c1.r[i] = t[1].c0.r(i) + t[1].c1.r(i) + t[1].c2.r(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    vout.c1.i[i] = t[1].c0.i(i) + t[1].c1.i(i) + t[1].c2.i(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    vout.c2.r[i] = t[2].c0.r(i) + t[2].c1.r(i) + t[2].c2.r(i);
  }for(unsigned int i = 0; i < vlength; ++i){
    vout.c2.i[i] = t[2].c0.i(i) + t[2].c1.i(i) + t[2].c2.i(i);
  }

//  for(unsigned int i = 0; i < vlength; ++i){
//    i_cXc( t[0].c0 , min.c00 , vin.c0 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    i_cXc( t[1].c0 , min.c10 , vin.c0 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    i_cXc( t[2].c0 , min.c20 , vin.c0 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    i_cXc( t[0].c1 , min.c01 , vin.c1 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    i_cXc( t[1].c1 , min.c11 , vin.c1 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    i_cXc( t[2].c1 , min.c21 , vin.c1 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    i_cXc( t[0].c2 , min.c02 , vin.c2 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    i_cXc( t[1].c2 , min.c12 , vin.c2 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    i_cXc( t[2].c2 , min.c22 , vin.c2 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    gather_3v( vout.c0, t[0].c0, t[0].c1, t[0].c2 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    gather_3v( vout.c1, t[1].c0, t[1].c1, t[1].c2 );
//  }for(unsigned int i = 0; i < vlength; ++i){
//    gather_3v( vout.c2, t[2].c0, t[2].c1, t[2].c2 );    
//  }
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



