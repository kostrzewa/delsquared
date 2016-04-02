#ifndef _SU3M_H
#define _SU3M_H

#include "dsComplex.h"

namespace delsquared {

template < typename T >
class su3m {
public:
  su3m() {}
  su3m( const su3m<T>& other ) {}
  su3m( const unsigned vlength ) : c00(vlength), c01(vlength), c02(vlength),
                                   c10(vlength), c11(vlength), c12(vlength),
                                   c20(vlength), c21(vlength), c22(vlength) {}
  dsComplex<T> c00, c01, c02;
  dsComplex<T> c10, c11, c12;
  dsComplex<T> c20, c21, c22;

  //void assign(su3matrix<T>&);
  //void assign_transpose(su3matrix<T>&);
  //int N;
};

//typedef enum SU3M_COMPONENT {
//  c00 = 0, 
//       c01, c02,
//  c10, c11, c12,
//  c20, c21, c22
//} SU3M_COMPONENT;

}

#endif
