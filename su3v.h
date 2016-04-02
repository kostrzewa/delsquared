#ifndef _SU3V_H
#define _SU3V_H

#include <vector>
#include "dsComplex.h"

namespace delsquared {

template < class T >
class su3v {
public:
  su3v() {}
  su3v( const su3v<T>& other ) {}
  su3v(unsigned int vlength) : c0(vlength), c1(vlength), c2(vlength) {}
  
  dsComplex<T> c0;
  dsComplex<T> c1;
  dsComplex<T> c2;
};

}

//typedef enum SU3V_COMPONENT {
//  c0 = 0, c1, c2
//} SU3V_COMPONENT;



#endif
