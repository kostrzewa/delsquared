#ifndef _DSCOMPLEX_H
#define _DSCOMPLEX_H

#include "vectorable.h"

namespace delsquared {

template < class T >
class dsComplex {
public:
  dsComplex() {}
  dsComplex( const dsComplex<T>& other ) {}
  dsComplex(unsigned int vlength) : r(vlength), i(vlength) {}

  vectorable<T> r,i;
};

}

#endif // _DSCOMPLEX_H
