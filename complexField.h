#ifndef _COMPLEXFIELD_H
#define _COMPLEXFIELD_H

#include <stddef.h>
#include <vector>
#include "dsComplex.h"

namespace delsquared {

template < class T >
class complexField {
public:
  complexField( const size_t vol, const size_t halo[4], const unsigned int vlength );
  ~complexField();

  std::vector< dsComplex<T> > field;
  T* rawmem;

};

}

#endif // _COMPLEXFIELD_H
