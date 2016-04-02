#ifndef _SU3VFIELD_H
#define _SU3VFIELD_H

#include <stddef.h>
#include <vector>
#include "su3v.h"

namespace delsquared {

template < class T >
class su3vField {
public:
  su3vField( const size_t vol, const size_t halo[4], const unsigned int vlength );
  ~su3vField();

  std::vector< su3v<T> > field;
  T* rawmem;

};

}

#endif
