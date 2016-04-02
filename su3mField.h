#ifndef _SU3MFIELD_H
#define _SU3MFIELD_H

#include <stddef.h>
#include <vector>
#include "su3m.h"

namespace delsquared {

template < class T >
class su3mField {
public:
  su3mField( const size_t vol, const size_t halo[4], const unsigned int vlength );
  ~su3mField();

  std::vector< su3m<T> > field;
  T* rawmem;

};

}

#endif
