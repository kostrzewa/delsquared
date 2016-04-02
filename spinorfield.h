#ifndef _SPINORFIELD_H
#define _SPINORFIELD_H

#include "spinor.h"

template< typename T >
class spinorfield {
public:
  spinor<T>* mem;
};

#endif // _SPINORFIELD_H
