#ifndef _GAUGEFIELD_H
#define _GAUGEFIELD_H

#include "gauge.h"

template < typename T >
class gaugefield {

public:
  gauge<T>* ptrmem;
  T*        rawmem;
};

#endif // _GAUGEFIELD_H

