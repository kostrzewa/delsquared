#ifndef _GAUGE_H
#define _GAUGE_H

#include "su3matrix.h"

template < typename T >
class gauge {
public:
  su3matrix<T> v[4];
};

#endif // _GAUGE_H

