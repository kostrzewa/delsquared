#ifndef _HALFSPINOR_H
#define _HALFSPINOR_H

#include "su3vector.h"

template< typename T >
class halfspinor {
public:
  su3vector<T> v[2];
};

typedef enum HALFSPINV_COMPONENT {
  s0 = 0, s1
} HALFSPINV_COMPONENT;


#endif // _HALFSPINOR_H
