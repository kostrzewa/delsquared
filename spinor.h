#ifndef _SPINOR_H
#define _SPINOR_H

#include "su3vector.h"

template< typename T >
class spinor {
public:
  su3vector<T> v[4];
};

typedef enum SPINV_COMPONENT {
  s0 = 0, s1, s2, s3 
} SPINV_COMPONENT;


#endif // _SPINOR_H
