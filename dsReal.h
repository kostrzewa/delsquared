#ifndef _DSREAL_H
#define _DSREAL_H

#include "vectorable.h"

namespace delsquared {

template < typename T >
class dsReal {
public:
  vectorable<T> r;
};

}

#endif // _DSREAL_H
