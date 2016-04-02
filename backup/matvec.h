#ifndef _MATVEC_H
#define _MATVEC_H

#include "../su3vector.h"
#include "../su3matrix.h"

template<typename T>
void Ax(su3vector<T>& y, su3matrix<T>& A, su3vector<T>& x);

#endif
