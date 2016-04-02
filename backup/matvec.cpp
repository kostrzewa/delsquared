#include <complex>
#include <assert.h>
#include "matvec.h"

// NOTE: this can probably be templated even more by using an abstract
// vector type 

template <typename T>
void c_Ax(su3vector<T>& Y, const su3matrix<T>& A, const su3vector<T>& X){
  assert( ((Y.N==A.N) && (A.N==X.N)) );
  // T tmp[3]; 
  /* once vectorization is really used, tmp[i] would be SIMD
   * types such that all operations below make use of the full
   * SIMD width */
  T* __restrict__ m1=A.mT[0].v;
  T* __restrict__ m2=A.mT[1].v;
  T* __restrict__ m3=A.mT[2].v;
  T* __restrict__ x=X.v.v;
  T* __restrict__ y=Y.v.v;

  for(int i = 0; i < A.N; ++i){
    //tmp[0] = A.mT[0][i] * x[i];
    //tmp[1] = A.mT[1][i] * x[i];
    //tmp[3] = A.mT[2][i] * x[i];
    // we multiply with the transpose of the matrix such that the colums are 
    // contiguous in memory
    y[i] = ( m1[i] + m2[i] + m3[i] ) * x[i];
  } 
}

template void c_Ax<double>(su3vector<double>&,const su3matrix<double>&,const su3vector<double>&);
template void c_Ax<float>(su3vector<float>&,const su3matrix<float>&,const su3vector<float>&);

