#include <stdlib.h>
#include <stdint.h>

#include "su3mField.h"

namespace delsquared {

  template < class T >
  su3mField<T>::su3mField( const size_t vol, const size_t halo[4], const unsigned int vlength ) {

    field.resize(vol);

    // TODO: call DS_MEMALLOC here
    void* temp = NULL;
    void** ptr = (void**) NULL;
    temp = malloc(9*2*vlength*vol*sizeof(T)+sizeof(void*)+32);
    ptr = (void**)(((uintptr_t)temp+(uintptr_t)32+sizeof(void*))&~(uintptr_t)32);
    ptr[-1] = temp;
    rawmem = (T*)ptr;
    
    for(size_t i = 0; i < vol; ++i){
      field[i].c00.r.v = rawmem + 18*i*vlength;
      field[i].c01.r.v = rawmem + 18*i*vlength + vlength;
      field[i].c02.r.v = rawmem + 18*i*vlength + 2*vlength;
      field[i].c10.r.v = rawmem + 18*i*vlength + 3*vlength;
      field[i].c11.r.v = rawmem + 18*i*vlength + 4*vlength;
      field[i].c12.r.v = rawmem + 18*i*vlength + 5*vlength;
      field[i].c20.r.v = rawmem + 18*i*vlength + 6*vlength;
      field[i].c21.r.v = rawmem + 18*i*vlength + 7*vlength;
      field[i].c22.r.v = rawmem + 18*i*vlength + 8*vlength;

      field[i].c00.i.v = rawmem + 18*i*vlength + 9*vlength;
      field[i].c01.i.v = rawmem + 18*i*vlength + 9*vlength + vlength;
      field[i].c02.i.v = rawmem + 18*i*vlength + 9*vlength + 2*vlength;
      field[i].c10.i.v = rawmem + 18*i*vlength + 9*vlength + 3*vlength;
      field[i].c11.i.v = rawmem + 18*i*vlength + 9*vlength + 4*vlength;
      field[i].c12.i.v = rawmem + 18*i*vlength + 9*vlength + 5*vlength;
      field[i].c20.i.v = rawmem + 18*i*vlength + 9*vlength + 6*vlength;
      field[i].c21.i.v = rawmem + 18*i*vlength + 9*vlength + 7*vlength;
      field[i].c22.i.v = rawmem + 18*i*vlength + 9*vlength + 8*vlength;
    }
  }

  template < class T >
  su3mField<T>::~su3mField(){
    if( (void*)rawmem != NULL )
      free(((void**)rawmem)[-1]);
  }

}

template class delsquared::su3mField<double>;
template class delsquared::su3mField<float>;
