#include <stdlib.h>
#include <stdint.h>
#include <cstring>
#include "su3vField.h"

namespace delsquared {

  template < class T >
  su3vField<T>::su3vField( const size_t vol, const size_t halo[4], const unsigned int vlength ) {

    field.resize(vol);

    // TODO: call DS_MEMALLOC here
    void* temp = NULL;
    void** ptr = (void**) NULL;
    temp = malloc(3*2*vlength*vol*sizeof(T)+sizeof(void*)+31);
    ptr = (void**)(((uintptr_t)temp+(uintptr_t)31+sizeof(void*))&~(uintptr_t)31);
    ptr[-1] = temp;
    rawmem = (T*)ptr;
    
    #pragma omp parallel for    
    for(size_t i = 0; i < vol; ++i){
      field[i].c0.r.v = rawmem + 6*i*vlength;
      memset( field[i].c0.r.v, 0, vlength*sizeof(T) );
      field[i].c1.r.v = rawmem + 6*i*vlength + vlength;
      memset( field[i].c1.r.v, 0, vlength*sizeof(T) );
      field[i].c2.r.v = rawmem + 6*i*vlength + 2*vlength;
      memset( field[i].c2.r.v, 0, vlength*sizeof(T) );
      
      field[i].c0.i.v = rawmem + 6*i*vlength + vlength*3;
      memset( field[i].c0.i.v, 0, vlength*sizeof(T) );
      field[i].c1.i.v = rawmem + 6*i*vlength + vlength*3 + vlength;
      memset( field[i].c1.i.v, 0, vlength*sizeof(T) );
      field[i].c2.i.v = rawmem + 6*i*vlength + vlength*3 + 2*vlength;
      memset( field[i].c2.i.v, 0, vlength*sizeof(T) );
    }
  }

  template < class T >
  su3vField<T>::~su3vField(){
    if( (void*)rawmem != NULL )
      free(((void**)rawmem)[-1]);
  }

}

template class delsquared::su3vField<double>;
template class delsquared::su3vField<float>;
