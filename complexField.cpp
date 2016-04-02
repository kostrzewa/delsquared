#include <stdlib.h>
#include <stdint.h>

#include "complexField.h"

namespace delsquared {

  template < class T >
  complexField<T>::complexField( const size_t vol, const size_t halo[4], const unsigned int vlength ) {

    field.resize(vol);

    // TODO: call DS_MEMALLOC here
    void* temp = NULL;
    void** ptr = (void**) NULL;
    temp = malloc(2*vlength*vol*sizeof(T)+sizeof(void*)+32);
    ptr = (void**)(((uintptr_t)temp+(uintptr_t)32+sizeof(void*))&~(uintptr_t)32);
    ptr[-1] = temp;
    rawmem = (T*)ptr;
    

    for(size_t i = 0; i < vol; ++i){
      field[i].r.v = rawmem + i*vlength;
      field[i].i.v = rawmem + i*vlength + vlength*vol;
    }
  }

  template < class T >
  complexField<T>::~complexField(){
    if( (void*)rawmem != NULL )
      free( ((void**)rawmem)[-1] );
  }

}

template class delsquared::complexField<double>;
template class delsquared::complexField<float>;
