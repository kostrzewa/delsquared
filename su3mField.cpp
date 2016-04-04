#include <stdlib.h>
#include <stdint.h>
#include <cstring>
#include <iostream>

#include "su3mField.h"

namespace delsquared {

  template < class T >
  su3mField<T>::su3mField( const size_t vol, const size_t halo[4], const unsigned int vlength ) {

    field.resize(vol);

    // TODO: call DS_MEMALLOC here
    void* temp = NULL;
    void** ptr = (void**) NULL;
    temp = malloc(9*2*vlength*vol*sizeof(T)+sizeof(void*)+31);
    ptr = (void**)(((uintptr_t)temp+(uintptr_t)31+sizeof(void*))&~(uintptr_t)31);
    ptr[-1] = temp;
    rawmem = (T*)ptr;
    
    // first touch allocation, hopefully this makes it NUMA-aware
    #pragma omp parallel for
    for(size_t i = 0; i < vol; ++i){
      // r and i separated by n*vlength
      //field[i].c00.r.v = rawmem + 18*i*vlength;
      //memset( field[i].c00.r.v, 0, vlength*sizeof(T) );
      //field[i].c01.r.v = rawmem + 18*i*vlength + vlength;
      //memset( field[i].c01.r.v, 0, vlength*sizeof(T) );
      //field[i].c02.r.v = rawmem + 18*i*vlength + 2*vlength;
      //memset( field[i].c02.r.v, 0, vlength*sizeof(T) );
      //field[i].c10.r.v = rawmem + 18*i*vlength + 3*vlength;
      //memset( field[i].c10.r.v, 0, vlength*sizeof(T) );
      //field[i].c11.r.v = rawmem + 18*i*vlength + 4*vlength;
      //memset( field[i].c11.r.v, 0, vlength*sizeof(T) );
      //field[i].c12.r.v = rawmem + 18*i*vlength + 5*vlength;
      //memset( field[i].c12.r.v, 0, vlength*sizeof(T) );
      //field[i].c20.r.v = rawmem + 18*i*vlength + 6*vlength;
      //memset( field[i].c20.r.v, 0, vlength*sizeof(T) );
      //field[i].c21.r.v = rawmem + 18*i*vlength + 7*vlength;
      //memset( field[i].c21.r.v, 0, vlength*sizeof(T) );
      //field[i].c22.r.v = rawmem + 18*i*vlength + 8*vlength;
      //memset( field[i].c22.r.v, 0, vlength*sizeof(T) );

      //field[i].c00.i.v = rawmem + 18*i*vlength + 9*vlength;
      //memset( field[i].c00.i.v, 0, vlength*sizeof(T) );
      //field[i].c01.i.v = rawmem + 18*i*vlength + 9*vlength + vlength;
      //memset( field[i].c01.i.v, 0, vlength*sizeof(T) );
      //field[i].c02.i.v = rawmem + 18*i*vlength + 9*vlength + 2*vlength;
      //memset( field[i].c02.i.v, 0, vlength*sizeof(T) );
      //field[i].c10.i.v = rawmem + 18*i*vlength + 9*vlength + 3*vlength;
      //memset( field[i].c10.i.v, 0, vlength*sizeof(T) );
      //field[i].c11.i.v = rawmem + 18*i*vlength + 9*vlength + 4*vlength;
      //memset( field[i].c11.i.v, 0, vlength*sizeof(T) );
      //field[i].c12.i.v = rawmem + 18*i*vlength + 9*vlength + 5*vlength;
      //memset( field[i].c12.i.v, 0, vlength*sizeof(T) );
      //field[i].c20.i.v = rawmem + 18*i*vlength + 9*vlength + 6*vlength;
      //memset( field[i].c20.i.v, 0, vlength*sizeof(T) );
      //field[i].c21.i.v = rawmem + 18*i*vlength + 9*vlength + 7*vlength;
      //memset( field[i].c21.i.v, 0, vlength*sizeof(T) );
      //field[i].c22.i.v = rawmem + 18*i*vlength + 9*vlength + 8*vlength;
      //memset( field[i].c22.i.v, 0, vlength*sizeof(T) );

      // r and i separated by n*vol
      field[i].c00.r.v = rawmem + 9*i*vlength;
      memset( field[i].c00.r.v, 0, vlength*sizeof(T) );
      field[i].c01.r.v = rawmem + 9*i*vlength + vlength;
      memset( field[i].c01.r.v, 0, vlength*sizeof(T) );
      field[i].c02.r.v = rawmem + 9*i*vlength + 2*vlength;
      memset( field[i].c02.r.v, 0, vlength*sizeof(T) );
      field[i].c10.r.v = rawmem + 9*i*vlength + 3*vlength;
      memset( field[i].c10.r.v, 0, vlength*sizeof(T) );
      field[i].c11.r.v = rawmem + 9*i*vlength + 4*vlength;
      memset( field[i].c11.r.v, 0, vlength*sizeof(T) );
      field[i].c12.r.v = rawmem + 9*i*vlength + 5*vlength;
      memset( field[i].c12.r.v, 0, vlength*sizeof(T) );
      field[i].c20.r.v = rawmem + 9*i*vlength + 6*vlength;
      memset( field[i].c20.r.v, 0, vlength*sizeof(T) );
      field[i].c21.r.v = rawmem + 9*i*vlength + 7*vlength;
      memset( field[i].c21.r.v, 0, vlength*sizeof(T) );
      field[i].c22.r.v = rawmem + 9*i*vlength + 8*vlength;
      memset( field[i].c22.r.v, 0, vlength*sizeof(T) );

      field[i].c00.i.v = rawmem + 9*i*vlength + 9*vlength*vol;
      memset( field[i].c00.i.v, 0, vlength*sizeof(T) );
      field[i].c01.i.v = rawmem + 9*i*vlength + vlength + 9*vlength*vol; 
      memset( field[i].c01.i.v, 0, vlength*sizeof(T) );
      field[i].c02.i.v = rawmem + 9*i*vlength + 2*vlength + 9*vlength*vol; 
      memset( field[i].c02.i.v, 0, vlength*sizeof(T) );
      field[i].c10.i.v = rawmem + 9*i*vlength + 3*vlength + 9*vlength*vol; 
      memset( field[i].c10.i.v, 0, vlength*sizeof(T) );
      field[i].c11.i.v = rawmem + 9*i*vlength + 4*vlength + 9*vlength*vol; 
      memset( field[i].c11.i.v, 0, vlength*sizeof(T) );
      field[i].c12.i.v = rawmem + 9*i*vlength + 5*vlength + 9*vlength*vol; 
      memset( field[i].c12.i.v, 0, vlength*sizeof(T) );
      field[i].c20.i.v = rawmem + 9*i*vlength + 6*vlength + 9*vlength*vol; 
      memset( field[i].c20.i.v, 0, vlength*sizeof(T) );
      field[i].c21.i.v = rawmem + 9*i*vlength + 7*vlength + 9*vlength*vol; 
      memset( field[i].c21.i.v, 0, vlength*sizeof(T) );
      field[i].c22.i.v = rawmem + 9*i*vlength + 8*vlength + 9*vlength*vol; 
      memset( field[i].c22.i.v, 0, vlength*sizeof(T) );
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
