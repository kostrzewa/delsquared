/* vectorable is a linear chunk of aligned memory for some elementary 
 * data type with padding at the end to fulfill the vetorization requirements
 *
 * in case of communication in the vectorization direction, the beginning
 * and end of this chunk of memory will contain the boundaries so that
 * shifted access becomes possible */

#ifndef _VECTORABLE_H
#define _VECTORABLE_H

#include <cstdlib>
#include <cstdint>

namespace delsquared {

template < class T >
class vectorable {
public:
  vectorable(){ 
    selfAlloc = false; 
    v = (T*) NULL;
  }
  // TODO: call DS_MEMALLOC here, probably use smartpointers
  vectorable(unsigned int vlength){
    void* temp = NULL;                                                                                                                                                              
    void** ptr = (void**) NULL;
    temp = malloc(vlength*sizeof(T)+sizeof(void*)+32);
    ptr = (void**)(((uintptr_t)temp+(uintptr_t)31+sizeof(void*))&~(uintptr_t)31);
    ptr[-1] = temp;
    v = (T*)ptr;
    selfAlloc = true;
  }
  ~vectorable(){
    // TODO: convert to a smartpointer free or a reference return to the memory manager
    if(selfAlloc){
      if((void*)v != NULL){
        free( ((void**)v)[-1] );
      }
    }
  }

  const T & operator ()(size_t i) __restrict__ const { return v[i]; }
  T & operator [](size_t i) __restrict__ { return v[i]; }
  T* v __attribute__ ((aligned (32)));

private:
  bool selfAlloc;
};

}

#endif

