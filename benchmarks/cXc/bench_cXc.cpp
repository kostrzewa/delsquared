#include "../../complexField.h"
#include "../../basicmath.h"

#include <iostream>
#include <random>
#include <chrono>
#include <cmath>
#include <omp.h>

using namespace delsquared;

int main(void){
  size_t halo[4] = {0,0,0,0};

  const size_t         cache = 3*pow(1024,2);
  const size_t          wset = 6*sizeof(float);
  const unsigned int vlength = 128;
  const size_t           vol = cache/wset/vlength;
  const size_t          reps = 20000;
  const size_t         psize = reps*vol*vlength;
  
  std::cout << wset*vol*vlength/pow(1024,2) << " MB working set" << std::endl;
  std::cout << "Problem size " << psize     << std::endl;
  std::cout << "3-Volume "     << vol       << std::endl;
  std::cout << "VLength "      << vlength   << std::endl;
  std::cout << "Repetitions "  << reps      << std::endl << std::endl;

  complexField<float> cf_a( vol, halo, vlength );
  complexField<float> cf_b( vol, halo, vlength );
  complexField<float> cf_c( vol, halo, vlength );

  std::random_device r;
  std::mt19937 gen(r());
  std::uniform_real_distribution<> dis(-0.5,0.5);
  for(size_t i = 0; i < vol; ++i){
    for(unsigned int j = 0; j < vlength; ++j){
      cf_a.field[i].r[j] = dis(gen);
      cf_a.field[i].i[j] = dis(gen);
      cf_b.field[i].r[j] = dis(gen);
      cf_b.field[i].i[j] = dis(gen);
      cf_c.field[i].r[j] = dis(gen);
      cf_c.field[i].i[j] = dis(gen);
    }
  }

  // TODO: look at modern C++ loops
  std::chrono::time_point<std::chrono::steady_clock> start;
  start = std::chrono::steady_clock::now();
  const size_t Vs = vol;
#ifdef DS_OMP
  omp_set_num_threads(4);
  #pragma omp parallel
  { 
#endif
  for(int j = 0; j < reps; ++j){
#ifdef DS_OMP
    #pragma omp for nowait
#endif
    for(size_t i = 0; i < Vs; ++i){
        cXc(cf_c.field[i], cf_a.field[i], cf_b.field[i], vlength);
        cXc(cf_a.field[i], cf_c.field[i], cf_b.field[i], vlength);
        cXc(cf_b.field[i], cf_a.field[i], cf_c.field[i], vlength);
    }
  }
#ifdef DS_OMP
  }
#endif
  std::chrono::duration<double> elapsed_seconds = std::chrono::steady_clock::now()-start;

  std::cout << cf_b.field[33].r[20] << std::endl;
  std::cout << elapsed_seconds.count() << " seconds" << std::endl;
  // 6 flops per complex multiply, we do it 3 times, 1e6 to get mflops
  std::cout << 3*6*psize/elapsed_seconds.count()/1e6 << " mflop/s" << std::endl << std::endl;

  return 0;
}
