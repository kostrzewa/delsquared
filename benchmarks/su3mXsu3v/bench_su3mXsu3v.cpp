#include "../../su3mField.h"
#include "../../su3vField.h"
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
  const size_t          wset = 3*(18+6)*sizeof(float);
  const unsigned int vlength = 32;
  const size_t           vol = cache/wset/vlength;
  const size_t          reps = 10000;
  const size_t         psize = reps*vol*vlength;
  
  std::cout << wset*vol*vlength/pow(1024,2) << " MB working set" << std::endl;
  std::cout << "Problem size " << psize     << std::endl;
  std::cout << "3-Volume "     << vol       << std::endl;
  std::cout << "VLength "      << vlength   << std::endl;
  std::cout << "Repetitions "  << reps      << std::endl << std::endl;

  su3mField<float> mf_a( vol, halo, vlength );
  su3mField<float> mf_b( vol, halo, vlength );
  su3mField<float> mf_c( vol, halo, vlength );
  
  su3vField<float> vf_a( vol, halo, vlength );
  su3vField<float> vf_b( vol, halo, vlength );
  su3vField<float> vf_c( vol, halo, vlength );

  std::random_device r;
  std::mt19937 gen(r());
  std::uniform_real_distribution<float> dis(-0.5,0.5);
  for(size_t i = 0; i < 18*vlength*vol; ++i){
    mf_a.rawmem[i] = dis(gen);
    mf_b.rawmem[i] = dis(gen);
    mf_c.rawmem[i] = dis(gen);
  }
  for(size_t i = 0; i < 6*vlength*vol; ++i){
    vf_a.rawmem[i] = dis(gen);
    vf_b.rawmem[i] = dis(gen);
    vf_c.rawmem[i] = dis(gen);
  }
  
  const size_t Vs = vol;
  omp_set_num_threads(4);
  // TODO: look at modern C++ loops
  std::chrono::time_point<std::chrono::steady_clock> start;
  std::chrono::duration<float> elapsed_seconds;
  
  
  start = std::chrono::steady_clock::now();
  #pragma omp parallel
  { 
  su3v<float> t[3] = { su3v<float>(vlength), su3v<float>(vlength), su3v<float>(vlength) };
  for(int j = 0; j < reps; ++j){
    #pragma omp for nowait
    for(size_t i = 0; i < Vs; ++i){
      su3mXsu3v(vf_c.field[i], mf_a.field[i], vf_b.field[i], vlength, t);
      su3mXsu3v(vf_a.field[i], mf_b.field[i], vf_c.field[i], vlength, t);
      su3mXsu3v(vf_b.field[i], mf_c.field[i], vf_a.field[i], vlength, t);
    }
  }
  }
  elapsed_seconds = std::chrono::steady_clock::now()-start;

  std::cout << vf_b.field[33].c0.r[20] << std::endl;
  std::cout << elapsed_seconds.count() << " seconds" << std::endl;
  // 9*3*2 multiplications, 3*2*2 additions, 3 times
  std::cout << (9*3*2 + 3*2*2)*psize/elapsed_seconds.count()/1e6 << " mflop/s" << std::endl << std::endl;

  start = std::chrono::steady_clock::now();
  #pragma omp parallel
  { 
  su3v<float> t[3] = { su3v<float>(vlength), su3v<float>(vlength), su3v<float>(vlength) };
  for(int j = 0; j < reps; ++j){
    #pragma omp for nowait
    for(size_t i = 0; i < Vs; ++i){
      su3mXsu3v_direct(vf_c.field[i], mf_a.field[i], vf_b.field[i], vlength, t);
      su3mXsu3v_direct(vf_a.field[i], mf_b.field[i], vf_c.field[i], vlength, t);
      su3mXsu3v_direct(vf_b.field[i], mf_c.field[i], vf_a.field[i], vlength, t);
    }
  }
  }
  elapsed_seconds = std::chrono::steady_clock::now()-start;

  std::cout << vf_b.field[33].c0.r[20] << std::endl;
  std::cout << elapsed_seconds.count() << " seconds" << std::endl;
  // 9*3*2 multiplications, 3*2*2 additions, 3 times
  std::cout << (9*3*2 + 3*2*2)*psize/elapsed_seconds.count()/1e6 << " mflop/s" << std::endl << std::endl;
  
  
  start = std::chrono::steady_clock::now();
  #pragma omp parallel
  { 
  su3v<float> t[3] = { su3v<float>(vlength), su3v<float>(vlength), su3v<float>(vlength) };
  for(int j = 0; j < reps; ++j){
    #pragma omp for nowait
    for(size_t i = 0; i < Vs; ++i){
      su3mXsu3v_intrin_float(vf_c.field[i], mf_a.field[i], vf_b.field[i], vlength);
      su3mXsu3v_intrin_float(vf_a.field[i], mf_b.field[i], vf_c.field[i], vlength);
      su3mXsu3v_intrin_float(vf_b.field[i], mf_c.field[i], vf_a.field[i], vlength);
    }
  }
  }
  elapsed_seconds = std::chrono::steady_clock::now()-start;

  std::cout << vf_b.field[33].c0.r[20] << std::endl;
  std::cout << elapsed_seconds.count() << " seconds" << std::endl;
  // 9*3*2 multiplications, 3*2*2 additions, 3 times
  std::cout << (9*3*2 + 3*2*2)*psize/elapsed_seconds.count()/1e6 << " mflop/s" << std::endl << std::endl;
  
  return 0;
}
