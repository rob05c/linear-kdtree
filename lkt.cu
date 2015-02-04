#include "lkt.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>    
#include <linux/cuda.h>

/// \returns the device totalGlobalMem
inline size_t GetDeviceMemory() {
  cudaDeviceProp properties;
  int deviceNum;
  cudaGetDevice(&deviceNum);
  cudaGetDeviceProperties(&properties, deviceNum);
  return properties.totalGlobalMem;
}

__global__ void create_mortoncodes_kernel(lkt_point* points, mortoncode_t* codes, const fixlentree<lkt_split_point>::node* splitpoints, size_t len) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i >= len)
    return; // skip the final block remainder

  mortoncode_t& code = codes[i];
  code = 0;

  const lkt_point& point = points[i];

  bool is_x = true;
  size_t code_i = 0;
  for(index_t j = 0; j != fixlentree<lkt_split_point>::tree_end;) {
    const lkt_split_point& splitpoint = splitpoints[j].value;

    const int left = is_x * (point.x < splitpoint.value) + !is_x * (point.y < splitpoint.value);

    code = code | (left << code_i);

    j = splitpoints[j].left * left + splitpoints[j].right * !left;
    is_x = !is_x;
    ++code_i;
  }
}

/// \return array of morton codes, of len length. Caller takes ownership.
mortoncode_t* lkt_create_mortoncodes_parallel(lkt_point* points, size_t len, const fixlentree<lkt_split_point>::node* splitpoints) {
  const size_t THREADS_PER_BLOCK = 512;

  lkt_point*                         cuda_points;
  mortoncode_t*                      cuda_codes;
  fixlentree<lkt_split_point>::node* cuda_splitpoints;

  cudaMalloc((void**)&cuda_points,      len * sizeof(lkt_point));
  cudaMalloc((void**)&cuda_codes,       len * sizeof(mortoncode_t));
  cudaMalloc((void**)&cuda_splitpoints, len * sizeof(fixlentree<lkt_split_point>::node));
  cudaMemcpy(cuda_points,      points,      len * sizeof(lkt_point),                         cudaMemcpyHostToDevice);
  cudaMemcpy(cuda_splitpoints, splitpoints, len * sizeof(fixlentree<lkt_split_point>::node), cudaMemcpyHostToDevice);

  create_mortoncodes_kernel<<<(len + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(cuda_points, cuda_codes, cuda_splitpoints, len);

  mortoncode_t* codes = new mortoncode_t[len];
  cudaMemcpy(codes, cuda_codes, len * sizeof(mortoncode_t), cudaMemcpyDeviceToHost);
  cudaFree(cuda_points);
  cudaFree(cuda_codes);
  cudaFree(cuda_splitpoints);
  return codes;
}
