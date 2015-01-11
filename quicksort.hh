#ifndef quicksort_hh
#define quicksort_hh

#include <cstdlib>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <atomic>
#include "tbb/tbb.h"
#include <iostream> // debug

using std::cout; // debug!
using std::endl; // debug!

const size_t BLOCK_SIZE = 10; ///< @todo find the optimum. Probably ~1/2 cache size
const uint_least64_t block_end = ~0;

template <typename T>
void print_block(T* block, const size_t block_end) {
  cout << "{";
  for(size_t i = 0, end = block_end; i != end; ++i) {
    cout << block[i] << ",\n";
  }
  cout << "}" << endl;
}

template <typename T>
bool test_neutralised_block(const T* block, const size_t block_len, const T& pivot, const uint_least64_t block_pos, const uint_least64_t left_block_pos, const uint_least64_t right_block_pos) {
  /// \return whether all elements of the block are less than or equal to the pivot
  auto block_is_less = [&](const T* block, const size_t len) -> bool {
    for(size_t i = 0, end = len; i != end; ++i) {
      if(!(block[i] < pivot))
        return false;
    }
    return true;
  };

  /// \return whether all elements of the block are greater than or equal to the pivot
  auto block_is_greater = [&](const T* block, const size_t len) -> bool {
    for(size_t i = 0, end = len; i != end; ++i) {
      if(block[i] < pivot)
        return false;
    }
    return true;
  };

  if(block_pos < left_block_pos && !block_is_less(block, block_len)) {
    cout << "FAILURE: block " << block_pos << " contains an element greater than the pivot" << endl;
    print_block(block, block_len);
    return false;
  } else if(block_pos > right_block_pos && !block_is_greater(block, block_len)) {
    cout << "FAILURE: block " << block_pos << " contains an element less than the pivot" << endl;
    print_block(block, block_len);
    return false;
  }
  return true;
}

template <typename T>
void test_partitioner(T* array, T* array_end, const T& pivot, const size_t block_size, const uint_least64_t left_block_pos, const uint_least64_t right_block_pos, 
                 const uint_least64_t remaining_blocks_len, const uint_least64_t* remaining_block_indices) {
  const size_t num_blocks = ceil((double)(array_end - array) / BLOCK_SIZE);
  const size_t len = array_end - array;

  cout << "left_block_pos:" << left_block_pos << endl;
  cout << "right_block_pos:" << right_block_pos << endl;
  cout << "array len:" << len << endl;
  cout << "block size:" << block_size << endl;
  cout << "block num:" << num_blocks << endl;
  cout << "remaining blocks: {";
  for(size_t i = 0, end = remaining_blocks_len; i != end; ++i) {
    cout << remaining_block_indices[i] << ", ";
  }
  cout << "}" << endl;

  auto block_is_remaining = [=](uint_least64_t blocki) -> bool {
    for(size_t i = 0, end = remaining_blocks_len; i != end; ++i) {
      if(remaining_block_indices[i] == blocki)
        return true;
    }
    return false;
  };

  bool failed = false;
  for(size_t i = 0, end = num_blocks; i != end; ++i) {
    if(block_is_remaining(i))
      continue; // skip unneutralized blocks
    const T*     block     = &array[i * block_size];
    const size_t block_len = i * block_size + block_size < len ? block_size : array_end - block; // this could be made more efficient by only checking if the block is the last block. I think.
    const bool block_good = test_neutralised_block(block, block_len, pivot, (uint_least64_t)i, left_block_pos, right_block_pos);
    if(!block_good)
      failed = true;
  }
  if(failed)
    cout << "FAILED" << endl;
  else
    cout << "SUCCESS" << endl;
}

/// \return index of the next left block, or block_end if there are no left blocks left.
uint_least64_t get_left_block(std::atomic_uint_least64_t* next_left_block, const std::atomic_uint_least64_t* next_right_block) {
  const uint_least64_t lblock = next_left_block->fetch_add(1);
  const uint_least64_t rblock = next_right_block->load() + 1; // +1 because load() is one greater than the last block someone has (they fetched-and-subbed).
  
  const uint_least64_t next = lblock >= rblock ? block_end : lblock;
  return next;
}

/// \return index of the next right block, or block_end if there are no left blocks left.
uint_least64_t get_right_block(std::atomic_uint_least64_t* next_right_block, const std::atomic_uint_least64_t* next_left_block) {
  const uint_least64_t rblock = next_right_block->fetch_sub(1);
  const uint_least64_t lblock = next_left_block->load() - 1; // -1 because load() is one less than the last block someone has (they fetched-and-added).
  const uint_least64_t next = rblock <= lblock ? block_end : rblock;
  return next;
}

enum Neutralised {
  NEUTRALISED_LEFT,
  NEUTRALISED_RIGHT,
  NEUTRALISED_BOTH,
};

template <typename T>
Neutralised neutralise(T* lblock, const uint_least64_t lend, T* rblock, const uint_least64_t rend, const T& pivot) {
  auto get_next_l = [&](uint_least64_t li) -> uint_least64_t {
    for(; li != lend && lblock[li] < pivot; ++li);
    return li;
  };

  auto get_next_r = [&](uint_least64_t ri) -> uint_least64_t {
    for(; ri != rend  && !(rblock[ri] < pivot); ++ri);
    return ri;
  };
  
  uint_least64_t li;
  uint_least64_t ri;
  for(li = get_next_l(0), ri = get_next_r(0); li != lend && ri != rend; li = get_next_l(li + 1), ri = get_next_r(ri + 1))
    std::swap(lblock[li], rblock[ri]);

  if(li == lend) {
    if(ri == rend) {
      return NEUTRALISED_BOTH;
    }
    else {
      return NEUTRALISED_LEFT;
    }
  } else {
    return NEUTRALISED_RIGHT;
  }
}

template <typename T>
void partitioner(T* array, T* end, const size_t block_size, std::atomic_uint_least64_t* next_left_block, std::atomic_uint_least64_t* next_right_block, 
                 uint_least64_t* remaining_block_indices, std::atomic_uint_least64_t* next_remaining, const T& pivot) {
  uint_least64_t lblock_i = get_left_block(next_left_block, next_right_block);
  uint_least64_t rblock_i = get_right_block(next_right_block, next_left_block);

  while(lblock_i != block_end && rblock_i != block_end) {
    const uint_least64_t lend = &array[lblock_i * block_size + block_size] >= end ? (end - array) - (lblock_i * block_size) : block_size;
    const uint_least64_t rend = &array[rblock_i * block_size + block_size] >= end ? (end - array) - (rblock_i * block_size) : block_size;
    const Neutralised neutralised = neutralise(&array[lblock_i * block_size], lend, &array[rblock_i * block_size], rend, pivot);
    switch(neutralised) {
    case NEUTRALISED_LEFT:
      lblock_i = get_left_block(next_left_block, next_right_block);
      break;
    case NEUTRALISED_RIGHT:
      rblock_i = get_right_block(next_right_block, next_left_block);
      break;
    case NEUTRALISED_BOTH:
      rblock_i = get_right_block(next_right_block, next_left_block);
      lblock_i = get_left_block(next_left_block, next_right_block);
      break;
    }
  }

  if(lblock_i == block_end) {
    next_left_block->fetch_sub(1); // undo the atomic add we just did
  }
  if(rblock_i == block_end) {
    next_right_block->fetch_add(1); // undo the atomic sub we just did
  }

  if(lblock_i != block_end) {
    const uint_least64_t remaining_i = next_remaining->fetch_add(1);
    remaining_block_indices[remaining_i] = lblock_i;
  } else if(rblock_i != block_end) {
    const uint_least64_t remaining_i = next_remaining->fetch_add(1);
    remaining_block_indices[remaining_i] = rblock_i;
  }
}

template <typename T>
void neutralise_remaining(T* array, T* end, const size_t block_size, uint_least64_t* remaining_block_indices, size_t remaining_blocks) {

}

template <typename T>
T* parallel_quicksort_partition(T* array, T* end, const T& pivot, size_t threads) {
  tbb::task_scheduler_init init(threads + 1); // +1 for the manager thread
  const size_t num_blocks = ceil((double)(end - array) / BLOCK_SIZE);

  std::unique_ptr<uint_least64_t> remaining_block_indices(new uint_least64_t[threads]);

  std::atomic_uint_least64_t next_left_block(0);
  std::atomic_uint_least64_t next_right_block(num_blocks - 1);
  std::atomic_uint_least64_t next_remaining(0);

  /// \todo change to spawn lots o' microthreads
  tbb::parallel_for(size_t(0), threads, size_t(1), [&](size_t) {
      partitioner(array, end, BLOCK_SIZE, &next_left_block, &next_right_block, remaining_block_indices.get(), &next_remaining, pivot);
    });

  const uint_least64_t left_block_pos = next_left_block.load();
  const uint_least64_t right_block_pos = next_right_block.load();
  const uint_least64_t remaining_blocks_len = next_remaining.load();

  // debug
  test_partitioner(array, end, pivot, BLOCK_SIZE, left_block_pos, right_block_pos, remaining_blocks_len, remaining_block_indices.get());

//  neutralise_remaining(array, end, block_size, remaining_block_indices.get(), next_remaining.load());

  return array;
}


#endif // quicksort_hh
