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

/// \return whether all elements of the block are less than or equal to the pivot
template <typename T>
bool block_is_less(const T* block, const size_t len, const T& pivot) {
  for(size_t i = 0, end = len; i != end; ++i) {
    if(!(block[i] < pivot))
      return false;
  }
  return true;
}

/// \return whether all elements of the block are greater than or equal to the pivot
template <typename T>
bool block_is_greater(const T* block, const size_t len, const T& pivot) {
  for(size_t i = 0, end = len; i != end; ++i) {
    if(block[i] < pivot)
      return false;
  }
  return true;
}

template <typename T>
uint_least64_t get_block_end(const T* array, const T* end, const uint_least64_t block_size, uint_least64_t block_i) {
  return  &array[block_i * block_size + block_size] >= end ? (end - array) - block_i * block_size : block_size;
}


template <typename T>
bool test_neutralised_block(const T* block, const size_t block_len, const T& pivot, const uint_least64_t block_pos, const uint_least64_t left_block_pos, const uint_least64_t right_block_pos) {
  if(block_pos < left_block_pos && !block_is_less(block, block_len, pivot)) {
    cout << "FAILURE: block " << block_pos << " contains an element greater than the pivot" << endl;
    print_block(block, block_len);
    return false;
  } else if(block_pos >= right_block_pos && !block_is_greater(block, block_len, pivot)) {
    cout << "FAILURE: block " << block_pos << " contains an element less than the pivot" << endl;
    print_block(block, block_len);
    return false;
  }
  return true;
}

template <typename T>
void test_partitioner(T* array, T* array_end, const T& pivot, const size_t block_size, const uint_least64_t left_block_pos, const uint_least64_t right_block_pos, 
                 const uint_least64_t remaining_blocks_len, const uint_least64_t* remaining_block_indices) {
  cout << "test_partitioner" << endl;

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

/// @todo change both to be left & right, to make ifs simpler
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

/// Note that this function changes the left_block_pos and right_block_pos. 
///      Thus, the caller must be aware their positions are no longer valid.
///      But after this, the array is partitioned, so the caller should never need them.
/// \return the index of the partition
template <typename T>
uint_least64_t neutralise_remaining(T* array, T* array_end, const size_t block_size, uint_least64_t* remaining_block_indices, const uint_least64_t remaining_blocks,
                          uint_least64_t left_block_pos, uint_least64_t right_block_pos, const T& pivot) {
//  cout << "neutralise_remaining" << endl;

  auto block_swap = [&](uint_least64_t li, uint_least64_t lend, uint_least64_t ri, uint_least64_t rend) {
//    cout << "swapping blocks " << li << " and " << ri << endl;
    /// \todo verify blocks are never the end, i.e. less than block_size
    std::swap_ranges(&array[li * block_size], &array[li * block_size + lend], &array[ri * block_size]);
//    cout << "swapped blocks " << endl;
  };

  uint_least64_t* remaining_end = &remaining_block_indices[remaining_blocks];

  std::sort(remaining_block_indices, remaining_block_indices + remaining_blocks);
  
  uint_least64_t remaining_left_i  = 0;
  uint_least64_t remaining_right_i = remaining_blocks - 1;
  uint_least64_t li = remaining_block_indices[remaining_left_i];
  uint_least64_t ri = remaining_block_indices[remaining_right_i];
  uint_least64_t lend = get_block_end(array, array_end, block_size, li);
  uint_least64_t rend = get_block_end(array, array_end, block_size, ri);


  while(remaining_left_i < remaining_right_i) {
//    cout << "neutralising " << li << " and " << ri << endl;
    const Neutralised neutralised = neutralise(&array[li * block_size], lend, &array[ri * block_size], rend, pivot);
    if(neutralised == NEUTRALISED_LEFT || neutralised == NEUTRALISED_BOTH) {
//      cout << "neutralised left " << li << endl;
      if(li >= left_block_pos) { // neutralised on the wrong side: swap this one with the last left block and decrease left_block_pos.
        if(li != left_block_pos) {
          block_swap(li, lend, left_block_pos, get_block_end(array, array_end, block_size, left_block_pos));

          uint_least64_t* swap_unneutralised_pos = std::find(remaining_block_indices, remaining_end, left_block_pos);
          if(swap_unneutralised_pos != remaining_end) {
//            cout << "swapping (left) unneutralised block in remaining_block_indices " << swap_unneutralised_pos - remaining_block_indices << endl;
            *swap_unneutralised_pos = li;
            if(&remaining_block_indices[remaining_right_i] == swap_unneutralised_pos) {
              ri = remaining_block_indices[remaining_right_i];
              rend = get_block_end(array, array_end, block_size, ri);
            }
          }
        }

//        remaining_block_indices[remaining_left_i] = block_end; // li is now neutralised (we don't want the swap_unneutralised_pos finding it)
        ++left_block_pos;
        ++right_block_pos;
      }

      remaining_block_indices[remaining_left_i] = block_end; // li is now neutralised (we don't want the swap_unneutralised_pos finding it)
      ++remaining_left_i;
      li = remaining_block_indices[remaining_left_i];
      lend = get_block_end(array, array_end, block_size, li);
    }
    if(neutralised == NEUTRALISED_RIGHT || neutralised == NEUTRALISED_BOTH) {
//      cout << "neutralised right " << ri << endl;
      if(ri < right_block_pos) { // neutralised on the wrong side: swap this one with the last left block and decrease left_block_pos.
        if(ri != right_block_pos - 1) {
          block_swap(ri, rend, right_block_pos - 1, get_block_end(array, array_end, block_size, right_block_pos - 1));

          uint_least64_t* swap_unneutralised_pos = std::find(remaining_block_indices, remaining_end, right_block_pos - 1);
          if(swap_unneutralised_pos != remaining_end) {
//            cout << "swapping (right) unneutralised block in remaining_block_indices " << swap_unneutralised_pos - remaining_block_indices << endl;
            *swap_unneutralised_pos = ri;
           if(&remaining_block_indices[remaining_left_i] == swap_unneutralised_pos) {
              li = remaining_block_indices[remaining_left_i];
              lend = get_block_end(array, array_end, block_size, li);
            }
          }
        }

//        remaining_block_indices[remaining_right_i] = block_end; // ri is now neutralised (we don't want the swap_unneutralised_pos finding it)
        --left_block_pos;
        --right_block_pos;
      }

      remaining_block_indices[remaining_right_i] = block_end; // ri is now neutralised (we don't want the swap_unneutralised_pos finding it)
      --remaining_right_i;
      ri = remaining_block_indices[remaining_right_i];
      rend = get_block_end(array, array_end, block_size, ri);
    }
  }
  /// \todo determine if left_block_pos and right_block_pos are always the same? Can they differ?


//  cout << "done swapping" << endl;

  // sort the last remaining block, and place it in the middle.

  lend = get_block_end(array, array_end, block_size, remaining_block_indices[remaining_left_i]);
  const uint_least64_t remaining_block = remaining_block_indices[remaining_left_i]; // note remaining_left_i == remaining_right_i at this point

//  cout << "got remaining block " << remaining_block << " (block_end " << block_end << ")" << endl;

  if(remaining_block != block_end) {
//    cout << "sorting from " << remaining_block * block_size << " to " << remaining_block * block_size + lend << endl;
    std::sort(&array[remaining_block * block_size], &array[remaining_block * block_size + lend]);
//    cout << "sorted remaining block" << endl;

    if(remaining_block < left_block_pos) {
//      cout << "remaining block goes left" << endl;
      const uint_least64_t lposend = get_block_end(array, array_end, block_size, left_block_pos - 1);
      block_swap(remaining_block, lend, left_block_pos - 1, lposend);
//      cout << "remaining block went left" << endl;
    } else {
//      cout << "remaining block goes right" << endl;
      const uint_least64_t rposend = get_block_end(array, array_end, block_size, right_block_pos - 1);
      block_swap(remaining_block, lend, right_block_pos, rposend);
//      cout << "remaining block went right" << endl;
    }
  }

/// \todo necessary??
//  if(ri != right_block_pos && ri != left_block_pos - 1) {
//    block_swap(ri, rend, right_block_pos - 1, get_block_end(array, end, block_size, right_block_pos - 1));
//  }

//  cout << "getting ends" << endl;
  uint_least64_t lminusend = get_block_end(array, array_end, block_size, left_block_pos - 1);
  lend = get_block_end(array, array_end, block_size, left_block_pos);
  rend = get_block_end(array, array_end, block_size, right_block_pos);
//  cout << "sorting again" << endl;
  std::sort(&array[left_block_pos * block_size - block_size], &array[left_block_pos * block_size - block_size + lminusend + lend + rend]); // sort from one before left_block_pos to one after (3 blocks).
//  cout << "sorted again, getting pivot" << endl;

//  std::sort(&array[right_block_pos], &array[right_block_pos + rend]);

  uint_least64_t pos = (left_block_pos - 1) * block_size;
  for(T* i = &array[pos]; i != array_end; ++i, ++pos) {
    if(pivot < *i)
      break;
  }

//  cout << "got pivot" << endl;

  return pos;
}

template <typename T>
void test_partitioned(T* array, T* end, const T& pivot, const uint_least64_t pos) {
  cout << "test_partitioned" << endl;
  cout << "pivot:\t" << pivot << endl;
  cout << "pos:\t" << pos << endl;
//  for(int i = pos - 10, end = pos + 10; i != end; ++i) {
//    cout << "array[" << i << "]: " << array[i] << endl;
//  }

  const uint_least64_t len = end - array;

  bool failed = false;
  
  for(int i = 0, end = pos; i != end; ++i) {
    if(pivot < array[i]) {
      cout << "failed index " << i << " value " << array[i] << " from block " << i / BLOCK_SIZE << endl;
      failed = true;
    }
  }

  for(int i = pos, end = len; i != end; ++i) {
    if(array[i] < pivot) {
      cout << "failed index " << i << " value " << array[i] << " from block " << i / BLOCK_SIZE << endl;
      failed = true;
    }
  }

  if(failed) {
    cout << "FAILED" << endl;
  } else {
    cout << "SUCCESS" << endl;
  }
}

/// \return the position of the partition
template <typename T>
uint_least64_t parallel_quicksort_partition(T* array, T* end, const T& pivot, size_t threads) {
  tbb::task_scheduler_init init(threads + 1); // +1 for the manager thread
  const size_t num_blocks = ceil((double)(end - array) / BLOCK_SIZE);

  std::unique_ptr<uint_least64_t> remaining_block_indices(new uint_least64_t[threads]);

  std::atomic_uint_least64_t next_left_block(0);
  std::atomic_uint_least64_t next_right_block(num_blocks - 1);
  std::atomic_uint_least64_t next_remaining(0);

  tbb::parallel_for(size_t(0), threads, size_t(1), [&](size_t) {
      partitioner(array, end, BLOCK_SIZE, &next_left_block, &next_right_block, remaining_block_indices.get(), &next_remaining, pivot);
    });

  uint_least64_t left_block_pos = next_left_block.load();
  uint_least64_t right_block_pos = next_right_block.load() + 1;
  const uint_least64_t remaining_blocks_len = next_remaining.load();
/*
  if(right_block_pos < left_block_pos) {
    const uint_least64_t* remaining_blocks_ptr = remaining_block_indices.get();
    const uint_least64_t* remaining_end = &remaining_block_indices.get()[remaining_blocks_len];
    const uint_least64_t* found = std::find(remaining_blocks_ptr, remaining_end, left_block_pos);
    if(found != remaining_end) { // if the middle block is unneutralised
      right_block_pos = left_block_pos;
    } else {
      const uint_least64_t* found = std::find(remaining_blocks_ptr, remaining_end, right_block_pos);
      if(found != remaining_end) { // if the middle block is unneutralised
        left_block_pos = right_block_pos;
      } else {
        if(block_is_greater(left_block_pos, get_block_end(array, end, block_size, left_block_pos), pivot))
          right_block_pos = left_block_pos;
        else
          left_block_pos = right_block_pos;
      }
    }
  }
*/

  // debug
//  test_partitioner(array, end, pivot, BLOCK_SIZE, left_block_pos, right_block_pos, remaining_blocks_len, remaining_block_indices.get());

  const uint_least64_t pos = neutralise_remaining(array, end, BLOCK_SIZE, remaining_block_indices.get(), remaining_blocks_len, left_block_pos, right_block_pos, pivot);

  // debug
//  test_partitioned(array, end, pivot, pos);

  return pos;
}


#endif // quicksort_hh
