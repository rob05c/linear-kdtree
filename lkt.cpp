#include "lkt.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <iostream>

using std::cout;
using std::endl;

void lkt_delete(linear_kdtree tree) {
  delete[] tree.points;
  delete[] tree.split_points;
  delete[] tree.morton_codes;
}

template <typename T>
inline void lkt_swap(T* a, T* b) {
  T old_a = *a;
  *a = *b;
  *b = old_a;
}

/// \todo make static, once tested
/// \todo replace with parallel algorithm
/// Unlike traditional quicksort partition, We don't actually use a pivot, since we only have the value to partion, not an element.
/// \param xaxis partition based on the x-axis. Else, the y-axis.
/// \return the index of the partition. Everything less than pivot_value is before the returned index. points[return] is the first element greater than pivot_value
size_t quicksort_partition(lkt_point* points, const size_t len, const ord_t pivot_value, const bool xaxis) {
//  const size_t pivot_index = points + len - 1;
//  const lkt_point* pivot_value = points[pivot_index];
//  lkt_swap(pivot_value, points[len - 1]);
//  if(len == 2)
//    return 1;
//  fprintf(stderr, "quicksort_partitioned for len %lu\n", len);

  if(len == 0)
    return 0;

  long i = 0;
  long j = len - 1;

//  fprintf(stderr, "quicksort_partitioned i %ld j %ld\n", i, j);

  /// we duplicate the loops rather than check inside, for efficiency
  if(xaxis) {
    while(i < j) {
      for(; points[i].x < pivot_value && i < (long)len; ++i);
      if(i >= (long)len)
        break;
      for(; points[j].x > pivot_value && j > -1; --j);
      if(j <= 0)
        break;
//      fprintf(stderr, "quicksort_partitioned swapping i %ld j %ld \n", i, j);
      lkt_swap(&points[i], &points[j]);
    }
//    fprintf(stderr, "quicksort_partitioned finished i %ld j %ld\n", i, j);

    // swap i,j such that j is the greater
    if(i > j) {
      int old_i = i;
      i = j;
      j = old_i;
    }

    if(j > 0) {
      if(i > j) {
        if(points[i].x < pivot_value || points[j].x > pivot_value) // less than, because i is now greater than j.
          lkt_swap(&points[i], &points[j]);
      } else {
        if(points[j].x < pivot_value || points[i].x > pivot_value) // less than, because i is now greater than j.
          lkt_swap(&points[i], &points[j]);
      }
    }
  } else {
    while(i < j) {
      for(; points[i].y < pivot_value && i < (long)len; ++i);
      if(i >= (long)len)
        break;
      for(; points[j].y > pivot_value && j > -1; --j);
      if(j <= 0)
        break;
//      fprintf(stderr, "quicksort_partitioned swapping i %ld j %ld\n", i, j);
      lkt_swap(&points[i], &points[j]);
    }
//    fprintf(stderr, "quicksort_partitioned finished i %ld j %ld\n", i, j);

    // swap i,j such that j is the greater
    if(i > j) {
      int old_i = i;
      i = j;
      j = old_i;
    }

    if(j > 0) {
      if(i > j) {
        if(points[i].y < pivot_value || points[j].y > pivot_value) // less than, because i is now greater than j.
          lkt_swap(&points[i], &points[j]);
      } else {
        if(points[j].y < pivot_value || points[i].y > pivot_value) // less than, because i is now greater than j.
          lkt_swap(&points[i], &points[j]);
      }
    }
  }

  if(j < 0)
    j = 0;


  // debug - sanity check
  if(len > 1) {
    for(long k = 0, kend = j; k != kend; ++k) {
      if(xaxis) {
        if(points[k].x > pivot_value) {
          fprintf(stderr, "quicksort_partitioned ERROR: SORT FAILED i %ld j %ld k %ld\n", i, j, k);
        }
      } else {
        if(points[k].y > pivot_value) {
          fprintf(stderr, "quicksort_partitioned ERROR: SORT FAILED i %ld j %ld k %ld\n", i, j, k);
        }
      }
    }
    for(long k = j, kend = len; k != kend; ++k) {
      if(xaxis) {
        if(points[k].x < pivot_value) {
          fprintf(stderr, "quicksort_partitioned ERROR: SORT FAILED i %ld j %ld k %ld\n", i, j, k);
        }
      } else {
        if(points[k].y < pivot_value) {
          fprintf(stderr, "quicksort_partitioned ERROR: SORT FAILED i %ld j %ld k %ld\n", i, j, k);
        }
      }
    }
  }
//  fprintf(stderr, "quicksort_partitioned on %ld\n", j);
  return j;
}

/// \return array of morton codes, of len length. Caller takes ownership.
mortoncode_t* lkt_create_mortoncodes_sisd(lkt_point* points, size_t len, const fixlentree<lkt_split_point>::node* splitpoints) {
  mortoncode_t* codes = new mortoncode_t[len];
  for(size_t i = 0, end = len; i != end; ++i) {
    const lkt_point& point = points[i];
    mortoncode_t& code = codes[i];
    code = 0;
    bool is_x = true;
    for(size_t code_i = 0, j = 0; j != fixlentree<lkt_split_point>::tree_end; ++code_i, is_x = !is_x) {
      const lkt_split_point& splitpoint = splitpoints[j].value;
      const int left = is_x * (point.x < splitpoint.value) + !is_x * (point.y < splitpoint.value);
      code = code | (left << code_i);

      j = splitpoints[j].left * left + splitpoints[j].right * !left;
    }

  }
  return codes;
}
