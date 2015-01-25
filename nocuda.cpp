/// this file exists to remove C++11 from CUDA, to support outdated nvcc compilers
#include "lkt.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>    
#include <assert.h>
#include <math.h>
//#include "mergesort.hh"
#include "tbb/tbb.h"
#include "quicksort.hh"

/*
std::ostream& operator<<(std::ostream& s, const lkt_point& p) {
  s << "{" << p.x << ", " << p.y << ", " << p.key << "}";
  return s;
}
*/

/// \todo change these to C++ and use templates. Or Macros. Something.

/// finds a heuristic value, using the given sample rate, splitting on the x-axis
/// \param retval[out] next splitpoint. Param instead of return, because TBB.
static void lkt_find_splitpoint_x(ord_t* retval, lkt_point* begin, lkt_point* end, size_t sample_rate) {
  double average = 0.0;
  size_t samples_taken = 0;
  for(lkt_point* i = begin; i < end; i += sample_rate, ++samples_taken) {
    average += i->x;
  }
  average /= samples_taken;
  *retval = average;
}
/*
/// finds a heuristic value, using the given sample rate, splitting on the y-axis
static ord_t lkt_find_splitpoint_y(lkt_point* begin, lkt_point* end, size_t sample_rate) {
  double average = 0.0;
  size_t samples_taken = 0;
  for(lkt_point* i = begin; i < end; i += sample_rate, ++samples_taken) {
    average += i->y;
  }
  average /= samples_taken;
  return average;
}
*/
/// find the next splitpoint on the y axis, ignoring values less than the given min
/// \param retval[out] next splitpoint. Param instead of return, because TBB.
static void lkt_find_next_splitpoint_y_l(ord_t* retval, lkt_point* begin, lkt_point* end, size_t sample_rate, ord_t min) {
  double average = 0.0;
  size_t samples_taken = 0;
  const size_t sample_distance = ((end - begin) / sample_rate) + 1;
  for(lkt_point* i = begin; i < end; i += sample_distance, ++samples_taken) {
    if(i->y < min)
      continue;
    average += i->y;
  }
  average /= samples_taken;
  *retval = average;
}

/// find the next splitpoint on the y axis, ignoring values greater than the given max
/// \param retval[out] next splitpoint. Param instead of return, because TBB.
static void lkt_find_next_splitpoint_y_r(ord_t* retval, lkt_point* begin, lkt_point* end, size_t sample_rate, ord_t max) {
  double average = 0.0;
  size_t samples_taken = 0;
  const size_t sample_distance = ((end - begin) / sample_rate) + 1;
  for(lkt_point* i = begin; i < end; i += sample_distance, ++samples_taken) {
    if(i->y > max)
      continue;
    average += i->y;
  }
  average /= samples_taken;
  *retval = average;
}
/// find the next splitpoint on the x axis, ignoring values less than the given min
/// \param retval[out] next splitpoint. Param instead of return, because TBB.
static void lkt_find_next_splitpoint_x_l(ord_t* retval, lkt_point* begin, lkt_point* end, size_t sample_rate, ord_t min) {
  double average = 0.0;
  size_t samples_taken = 0;
  const size_t sample_distance = ((end - begin) / sample_rate) + 1;
  for(lkt_point* i = begin; i < end; i += sample_distance, ++samples_taken) {
    if(i->x < min)
      continue;
    average += i->x;
  }
  average /= samples_taken;
  *retval = average;
}
/// find the next splitpoint on the x axis, ignoring values greater than the given max
/// \param retval[out] next splitpoint. Param instead of return, because TBB.
static void lkt_find_next_splitpoint_x_r(ord_t* retval, lkt_point* begin, lkt_point* end, size_t sample_rate, ord_t max) {
  double average = 0.0;
  size_t samples_taken = 0;
  const size_t sample_distance = ((end - begin) / sample_rate) + 1;
  for(lkt_point* i = begin; i < end; i += sample_distance, ++samples_taken) {
    if(i->x > max)
      continue;
    average += i->x;
  }
  average /= samples_taken;
  *retval = average;
}

/// DO NOT change this to use the XOR method. It is slow.
static inline void lkt_swap(lkt_point* a, lkt_point* b) {
  lkt_point old_a = *a;
  *a = *b;
  *b = old_a;
}

static inline size_t get_heap_child_l(const size_t i) {return i * 2 + 1;}
static inline size_t get_heap_child_r(const size_t i) {return i * 2 + 2;}
static inline size_t get_heap_parent(const size_t i) {return (i - 1) / 2;}

static bool point_comparator_x(const lkt_point& a, const lkt_point& b) {
  return a.x < b.x;
}

static bool point_comparator_y(const lkt_point& a, const lkt_point& b) {
  return a.y < b.y;
}

/// \todo change splits to use copied double buffer, for parallel execution
/// \param sample_rate the rate to sample when finding the split point
static void lkt_sort_parallel(lkt_point* points, size_t len, /* lkt_point* buffer, */
                     size_t sample_rate, lkt_split_point* splitpoints, ord_t splitpoint, size_t splitpoint_i, 
                     bool xaxis, 
                     const unsigned short current_depth, const unsigned short max_depth) {
  const size_t PARALLEL_QUICKSORT_THREADS = 8;
/*
  fprintf(stderr, "lkt_sort called for splitpoint %f\n", splitpoint);
  fprintf(stderr, "lkt_sort called for points %p, len %lu\n", (void*)points, len);
  fprintf(stderr, "lkt_sort called for splitpoint_i %lu\n", splitpoint_i);
  fprintf(stderr, "lkt_sort next splitpoint_is: %lu, %lu\n", get_heap_child_l(splitpoint_i), get_heap_child_r(splitpoint_i));

  fprintf(stderr, "sort at splitpoint_i = %lu\n", splitpoint_i);
  fflush(stdout);
*/
  if(len < 2 || current_depth == max_depth || splitpoint_i > len)
    return;

//  memcpy(buffer, points, sizeof(lkt_point) * len);

  // splitpoint is the value in the points array, by which the points will be partitioned
  // splitpoint_i is the index into the splitpoints array, of the current split.
  // splitpoint_val is the (local) index in the points array, before which values are less than splitpoint.

  const lkt_point splitpoint_point = {splitpoint, splitpoint, 0}; // cheaper to just assign both, than a conditional

  typedef void (*splitpoint_finder_func_t)(ord_t* retval, lkt_point* begin, lkt_point* end, size_t sample_rate, ord_t min);
  typedef bool (*comparator_func_t)(const lkt_point&, const lkt_point&);

  ord_t          next_split_l;
  ord_t          next_split_r;
  uint_least64_t splitpoint_val;

  splitpoint_finder_func_t next_split_l_func;
  splitpoint_finder_func_t next_split_r_func;
  comparator_func_t        comparator_func;

  // avoids conditionals.
  next_split_l_func = (splitpoint_finder_func_t)((intptr_t)lkt_find_next_splitpoint_x_l * xaxis + (intptr_t)lkt_find_next_splitpoint_y_l * !xaxis);
  next_split_r_func = (splitpoint_finder_func_t)((intptr_t)lkt_find_next_splitpoint_x_r * xaxis + (intptr_t)lkt_find_next_splitpoint_y_r * !xaxis);
  comparator_func   = (comparator_func_t)((intptr_t)point_comparator_x * xaxis + (intptr_t)point_comparator_y * !xaxis);

  next_split_l_func(&next_split_l, points, points + len, sample_rate, splitpoint);
  next_split_r_func(&next_split_r, points, points + len, sample_rate, splitpoint);
  parallel_quicksort_partition(&splitpoint_val, points, &points[len], splitpoint_point, PARALLEL_QUICKSORT_THREADS, comparator_func);

/* 
  tbb::parallel_invoke([&]() {next_split_l_func(&next_split_l, buffer, buffer + len, sample_rate, splitpoint);},
                       [&]() {next_split_r_func(&next_split_r, buffer, buffer + len, sample_rate, splitpoint);},
                       [&]() {parallel_quicksort_partition(&splitpoint_val, points, &points[len], splitpoint_point, PARALLEL_QUICKSORT_THREADS, comparator_func);});
*/

/*
  if(splitpoint_val == 0 || splitpoint_val == len) {
    fprintf(stderr, "not setting splitpoints: splitpoint_val is %lu", splitpoint_val);
    return;
  }
*/
//  fprintf(stderr, "DEBUG: setting splitpoints[%lu] = {%lu,%f}\n", splitpoint_i, splitpoint_val, splitpoint);
//  fprintf(stderr, "lkt_sort [%lu:%lu) recursing [%lu:%lu) and [%lu:%lu)\n", 0ul, len,  0ul, splitpoint_val, splitpoint_val, splitpoint_val + (len - splitpoint_val));

  splitpoints[splitpoint_i].value = splitpoint;
  splitpoints[splitpoint_i].index = splitpoint_val;

  tbb::parallel_invoke([&]() {lkt_sort_parallel(points, splitpoint_val, /* buffer, */ sample_rate, splitpoints, 
                                                next_split_l, get_heap_child_l(splitpoint_i), 
                                                !xaxis, current_depth + 1, max_depth);},
                       [&]() {lkt_sort_parallel(points + splitpoint_val, len - splitpoint_val, /* buffer + splitpoint_val, */ sample_rate, splitpoints, 
                                                next_split_r, get_heap_child_r(splitpoint_i), 
                                                !xaxis, current_depth + 1, max_depth);});
}

/// returns a heap
linear_kdtree lkt_create_parallel(lkt_point* points, size_t len) {
//  fprintf(stderr, "lkt_create called for points %p, true end %p \n", (void*)points, (void*)(points + len));

  if(sizeof(mortoncode_t) != 4) {
    fprintf(stderr, "mortoncode_t NOT 32 BITS! ERROR!ERROR!ERROR!"); /// \todo fix to be static_assert
    exit(1);
  }

  const unsigned short depth = sizeof(mortoncode_t) * CHAR_BIT;

  linear_kdtree tree;
  tree.points = points;
  tree.len = len;
  tree.split_points_len = len; // min(len, pow(2, depth)), but the length will always be less. We'd run out of memory before it wasn't.

//  fprintf(stderr, "lkt_create depth %lu\n", (size_t)depth);
//  fprintf(stderr, "lkt_create split_points_len %lu\n", (size_t)tree.split_points_len);
//  fprintf(stderr, "lkt_create newing split_points size %lu\n", sizeof(lkt_split_point) * tree.split_points_len);
  tree.split_points = new lkt_split_point[tree.split_points_len];
//  fprintf(stderr, "lkt_create newed split_points\n");
  tree.split_depth = depth;
  memset(tree.split_points, '\0', sizeof(lkt_split_point) * tree.split_points_len); // debug

  const size_t sample_rate = 100;
  ord_t initial_splitpoint;
  lkt_find_splitpoint_x(&initial_splitpoint, points, points + len, sample_rate);
  const size_t initial_splitpoint_i = 0;

//  fprintf(stderr, "lkt_create sorting\n");

//  lkt_point* buffer = new lkt_point[len];

  lkt_sort_parallel(points, len, /* buffer, */ sample_rate, tree.split_points, initial_splitpoint, initial_splitpoint_i, true, 0, depth);

//  delete[] buffer;

//  fprintf(stderr, "lkt_create coding\n");

  tree.morton_codes = lkt_create_mortoncodes_parallel(tree.points, tree.len, tree.split_points, tree.split_points_len, depth);

//  fprintf(stderr, "lkt_create returning\n");
  return tree;
}

/// \return array of morton codes, of len length. Caller takes ownership.
mortoncode_t* lkt_create_mortoncodes_parallel(lkt_point* points, size_t len, lkt_split_point* split_points, size_t split_points_len, size_t split_depth) {
  if(sizeof(mortoncode_t) * CHAR_BIT < split_depth) {
    fprintf(stderr, "mortoncode_t LESS THAN split_depth! ERROR!ERROR!ERROR!"); /// \todo fix to be static_assert
    exit(1);
  }

//  fprintf(stderr, "lkt_create_mortoncodes len %lu, split_points_len %lu, split_depth %lu\n", len, split_points_len, split_depth);

  mortoncode_t* codes = new mortoncode_t[len];

  return codes; // debug !! undo !!

/*
  for(size_t i = 0, end = len; i != end; ++i) { /// \todo vectorize
    mortoncode_t code = 0;
    const lkt_point point = points[i];

//    fprintf(stderr, "lkt_create_mortoncodes point {%f, %f}\n", point.x, point.y);

    int j = 0;
    // in practice, this can be optimised to remove 'j < jend' because we'd run out of memory before using 2^h instead of len(points)
    for(long jend = split_depth, xaxis = true, split_pos = 0; j < jend && split_pos < (long) split_points_len; ++j, xaxis = !xaxis) {
      const lkt_split_point splitpoint = split_points[split_pos];
      const ord_t split_point_val      = splitpoint.value;
      const ord_t point_ord            = point.x * xaxis + point.y * !xaxis; // xaxis ? point.x : point.y;
      const unsigned int left          = point_ord < split_point_val;

      code = code | (left << j);

//      fprintf(stderr, "\tj %d,\tsplitpoint_pos %ld,\tsplitpoint_val %f,\t%s,\tp_ord %f,\t%s, code %u\n", j, split_pos, split_point_val, xaxis ? "x" : "y", point_ord, left ? "left" : "right", code);

      split_pos = 2 * split_pos + (1 + !left); // left ? 1 : 2 (heap left child is 2*i+1, right child is 2*i+2)
    }
    
//    fprintf(stderr, "lkt_create_mortoncodes j %d code_size %lu \n", j, (size_t)(sizeof(mortoncode_t) * CHAR_BIT));
//    fprintf(stderr, "lkt_create_mortoncodes shifting the remaining %lu\n", (size_t)(sizeof(mortoncode_t) * CHAR_BIT - j - 1));
//    code = code << (sizeof(mortoncode_t) * CHAR_BIT - j - 1); // shift the rest of the int, so the most-significant-bit is the first split pos

//    fprintf(stderr, "lkt_create_mortoncodes code %lu: %u\n", i, code);
    codes[i] = code;
  }
*/
  return codes;
}


/*
// x value ALONE is used for comparison, to create an xpack
bool operator<(const lkt_point& rhs, const lkt_point& lhs) {
  return rhs.x < lhs.x;
}

bool operator<(const lqt_unified_node& rhs, const lqt_unified_node& lhs) {
  return rhs.location < lhs.location;
}

linear_quadtree_unified tbb_sortify_unified(linear_quadtree_unified lqt, const size_t threads) {
//  auto lowxpack = [](const rtree_point& rhs, const rtree_point& lhs) {
//    return rhs.x < rhs.y;
//  };
  tbb::task_scheduler_init init(threads);
  tbb::parallel_sort(lqt.nodes, lqt.nodes + lqt.length);
  return lqt;
}

/// does not block for GPU memory. Will fail, if GPU memory is insufficient.
linear_quadtree_unified lqt_create_heterogeneous(lkt_point* points, size_t len, 
                                                       ord_t xstart, ord_t xend, 
                                                       ord_t ystart, ord_t yend,
                                                       size_t* depth, const size_t threads) {
  return tbb_sortify_unified(lqt_nodify_cuda_unified(points, len, xstart, xend, ystart, yend, depth), threads);
}

/// \param threads the number of threads to use when sorting. ONLY used in the 'sort' part of the algorithm
rtree cuda_create_rtree_heterogeneously_mergesort(rtree_point* points, const size_t len, const size_t threads) {
  rtree_leaf* leaves = cuda_create_leaves_together(parallel_mergesort(points, points + len, threads), len);
  const size_t leaves_len = DIV_CEIL(len, RTREE_NODE_SIZE);

  rtree_node* previous_level = (rtree_node*) leaves;
  size_t      previous_len = leaves_len;
  size_t      depth = 1; // leaf level is 0
  while(previous_len > RTREE_NODE_SIZE) {
    previous_level = cuda_create_level(previous_level, previous_len);
    previous_len = DIV_CEIL(previous_len, RTREE_NODE_SIZE);
    ++depth;
  }

  rtree_node* root = (rtree_node*) malloc(sizeof(rtree_node));
  init_boundary(&root->bounding_box);
  root->num = previous_len;
  root->children = previous_level;
  for(size_t i = 0, end = previous_len; i != end; ++i)
    update_boundary(&root->bounding_box, &root->children[i].bounding_box);
  ++depth;

  rtree tree = {depth, root};
  return tree;
}
*/
/*
/// SISD sort via single CPU core (for benchmarks)
rtree cuda_create_rtree_sisd(rtree_point* points, const size_t len) {
  std::sort(points, points + len);
  rtree_leaf* leaves = cuda_create_leaves_together(points, len);
  const size_t leaves_len = DIV_CEIL(len, RTREE_NODE_SIZE);

  rtree_node* previous_level = (rtree_node*) leaves;
  size_t      previous_len = leaves_len;
  size_t      depth = 1; // leaf level is 0
  while(previous_len > RTREE_NODE_SIZE) {
    previous_level = cuda_create_level(previous_level, previous_len);
    previous_len = DIV_CEIL(previous_len, RTREE_NODE_SIZE);
    ++depth;
  }

  rtree_node* root = (rtree_node*) malloc(sizeof(rtree_node));
  init_boundary(&root->bounding_box);
  root->num = previous_len;
  root->children = previous_level;
  for(size_t i = 0, end = previous_len; i != end; ++i)
    update_boundary(&root->bounding_box, &root->children[i].bounding_box);
  ++depth;

  rtree tree = {depth, root};
  return tree;
}
*/
