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
static ord_t lkt_find_splitpoint_x(lkt_point* begin, lkt_point* end, size_t sample_rate) {
  double average = 0.0;
  size_t samples_taken = 0;
  sample_rate = (end - begin) / sample_rate + 1;
  for(lkt_point* i = begin; i < end; i += sample_rate, ++samples_taken) {
    average += i->x;
  }
  average /= samples_taken;
  return average;
}

/// finds a heuristic value, using the given sample rate, splitting on the y-axis
static ord_t lkt_find_splitpoint_y(lkt_point* begin, lkt_point* end, size_t sample_rate) {
  double average = 0.0;
  size_t samples_taken = 0;
  sample_rate = (end - begin) / sample_rate + 1;
  for(lkt_point* i = begin; i < end; i += sample_rate, ++samples_taken) {
    average += i->y;
  }
  average /= samples_taken;
  return average;
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

/// \param sample_rate the rate to sample when finding the split point
static void lkt_sort_parallel(lkt_point* points, size_t len,
                              size_t sample_rate, fixlentree<lkt_split_point>& splitpoints, index_t splitpoint_parent, bool splitpoint_thisisleft,
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
  if(len < 2 || current_depth == max_depth)
    return;

  // splitpoint is the value in the points array, by which the points will be partitioned
  // splitpoint_val is the (local) index in the points array, before which values are less than splitpoint.

  typedef ord_t (*splitpoint_finder_func_t)(lkt_point* begin, lkt_point* end, size_t sample_rate);
  typedef bool (*comparator_func_t)(const lkt_point&, const lkt_point&);


  // avoids conditionals.
  const splitpoint_finder_func_t find_split_func = (splitpoint_finder_func_t)((intptr_t)lkt_find_splitpoint_x * xaxis + (intptr_t)lkt_find_splitpoint_y * !xaxis);
  const comparator_func_t        comparator_func = (comparator_func_t)((intptr_t)point_comparator_x * xaxis + (intptr_t)point_comparator_y * !xaxis);

  const ord_t          splitpoint       = find_split_func(points, points + len, sample_rate);
  const lkt_point      splitpoint_point = {splitpoint, splitpoint, 0}; // cheaper to just assign both, than a conditional

//  cout << "debug len:" << len << " splitpoint: " << splitpoint << endl;

//  cout << "debug partition started" << endl;
  const uint_least64_t splitpoint_val   = parallel_quicksort_partition(points, points + len, splitpoint_point, PARALLEL_QUICKSORT_THREADS, comparator_func); ///< \todo fix last (?) block bug
//  cout << "debug partition finished" << endl;

  const lkt_split_point lkt_splitpoint         = {splitpoint, splitpoint_val};

  const index_t         splitpoint_next_parent = splitpoints.insert(splitpoint_parent, splitpoint_thisisleft, lkt_splitpoint);

  if(splitpoint_next_parent == splitpoints.tree_end)
    return;

  if(splitpoint_val == 0 || splitpoint_val == len - 1)
    return;

//  cout << "splitpoint " << (xaxis ? "x" : "y") << " val: " << splitpoint << " location: " << splitpoint_val << endl;

/*
lkt_sort_parallel(points, splitpoint_val, sample_rate, 
                  splitpoints, splitpoint_next_parent, true,
                  !xaxis, current_depth + 1, max_depth);
lkt_sort_parallel(points + splitpoint_val, len - splitpoint_val, sample_rate, 
                  splitpoints, splitpoint_next_parent, false,
                  !xaxis, current_depth + 1, max_depth);
*/

  tbb::parallel_invoke([&]() {lkt_sort_parallel(points, splitpoint_val, sample_rate, 
                                                splitpoints, splitpoint_next_parent, true,
                                                !xaxis, current_depth + 1, max_depth);},
                       [&]() {lkt_sort_parallel(points + splitpoint_val, len - splitpoint_val, sample_rate, 
                                                splitpoints, splitpoint_next_parent, false,
                                              !xaxis, current_depth + 1, max_depth);});

}

/// returns a heap
linear_kdtree lkt_create_parallel(lkt_point* points, size_t len) {
  const size_t PARALLEL_QUICKSORT_THREADS = 8;

  fprintf(stderr, "lkt_create called for points %p, true end %p \n", (void*)points, (void*)(points + len));

  if(sizeof(mortoncode_t) != 4) {
    fprintf(stderr, "mortoncode_t NOT 32 BITS! ERROR!ERROR!ERROR!"); /// \todo fix to be static_assert
    exit(1);
  }

  const unsigned short max_depth = sizeof(mortoncode_t) * CHAR_BIT;

  linear_kdtree tree;
  tree.points = points;
  tree.len    = len;

//  fprintf(stderr, "lkt_create depth %lu\n", (size_t)depth);
//  fprintf(stderr, "lkt_create split_points_len %lu\n", (size_t)tree.split_points_len);
//  fprintf(stderr, "lkt_create newing split_points size %lu\n", sizeof(lkt_split_point) * tree.split_points_len);
//  fprintf(stderr, "lkt_create newed split_points\n");

  const size_t sample_rate = 100;

//  fprintf(stderr, "lkt_create sorting\n");

  tree.split_points_len = len;

  fixlentree<lkt_split_point> splitpoints(tree.split_points_len); ///< \todo scope

  const ord_t           splitpoint       = lkt_find_splitpoint_x(points, points + len, sample_rate);
  const lkt_point       splitpoint_point = {splitpoint, splitpoint, 0}; // cheaper to just assign both, than a conditional
  const uint_least64_t  splitpoint_val   = parallel_quicksort_partition(points, &points[len], splitpoint_point, PARALLEL_QUICKSORT_THREADS, point_comparator_x);

  const lkt_split_point lkt_splitpoint = {splitpoint, splitpoint_val};
  const index_t         root           = splitpoints.insert_root(lkt_splitpoint);

//  cout << "first splitpoint x val: " << splitpoint << " location: " << splitpoint_val << endl;
  
  cout << "debug 0 parallel invoking" << endl;


  tbb::parallel_invoke([&]() {lkt_sort_parallel(points, splitpoint_val, sample_rate, 
                                                splitpoints, root, true,
                                                false, 1, max_depth);},
                       [&]() {lkt_sort_parallel(points + splitpoint_val, len - splitpoint_val, sample_rate, 
                                                splitpoints, root, false,
                                                false, 1, max_depth);});

/*
  lkt_sort_parallel(points, splitpoint_val, sample_rate, 
                    splitpoints, root, true,
                    false, 1, max_depth);
  lkt_sort_parallel(points + splitpoint_val, len - splitpoint_val, sample_rate, 
                    splitpoints, root, false,
                    false, 1, max_depth);
*/

  cout << "debug 1 splitpoints releasing" << endl;

  tree.split_points = splitpoints.release();

  fprintf(stderr, "lkt_create coding\n");

  tree.morton_codes = lkt_create_mortoncodes_parallel(tree.points, tree.len, tree.split_points, tree.split_points_len, max_depth);

//  fprintf(stderr, "lkt_create returning\n");
  return tree;
}

/// \return array of morton codes, of len length. Caller takes ownership.
mortoncode_t* lkt_create_mortoncodes_parallel(lkt_point* points, size_t len, fixlentree<lkt_split_point>::node* splitpoints, size_t splitpoints_len, size_t split_depth) {
  if(sizeof(mortoncode_t) * CHAR_BIT < split_depth) {
    fprintf(stderr, "mortoncode_t LESS THAN split_depth! ERROR!ERROR!ERROR!"); /// \todo fix to be static_assert
    exit(1);
  }

//  fprintf(stderr, "lkt_create_mortoncodes len %lu, split_points_len %lu, split_depth %lu\n", len, split_points_len, split_depth);

  mortoncode_t* codes = new mortoncode_t[len];

//  return codes; // debug !! undo !!

  for(size_t i = 0, end = len; i != end; ++i) { /// \todo vectorize
    mortoncode_t code = 0;

    const lkt_point& point = points[i];
//    fprintf(stderr, "lkt_create_mortoncodes point {%f, %f}\n", point.x, point.y);

    bool is_x = true;
    size_t code_i = 0;
    for(index_t j = 0; j != fixlentree<lkt_split_point>::tree_end;) {
//      fprintf(stderr, "j: %lu", j);
      const lkt_split_point& splitpoint = splitpoints[j].value;
      const bool             left       = is_x ? point.x < splitpoint.value : point.y < splitpoint.value;

      code = code | (left << code_i);

      j = splitpoints[j].left * left + splitpoints[j].right * !left;
      is_x = !is_x;
      ++code_i;
    }

//      fprintf(stderr, "\tj %d,\tsplitpoint_pos %ld,\tsplitpoint_val %f,\t%s,\tp_ord %f,\t%s, code %u\n", j, split_pos, split_point_val, xaxis ? "x" : "y", point_ord, left ? "left" : "right", code);
    
//    fprintf(stderr, "lkt_create_mortoncodes j %d code_size %lu \n", j, (size_t)(sizeof(mortoncode_t) * CHAR_BIT));
//    fprintf(stderr, "lkt_create_mortoncodes shifting the remaining %lu\n", (size_t)(sizeof(mortoncode_t) * CHAR_BIT - j - 1));
//    code = code << (sizeof(mortoncode_t) * CHAR_BIT - j - 1); // shift the rest of the int, so the most-significant-bit is the first split pos

//    fprintf(stderr, "lkt_create_mortoncodes code %lu: %u\n", i, code);
    codes[i] = code;
  }

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

p  rtree_node* root = (rtree_node*) malloc(sizeof(rtree_node));
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
