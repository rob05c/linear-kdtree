/// this file exists to remove C++11 from CUDA, to support outdated nvcc compilers
#include "lkt.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>    
#include <assert.h>
#include <math.h>
#include <future>
#include "tbb/tbb.h"
#include "quicksort.hh"

using std::vector;
using std::pair;
using std::promise;
using std::future;

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
static void lkt_sort_mimd(lkt_point* points, size_t len,
                              size_t sample_rate, fixlentree<lkt_split_point>& splitpoints, index_t splitpoint_parent, bool splitpoint_thisisleft,
                              bool xaxis, 
                              const unsigned short current_depth, const unsigned short max_depth) {
  const size_t PARALLEL_QUICKSORT_THREADS = 8;
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

  const uint_least64_t splitpoint_val   = parallel_quicksort_partition(points, points + len, splitpoint_point, PARALLEL_QUICKSORT_THREADS, comparator_func); ///< \todo fix last (?) block bug

  const lkt_split_point lkt_splitpoint         = {splitpoint, splitpoint_val};

  const index_t         splitpoint_next_parent = splitpoints.insert(splitpoint_parent, splitpoint_thisisleft, lkt_splitpoint);

  if(splitpoint_next_parent == splitpoints.tree_end)
    return;
  if(splitpoint_val == 0 || splitpoint_val == len - 1)
    return;

  tbb::parallel_invoke([&]() {lkt_sort_mimd(points, splitpoint_val, sample_rate, 
                                                splitpoints, splitpoint_next_parent, true,
                                                !xaxis, current_depth + 1, max_depth);},
                       [&]() {lkt_sort_mimd(points + splitpoint_val, len - splitpoint_val, sample_rate, 
                                                splitpoints, splitpoint_next_parent, false,
                                              !xaxis, current_depth + 1, max_depth);});
}

static linear_kdtree lkt_create_mimd_codeless(lkt_point* points, size_t len) {
  static_assert(sizeof(mortoncode_t) == 4, "mortoncode_t must be 32 bits");
  const size_t         PARALLEL_QUICKSORT_THREADS = 8;
  const unsigned short max_depth                  = sizeof(mortoncode_t) * CHAR_BIT;
  const size_t         sample_rate                = 100;

  linear_kdtree tree;
  tree.points           = points;
  tree.len              = len;
  tree.split_points_len = len;

  fixlentree<lkt_split_point> splitpoints(tree.split_points_len); ///< \todo scope

  const ord_t           splitpoint       = lkt_find_splitpoint_x(points, points + len, sample_rate);
  const lkt_point       splitpoint_point = {splitpoint, splitpoint, 0}; // cheaper to just assign both, than a conditional
  const uint_least64_t  splitpoint_val   = parallel_quicksort_partition(points, &points[len], splitpoint_point, PARALLEL_QUICKSORT_THREADS, point_comparator_x);

  const lkt_split_point lkt_splitpoint = {splitpoint, splitpoint_val};
  const index_t         root           = splitpoints.insert_root(lkt_splitpoint);

  tbb::parallel_invoke([&]() {lkt_sort_mimd(points, splitpoint_val, sample_rate, 
                                                splitpoints, root, true,
                                                false, 1, max_depth);},
                       [&]() {lkt_sort_mimd(points + splitpoint_val, len - splitpoint_val, sample_rate, 
                                                splitpoints, root, false,
                                                false, 1, max_depth);});
  tree.split_points = splitpoints.release();
  return tree;
}

/// \return array of morton codes, of len length. Caller takes ownership.
mortoncode_t* lkt_create_mortoncodes_mimd(lkt_point* points, size_t len, const fixlentree<lkt_split_point>::node* splitpoints) {
  mortoncode_t* codes = new mortoncode_t[len];

  tbb::parallel_for(size_t(0), len, size_t(1), [&](size_t i) {
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
    });

  return codes;
}

linear_kdtree lkt_create_heterogeneous(lkt_point* points, size_t len) {
  linear_kdtree tree = lkt_create_mimd_codeless(points, len);
  tree.morton_codes = lkt_create_mortoncodes_simd(tree.points, tree.len, tree.split_points);
  return tree;
}

linear_kdtree lkt_create_mimd(lkt_point* points, size_t len) {
  linear_kdtree tree = lkt_create_mimd_codeless(points, len);
  tree.morton_codes = lkt_create_mortoncodes_mimd(tree.points, tree.len, tree.split_points);
  return tree;
}


/// pipeline

vector<linear_kdtree> lkt_create_pipelined(vector<pair<lkt_point*, size_t>> pointses) {
  vector<promise<bool>> promises(pointses.size());
  vector<future<bool>> futures;
  for(vector<promise<bool>>::iterator i = promises.begin(), end = promises.end(); i != end; ++i)
    futures.push_back(i->get_future());

  vector<linear_kdtree> trees;

  tbb::task_group tasks;
  tasks.run([&]{
      for(size_t i = 0, end = futures.size(); i != end; ++i) {
        futures[i].wait();
        linear_kdtree& tree = trees[i];
        tree.morton_codes = lkt_create_mortoncodes_simd(tree.points, tree.len, tree.split_points);
      }
    });

  for(size_t i = 0, end =  pointses.size(); i != end; ++i) {
    trees.push_back(lkt_create_mimd_codeless(pointses[i].first, pointses[i].second));
    promises[i].set_value(true);
  }
  tasks.wait();
  return trees;
}
