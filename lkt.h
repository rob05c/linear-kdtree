#ifndef lqtH
#define lqtH
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <limits.h>

// nvcc is C++, not C
#ifdef __cplusplus
extern "C" {
#endif

typedef int key_t;
typedef float ord_t; // ordinate. There's only one, so it's not a coordinate.
struct lqt_point {
  ord_t x;
  ord_t y;
  key_t key;
};

typedef unsigned int mortoncode_t;
struct lkt_split_point {
  ord_t  value;
  size_t index;
};
struct linear_kdtree {
  lqt_point*        points;
  size_t            len;
  lkt_split_point*  split_points;
  size_t            split_points_len;
  size_t            split_depth;
  mortoncode_t*     morton_codes;
};

void lkt_delete(linear_kdtree tree);
linear_kdtree lkt_create(lqt_point* points, size_t len);
size_t quicksort_partition(lqt_point* points, const size_t len, const ord_t pivot_value, const bool xaxis);
mortoncode_t* lkt_create_mortoncodes(lqt_point* points, size_t len, lkt_split_point* split_points, size_t split_points_len, size_t split_depth);

mortoncode_t* lkt_create_mortoncodes_parallel(lqt_point* points, size_t len, lkt_split_point* split_points, size_t split_points_len, size_t split_depth);

linear_kdtree lkt_create_parallel(lqt_point* points, size_t len);

#ifdef __cplusplus
}
#endif
#endif
