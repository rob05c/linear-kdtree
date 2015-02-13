#ifndef lqtH
#define lqtH
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include <vector>
#include <utility>

#include "fixlentree.hh"

typedef int key_t;
typedef float ord_t; // ordinate. There's only one, so it's not a coordinate.
struct lkt_point {
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
  lkt_point*        points;
  size_t            len;
  fixlentree<lkt_split_point>::node* split_points;
  size_t            split_points_len;
  size_t            split_depth;
  mortoncode_t*     morton_codes;
};

void lkt_delete(linear_kdtree tree);
linear_kdtree lkt_create(lkt_point* points, size_t len);
size_t quicksort_partition(lkt_point* points, const size_t len, const ord_t pivot_value, const bool xaxis);
mortoncode_t* lkt_create_mortoncodes_sisd(lkt_point* points, size_t len, const fixlentree<lkt_split_point>::node* splitpoints);
mortoncode_t* lkt_create_mortoncodes_simd(lkt_point* points, size_t len, const fixlentree<lkt_split_point>::node* split_points);
mortoncode_t* lkt_create_mortoncodes_mimd(lkt_point* points, size_t len, const fixlentree<lkt_split_point>::node* split_points);

linear_kdtree lkt_create_mimd(lkt_point* points, size_t len);
linear_kdtree lkt_create_heterogeneous(lkt_point* points, size_t len);
linear_kdtree lkt_create_mimd_codeless(lkt_point* points, size_t len);

std::vector<linear_kdtree> lkt_create_pipelined(std::vector<std::pair<lkt_point*, size_t>> pointses);

#endif
