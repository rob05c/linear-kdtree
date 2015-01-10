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

typedef uint64_t location_t;
extern const location_t location_t_max;

struct linear_quadtree {
  location_t*       locations;
  lqt_point* points;
  size_t            length;
};

#define LINEAR_QUADTREE_DEPTH (sizeof(location_t) * CHAR_BIT / 2)

linear_quadtree lqt_create(lqt_point* points, size_t len, 
                                  ord_t xstart, ord_t xend, 
                                  ord_t ystart, ord_t yend,
                                  size_t* depth);
linear_quadtree lqt_nodify(lqt_point* points, size_t len, 
                                  ord_t xstart, ord_t xend, 
                                  ord_t ystart, ord_t yend,
                                  size_t* depth);
linear_quadtree lqt_sortify(linear_quadtree);

linear_quadtree lqt_create_cuda(lqt_point* points, size_t len, 
                                       ord_t xstart, ord_t xend, 
                                       ord_t ystart, ord_t yend,
                                       size_t* depth);
linear_quadtree lqt_create_cuda_slow(lqt_point* points, size_t len, 
                                            ord_t xstart, ord_t xend, 
                                            ord_t ystart, ord_t yend,
                                            size_t* depth);
linear_quadtree lqt_nodify_cuda(lqt_point* points, size_t len, 
                                       ord_t xstart, ord_t xend, 
                                       ord_t ystart, ord_t yend,
                                       size_t* depth);
linear_quadtree lqt_sortify_cuda(linear_quadtree);

void lqt_copy(linear_quadtree* destination, linear_quadtree* source);
void lqt_delete(linear_quadtree);
void lqt_print_node(const location_t* location, const lqt_point* point, const bool verbose);
void lqt_print_nodes(linear_quadtree lqt, const bool verbose);

struct linear_quadtree_cuda {
  lqt_point*  points;
  location_t* cuda_locations;
  lqt_point*  cuda_points;
  size_t      length;
};
linear_quadtree_cuda lqt_nodify_cuda_mem(lqt_point* points, size_t len, 
                                                ord_t xstart, ord_t xend, 
                                                ord_t ystart, ord_t yend,
                                                size_t* depth);
linear_quadtree lqt_sortify_cuda_mem(linear_quadtree_cuda);


///
/// unified / heterogeneous
///
struct lqt_unified_node {
  location_t location;
  lqt_point  point;
};

struct linear_quadtree_unified {
  lqt_unified_node* nodes;
  size_t            length;
};
void lqt_delete_unified(linear_quadtree_unified);

linear_quadtree_unified lqt_nodify_cuda_unified(lqt_point* points, size_t len, 
                                                       ord_t xstart, ord_t xend, 
                                                       ord_t ystart, ord_t yend,
                                                       size_t* depth);
linear_quadtree_unified tbb_sortify_unified(linear_quadtree_unified lqt, const size_t threads);

linear_quadtree_unified lqt_create_heterogeneous(lqt_point* points, size_t len, 
                                                        ord_t xstart, ord_t xend, 
                                                        ord_t ystart, ord_t yend,
                                                        size_t* depth, const size_t threads);

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

linear_kdtree lkt_create(lqt_point* points, size_t len);
size_t quicksort_partition(lqt_point* points, const size_t len, const ord_t pivot_value, const bool xaxis);
mortoncode_t* lkt_create_mortoncodes(lqt_point* points, size_t len, lkt_split_point* split_points, size_t split_points_len, size_t split_depth);

#ifdef __cplusplus
}
#endif
#endif
