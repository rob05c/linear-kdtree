#include "lkt.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

const location_t location_t_max = ~0ULL;

void lqt_delete(linear_quadtree q) {
  delete[] q.locations;
  delete[] q.points;
}

/// @param points points to cona quadtree from. Takes ownership. MUST be dynamically allocated
/// @return linear quadtree. Caller takes ownership and must call lqt_delete()
linear_quadtree lqt_create(lqt_point* points, size_t len, 
             ord_t xstart, ord_t xend, 
             ord_t ystart, ord_t yend,
             size_t* depth) {
  return lqt_sortify(lqt_nodify(points, len, xstart, xend, ystart, yend, depth));
}
/* 
 * Turn an array of points into an unsorted quadtree of nodes.
 * You'll probably want to call sortify() to sort the list into a
 * useful quadtree.
 *
 * @param points points to create a quadtree from. Takes ownership. MUST be dynamically allocated
 *
 * @param[out] depth the depth of the quadtree. This is important for
 *             a linear quadtree, as it signifies the number of
 *             identifying bit-pairs preceding the node
 *
 * @return a new unsorted linear_quadtree. caller takes ownership, and must call lqt_delete()
 */
linear_quadtree lqt_nodify(lqt_point* points, size_t len, 
             ord_t xstart, ord_t xend, 
             ord_t ystart, ord_t yend,
             size_t* depth) {
  *depth = LINEAR_QUADTREE_DEPTH;

  linear_quadtree lqt;
  lqt.locations = new location_t[len];
  memset(lqt.locations, 0, sizeof(location_t) * len);
  lqt.points = points;
  lqt.length = len;

  for(size_t i = 0, end = len; i != end; ++i) {
    lqt_point* thisPoint = &lqt.points[i];

    ord_t currentXStart = xstart;
    ord_t currentXEnd = xend;
    ord_t currentYStart = ystart;
    ord_t currentYEnd = yend;
    for(size_t j = 0, jend = *depth; j != jend; ++j) {
      const location_t bit1 = thisPoint->y > (currentYStart + (currentYEnd - currentYStart) / 2);
      const location_t bit2 = thisPoint->x > (currentXStart + (currentXEnd - currentXStart) / 2);
      const location_t currentPosBits = (bit1 << 1) | bit2;
      lqt.locations[i] = (lqt.locations[i] << 2) | currentPosBits;

      const ord_t newWidth = (currentXEnd - currentXStart) / 2;
      currentXStart = floor((thisPoint->x - currentXStart) / newWidth) * newWidth + currentXStart;
      currentXEnd = currentXStart + newWidth;
      const ord_t newHeight = (currentYEnd - currentYStart) / 2;
      currentYStart = floor((thisPoint->y - currentYStart) / newHeight) * newHeight + currentYStart;
      currentYEnd = currentYStart + newHeight;
    }
  }
  return lqt;
}

struct rs_list_node {
  location_t           location;
  lqt_point     point;
  rs_list_node* next;
};
struct rs_list {
  rs_list_node* head;
  rs_list_node* tail;
};
/// @todo determine if a location pointer is faster
void rs_list_insert(rs_list* l, const location_t location, const lqt_point* point) {
  rs_list_node* n = new rs_list_node;
  n->location = location;
  n->point    = *point;
  n->next     = NULL;

  if(l->head == NULL) {
    l->head = n;
    l->tail = n;
    return;
  }
  l->tail->next = n;
  l->tail = n;
}
void rs_list_init(rs_list* l) {
  l->head = NULL;
  l->tail = NULL;
}
void rs_list_clear(rs_list* l) {
  for(rs_list_node* node = l->head; node;) {
    rs_list_node* toDelete = node;
    node = node->next;
    delete toDelete;
  }
  l->head = NULL;
  l->tail = NULL;
}

/// @todo fix this to not be global
#define BASE 10 
#define MULT_WILL_OVERFLOW(a, b, typemax) ((b) > (typemax) / (a))

// radix sort an unsorted quadtree
linear_quadtree lqt_sortify(linear_quadtree lqt) {
  rs_list buckets[BASE];
  for(int i = 0, end = BASE; i != end; ++i) 
    rs_list_init(&buckets[i]);

  const location_t max = location_t_max; ///< @todo pass max? iterate to find?

  size_t i;
  for(location_t n = 1; max / n > 0; n *= BASE) {
    // sort list of numbers into buckets
    for(i = 0; i < lqt.length; ++i) {
      const location_t location = lqt.locations[i];
      // replace array[i] in bucket_index with position code
      const size_t bucket_index = (location / n) % BASE;
      rs_list_insert(&buckets[bucket_index], location, &lqt.points[i]);
    }

    // merge buckets back into list
    for(int k = i = 0; i < BASE; rs_list_clear(&buckets[i++])) {
      for(rs_list_node* j = buckets[i].head; j != NULL; j = j->next) {
        lqt.locations[k] = j->location;
        lqt.points[k]    = j->point;
        ++k;
      }
    }
    if(MULT_WILL_OVERFLOW(n, BASE, location_t_max))
      break;
  }
  for(int i = 0, end = BASE; i != end; ++i) 
    rs_list_clear(&buckets[i]);
  return lqt;
}

/*
 * print out a quadtree node
 * @param depth the quadtree depth. Necessary, because it indicates
 *              the number of position bit-pairs
 */
void lqt_print_node(const location_t* location, const lqt_point* point, const bool verbose) {
  if(verbose)
  {
    for(int j = sizeof(location_t) * CHAR_BIT - 1, jend = 0; j >= jend; j -= 2)
      printf("%lu%lu ", (*location >> j) & 0x01, (*location >> (j - 1)) & 0x01);
    printf("%lu ", *location);
  }
  printf("%.15f\t%.15f\t%d\n", point->x, point->y, point->key);
}

/* 
 * print out all the nodes in a linear quadtree
 * @param array the linear quadtree
 * @param len the number of nodes in the quadtree
 * @param depth the depth of the quadtree.
 */
void lqt_print_nodes(linear_quadtree lqt, const bool verbose) {
  printf("linear quadtree: \n");
  if(verbose) {
    for(size_t i = 0, end = sizeof(location_t); i != end; ++i)
      printf("            ");
  }

  printf("x\ty\tkey\n");
  for(size_t i = 0, end = lqt.length; i != end; ++i) {
    lqt_print_node(&lqt.locations[i], &lqt.points[i], verbose);
  }
  printf("\n");
}

/// copies the tree from the source into destination.
/// caller takes ownership of destination, and must call delete_linear_quadtree()
/// does not delete destination, if destination is an allocated quadtree. Call delete_linear_quadtree(destination) first.
void lqt_copy(linear_quadtree* destination, linear_quadtree* source) {
  destination->length = source->length;
  destination->locations = new location_t[destination->length];
  memcpy(destination->locations, source->locations, source->length * sizeof(location_t));
  destination->points = new lqt_point[destination->length];
  memcpy(destination->points, source->points, source->length * sizeof(lqt_point));
}

///
/// unified
///

void lqt_delete_unified(linear_quadtree_unified q) {
  delete[] q.nodes;
}

#undef ENDIANSWAP

////////////////////////////////////////////////////////////////////////////////
//                              lkt                                           //
////////////////////////////////////////////////////////////////////////////////


/// \todo change these to C++ and use templates. Or Macros. Something.

/// finds a heuristic value, using the given sample rate, splitting on the x-axis
static ord_t lkt_find_splitpoint_x(lqt_point* begin, lqt_point* end, size_t sample_rate) {
  double average = 0.0;
  size_t samples_taken = 0;
  for(lqt_point* i = begin; i < end; i += sample_rate, ++samples_taken) {
    average += i->x;
  }
  average /= samples_taken;
  return average;
}
/*
/// finds a heuristic value, using the given sample rate, splitting on the y-axis
static ord_t lkt_find_splitpoint_y(lqt_point* begin, lqt_point* end, size_t sample_rate) {
  double average = 0.0;
  size_t samples_taken = 0;
  for(lqt_point* i = begin; i < end; i += sample_rate, ++samples_taken) {
    average += i->y;
  }
  average /= samples_taken;
  return average;
}
*/
/// find the next splitpoint on the y axis, ignoring values less than the given min
static ord_t lkt_find_next_splitpoint_y_l(lqt_point* begin, lqt_point* end, size_t sample_rate, ord_t min) {
  double average = 0.0;
  size_t samples_taken = 0;
  const size_t sample_distance = ((end - begin) / sample_rate) + 1;
  for(lqt_point* i = begin; i < end; i += sample_distance, ++samples_taken) {
    if(i->y < min)
      continue;
    average += i->y;
  }
  average /= samples_taken;
  return average;
}
/// find the next splitpoint on the y axis, ignoring values greater than the given max
static ord_t lkt_find_next_splitpoint_y_r(lqt_point* begin, lqt_point* end, size_t sample_rate, ord_t max) {
  double average = 0.0;
  size_t samples_taken = 0;
  const size_t sample_distance = ((end - begin) / sample_rate) + 1;
  for(lqt_point* i = begin; i < end; i += sample_distance, ++samples_taken) {
    if(i->y > max)
      continue;
    average += i->y;
  }
  average /= samples_taken;
  return average;
}
/// find the next splitpoint on the x axis, ignoring values less than the given min
static ord_t lkt_find_next_splitpoint_x_l(lqt_point* begin, lqt_point* end, size_t sample_rate, ord_t min) {
  double average = 0.0;
  size_t samples_taken = 0;
  const size_t sample_distance = ((end - begin) / sample_rate) + 1;
  for(lqt_point* i = begin; i < end; i += sample_distance, ++samples_taken) {
    if(i->x < min)
      continue;
    average += i->x;
  }
  average /= samples_taken;
  return average;
}
/// find the next splitpoint on the x axis, ignoring values greater than the given max
static ord_t lkt_find_next_splitpoint_x_r(lqt_point* begin, lqt_point* end, size_t sample_rate, ord_t max) {
  double average = 0.0;
  size_t samples_taken = 0;
  const size_t sample_distance = ((end - begin) / sample_rate) + 1;
  for(lqt_point* i = begin; i < end; i += sample_distance, ++samples_taken) {
    if(i->x > max)
      continue;
    average += i->x;
  }
  average /= samples_taken;
  return average;
}

/// DO NOT change this to use the XOR method. It is slow.
static inline void lkt_swap(lqt_point* a, lqt_point* b) {
  lqt_point old_a = *a;
  *a = *b;
  *b = old_a;
}

/// \todo make static, once tested
/// \todo replace with parallel algorithm
/// Unlike traditional quicksort partition, We don't actually use a pivot, since we only have the value to partion, not an element.
/// \param xaxis partition based on the x-axis. Else, the y-axis.
/// \return the index of the partition. Everything less than pivot_value is before the returned index. points[return] is the first element greater than pivot_value
size_t quicksort_partition(lqt_point* points, const size_t len, const ord_t pivot_value, const bool xaxis) {
//  const size_t pivot_index = points + len - 1;
//  const lqt_point* pivot_value = points[pivot_index];
//  lkt_swap(pivot_value, points[len - 1]);
//  if(len == 2)
//    return 1;
  fprintf(stderr, "quicksort_partitioned for len %lu\n", len);

  if(len == 0)
    return 0;

  long i = 0;
  long j = len - 1;

  fprintf(stderr, "quicksort_partitioned i %ld j %ld\n", i, j);

  /// we duplicate the loops rather than check inside, for efficiency
  if(xaxis) {
    while(i < j) {
      for(; points[i].x < pivot_value && i < (long)len; ++i);
      if(i >= (long)len)
        break;
      for(; points[j].x > pivot_value && j > -1; --j);
      if(j <= 0)
        break;
      fprintf(stderr, "quicksort_partitioned swapping i %ld j %ld \n", i, j);
      lkt_swap(&points[i], &points[j]);
    }
    fprintf(stderr, "quicksort_partitioned finished i %ld j %ld\n", i, j);

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
      fprintf(stderr, "quicksort_partitioned swapping i %ld j %ld\n", i, j);
      lkt_swap(&points[i], &points[j]);
    }
    fprintf(stderr, "quicksort_partitioned finished i %ld j %ld\n", i, j);

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
  fprintf(stderr, "quicksort_partitioned on %ld\n", j);
  return j;
}

static inline size_t get_heap_child_l(const size_t i) {return i * 2 + 1;}
static inline size_t get_heap_child_r(const size_t i) {return i * 2 + 2;}
static inline size_t get_heap_parent(const size_t i) {return (i - 1) / 2;}

/// \todo change splits to use copied double buffer, for parallel execution
/// \param sample_rate the rate to sample when finding the split point
static void lkt_sort(lqt_point* points, size_t len, 
                     size_t sample_rate, lkt_split_point* splitpoints, ord_t splitpoint, size_t splitpoint_i, 
                     bool xaxis, 
                     const unsigned short current_depth, const unsigned short max_depth) {
  fprintf(stderr, "lkt_sort called for splitpoint %f\n", splitpoint);
  fprintf(stderr, "lkt_sort called for points %p, len %lu\n", (void*)points, len);
  fprintf(stderr, "lkt_sort called for splitpoint_i %lu\n", splitpoint_i);
  fprintf(stderr, "lkt_sort next splitpoint_is: %lu, %lu\n", get_heap_child_l(splitpoint_i), get_heap_child_r(splitpoint_i));

  fprintf(stderr, "sort at splitpoint_i = %lu\n", splitpoint_i);
  fflush(stdout);

  if(len < 2 || current_depth == max_depth || splitpoint_i > len)
    return;

  /// \todo fix conditionals

  // splitpoint is the value in the points array, by which the points will be partitioned
  // splitpoint_i is the index into the splitpoints array, of the current split.
  // splitpoint_val is the (local) index in the points array, before which values are less than splitpoint.

  const ord_t next_split_l = xaxis ? lkt_find_next_splitpoint_x_l(points, points + len, sample_rate, splitpoint)
    : lkt_find_next_splitpoint_y_l(points, points + len, sample_rate, splitpoint); // parallel
  const ord_t next_split_r = xaxis ? lkt_find_next_splitpoint_x_r(points, points + len, sample_rate, splitpoint)
    : lkt_find_next_splitpoint_y_r(points, points + len, sample_rate, splitpoint); // parallel
  const size_t splitpoint_val = quicksort_partition(points, len, splitpoint, xaxis); // parallel

/*
  if(splitpoint_val == 0 || splitpoint_val == len) {
    fprintf(stderr, "not setting splitpoints: splitpoint_val is %lu", splitpoint_val);
    return;
  }
*/
//  fprintf(stderr, "DEBUG: setting splitpoints[%lu] = {%lu,%f}\n", splitpoint_i, splitpoint_val, splitpoint);

  fprintf(stderr, "lkt_sort [%lu:%lu) recursing [%lu:%lu) and [%lu:%lu)\n", 0ul, len,  0ul, splitpoint_val, splitpoint_val, splitpoint_val + (len - splitpoint_val));

  splitpoints[splitpoint_i].value = splitpoint;
  splitpoints[splitpoint_i].index = splitpoint_val;

  // parallel_wait

  lkt_sort(points,             splitpoint_val, sample_rate, 
           splitpoints, next_split_l, get_heap_child_l(splitpoint_i), 
           !xaxis, current_depth + 1, max_depth); // parallel
  lkt_sort(points + splitpoint_val, len - splitpoint_val, sample_rate, 
           splitpoints, next_split_r, get_heap_child_r(splitpoint_i), 
           !xaxis, current_depth + 1, max_depth); // parallel
}

/// returns a heap
linear_kdtree lkt_create(lqt_point* points, size_t len) {
  fprintf(stderr, "lkt_create called for points %p, true end %p \n", (void*)points, (void*)(points + len));

  if(sizeof(mortoncode_t) != 4) {
    fprintf(stderr, "mortoncode_t NOT 32 BITS! ERROR!ERROR!ERROR!"); /// \todo fix to be static_assert
    exit(1);
  }

  const unsigned short depth = sizeof(mortoncode_t) * CHAR_BIT;

  linear_kdtree tree;
  tree.points = points;
  tree.len = len;
  tree.split_points_len = len; // min(len, pow(2, depth)), but the length will always be less. We'd run out of memory before it wasn't.

  fprintf(stderr, "lkt_create depth %lu\n", (size_t)depth);
  fprintf(stderr, "lkt_create split_points_len %lu\n", (size_t)tree.split_points_len);
  fprintf(stderr, "lkt_create newing split_points size %lu\n", sizeof(lkt_split_point) * tree.split_points_len);
  tree.split_points = new lkt_split_point[tree.split_points_len];
  fprintf(stderr, "lkt_create newed split_points\n");
  tree.split_depth = depth;
  memset(tree.split_points, '\0', sizeof(lkt_split_point) * tree.split_points_len); // debug

  const size_t sample_rate = 100;
  const ord_t initial_splitpoint = lkt_find_splitpoint_x(points, points + len, sample_rate);
  const size_t initial_splitpoint_i = 0;

  fprintf(stderr, "lkt_create sorting\n");

  lkt_sort(points, len, sample_rate, tree.split_points, initial_splitpoint, initial_splitpoint_i, true, 0, depth);

  fprintf(stderr, "lkt_create coding\n");

  tree.morton_codes = lkt_create_mortoncodes(tree.points, tree.len, tree.split_points, tree.split_points_len, depth);
  fprintf(stderr, "lkt_create returning\n");
  return tree;
}

/// \return array of morton codes, of len length. Caller takes ownership.
mortoncode_t* lkt_create_mortoncodes(lqt_point* points, size_t len, lkt_split_point* split_points, size_t split_points_len, size_t split_depth) {
  if(sizeof(mortoncode_t) * CHAR_BIT < split_depth) {
    fprintf(stderr, "mortoncode_t LESS THAN split_depth! ERROR!ERROR!ERROR!"); /// \todo fix to be static_assert
    exit(1);
  }

  fprintf(stderr, "lkt_create_mortoncodes len %lu, split_points_len %lu, split_depth %lu\n", len, split_points_len, split_depth);

  mortoncode_t* codes = new mortoncode_t[len];

  for(size_t i = 0, end = len; i != end; ++i) { /// \todo vectorize
    mortoncode_t code = 0;
    const lqt_point point = points[i];

    fprintf(stderr, "lkt_create_mortoncodes point {%f, %f}\n", point.x, point.y);

    int j = 0;
    // in practice, this can be optimised to remove 'j < jend' because we'd run out of memory before using 2^h instead of len(points)
    for(long jend = split_depth, xaxis = true, split_pos = 0; j < jend && split_pos < (long) split_points_len; ++j, xaxis = !xaxis) {
      const lkt_split_point splitpoint = split_points[split_pos];
      const ord_t split_point_val      = splitpoint.value;
      const ord_t point_ord            = point.x * xaxis + point.y * !xaxis; // xaxis ? point.x : point.y;
      const unsigned int left          = point_ord < split_point_val;

      code = code | (left << j);

      fprintf(stderr, "\tj %d,\tsplitpoint_pos %ld,\tsplitpoint_val %f,\t%s,\tp_ord %f,\t%s, code %u\n", j, split_pos, split_point_val, xaxis ? "x" : "y", point_ord, left ? "left" : "right", code);

      split_pos = 2 * split_pos + (1 + !left); // left ? 1 : 2 (heap left child is 2*i+1, right child is 2*i+2)
    }
    
//    fprintf(stderr, "lkt_create_mortoncodes j %d code_size %lu \n", j, (size_t)(sizeof(mortoncode_t) * CHAR_BIT));
//    fprintf(stderr, "lkt_create_mortoncodes shifting the remaining %lu\n", (size_t)(sizeof(mortoncode_t) * CHAR_BIT - j - 1));
//    code = code << (sizeof(mortoncode_t) * CHAR_BIT - j - 1); // shift the rest of the int, so the most-significant-bit is the first split pos

    fprintf(stderr, "lkt_create_mortoncodes code %lu: %u\n", i, code);
    codes[i] = code;
  }

  return codes;
}
