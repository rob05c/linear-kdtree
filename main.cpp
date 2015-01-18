#include "lkt.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <chrono>
#include <tbb/tbb.h>
#include "quicksort.hh"

using std::cout;
using std::endl;

//#include "quicksort.hh"

// x value ALONE is used for comparison, to create an xpack
bool operator<(const lqt_point& rhs, const lqt_point& lhs) {
  return rhs.x < lhs.x;
}
/*
std::ostream& operator<<(std::ostream& s, const lqt_point& p) {
  s << "{" << p.x << ", " << p.y << ", " << p.key << "}";
  return s;
}
*/
// generate a uniform random between min and max exclusive
static inline ord_t uniformFrand(const ord_t min, const ord_t max) {
  const double r = (double)rand() / RAND_MAX;
  return min + r * (max - min);
}

static inline void print_points(lqt_point* points, const size_t len) {
  for(int i = 0, end = len; i != end; ++i) {
    printf("(%f,%f,%d)\n", points[i].x, points[i].y, points[i].key);
  }
}

static inline void print_code(mortoncode_t c) {
  mortoncode_t one = 1;
  printf("[");
  for(mortoncode_t i = 0, end = sizeof(c) * CHAR_BIT; i != end; ++i) {
    printf("%d", (c & (one << i)) != 0);
  }
  printf("] (%u)", c);
}

static inline void print_points_codes(lqt_point* points, mortoncode_t* codes, const size_t len) {
  for(int i = 0, end = len; i != end; ++i) {
    printf("(%f,%f,%d,\t", points[i].x, points[i].y, points[i].key);
    print_code(codes[i]);
    printf(")\n");
  }
}

static inline void print_splitpoints(lkt_split_point* points, const size_t len) {
  printf("[\n");
  for(int i = 0, end = len; i != end; ++i)
    printf("{%d, %lu, %f}\n", i, points[i].index, points[i].value);
  printf("]\n");
}

/// \return an array of len points. Caller takes ownership, and must call delete
static inline lqt_point* create_points(const size_t len, const ord_t min, const ord_t max) {
  lqt_point* points = new lqt_point[len];
  for(int i = 0, end = len; i != end; ++i) {
    points[i].x = uniformFrand(min, max);
    points[i].y = uniformFrand(min, max);
    points[i].key = i;
  }
  return points;
}


static inline void test_quicksort_partition(const size_t numPoints, const size_t threads) {
  printf("test_quicksort_partition\n");
  const size_t len = 25;
  const ord_t min = 0.0;
  const ord_t max = 10.0;

  lqt_point* points = create_points(len, min, max);

  printf("created points:\n");
  print_points(points, len);
  printf("\n\nsorting...\n");

  const ord_t pivot_val = 4;
  const size_t pivot_i = quicksort_partition(points, len, pivot_val, true);

  printf("points\n");
  print_points(points, len);
  printf("\npivot index: %lu\n", pivot_i);

  delete[] points;
}

static inline void test_lkt(const size_t len, const size_t threads) {
  printf("test_lkt\n");
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  lqt_point* points = create_points(len, min, max);

  printf("created points:\n");
  print_points(points, len);
  printf("creating lkt:\n");

  linear_kdtree lkt = lkt_create(points, len);

  printf("\n\nlkt points:\n");
  print_points_codes(lkt.points, lkt.morton_codes, lkt.len);

  printf("\n\nlkt splitpoints:\n");
  print_splitpoints(lkt.split_points, lkt.split_points_len);
  lkt_delete(lkt);
}

static inline void test_lkt_parallel(const size_t len, const size_t threads) {
  printf("test_lkt_parallel\n");
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  lqt_point* points = create_points(len, min, max);

  printf("created points:\n");
  print_points(points, len);
  printf("creating lkt:\n");

  linear_kdtree lkt = lkt_create_parallel(points, len);

  printf("\n\nlkt points:\n");
  print_points_codes(lkt.points, lkt.morton_codes, lkt.len);

  fflush(stdout);
  printf("\n\nlkt splitpoints:\n");
  print_splitpoints(lkt.split_points, lkt.split_points_len);
  printf("\n\nlkt splitpoints:\n");
  print_splitpoints(lkt.split_points, lkt.split_points_len);
  printf("FIN\n");
  lkt_delete(lkt);
}

static inline void test_lkt_compare(const size_t len, const size_t threads) {
  cout << "test_lkt_compare\n" << endl;
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  lqt_point* points = create_points(len, min, max);
  lqt_point* points2 = new lqt_point[len];
  memcpy(points2, points, sizeof(lqt_point) * len);

  cout << "created points:" << endl;

  cout << "creating lkt:" << endl;
  const auto start = std::chrono::high_resolution_clock::now();
  linear_kdtree lkt = lkt_create(points, len);
  const auto end = std::chrono::high_resolution_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  cout << "sisd time (ms): " << elapsed_ms << endl;
  lkt_delete(lkt);

  cout << "creating parallel lkt:" << endl;
  const auto pstart = std::chrono::high_resolution_clock::now();
  linear_kdtree plkt = lkt_create_parallel(points2, len);
  const auto pend = std::chrono::high_resolution_clock::now();
  const auto pelapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(pend - pstart).count();
  lkt_delete(plkt);

  cout << "mimd time (ms): " << pelapsed_ms << endl;

  printf("\n");

}



static bool point_comparator_x(const lqt_point& a, const lqt_point& b) {
  return a.x < b.x;
}

static inline void test_quicksort(const size_t len, const size_t threads) {
  printf("test_lkt\n");
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  lqt_point* points = create_points(len, min, max);

  printf("created points:\n");
//  print_points(points, len);
  printf("quicksorting\n");

//  lqt_point pivot = {50.0, 50.0, ~0};
//  parallel_quicksort_partition(points, points + len, pivot, threads, point_comparator_x);

  /// \todo move quicksort.hh validity tests here

  delete[] points;
}

static inline void test_quicksort_compare(const size_t len, const size_t threads) {
  printf("test_quicksort_compare\n");
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  printf("creating points:\n");
  lqt_point* points = create_points(len, min, max);
  lqt_point* points2 = new lqt_point[len];
  lqt_point pivot = {50.0, 50.0, ~0};
  uint_least64_t splitpoint_val;
  memcpy(points2, points, sizeof(lqt_point) * len);
  printf("quicksorting\n");

  const auto start = std::chrono::high_resolution_clock::now();
  quicksort_partition(points, len, pivot.x, true);
  const auto end = std::chrono::high_resolution_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  delete[] points;
  cout << "sisd time (ms): " << elapsed_ms << endl;

  printf("parallel quicksorting\n");
  const auto pstart = std::chrono::high_resolution_clock::now();
  parallel_quicksort_partition(&splitpoint_val, points2, &points2[len], pivot, threads, &point_comparator_x);
  const auto pend = std::chrono::high_resolution_clock::now();
  const auto pelapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(pend - pstart).count();
  delete[] points2;
  cout << "mimd time (ms): " << pelapsed_ms << endl;
}

void(*test_funcs[])(const size_t, const size_t threads) = {
  test_quicksort_partition,
  test_quicksort,
  test_quicksort_compare,
  test_lkt,
  test_lkt_parallel,
  test_lkt_compare,
};

static const char* default_app_name = "lkt";

const char* tests[][2] = {
  {"test_quicksort_partition", "test quicksort_partition()"},
  {"test_quicksort"          , "test parallel quicksort partitioner"},
  {"test_quicksort_compare"  , "benchmark quicksort_partition() vs parallel_quicksort_partition()"},  
  {"test_lkt"                , "test lkt_create()"},
  {"test_lkt_parallel"       , "test lkt_create_parallel()"},
  {"test_lkt_compare"        , "benchmark lkt_create() vs lkt_create_parallel()"},  
};

const size_t test_num = sizeof(tests) / (sizeof(const char*) * 2);

struct app_arguments {
  bool        success;
  const char* app_name;
  size_t      test_num;
  size_t      array_size;
  size_t      threads;
};

static app_arguments parseArgs(const int argc, const char** argv) {
  app_arguments args;
  args.success = false;

  if(argc < 1)
    return args;
  args.app_name = argv[0];

  if(argc < 2)
    return args;
  args.test_num = strtol(argv[1], NULL, 10);

  if(argc < 3)
    return args;
  args.array_size = strtol(argv[2], NULL, 10);

  if(argc < 4)
    return args;
  args.threads = strtol(argv[3], NULL, 10);

  args.success = true;
  return args;
}

/// \param[out] msg
/// \param[out] msg_len
static void print_usage(const char* app_name) {
  printf("usage: %s test_num array_size threads\n", strlen(app_name) == 0 ? default_app_name : app_name);
  printf(" (threads is only used for heterogeneous test(s)\n");
  printf("\n");
  printf("       num test            description\n");
  for(size_t i = 0, end = test_num; i != end; ++i) {
    printf("       %-3.1lu %-15.15s %s\n", i, tests[i][0], tests[i][1]);
  }
  printf("\n");
}

int main(const int argc, const char** argv) {
  const time_t now = time(NULL);
//  const time_t now = 1420954039; // debug fail with 100 elements!!
  cout << "seed: " << now << endl;
  srand(now); // constant seed for debugging

  const app_arguments args = parseArgs(argc, argv);
  if(!args.success) {
    print_usage(args.app_name);
    return 0;
  }

  tbb::task_scheduler_init init(args.threads); // +1 for the manager thread

  test_funcs[args.test_num](args.array_size, args.threads);
  printf("\n");
  return 0;
}
