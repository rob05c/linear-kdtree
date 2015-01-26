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
#include "fixlentree.hh"

using std::cout;
using std::endl;

//#include "quicksort.hh"

// x value ALONE is used for comparison, to create an xpack
bool operator<(const lkt_point& rhs, const lkt_point& lhs) {
  return rhs.x < lhs.x;
}

std::ostream& operator<<(std::ostream& s, const lkt_point& p) {
  s << "{" << p.x << ", " << p.y << ", " << p.key << "}";
  return s;
}

std::ostream& operator<<(std::ostream& s, const lkt_split_point& p) {
  s << "{" << p.value << ", " << p.index << "}";
  return s;
}

template <typename T>
static void print_fixlentree(const fixlentree<T>& tree) {
  cout << "{" << endl;
  const typename fixlentree<T>::node* nodes = tree.get_array();
  for(int i = 0, end = tree.size(); i != end; ++i) {
    const typename fixlentree<T>::node& n = nodes[i];
    cout << "[" << n.value << ", ";

    if(n.right != fixlentree<T>::tree_end)
      cout << n.right;
    else
      cout << "-";
    cout << ", ";

    if(n.left != fixlentree<T>::tree_end)
      cout << n.left;
    else
      cout << "-";

    cout << "]";

    cout << endl;
  }
  cout << "}" << endl;
}

// generate a uniform random between min and max exclusive
static inline ord_t uniformFrand(const ord_t min, const ord_t max) {
  const double r = (double)rand() / RAND_MAX;
  return min + r * (max - min);
}

static inline void print_points(lkt_point* points, const size_t len) {
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

static inline void print_points_codes(lkt_point* points, mortoncode_t* codes, const size_t len) {
  for(int i = 0, end = len; i != end; ++i) {
    printf("(%f,%f,%d,\t", points[i].x, points[i].y, points[i].key);
    print_code(codes[i]);
    printf(")\n");
  }
}

void print_splitpoints(fixlentree<lkt_split_point>::node* nodes, size_t len) {
  cout << "{" << endl;
  for(int i = 0, end = len; i != end; ++i) {
    const typename fixlentree<lkt_split_point>::node& n = nodes[i];
    cout << "[" << n.value << ", ";
    if(n.right != fixlentree<lkt_split_point>::tree_end)
      cout << n.right;
    else
      cout << "-";
    cout << ", ";
    if(n.left != fixlentree<lkt_split_point>::tree_end)
      cout << n.left;
    else
      cout << "-";
    cout << "]";
    cout << endl;
  }
  cout << "}" << endl;
}

/// \return an array of len points. Caller takes ownership, and must call delete
static inline lkt_point* create_points(const size_t len, const ord_t min, const ord_t max) {
  lkt_point* points = new lkt_point[len];
  for(int i = 0, end = len; i != end; ++i) {
    points[i].x = uniformFrand(min, max);
    points[i].y = uniformFrand(min, max);
    points[i].key = i;
  }
  return points;
}


static inline void test_quicksort_partition(const size_t numPoints, const size_t threads) {
  printf("test_quicksort_partition\n");
  cout << "test_quicksort_partition undergoing maintenance. Be back soon!" << endl;
/*
  const size_t len = 25;
  const ord_t min = 0.0;
  const ord_t max = 10.0;

  lkt_point* points = create_points(len, min, max);

  printf("created points:\n");
  print_points(points, len);
  printf("\n\nsorting...\n");

  const ord_t pivot_val = 4;
  const size_t pivot_i = quicksort_partition(points, len, pivot_val, true);

  printf("points\n");
  print_points(points, len);
  printf("\npivot index: %lu\n", pivot_i);

  delete[] points;
*/
}

static inline void test_lkt(const size_t len, const size_t threads) {
  printf("test_lkt\n");
  cout << "test_lkt is down for maintenance" << endl;
/*
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  lkt_point* points = create_points(len, min, max);

  printf("created points:\n");
  print_points(points, len);
  printf("creating lkt:\n");

  linear_kdtree lkt = lkt_create(points, len);

  printf("\n\nlkt points:\n");
  print_points_codes(lkt.points, lkt.morton_codes, lkt.len);

  printf("\n\nlkt splitpoints:\n");
  print_splitpoints(lkt.split_points, lkt.split_points_len);
  lkt_delete(lkt);
*/
}

static inline void test_lkt_parallel(const size_t len, const size_t threads) {
  printf("test_lkt_parallel\n");
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  lkt_point* points = create_points(len, min, max);

  printf("created points:\n");
//  print_points(points, len);
  printf("creating lkt:\n");

  linear_kdtree lkt = lkt_create_parallel(points, len);

  printf("\n\nlkt points:\n");
  print_points_codes(lkt.points, lkt.morton_codes, lkt.len);

  printf("\n\nlkt splitpoints:\n");
  print_splitpoints(lkt.split_points, lkt.split_points_len);

  lkt_delete(lkt);
}

static inline void test_lkt_compare(const size_t len, const size_t threads) {
  cout << "test_lkt_compare" << endl;
  cout << "What are you looking at?" << endl;
/*
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  lkt_point* points = create_points(len, min, max);
  lkt_point* points2 = new lkt_point[len];
  memcpy(points2, points, sizeof(lkt_point) * len);

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
*/
}

static bool point_comparator_x(const lkt_point& a, const lkt_point& b) {
  return a.x < b.x;
}

static inline void test_quicksort(const size_t len, const size_t threads) {
  printf("test_lkt\n");
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  lkt_point* points = create_points(len, min, max);

  printf("created points:\n");
//  print_points(points, len);
  printf("quicksorting\n");

//  lkt_point pivot = {50.0, 50.0, ~0};
//  parallel_quicksort_partition(points, points + len, pivot, threads, point_comparator_x);

  /// \todo move quicksort.hh validity tests here

  delete[] points;
}

static inline void test_quicksort_compare(const size_t len, const size_t threads) {
  printf("test_quicksort_compare\n");
  const ord_t min = 0.0;
  const ord_t max = 100.0;

  printf("creating points:\n");
  lkt_point* points = create_points(len, min, max);
  lkt_point* points2 = new lkt_point[len];
  lkt_point pivot = {50.0, 50.0, ~0};

  memcpy(points2, points, sizeof(lkt_point) * len);
  printf("quicksorting\n");

  const auto start = std::chrono::high_resolution_clock::now();
  quicksort_partition(points, len, pivot.x, true);
  const auto end = std::chrono::high_resolution_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  delete[] points;
  cout << "sisd time (ms): " << elapsed_ms << endl;

  printf("parallel quicksorting\n");
  const auto pstart = std::chrono::high_resolution_clock::now();
  parallel_quicksort_partition(points2, &points2[len], pivot, threads, &point_comparator_x);
  const auto pend = std::chrono::high_resolution_clock::now();
  const auto pelapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(pend - pstart).count();
  delete[] points2;
  cout << "mimd time (ms): " << pelapsed_ms << endl;
}

static inline void test_fixlentree(const size_t len, const size_t threads) {
  cout << "test_fixlentree" << endl;
  if(len == 0) {
    cout << "Zero-length tree? really?" << endl;
    return;
  }

  const ord_t min = 0.0;
  const ord_t max = 100.0;

  cout << "creating points" << endl;
  lkt_point* points = create_points(len, min, max);
  cout << "created points" << endl;

  cout << "creating tree" << endl;
  fixlentree<lkt_point> tree(len);
  cout << "size: " << tree.size() << endl;

  index_t pos = tree.insert_root(points[0]);
  for(int i = 1, end = len; i != end; ++i) {
    const bool left = i % 2 == 0;
    pos = tree.insert(pos, left, points[i]);
  }

  cout << "created tree" << endl;

  print_fixlentree(tree);
}

void(*test_funcs[])(const size_t, const size_t threads) = {
  test_quicksort_partition,
  test_quicksort,
  test_quicksort_compare,
  test_lkt,
  test_lkt_parallel,
  test_lkt_compare,
  test_fixlentree,
};

static const char* default_app_name = "lkt";

const char* tests[][2] = {
  {"test_quicksort_partition", "test quicksort_partition()"},
  {"test_quicksort"          , "test parallel quicksort partitioner"},
  {"test_quicksort_compare"  , "benchmark quicksort_partition() vs parallel_quicksort_partition()"},  
  {"test_lkt"                , "test lkt_create()"},
  {"test_lkt_parallel"       , "test lkt_create_parallel()"},
  {"test_lkt_compare"        , "benchmark lkt_create() vs lkt_create_parallel()"},
  {"test_fixlentree"        , "test fixlentree structure"},
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
  //  const time_t now = time(NULL);
  const time_t now = 1422251841; // debug freeze at 10 000 000
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
