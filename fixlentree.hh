#ifndef fixlentree_hh
#define fixlentree_hh

#include <cstdlib>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <atomic>
#include "tbb/tbb.h"
#include <iostream> // debug
#include <chrono> // debug

using std::cout; // debug!
using std::endl; // debug!

typedef size_t index_t;

template <typename T>
class fixlentree {
public:
  static const index_t tree_end = std::numeric_limits<index_t>::max();

  struct node {
//    static const index_t tree_end = fixlentree<T>::tree_end;
    T       value;
    index_t right;
    index_t left;
  };

  fixlentree(const size_t len) 
    : len_(len)
    , array_(new node[len])
    , node_acquirer_(new std::atomic<index_t>(0)) 
  {}

  /// result is undefined if the root already exists.
  /// result is undefined if release() has been called.
  /// \return the index of the root.
  index_t insert_root(const T& val) {
    const index_t pos = node_acquirer_->fetch_add(1);
    node* node = &array_[pos];
    node->value = val;
    node->right = tree_end;
    node->left  = tree_end;
    return pos;
  }

  /// result is undefined if parent has not been previously returned by insert_root or insert.
  /// \return the index of the newly created child, or tree_end if fixlentree is full.
  index_t insert(const index_t parent, const bool left, const T& val) {
    const index_t pos = node_acquirer_->fetch_add(1);
    if(pos >= len_)
      return tree_end;
    
    node* node = &array_[pos];
    node->value = val;
    node->right = tree_end;
    node->left  = tree_end;

    /// \todo replace conditional with arithmetic?
    if(left)
      array_[parent].left = pos;
    else
      array_[parent].right = pos;

    return pos;
  }
  
  const node& get(const index_t i) const   {return &array_[i];} ///< result is undefined if i has not been previously returned by insert_root or insert.
  node        get_mutable(const index_t i) {return array_[i];}
  const node* get_array() const            {return array_.get();}
  index_t     size()      const            {return len_;}
  node*       release()                    {return array_.release();} ///< Release and return the array pointer. Caller takes ownership. Tree is rendered useless. The result of all subsequent member function calls is undefined.
  
private:
  const index_t                         len_;
  std::unique_ptr<node[]>               array_;
  std::unique_ptr<std::atomic<index_t>> node_acquirer_; ///< this is a pointer, so the tree is copyable.
};


#endif // fixlentree_hh
