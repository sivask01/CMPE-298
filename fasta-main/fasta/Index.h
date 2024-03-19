/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once

#include <fasta/types.h>
#include <fasta/Matrix.h>
#include <fasta/utils/Timers.h>
#include <fasta/utils/Counters.h>

namespace fasta {

  using namespace std;
  using namespace timer;
  using namespace counter;

  /**
   * Base Index class for database of vectors to be searched.
   */
  template <typename T>
  struct Index {
    idx_t nsamples;                  // total number of samples in the index (may be in different structures)
    idx_t ncols;                     // the dimensionality of the samples
    idx_t* ids;                      // IDs associated with each sample (if null, IDs are 0-->n-1).
    bool is_trained;                 // whether indexed samples have been processed
    Timer timers[NTIMERS];           // timers
    Counter counters[NCOUNTERS];     // counters

    explicit Index(const idx_t ncols)
      : ncols(ncols),
        nsamples(0),
        is_trained(false)
    {
      ids = nullptr;
    }

    explicit Index(
      const idx_t ncols,
      const idx_t n,
      const T* x
    ) : ncols(ncols),
        nsamples(0),
        is_trained(false)
    {
      ids = nullptr;
      add(n, x);
    }

    explicit Index(const Matrix<T>& x)
      : ncols(x.ncols),
        nsamples(0),
        is_trained(false)
    {
      ids = nullptr;
      add(x);
    }

    virtual ~Index(){
      if(ids){
        free(ids);
      }
    }

    /**
     * Add a set of vectors to the index. IDs are assumed 0 -> n-1.
     *
     * @param n      number of vectors being added
     * @param x      sample vecors, size n * ncols
     */
    virtual void add(const idx_t n, const T* x) = 0;

    /**
     * Add a set of vectors to the index. IDs are assumed 0 -> n-1.
     *
     * @param x      sample vecors
     */
    void add(const Matrix<T>& x);

    /**
     * Add a set of vectors to the index.
     *
     * @param n      number of vectors being added
     * @param x      sample vecors, size n * ncols
     * @param ids    IDs of sample vecors, size n
     */
    virtual void add(const idx_t n, const T* x, const idx_t* ids) = 0;

    /**
     * Perform training on vectors that have already been added.
     */
    virtual void train() = 0;

    /**
     * Write the index to a file.
     * @param fpath Path where index should be written
     */
    virtual void write(const string fpath) = 0;

    /**
     * Load index data from file.
     * @param fpath Path where index resides on disk.
     */
    virtual void read(const string fpath) = 0;

    /**
     * Find neighbors with min similarity `epsilon` within the index for `n` vectors of dimension `ncols`.
     *
     * Return all vectors with similality >= epsilon.
     *
     * @param n           number of vectors to search
     * @param x           input vectors to search, size n * d
     * @param epsilon     search radius / minimum similarity
     * @param result      result table (sparse matrix)
     */
    virtual void search(
            const idx_t n,
            const T* x,
            const float epsilon,
            csr_t& result) = 0;

    /**
     * Find `k` nearest neighbors within the index for `n` vectors of dimension `ncols`.
     *
     * @param n           number of vectors to search
     * @param x           input vectors to search, size n * d
     * @param k           number of nearest neighbors to retrieve for each sample
     * @param nbrs        nearest neighbor ids table, of size n * k
     * @param sims        nearest neighbor similarities table, of size n * k
     */
    virtual void knn_search(
            const idx_t n,
            const T* x,
            const idx_t k,
            idx_t* nbrs,
            val_t* sims) = 0;

    /**
     * Find `k` nearest neighbors within the index for `n` vectors of dimension `ncols`.
     *
     * @param n           number of vectors to search
     * @param x           input vectors to search, size n * d
     * @param k           number of nearest neighbors to retrieve for each sample
     * @param nbrs        nearest neighbor table, of size n * k
     */
    void knn_search(
        const idx_t n,
        const T* x,
        const idx_t k,
        Matrix<sim_t>& nbrs);

  };  // struct Index

  /**
   * Replace row indexes with IDs in the result if `ids` are present in the index.
   * @param result Result structure
   */
  void set_result_ids(const idx_t* ids, csr_t& result);

  /**
   * Replace row indexes with IDs in the result if `ids` are present in the index.
   * @param n        number of vectors to search
   * @param k        number of nearest neighbors to retrieve for each sample
   * @param result   nearest neighbor ids table, of size n * k
   */
  void set_result_ids(const idx_t* ids, const idx_t n, const idx_t k, idx_t* result);

} // namespace fasta

#include <fasta/impl/Index.h>
