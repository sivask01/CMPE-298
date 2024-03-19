/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once

#include <vector>

#include <fmt/format.h>

#include <fasta/types.h>
#include <fasta/Matrix.h>
#include <fasta/Index.h>

namespace fasta {

  using namespace std;

  /**
   * GEMM Index class for database of vectors to be searched.
   */
  template <typename T>
  struct GemmIndex : public Index<T> {
    Matrix<T> * data;

    explicit GemmIndex(const idx_t ncols)
     : Index<T>(ncols)
    {
      data = new Matrix<T>(0, ncols);
    }

    explicit GemmIndex(
      const idx_t ncols,
      const idx_t n,
      const T* x
    )
     : Index<T>(ncols)
    {
      data = new Matrix<T>(0, ncols);
      add(n, x);
    }

    explicit GemmIndex(const Matrix<T>& x)
     : Index<T>(x.ncols)
    {
      data = new Matrix<T>(x);
      this->nsamples = x.nrows;
    }

    ~GemmIndex()
    {
      if(data){
        delete data;
      }
    }

    /**
     * Add a set of vectors to the index. IDs are assumed 0 -> n-1.
     *
     * @param n      number of vectors being added
     * @param x      sample vecors, size n * ncols
     */
    void add(
      const idx_t n,
      const T* x
    ) override
    {
      if(this->is_trained){
        throw runtime_error("Index has already been trained. No new data may be added after training.");
      }
      data->add(n, x);
      this->nsamples = data->nrows;
    }

    /**
     * Add a set of vectors to the index. IDs are assumed 0 -> n-1.
     *
     * @param x      sample vecors
     */
    void add(const Matrix<T>& x)
    {
      if(this->is_trained){
        throw runtime_error("Index has already been trained. No new data may be added after training.");
      }
      assert(this->ncols == x.ncols);
      data->add(x);
      this->nsamples = data->nrows;
    }

    /**
     * Add a set of vectors to the index.
     *
     * @param n      number of vectors being added
     * @param x      sample vecors, size n * ncols
     * @param add_ids    IDs of sample vecors, size n
     */
    void add(
      const idx_t n,
      const T* x,
      const idx_t* add_ids
    ) override
    {
      if(this->is_trained){
        throw runtime_error("Index has already been trained. No new data may be added after training.");
      }
      idx_t nr = data->nrows;
      idx_t mr = data->max_nrows;
      data->add(n, x);
      if(!this->ids){
        this->ids = (idx_t*) malloc(sizeof(idx_t) * data->get_max_nnz());
        if(!this->ids){
          throw runtime_error(fmt::format("Could not allocate {} elements for storing ids.", data->get_max_nnz()));
        }
      } else if(data->max_nrows > nr){
        this->ids = (idx_t*) realloc(this->ids, sizeof(idx_t) * data->get_max_nnz());
        if(!this->ids){
          throw runtime_error(fmt::format("Could not allocate {} elements for storing ids.", data->get_max_nnz()));
        }
      }
      memcpy((void*) (this->ids + nr), (void*) add_ids, sizeof(idx_t) * n);
      this->nsamples = data->nrows;
    }

    /**
     * Perform training on vectors that have already been added.
     */
    void train() override
    {
      if(this->is_trained){
        return;
      }
      if(!data->is_normalized()){
        data->normalize();
      }
      this->is_trained = true;
    }

    /**
     * Write the index to a file.
     * @param fpath Path where index should be written
     */
    void write(const string fpath) override
    {
      // write the data matrix to file
      // TODO: should also store is_trained status
      data->write(fpath, 'b');
    }

    /**
     * Load index data from file.
     * @param fpath Path where index resides on disk.
     */
    void read(const string fpath) override
    {
      if(data){
        delete data;
      }
      data = Matrix<T>::read(fpath, 'b');
      this->is_trained = false;
    }

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
    void search(
        const idx_t n,
        const T* x,
        const float epsilon,
        csr_t& result) override;

    /**
     * Find `k` nearest neighbors within the index for `n` vectors of dimension `ncols`.
     *
     * @param n           number of vectors to search
     * @param x           input vectors to search, size n * d
     * @param k           number of nearest neighbors to retrieve for each sample
     * @param nbrs        nearest neighbor ids table, of size n * k
     * @param sims        nearest neighbor similarities table, of size n * k
     */
    void knn_search(
        const idx_t n,
        const T* x,
        const idx_t k,
        idx_t* nbrs,
        val_t* sims) override;

  }; // class GemmIndex



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
  template <typename T>
  void GemmIndex<T>::search(
      const idx_t n,
      const T* x,
      const float epsilon,
      csr_t& result
  ){
    if(!this->is_trained){
      train();
    }

    this->timers[QUERY_PREPROCESS].start();
    Matrix<T> qmat(n, this->ncols, x);
    if(!qmat.is_normalized()){
      qmat.normalize();
    }
    this->timers[QUERY_PREPROCESS].stop();

    vector<vector<sim_t>> results(n, vector<sim_t>());
    ptr_t nnz = 0;
    #pragma omp parallel for reduction(+: nnz)
    for(idx_t i=0; i < n; ++i){
      // compute similarities
      vector<sim_t> r;
      r.reserve(this->nsamples);
      for(idx_t j=0; j < this->nsamples; ++j){
        acm_t v = 0.0;
        for(idx_t k=0; k < this->ncols; ++k){
          v += qmat(i,k) * data->operator()(j,k);
        }
        r.push_back((sim_t){j, v});
      }
      // filter scores below threshold
      idx_t k = 0;
      for(idx_t j=0; j < this->nsamples; ++j){
        if(r[j].score >= epsilon){
          r[k++] = r[j];
        }
      }
      r.resize(k);
      nnz += k;
      results[i] = r;
    }

    // transfer results to results structure
    result.reserve(n, data->nrows, nnz);
    for(idx_t i=0; i < n; ++i){
      result.ptr[i+1] = result.ptr[i] + results[i].size();
    }
    #pragma omp parallel for schedule(static, 20)
    for(idx_t i=0; i < n; ++i){
      ptr_t p = result.ptr[i];
      for(idx_t j=0; j < results[i].size(); ++j, ++p){
        result.ind[p] = results[i][j].cid;
        result.val[p] = results[i][j].score;
      }
    }
    set_result_ids(this->ids, result);
    this->counters[NUM_CANDIDATES] = (uint64_t) qmat.nrows * data->nrows;
    this->counters[DOT_PRODUCTS] = (uint64_t) qmat.nrows * data->nrows;
    this->counters[RESULTS] = (uint64_t) result.ptr[result.nrows];
  }


  /**
   * Find `k` nearest neighbors within the index for `n` vectors of dimension `ncols`.
   *
   * @param n           number of vectors to search
   * @param x           input vectors to search, size n * d
   * @param k           number of nearest neighbors to retrieve for each sample
   * @param nbrs        nearest neighbor ids table, of size n * k
   * @param sims        nearest neighbor similarities table, of size n * k
   */
  template <typename T>
  void GemmIndex<T>::knn_search(
      const idx_t n,
      const T* x,
      const idx_t k,
      idx_t* nbrs,
      val_t* sims
  ){
    assert(k <= this->nsamples);
    if(!this->is_trained){
      train();
    }

    this->timers[QUERY_PREPROCESS].start();
    Matrix<T> qmat(n, this->ncols, x);
    if(!qmat.is_normalized()){
      qmat.normalize();
    }
    this->timers[QUERY_PREPROCESS].stop();

    #pragma omp parallel
    {
      vector<sim_t> r;
      r.reserve(this->nsamples);
      #pragma omp for schedule(static, 20)
      for(idx_t i=0; i < n; ++i){
        // compute similarities
        for(idx_t j=0; j < this->nsamples; ++j){
          acm_t v = 0.0;
          for(idx_t l=0; l < this->ncols; ++l){
            v += qmat(i,l) * data->operator()(j,l);
          }
          r.push_back((sim_t){j, v});
        }
        // sort
        sort(
          begin(r),
          end(r),
          [](sim_t a, sim_t b) {
            return a.score > b.score;
          }
        );
        // transfer top k
        if(nbrs){
          for(idx_t j=0; j < k; ++j){
            nbrs[i*k + j] = r[j].cid;
          }
        }
        if(sims){
          for(idx_t j=0; j < k; ++j){
            sims[i*k + j] = r[j].score;
          }
        }
        // prepare for next query
        r.clear();
      }
    }
    set_result_ids(this->ids, n, k, nbrs);
    this->counters[NUM_CANDIDATES] = (uint64_t) qmat.nrows * data->nrows;
    this->counters[DOT_PRODUCTS] = (uint64_t) qmat.nrows * data->nrows;
    this->counters[RESULTS] = (uint64_t) qmat.nrows * k;
  }


} // namespace fasta
