/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once

#include <vector>
#include <assert.h>

#include <fmt/format.h>

#include <fasta/types.h>
#include <fasta/Matrix.h>
#include <fasta/Index.h>
#include <fasta/utils/BitVector.h>
#include <fasta/utils/Counters.h>
#include <fasta/utils/Timers.h>


namespace fasta {

  using namespace std;
  using namespace counter;
  using namespace timer;

  /**
   * Extended CSR structure to store sparse inverted indexes
   */
  typedef struct l2_index_t : public csr_t {
    val_t * nrm = nullptr;    // prefix or suffix norms, size nnz
    val_t * max = nullptr;    // row (columns in data) maximum absolute values, size nrows
    idx_t start = 0;          // first data row in this index
    idx_t end = 0;            // row after last data row in this index

    l2_index_t() : csr_t(){}

    /**
     * Reserve space for more rows or non-zeros. Structure may only grow, not shrink.
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param nnz   Number of non-zeros
     */
    void reserve(const idx_t nrows, const idx_t ncols, const ptr_t nnz)
    {
      ptr_t cnnz = ptr ? ptr[this->nrows] : 0; // current number of non-zeros
      csr_t::reserve(nrows, ncols, nnz);
      if(nnz > cnnz){
        if(nrm){
          nrm = (val_t*) realloc(nrm, sizeof(val_t) * nnz);
        } else {
          nrm = (val_t*) malloc(sizeof(val_t) * nnz);
        }
        if(!nrm){
          throw std::runtime_error("Could not allocate nrm array.");
        }
      }
      if(!max){
        max = (val_t*) malloc(sizeof(val_t) * ncols);
      }
    }

    /**
     * Reserve space for more rows or non-zeros. Structure may only grow, not shrink.
     * @param nrows Number of rows
     * @param nnz   Number of non-zeros
     */
    void reserve(const idx_t nrows, const ptr_t nnz)
    {
      throw std::runtime_error("The reserve method requires ncols for l2_index_t. Use l2_index_t::reserve(nrows, ncols, nnz).");
    }

    /**
     * Write matrix to a binary file
     * @param fout [description]
     */
    void write_binary(FILE *fout) const
    {
      csr_t::write_binary(fout);
      /** write out the nrm and cmax arrays. */
      // TODO
    }

    /**
     * Read matrix from a binary file
     * @param fin [description]
     */
    void read_binary(FILE *fin)
    {
      csr_t::read_binary(fin);
      /** read in the nrm and cmax arrays */
      // TODO
    }

    ~l2_index_t()
    {
      if(nrm){
        free(nrm);
      }
      if(max){
        free(max);
      }
    }
  } l2_index_t;

  /**
   * Structure for a thread to do work in the L2Index.
   */
  typedef struct l2_workspace_t {
    BitVector* status = nullptr;
    idx_t* cands = nullptr;
    acm_t* accum = nullptr;
    Timer timers[NTIMERS];
    Counter counters[NCOUNTERS];

    l2_workspace_t(uint64_t n)
    {
      status = new BitVector(n); // true means it has been prunned
      cands = (idx_t *) malloc(sizeof(idx_t) * n);
      accum = (acm_t *) calloc(n, sizeof(acm_t));
      if(!status || !cands || !accum){
        throw runtime_error("Could not allocate workspace arrays.");
      }
    }

    ~l2_workspace_t()
    {
      delete status;
      free(cands);
      free(accum);
    }
  } l2_workspace_t;


  /**
   * L2 Index class for database of vectors to be searched.
   */
  template <typename T>
  struct L2Index : public Index<T> {
    Matrix<T> * data = nullptr;    // data in index
    Matrix<T> * norms = nullptr;   // suffix norms for data
    idx_t * cptr = nullptr;        // column start index for the un-indexed part of each row
    vector<l2_index_t *> indexes;  // set of indexes
    bool permute_columns = false;  // whether to permute the columns
    idx_t * perm = nullptr;        // column permutation
    val_t min_epsilon = 0.75;      // minimum epsilon value
    ptr_t ninnz = 1e+6;            // max number of non-zeros per index
    idx_t ninrows = 1e+5;          // max number of rows per index
    idx_t mncands = 0;             // max number of candidates, i.e., actual max number of rows per index, after training

    explicit L2Index(
      const idx_t ncols,
      const val_t min_epsilon=0.75
    )
     : Index<T>(ncols),
       min_epsilon(min_epsilon),
       norms(nullptr)
    {
      data = new Matrix<T>(0, ncols);
      this->is_trained = false;
    }

    explicit L2Index(
      const idx_t ncols,
      const idx_t n,
      const T* x,
      const val_t min_epsilon=0.75
    )
     : Index<T>(ncols),
       min_epsilon(min_epsilon),
       norms(nullptr)
    {
      data = new Matrix<T>(n, ncols, x, true); // copy data
      this->nsamples = n;
      this->is_trained = false;
    }

    explicit L2Index(
      const Matrix<T>& x,
      const val_t min_epsilon=0.75
    )
     : Index<T>(x.ncols),
       min_epsilon(min_epsilon),
       norms(nullptr)
    {
      data = new Matrix<T>(x);
      this->nsamples = x.nrows;
      this->is_trained = false;
    }

    ~L2Index()
    {
      if(data){
        delete data;
      }
      if(norms){
        delete norms;
      }
      if(cptr){
        free(cptr);
      }
      if(indexes.size()){
        for(auto idx: indexes){
          delete idx;
        }
      }
      if(perm){
        free(perm);
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
      idx_t mr = data->get_max_nrows();
      data->add(n, x);
      if(!this->ids){
        this->ids = (idx_t*) malloc(sizeof(idx_t) * data->get_max_nnz());
        if(!this->ids){
          throw runtime_error(fmt::format("Could not allocate {} elements for storing ids.", data->get_max_nnz()));
        }
      } else if(data->get_max_nrows() > nr){
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
    void train() override;

    /**
     * Re-train index for a different epsilon value.
     * @param min_epsilon Minimum epsilon that can be searched with after training.
     */
    void retrain(const val_t min_epsilon)
    {
        this->min_epsilon = min_epsilon;
        this->is_trained = false;
        this->train();
    }

    /**
     * Write the index to a file.
     * @param fpath Path where index should be written
     */
    void write(const string fpath) override;

    /**
     * Load index data from file.
     * @param fpath Path where index resides on disk.
     */
    void read(const string fpath) override;

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
        val_t* sims) override
    {
      throw runtime_error("L2Index does not implement knn_search. Use the L2KnnIndex instead.");
    }

    /** private methods -- should not be called directly */
    void _index();
    idx_t _generate_candidates(
      const idx_t qid,
      const Matrix<T>* qmat,
      const Matrix<T>* qnorms,
      const float epsilon,
      const l2_index_t* idx,
      l2_workspace_t * work
    );
    void _verify_candidates(
      const idx_t qid,
      const Matrix<T>* qmat,
      const Matrix<T>* qnorms,
      const float epsilon,
      const idx_t start,
      const idx_t ncands,
      l2_workspace_t * work,
      vector<vector<sim_t>>& results
    );

  }; // class L2Index

} // namespace fasta

#include <fasta/impl/L2Index.h>
