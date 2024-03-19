/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <omp.h>

#include <fasta/Index.h>

namespace fasta {

  using namespace std;

  /**
   * Add a set of vectors to the index. IDs are assumed 0 -> n-1.
   *
   * @param x      sample vecors
   */
  template <typename T>
  void Index<T>::add(const Matrix<T>& x)
  {
    assert(x.ncols == ncols);
    add(x.nrows, x.data);
  }

  /**
   * Find `k` nearest neighbors within the index for `n` vectors of dimension `ncols`.
   *
   * @param n           number of vectors to search
   * @param x           input vectors to search, size n * d
   * @param k           number of nearest neighbors to retrieve for each sample
   * @param nbrs        nearest neighbor table, of size n * k
   */
  template <typename T>
  void Index<T>::knn_search(
      const idx_t n,
      const T* x,
      const idx_t k,
      Matrix<sim_t>& nbrs)
  {
    assert(n == nbrs.nrows && k == nbrs.ncols);
    idx_t* nids = (idx_t*) malloc(sizeof(idx_t) * n * k);
    val_t* sims = (val_t*) malloc(sizeof(val_t) * n * k);
    if(!nids || !sims){
      throw runtime_error(fmt::format("Could not allocate {} similarity elements.", (ptr_t) n * k));
    }
    knn_search(n, x, k, nids, sims);
    // transfer results to Matrix
    #pragma omp parallel for schedule(static, 20)
    for(idx_t i=0; i < n; ++i){
      for(idx_t j=0; j < k; ++j){
        nbrs(i,j) = (sim_t){nids[i*k+j], sims[i*k+j]};
      }
    }
    free(nids);
    free(sims);
  }

  /**
   * Replace row indexes with IDs in the result if `ids` are present in the index.
   * @param result Result structure
   */
  void set_result_ids(const idx_t* ids, csr_t& result)
  {
    if(!ids){
      return;
    }
    #pragma omp parallel for schedule(static, 20)
    for(ptr_t j=0; j < result.ptr[result.nrows]; ++j){
      result.ind[j] = ids[result.ind[j]];
    }
  }

  /**
   * Replace row indexes with IDs in the result if `ids` are present in the index.
   * @param n        number of vectors to search
   * @param k        number of nearest neighbors to retrieve for each sample
   * @param nbrs     nearest neighbor ids table, of size n * k
   */
  void set_result_ids(const idx_t* ids, const idx_t n, const idx_t k, idx_t* nbrs)
  {
    if(!ids){
      return;
    }
    ptr_t nnz = (ptr_t) n * k;
    #pragma omp parallel for schedule(static, 20)
    for(ptr_t j=0; j < nnz; ++j){
      nbrs[j] = ids[nbrs[j]];
    }
  }

}  // namespace fasta
