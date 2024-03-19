/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <omp.h>

#include <fasta/L2Index.h>
#include <fasta/utils/sort.h>
#include <fasta/utils/KVVectorView.h>

namespace fasta {

  using namespace std;


  /**
   * Compute the median of the values in the `varr` array
   * @param  n    Number of elements in the array
   * @param  varr The array to compute the median of
   * @return      Median of the array values.
   */
  template <typename T>
  T compute_median(
    const size_t n,
    T* varr
  ){
    // edge cases
    if(n == 0){
      return 0;
    } else if(n == 1){
      return varr[0];
    } else if(n == 2){
      return (varr[0] + varr[1])/2.0 ;
    }

    // move the lowest values to the first half of the vector
    size_t mid = n/2 + 1;
    selecti(n, varr, mid);

    T val, max1;

    // single middle point if size is odd
    if(n % 2 == 1){
      max1 = varr[0];
      for(std::size_t i=1; i < mid; ++i){
        if(varr[i] > max1){
          max1 = varr[i];
        }
      }
      return max1;
    }

    // 2 highest values: max1 is the low max, and max2 is the high max
    T max2;
    if(varr[0] < varr[1]){
      max1 = varr[0];
      max2 = varr[1];
    } else {
      max1 = varr[1];
      max2 = varr[0];
    }
    for(std::size_t i=1; i < mid; ++i){
      if(varr[i] > max2){
        max1 = max2;
        max2 = varr[i];
      } else if(varr[i] > max1){
        max1 = varr[i];
      }
    }
    return (max1 + max2) / 2.0;
  }
  /**
   * Extract absolute median values for each column.
   * @param  data Matrix to analyze
   * @return      Array of absolute median values for each column.
   */
  template <typename T>
  val_t * get_column_absolute_medians(const Matrix<T> * data)
  {
    val_t * medians = (val_t*) malloc(sizeof(val_t) * data->ncols);
    if(!medians){
      throw runtime_error(fmt::format("Could not allocate {} elements for storing medians.", data->ncols));
    }

    #pragma omp parallel
    {
      val_t * col = (val_t*) malloc(sizeof(val_t) * data->nrows);
      if(!col){
        throw runtime_error(fmt::format("Could not allocate {} elements for storing a column of data.", data->nrows));
      }
      #pragma omp for schedule(static, 1)
      for(idx_t j=0; j < data->ncols; ++j){
        for(idx_t i=0; i < data->nrows; ++i){
          col[i] = (*data)(i,j) >= 0 ? (*data)(i,j) : -(*data)(i,j);
        }
        medians[j] = compute_median<T>(data->nrows, col);
      }
      free(col);
    }

    return medians;
  }


  /**
   * Create partial indexes
   */
  template <typename T>
  void L2Index<T>::_index()
  {
    idx_t start = 0;
    while(start < data->nrows){
      // figure out nrows and nnz for current index
      idx_t nr = 0;  // number of rows that will be indexed
      idx_t mnc = 0; // maximum number of prefix columns participating in index
      ptr_t nz = 0;  // index number of non-zeros
      for(idx_t i=start; i < data->nrows && nr < ninrows; ++i, ++nr){
        if(nz + cptr[i] > ninnz){
          break;
        }
        nz += cptr[i];
        if(cptr[i] > mnc){
          mnc = cptr[i];
        }
      }
      // allocate index and transfer non-zeros
      auto idx = new l2_index_t();
      idx->start = start;
      idx->reserve(mnc, data->nrows, nz);
      // first build pointer array
      memset((void*) idx->ptr, 0, sizeof(ptr_t) * (mnc+1));
      for(idx_t i=start; i < start+nr; ++i){
        for(idx_t j=0; j < cptr[i]; ++j){
          idx->ptr[j+1]++;
        }
      }
      for(idx_t j=0; j < mnc; ++j){
        idx->ptr[j+1] += idx->ptr[j]; // prefix sum
      }
      // now transfer data
      ptr_t p;
      for(idx_t i=start; i < start+nr; ++i){
        for(idx_t j=0; j < cptr[i]; ++j){
          p = idx->ptr[j];
          idx->ind[p] = i - start;     // row id in data, localized to this index
          idx->val[p] = (*data)(i,j);  // value
          idx->nrm[p] = (*norms)(i,j); // prefix norm
          idx->ptr[j]++;
        }
      }
      // finally, fix pointer array
      for(idx_t j=mnc; j > 0; j--){
        idx->ptr[j] = idx->ptr[j-1]; // shift right by 1
      }
      idx->ptr[0] = 0;
      // add index to indexes lists and move on
      idx->end = start + nr;
      indexes.push_back(idx);
      start += nr;
    }
    for(auto idx : indexes){
      if(idx->end - idx->start > mncands){
        mncands = idx->end - idx->start;
        this->counters[INDEX_SIZE] += idx->ptr[idx->nrows];
      }
    }
  }


  /**
   * Perform training on vectors that have already been added.
   */
  template <typename T>
  void L2Index<T>::train()
  {
    if(this->is_trained){
      return;
    }

    // normalize if necessary
    if(!data->is_normalized()){
      data->normalize();
    }

    // figure out column order and permute
    if(permute_columns){
      auto medians = get_column_absolute_medians(data);
      auto corder_kvv = KVVectorView<idx_t, val_t>(data->ncols, medians);
      corder_kvv.sort_vkd(); // sort medians in non-increasing order
      perm = corder_kvv.keys();
      data->permute_columns(perm);
      free(medians);
    }

    // compute suffix norms for all vectors and track where the forward index starts
    norms = new Matrix<T>(data->nrows, data->ncols);
    cptr = (idx_t *) malloc(sizeof(idx_t) * data->nrows);
    #pragma omp parallel for schedule(static, 32)
    for(idx_t i=0; i < data->nrows; ++i){
      long double v = 1.0, n;
      cptr[i] = 0;
      for(idx_t j=0; j < data->ncols; ++j){
        v -= (*data)(i,j) * (*data)(i,j);
        (*norms)(i,j) = n = v > 0 ? sqrt(v) : 0.0;
        if(n < min_epsilon && cptr[i] == 0){
          cptr[i] = j + 1;
        }
      }
    }

    // index prefixes
    this->_index();
    this->is_trained = true;
  }

  /**
   * Write the index to a file.
   * @param fpath Path where index should be written
   */
  template <typename T>
  void L2Index<T>::write(const string fpath)
  {
    throw runtime_error("Not yet implemented.");
  }

  /**
   * Load index data from file.
   * @param fpath Path where index resides on disk.
   */
  template <typename T>
  void L2Index<T>::read(const string fpath)
  {
    throw runtime_error("Not yet implemented.");
  }

  /**
   * Generate candidates from given partial inverted index.
   *
   * @param  qid      Query ID
   * @param  qmat     Query matrix
   * @param  qnorms   Suffix norms of the query matrix
   * @param  epsilon  Search radius / minimum similarity
   * @param  idx      Index to search
   * @param  cands    Buffer to store identified candidates
   * @param  accum    Accumulator
   * @param  status   Prunning status of identified candidates
   * @return          Number of candidates found
   */
  template <typename T>
  idx_t L2Index<T>::_generate_candidates(
    const idx_t qid,
    const Matrix<T>* qmat,
    const Matrix<T>* qnorms,
    const float epsilon,
    const l2_index_t* idx,
    l2_workspace_t * work
  ){
    idx_t * cands = work->cands;
    acm_t * accum = work->accum;
    BitVector* status = work->status;

    idx_t nc = 0;

    // column 0 should establish all candidates
    acm_t qv, qn, v;
    T cv, cn;
    idx_t ci;
    qv = (*qmat)(qid, 0);
    qn = (*qnorms)(qid, 0);
    for(ptr_t j=0; j < idx->ptr[1]; ++j){
      ci = idx->ind[j];
      cv = idx->val[j];
      cn = idx->nrm[j];
      v = qv * cv;
      if(v + qn * cn >= epsilon){
        status->zero(ci);
        accum[ci] = v;
        cands[nc++] = ci;
      } else {
        status->set(ci);  // no chance, prune
      }
    }
    work->counters[NUM_CANDIDATES] += nc;

    // process the rest of the idx columns
    for(idx_t i=1; i < idx->nrows; ++i){
      qv = (*qmat)(qid, i);
      qn = (*qnorms)(qid, i);
      for(ptr_t j=idx->ptr[i]; j < idx->ptr[i+1]; ++j){
        ci = idx->ind[j];
        if(status->isset(ci)){
          continue;
        }
        cv = idx->val[j];
        cn = idx->nrm[j];
        accum[ci] += qv * cv;
        if(accum[ci] + qn * cn < epsilon){
          status->set(idx->ind[j]);  // no chance, prune
          accum[ci] = 0.0;  // reset accumulator
          #ifdef EXTRA_COUNTERS
          work->counters[PRUNE_CG_L2]++;
          #endif
        }
      }
    }

    // clean out pruned condidates
    idx_t ncands = 0;
    for(idx_t i=0; i < nc; ++i){
      if(!status->isset(cands[i])){
        cands[ncands++] = cands[i];
      }
    }
    return ncands;
  }

  /**
   * Verify candidates, add results to result vector, and clear accum for next search.
   *
   * @param  qid      Query ID
   * @param  qmat     Query matrix
   * @param  qnorms   Suffix norms of the query matrix
   * @param  epsilon  Search radius / minimum similarity
   * @param  start    Index start ID
   * @param  ncands   Number of identified candidates
   * @param  work     Thread workspace (contains data needed for this search)
   * @param  results  Vector to store results in
   */
  template <typename T>
  void L2Index<T>::_verify_candidates(
    const idx_t qid,
    const Matrix<T>* qmat,
    const Matrix<T>* qnorms,
    const float epsilon,
    const idx_t start,
    const idx_t ncands,
    l2_workspace_t * work,
    vector<vector<sim_t>>& results
  ){
    idx_t * cands = work->cands;
    acm_t * accum = work->accum;

    acm_t v;
    idx_t ci;
    for(idx_t i=0; i < ncands; ++i){
      ci = start + cands[i];
      v = accum[cands[i]];
      accum[cands[i]] = 0; // reset accumulator
      bool prunned = false;
      // process the rest of the dot-product using the forward index.
      for(idx_t j=cptr[ci]; j < this->ncols; ++j){
        v += (*data)(ci, j) * (*qmat)(qid, j); // accumulate
        // check CS bound
        if(v + (*norms)(ci, j) * (*qnorms)(qid, j) < epsilon){
          #ifdef EXTRA_COUNTERS
          work->counters[PRUNE_CV_L2]++;
          #endif
          prunned = true;
          break;
        }
      }
      if(!prunned){
        work->counters[DOT_PRODUCTS]++;
        if(v >= epsilon){
          // add result
          results[qid].push_back((sim_t){ci, v});
          #ifdef EXTRA_COUNTERS
          work->counters[RESULTS]++;
          #endif
        }
      }
    }
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
  template <typename T>
  void L2Index<T>::search(
      const idx_t n,
      const T* x,
      const float epsilon,
      csr_t& result
  ){
    if(!this->is_trained){
      train();
    } else if(epsilon < min_epsilon){
      throw runtime_error(fmt::format("Current index was trained using min_epsilon = {}. Must retrain to use lower values.",
        min_epsilon));
      // TODO: add a method to retrain using a lower (or higher) epsilon value.
    }

    // pre-process query vectors
    this->timers[QUERY_PREPROCESS].start();
    Matrix<T> qmat(n, this->ncols, x);
    // normalize
    if(!qmat.is_normalized()){
      qmat.normalize();
    }
    // permute columns
    if(perm){
      qmat.permute_columns(perm);
    }
    // compute suffix norms for the query matrix
    auto qnorms = new Matrix<T>(qmat.nrows, qmat.ncols);
    #pragma omp parallel for schedule(static, 32)
    for(idx_t i=0; i < qmat.nrows; ++i){
      long double v = 1.0;
      for(idx_t j=0; j < qmat.ncols; ++j){
        v -= qmat(i,j) * qmat(i,j);
        (*qnorms)(i,j) = v > 0 ? sqrt(v) : 0.0;
      }
    }
    this->timers[QUERY_PREPROCESS].stop();

    // allocate workspace for each thread
    uint32_t nthreads;
    l2_workspace_t** work = nullptr;
    #pragma omp parallel
    {
      uint32_t tid = omp_get_thread_num();
      #pragma omp single
      {
        nthreads = omp_get_num_threads();
        work = (l2_workspace_t**) malloc(sizeof(l2_workspace_t *) * nthreads);
        if(!work){
          throw runtime_error("Could not allocate thread workspace arrays.");
        }
      }
      work[tid] = new l2_workspace_t(mncands);
    }

    vector<vector<sim_t>> results_vec(n, vector<sim_t>());
    for(auto idx : indexes){
      #pragma omp parallel
      {
        idx_t tid = omp_get_thread_num();
        l2_workspace_t* w = work[tid];

        #pragma omp for schedule(static, 32)
        for(idx_t i=0; i < n; ++i){
          // generate candidates
          #ifdef EXTRA_TIMES
          w->timers[CANDIDATE_GENERATION].start();
          #endif
          idx_t ncands = this->_generate_candidates(i, &qmat, qnorms, epsilon, idx, w);
          #ifdef EXTRA_TIMES
          w->timers[CANDIDATE_GENERATION].stop();
          w->timers[CANDIDATE_VERIFICATION].start();
          #endif
          // verify candidates and store results
          this->_verify_candidates(i, &qmat, qnorms, epsilon, idx->start, ncands, w, results_vec);
          #ifdef EXTRA_TIMES
          w->timers[CANDIDATE_VERIFICATION].stop();
          #endif
        }
      }
    }

    // transfer results to results structure
    ptr_t nnz = 0;
    for(idx_t i=0; i < n; ++i){
      nnz += results_vec[i].size();
    }
    result.reserve(n, data->nrows, nnz);
    for(idx_t i=0; i < n; ++i){
      result.ptr[i+1] = result.ptr[i] + results_vec[i].size();
    }
    #pragma omp parallel for schedule(static, 20)
    for(idx_t i=0; i < n; ++i){
      ptr_t p = result.ptr[i];
      for(auto item : results_vec[i]){
        result.ind[p] = item.cid;
        result.val[p++] = item.score;
      }
    }
    set_result_ids(this->ids, result);

    // consolidate timers and counters and free workspace memory
    for(uint32_t tid=0; tid < nthreads; ++tid){
      this->timers[CANDIDATE_GENERATION] += work[tid]->timers[CANDIDATE_GENERATION];
      this->timers[CANDIDATE_VERIFICATION] += work[tid]->timers[CANDIDATE_VERIFICATION];
      this->counters[NUM_CANDIDATES] += work[tid]->counters[NUM_CANDIDATES];
      this->counters[PRUNE_CG_L2] += work[tid]->counters[PRUNE_CG_L2];
      this->counters[PRUNE_CV_L2] += work[tid]->counters[PRUNE_CV_L2];
      this->counters[DOT_PRODUCTS] += work[tid]->counters[DOT_PRODUCTS];
      this->counters[RESULTS] += work[tid]->counters[RESULTS];
      delete (work[tid]);
    }
    this->timers[CANDIDATE_GENERATION] /= nthreads;
    this->timers[CANDIDATE_VERIFICATION] /= nthreads;
    free(work);
    delete qnorms;
  }


}  // namespace fasta
