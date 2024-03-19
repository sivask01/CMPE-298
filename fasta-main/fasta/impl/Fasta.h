/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <omp.h>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <fasta/Fasta.h>
#include <fasta/GemmIndex.h>
#include <fasta/L2Index.h>

namespace fasta {

  using namespace std;


  /**
   * @brief Creates a new searcher
   * @param config Config parameters
   * @param dmat Database matrix to search against (samples are rows)
   */
  template <typename T>
  Fasta<T>::Fasta(
    Config& config,
    Matrix<T>& dmat
  ) :
    config(config),
    data(dmat)
  {
    // set number of threads
    spdlog::info("Setting omp_num_treads to {}.", config.nthreads);
    omp_set_num_threads(config.nthreads);
    create_index();
  }


  /**
   * @brief Create (partial) index needed for further search by the chosen search mode.
   * If index does not exist, method is called automatically when calling search.
   */
  template <typename T>
  void Fasta<T>::create_index()
  {
    spdlog::info("Indexing...");
    switch(config.method){
      case FASTA:
        if(config.mode == SIM){
          throw runtime_error("Not yet implemented.");
        } else if(config.mode == KNN){
          throw runtime_error("Not yet implemented.");
        } else {
          throw invalid_argument("Mode can either be knn or sim.");
        }
        break;
      case GEMM:
        index = new GemmIndex<T>(data);
        break;
      case BSI:
        throw runtime_error("Not yet implemented.");
        break;
      case L2:
        if(config.mode == SIM){
          index = new L2Index<T>(data, config.epsilon);
          auto idx = static_cast<L2Index<T>*>(index);
          idx->permute_columns = config.permute_columns;
          idx->ninnz = config.ninnz;
          idx->ninrows = config.ninrows;
          // TODO: pass other config options to index
        } else if(config.mode == KNN){
          throw runtime_error("Not yet implemented.");
        } else {
          throw invalid_argument("Mode can either be knn or sim.");
        }
        break;
      default:
        throw invalid_argument("Invalid method.");
        break;
    }
    // process the index data
    index->train();
  }

  /**
   * @brief Find dmat nearest neighbors for rows in qmat that have cosine similarity >= eps.
   * @param qmat Query matrix containing objects as rows.
   * @param result Vector of vectors of sim structures holding search results.
   */
  template <typename T>
  void Fasta<T>::search(
    Matrix<T>& qmat,
    csr_t& result
  ){
    spdlog::info("Searching...");
    index->search(qmat.get_nrows(), qmat.get_data(), config.epsilon, result);
  }

  /**
   * @brief Find the nnbrs dmat nearest neighbors for rows in qmat using cosine similarity.
   * @param qmat Query matrix containing objects as rows.
   * @param nnbrs The number of nearest neighbors to find.
   * @param result Matrix of size (qmat.nrows x nnbrs) storing results as sim structures.
   */
  template <typename T>
  void Fasta<T>::knn_search(
    Matrix<T>& qmat,
    const idx_t nnbrs,
    Matrix<sim_t>& result
  ){
    spdlog::info("Searching...");
    index->knn_search(qmat.get_nrows(), qmat.get_data(), nnbrs, result);
  }

  /**
   * @brief Find the nnbrs dmat nearest neighbors for rows in qmat using cosine similarity.
   * @param qmat Query matrix containing objects as rows.
   * @param nnbrs The number of nearest neighbors to find.
   * @param nbrs Array of size (qmat.nrows x nnbrs) storing neighbor IDs [optional, may be null]
   * @param sims Array of size (qmat.nrows x nnbrs) storing neighbor similarities [optional, may be null]
   */
  template <typename T>
  void Fasta<T>::knn_search(
    Matrix<T>& qmat,
    const idx_t nnbrs,
    idx_t * nbrs,
    val_t * sims
  ){
    spdlog::info("Searching...");
    index->knn_search(qmat.get_nrows(), qmat.get_data(), nnbrs, nbrs, sims);
  }

}  // namespace fasta
