/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once

#include <fasta/Matrix.h>
#include <fasta/GemmIndex.h>

namespace fasta {

  /**
   * @brief FASTA (Fast Angle-Based Search Target Aquisition) is an exact
   * Nearest Neighbor Search algorithm for dense data using
   * cosine similarity as the proximity measure of choice.
   */
  template <typename T>
  class Fasta {

  public:

    /**
     * @brief Creates a new searcher
     * @param config Config parameters
     * @param dmat Database matrix (samples are rows)
     */
    Fasta(
      Config& config,
      Matrix<T>& dmat
    );


    ~Fasta()
    {
      if(index){
        delete index;
      }
    }

    /**
     * @brief Find dmat nearest neighbors for rows in qmat that have cosine similarity >= eps.
     * @param qmat Query matrix containing objects as rows.
     * @param result Vector of vectors of sim structures holding search results.
     */
    void search(
      Matrix<T>& qmat,
      csr_t& result
    );

    /**
     * @brief Find the nnbrs dmat nearest neighbors for rows in qmat using cosine similarity.
     * @param qmat Query matrix containing objects as rows.
     * @param nnbrs The number of nearest neighbors to find.
     * @param result Matrix of size (qmat.nrows x nnbrs) storing results as sim structures.
     */
    void knn_search(
      Matrix<T>& qmat,
      const idx_t nnbrs,
      Matrix<sim_t>& result
    );

    /**
     * @brief Find the nnbrs dmat nearest neighbors for rows in qmat using cosine similarity.
     * @param qmat Query matrix containing objects as rows.
     * @param nnbrs The number of nearest neighbors to find.
     * @param nbrs Array of size (qmat.nrows x nnbrs) storing neighbor IDs [optional, may be null]
     * @param sims Array of size (qmat.nrows x nnbrs) storing neighbor similarities [optional, may be null]
     */
    void knn_search(
      Matrix<T>& qmat,
      const idx_t nnbrs,
      idx_t * nbrs,
      val_t * sims
    );

    inline Index<T> * get_index()
    {
      return index;
    }

  protected:

    /**
     * @brief Create (partial) index needed for further search by the chosen search mode.
     * If index does not exist, method is called automatically when calling search.
     */
    void create_index();

  private:
    Config& config;
    Matrix<T>& data;
    Index<T> * index;
  };

}  // namespace fasta

/** Implemenation */
#include <fasta/impl/Fasta.h>
