/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once

#include <assert.h>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

namespace fasta {

  /**
   * Default types
   */

  #ifdef IDX_TYPE_64
  using idx_t = std::uint64_t; ///< all indices are this type
  #else
  using idx_t = std::uint32_t; ///< all indices are this type
  #endif

  #ifdef VAL_TYPE_DOUBLE
  using val_t = double; ///< all values are this type
  #else
  using val_t = float; ///< all values are this type
  #endif

  using ptr_t = std::uintptr_t; ///< all pointers are this type

  #ifdef ACCUM_TYPE_FLOAT
  using acm_t = float;  ///< all similarity accumulations are this type
  #else
  using acm_t = double;  ///< all similarity accumulations are this type
  #endif

  /** temporary defines -- TODO: remove for production */
  #define EXTRA_TIMES
  #define EXTRA_COUNTERS

  /**
   * Structure to store a similarity for a given object.
   */
  typedef struct sim_t {
    idx_t cid;
    acm_t score;
  } sim_t;

  std::ostream& operator<< (std::ostream& out, const sim_t& s)
  {
      out << s.cid << " " << s.score;
      return out;
  }

  /**
   * CSR structure to store search results
   */
  typedef struct csr_t {
    idx_t nrows; // number of rows
    idx_t ncols; // number of rows
    idx_t * ind; // column ids
    val_t * val; // similarities
    ptr_t * ptr; // pointers (start of row in ind/sim)

    csr_t()
    {
      nrows = ncols = 0;
      ind = nullptr;
      val = nullptr;
      ptr = nullptr;
    }

    /**
     * Reserve space for more rows or non-zeros. Structure may only grow, not shrink.
     * @param nrows Number of rows
     * @param nnz   Number of non-zeros
     */
    void reserve(const idx_t nrows, const ptr_t nnz)
    {
      ptr_t cnnz = ptr ? ptr[this->nrows] : 0;
      if(nrows > this->nrows){
        if(ptr){
          ptr = (ptr_t*) realloc(ptr, sizeof(ptr_t) * (nrows+1));
        } else {
          ptr = (ptr_t*) malloc(sizeof(ptr_t) * (nrows+1));
          ptr[0] = 0;
        }
        if(!ptr){
          throw std::runtime_error("Could not allocate ptr array.");
        }
        this->nrows = nrows;
      }
      if(nnz > cnnz){
        if(ind){
          ind = (idx_t*) realloc(ind, sizeof(idx_t) * nnz);
        } else {
          ind = (idx_t*) malloc(sizeof(idx_t) * nnz);
        }
        if(!ind){
          throw std::runtime_error("Could not allocate ind array.");
        }
        if(val){
          val = (val_t*) realloc(val, sizeof(val_t) * nnz);
        } else {
          val = (val_t*) malloc(sizeof(val_t) * nnz);
        }
        if(!val){
          throw std::runtime_error("Could not allocate val array.");
        }
      }
    }

    /**
     * Reserve space for more rows or non-zeros. Structure may only grow, not shrink.
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param nnz   Number of non-zeros
     */
    void reserve(const idx_t nrows, const idx_t ncols, const ptr_t nnz)
    {
      this->ncols = ncols;
      reserve(nrows, nnz);
    }

    /**
     * Reduce the allocated space and size of matrix to `nrows` rows.
     * @param nrows Number of rows.
     */
    void trim(const idx_t nrows){
      if(nrows > this->nrows){
        throw std::runtime_error("Cannot trim to more rows than exist in the matrix. Try reserve.");
      }
      ptr_t nnz = ptr[nrows];
      ptr = (ptr_t*) realloc(ptr, sizeof(ptr_t) * (nrows+1));
      if(!ptr){
        throw std::runtime_error("Could not trim ptr array.");
      }
      ind = (idx_t*) realloc(ind, sizeof(idx_t) * nnz);
      if(!ind){
        throw std::runtime_error("Could not trim ind array.");
      }
      val = (val_t*) realloc(val, sizeof(val_t) * nnz);
      if(!val){
        throw std::runtime_error("Could not trim val array.");
      }
    }

    virtual ~csr_t()
    {
      if(ind){
        free(ind);
      }
      if(val){
        free(val);
      }
      if(ptr){
        free(ptr);
      }
    }

    /**
     * Write matrix to text file
     * @param output_fpath File to write to
     */
    void write(const std::string output_fpath, const bool header=false)
    {
      std::fstream resfile;
      resfile.open(output_fpath, std::ios::out);
      if(!resfile){
        throw std::runtime_error("Could not open output file for writing.");
      }
      if(header){
        resfile << nrows << " " << ncols << " " << ptr[nrows] << std::endl;
      }
      for(idx_t i=0; i < nrows; ++i){
        for(ptr_t j=ptr[i]; j < ptr[i+1]; ++j){
          resfile << ind[j] << " " << val[j];
          if(j+1 < ptr[i+1]){
            resfile << " ";
          }
        }
        resfile << std::endl;
      }
      resfile.close();
    }

    /**
     * Write matrix to a binary file
     * @param fout [description]
     */
    void write_binary(FILE *fout) const
    {
      /* write type ids for idx_t and val_t; ptr_t is always same type. */
      char c = std::is_same<idx_t, uint64_t>::value ? 'l' : 'i';
      if(fwrite(&c, sizeof(char), 1, fout) != 1){
        fclose(fout);
        throw std::runtime_error("Could not write idx_t type to file.");
      }
      c = std::is_same<val_t, double>::value ? 'd' : 'f';
      if(fwrite(&c, sizeof(char), 1, fout) != 1){
        fclose(fout);
        throw std::runtime_error("Could not write val_t type to file.");
      }
      /* write out nrows, ncols */
      if(fwrite(&nrows, sizeof(idx_t), 1, fout) != 1){
        fclose(fout);
        throw std::runtime_error("Could not write nrows to file.");
      }
      if(fwrite(&ncols, sizeof(idx_t), 1, fout) != 1){
        fclose(fout);
        throw std::runtime_error("Could not write ncols to file.");
      }
      /* write out the pointer array, followed by the indeces and values arrays */
      if(fwrite(ptr, sizeof(ptr_t), nrows+1, fout) != nrows+1){
        fclose(fout);
        throw std::runtime_error("Could not write ptr to file.");
      }
      if(nrows == 0){
        return;
      }
      if(fwrite(ind, sizeof(idx_t), ptr[nrows], fout) != ptr[nrows]){
        fclose(fout);
        throw std::runtime_error("Could not write ind to file.");
      }
      if(fwrite(val, sizeof(val_t), ptr[nrows], fout) != ptr[nrows]){
        fclose(fout);
        throw std::runtime_error("Could not write ind to file.");
      }
    }

    /**
     * Read matrix from a binary file
     * @param fin [description]
     */
    void read_binary(FILE *fin)
    {
      char c;
      /* read and check idx_t and val_t types */
      if(fread(&c, sizeof(char), 1, fin) != 1){
        fclose(fin);
        throw std::runtime_error("Could not read idx_t type from file.");
      }
      if( (c == 'l' && !std::is_same<idx_t, uint64_t>::value) ||
        (c == 'i' && !std::is_same<idx_t, uint32_t>::value) ||
        (c != 'l' && c!= 'i'))
      {
        fclose(fin);
        throw std::runtime_error("Mismatched idx_t type in binary CSR file.");
      }
      if(fread(&c, sizeof(char), 1, fin) != 1){
        fclose(fin);
        throw std::runtime_error("Could not read val_t type from file.");
      }
      if( (c == 'd' && !std::is_same<val_t, double>::value) ||
        (c == 'f' && !std::is_same<val_t, float>::value) ||
        (c != 'f' && c!= 'd'))
      {
        fclose(fin);
        throw std::runtime_error("Mismatched val_t type in binary CSR file.");
      }
      /* read nrows and ncols */
      idx_t n, m;
      if(fread(&n, sizeof(idx_t), 1, fin) != 1){
        fclose(fin);
        throw std::runtime_error("Could not read nrows from file.");
      }
      if(fread(&m, sizeof(idx_t), 1, fin) != 1){
        fclose(fin);
        throw std::runtime_error("Could not read nrows from file.");
      }
      /* allocate and read pointers array */
      if(ptr){
        ptr = (ptr_t*) realloc(ptr, sizeof(ptr_t) * (n+1));
      } else {
        ptr = (ptr_t*) malloc(sizeof(ptr_t) * (n+1));
      }
      if(fread(ptr, sizeof(ptr_t), n+1, fin) != n+1){
        fclose(fin);
        throw std::runtime_error("Could not read ptr array from file.");
      }
      /* allocate and read indeces and values arrays */
      if(ind){
        ind = (idx_t*) realloc(ind, sizeof(idx_t) * ptr[n]);
      } else {
        ind = (idx_t*) malloc(sizeof(idx_t) * ptr[n]);
      }
      if(fread(ind, sizeof(idx_t), ptr[n], fin) != ptr[n]){
        fclose(fin);
        throw std::runtime_error("Could not read ind array from file.");
      }
      if(val){
        val = (val_t*) realloc(val, sizeof(val_t) * ptr[n]);
      } else {
        val = (val_t*) malloc(sizeof(val_t) * ptr[n]);
      }
      if(fread(val, sizeof(val_t), ptr[n], fin) != ptr[n]){
        fclose(fin);
        throw std::runtime_error("Could not read val array from file.");
      }
      nrows = n;
      ncols = m;
    }
  } csr_t;

  /** Search mode */
  enum Modes {
      SIM, // Similarity Search
      KNN, // K-Nearest Neighbor Search
      NMODES
  };
  const std::string ModeNames[] = {"sim", "knn"};

  /** Methods */
  enum Methods {
      FASTA,
      GEMM,
      L2,
      BSI,
      NMETHODS
  };
  const std::string MethodNames[] = {"fasta", "gemm", "l2", "bsi"};

  /**
   * Configuration options for the search
   */
  typedef struct Config {
    /** General settings */
    uint32_t nthreads = 4;          // number of threads used to index and search
    uint8_t log_level = 0;          // no logging
    bool normalize = false;         // whether to normalize the rows of the matrix with L2 norm
    bool permute_columns = false;   // whether to permute columns (for methods that have that feature)
    double epsilon = 0.9;           // minimum similarity
    idx_t max_k = 100;              // number of neighbors in knn methods
    Modes mode = SIM;               // search mode
    Methods method = FASTA;         // search method

    uint64_t seed = 191876547;
    // which counters to enable

    /** Meta-parameters for search modes */
    // gemm
    // gemm_k
    // l2
    ptr_t ninnz = 1e+6;             // number of non-zeros for each inverted index
    idx_t ninrows = 1e+5;           // number of rows for each inverted index
    // l2_k
    // bsi
    // bsi_k
    // fasta
    // fasta_k
  } Config;

  /**
   * General exception for Fasta.
   */
  class FastaError : public std::logic_error {
    public:
    FastaError(const char* msg) : logic_error(msg) {}
  };




}
