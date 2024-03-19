/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once

#include <fmt/format.h>

#include <fasta/types.h>

namespace fasta {

  using namespace std;

  /**
   * Base Matrix class for storing dense matrices and vectors.
   * Once `ncols` is set, it cannot be changed. The matrix can
   * grow by adding rows. Use `reserve` to allocate space for
   * a given number of rows before adding them. Otherwise, the
   * matrix will use a 2.0 growth factor, doubling its available
   * data space.
   */
  template <typename T>
  struct Matrix {

    /**
     * Generic constructor
     */
    Matrix(
      const idx_t nrows,
      const idx_t ncols
    );

    /**
     * Generic constructor with reserve space
     */
    Matrix(
      const idx_t nrows,
      const idx_t ncols,
      const idx_t max_nrows
    );

    /**
     * Constructor with data payload
     */
    Matrix(
      const idx_t nrows,
      const idx_t ncols,
      T* data,
      const bool local_data
    );

    /**
     * Constructor for matrix with set default value
     */
    Matrix(
      const idx_t nrows,
      const idx_t ncols,
      const T val
    );

    /**
     * Copy constructor
     */
    Matrix(const Matrix<T> &other);

    /**
     * Copy constructor from parts
     */
    Matrix(
      const idx_t nrows,
      const idx_t ncols,
      const T* data
    );

    /**
     * Destructror
     */
    virtual ~Matrix() {
      if(data && local_data){
        free(data);
      }
    }

    // Growth methods
    /**
     * Add rows to the matrix. If not enough space exists,
     * space will double to allow the add.
     * @param n Number of rows being added
     * @param x Vectors being added, size n * ncols
     */
    void add(const idx_t n, const T* x);

    /**
     * Add rows of `other` to this matrix. If not enough space exists,
     * space will grow by `other.nrows * ncols`
     * @param other [description]
     */
    void add(const Matrix<T> &other);

    /**
     * Reserve enough space in the matrix to store `n` rows. If space
     * already exists to store `n` rows, do nothing.
     * @param n Number of rows to potentially store in the matrix.
     */
    void reserve(const idx_t n);

    /**
     * Remove all extra space not currently used by rows.
     */
    void trim();

    // Assignment operators
    /**
     * Assignment from a different type matrix
     */
    template <typename T2>
    Matrix<T>& operator=(const Matrix<T2>& rhs);

    /**
     * Assignment from same type matrix
     */
    Matrix<T>& operator=(const Matrix<T>& rhs);

    // Setters
    void resize(const idx_t new_nrows, const idx_t new_ncols);
    void set(const T val=0.0);
    void set_random(const T min=0.0, const T max=1.0, const float fill=1.0, const uint64_t seed=0);

    // I/O methods
    static Matrix<T>* read_binary(FILE* fin);
    static Matrix<T>* read(const string fpath, const char format='d');
    static void write_binary(const Matrix<T>& mat, FILE* fout);
    static void write(const Matrix<T>& mat, const string fpath, const char format='d', const char delimiter=',');
    void write(const string fpath, const char format='d', const char delimiter=',') const;
    void print(ostream& stream=std::cout) const;
    string get_info() const;

    // Normalization and Scaling
    bool is_normalized(const char norm='2', const long double error=1e-7) const;
    Matrix<T> normalized(const char norm='2');
    Matrix<T>& normalize(const char norm='2');
    void compute_scaling_factors(const string mode, double** factors=nullptr);
    void scale(const idx_t n, T* x, const string mode, const double* factors) const;
    void unscale(const idx_t n, T* x, const string mode, const double* factors) const;
    void scale(const string mode, double** factors=nullptr);
    void scale(const string mode, const double* factors);

    // Permutations
    /**
     * Permute the columns of the matrix given the permutation array perm
     * @param perm Permutation array containing new column id for each column.
     */
    void permute_columns(const idx_t * perm);

    // Access operators
    /**
     * Subscript operators
     * @param  i Row ID
     * @param  j Column ID
     * @return   Value at MAT[i,j]
     * @thows out_of_range
     */
    T& operator() (idx_t i, idx_t j);
    T  operator() (idx_t i, idx_t j) const;

    // Mathematical operations
    Matrix<T> operator+(const Matrix<T>& rhs);
    Matrix<T>& operator+=(const Matrix<T>& rhs);
    Matrix<T> operator-(const Matrix<T>& rhs);
    Matrix<T>& operator-=(const Matrix<T>& rhs);
    Matrix<T> operator*(const Matrix<T>& rhs);
    Matrix<T>& operator*=(const Matrix<T>& rhs);
    Matrix<T>* transposed();
    Matrix<T>& transpose();

    // Matrix/scalar operations
    Matrix<T> operator+(const T val);
    Matrix<T>& operator+=(const T val);
    Matrix<T> operator-(const T val);
    Matrix<T>& operator-=(const T val);
    Matrix<T> operator*(const T val);
    Matrix<T>& operator*=(const T val);
    Matrix<T> operator/(const T val);
    Matrix<T>& operator/=(const T val);

    // Equality operators
    #ifndef FASTA_MATRIX_EPSILON
    #define FASTA_MATRIX_EPSILON 1e-5
    #endif
    bool operator==(const Matrix<T>& rhs) const;

    /**
     * Returns the number of rows in the matrix.
     */
    inline idx_t get_nrows() const
    {
      return nrows;
    }

    /**
     * Returns the maximum number of rows in the matrix.
     */
    inline idx_t get_max_nrows() const
    {
      return max_nrows;
    }

    /**
     * Returns the number of columns in the matrix.
     */
    inline idx_t get_ncols() const
    {
      return ncols;
    }

    /**
     * Returns the number of non-zeros in the matrix.
     */
    inline ptr_t get_nnz() const
    {
      return (ptr_t) nrows * ncols;
    }

    /**
     * Returns the number of non-zeros in the matrix.
     */
    inline ptr_t get_max_nnz() const
    {
      return (ptr_t) max_nrows * ncols;
    }

    /**
     * Returns a pointer to our data array. Potentially unsafe.
     */
    inline T * get_data()
    {
      return data;
    }

    /**
     * Returns a const pointer to our data array. Potentially unsafe.
     */
    inline const T * get_data() const
    {
      return data;
    }

    idx_t nrows;
    idx_t max_nrows;
    idx_t ncols;
    bool local_data;
    T* data;

    static uint8_t get_type_id();
    static string get_type();
    static T* read_binary_data(FILE * fin, uint8_t type_id, ptr_t nnz);

  };

} // namespace fasta

/** Implementation */
#include <fasta/impl/Matrix.h>
