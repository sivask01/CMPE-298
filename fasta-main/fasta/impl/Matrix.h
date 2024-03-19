/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <omp.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <cmath>
#include <cassert>
#include <ctype.h>
#include <string_view>
#include <stdexcept>
#include <algorithm>
#include <omp.h>
#include <fmt/format.h>

#include <fasta/Matrix.h>

namespace fasta {

  /**
   * @brief Generic constructor
   * @param nrows Number of rows
   * @param ncols Number of columns
   */
  template <typename T>
  Matrix<T>::Matrix(
    const idx_t nrows,
    const idx_t ncols
  ) : nrows(nrows),
    max_nrows(nrows),
    ncols(ncols),
    local_data(true)
  {
    if(max_nrows < 1){
      max_nrows = 1; // make sure allocation occurs even for empty matrix.
    }
    ptr_t max_nnz = max_nrows * ncols;
    data = (T *) malloc(sizeof(T) * max_nnz);
    if(!data){
      throw runtime_error(fmt::format("Could not allocate {} elements.", max_nnz));
    }
    if(nrows > 0){
      memset((void*) data, 0, sizeof(T) * nrows * ncols);
    }
  }

  /**
   * @brief Generic constructor with reserve space
   * @param nrows Number of rows
   * @param ncols Number of columns
   * @param max_nrows Maximum number of rows
   */
  template <typename T>
  Matrix<T>::Matrix(
    const idx_t nrows,
    const idx_t ncols,
    const idx_t max_nrows
  ) : nrows(nrows),
    max_nrows(max_nrows),
    ncols(ncols),
    local_data(true)
  {
    assert(max_nrows >= nrows);
    ptr_t max_nnz = max_nrows * ncols;
    data = (T *) malloc(sizeof(T) * max_nnz);
    if(!data){
      throw runtime_error(fmt::format("Could not allocate {} elements.", max_nnz));
    }
    if(nrows > 0){
      memset((void*) data, 0, sizeof(T) * nrows * ncols);
    }
  }

  /**
   * @brief Constructor with data payload.
   * @param nrows Number of rows
   * @param ncols Number of columns
   * @param data Array containing the data for this matrix
   * @param local_data Whether data should be deleted by matrix destructor
   */
  template <typename T>
  Matrix<T>::Matrix(
    const idx_t nrows,
    const idx_t ncols,
    T* data,
    const bool local_data
  ) : nrows(nrows),
    max_nrows(nrows),
    ncols(ncols),
    data(data),
    local_data(local_data){}

  /**
   * @brief Constructor for matrix with set default value
   * @param nrows Number or rows
   * @param ncols Number of columns
   * @param val Value that the matrix should be initialized with
   */
  template <typename T>
  Matrix<T>::Matrix(
    const idx_t nrows,
    const idx_t ncols,
    const T val
  ) : nrows(nrows),
    max_nrows(nrows),
    ncols(ncols),
    local_data(true)
  {
    if(max_nrows < 1){
      max_nrows = 1;
    }
    ptr_t max_nnz = max_nrows * ncols;
    data = (T *) malloc(sizeof(T) * max_nnz);
    if(!data){
      throw runtime_error(fmt::format("Could not allocate {} elements.", max_nnz));
    }
    if(nrows > 0){
      fill_n(data, get_nnz(), val);
    }
  }

  /**
   * @brief Copy constructor
   * @param rhs The matrix that should be copied.
   */
  template <typename T>
  Matrix<T>::Matrix(const Matrix<T> &rhs)
  : nrows(rhs.nrows),
    max_nrows(rhs.max_nrows),
    ncols(rhs.ncols),
    local_data(true)
  {
    data = (T *) malloc(sizeof(T) * get_max_nnz());
    if(!data){
      throw runtime_error(fmt::format("Could not allocate {} elements.", get_max_nnz()));
    }
    if(nrows > 0){
      memcpy((void*) data, (void*) rhs.data, sizeof(T) * get_nnz());
    }
  }

  /**
   * @brief Copy constructor from parts
   * @param nrows Number of rows
   * @param ncols Number of columns
   * @param data Array containing the data for this matrix
   */
  template <typename T>
  Matrix<T>::Matrix(
    const idx_t nrows,
    const idx_t ncols,
    const T* data
  ) : nrows(nrows),
    max_nrows(nrows),
    ncols(ncols),
    local_data(true)
  {
    this->data = (T *) malloc(sizeof(T) * get_nnz());
    if(!this->data){
      throw runtime_error(fmt::format("Could not allocate {} elements.", get_nnz()));
    }
    if(nrows > 0){
      memcpy((void*) this->data, (void*) data, sizeof(T) * get_nnz());
    }
  }

  // Growth methods
  /**
   * Add rows to the matrix. If not enough space exists,
   * space will doubled to allow the add.
   * @param n Number of rows being added
   * @param x Vectors being added, size n * ncols
   */
  template <typename T>
  void Matrix<T>::add(const idx_t n, const T* x)
  {
    assert(n > 0);
    if(nrows + n > max_nrows){
      idx_t mnr = nrows + n > max_nrows * 2.0 ? nrows + n : max_nrows * 2.0;
      reserve(mnr);
    }
    memcpy((void *)(data + get_nnz()), (void *) x, sizeof(T) * n * ncols);
    nrows += n;
  }

  /**
   * Add rows of `other` to this matrix. If not enough space exists,
   * space will grow by `other.nrows` rows.
   * @param other [description]
   */
  template <typename T>
  void Matrix<T>::add(const Matrix<T> &other)
  {
    assert(other.nrows > 0);
    reserve(nrows + other.nrows);
    add(other.nrows, other.data);
  }

  /**
   * Reserve enough space in the matrix to store `n` rows. If space
   * already exists to store `n` rows, do nothing.
   * @param n Number of rows to potentially store in the matrix.
   */
  template <typename T>
  void Matrix<T>::reserve(const idx_t n)
  {
    if(n > max_nrows){
      max_nrows = n;
      data = (T *) realloc(data, sizeof(T) * max_nrows * ncols);
      if(!data){
        throw runtime_error(fmt::format("Could not re-allocate {} elements.", max_nrows * ncols));
      }
    }
  }

  /**
   * Remove all extra space not currently used by rows in the matrix.
   * This operation is the oposite of `reserve`.
   */
  template <typename T>
  void Matrix<T>::trim()
  {
    if(max_nrows > nrows){
      data = (T *) realloc(data, sizeof(T) * get_nnz());
      if(!data){
        throw runtime_error(fmt::format("Could not re-allocate {} elements.", nrows * ncols));
      }
      max_nrows = nrows;
    }
  }

  // Assignment methods
  /**
   * Assignment from a different type matrix
   */
  template <typename T>
  template <typename T2>
  Matrix<T>& Matrix<T>::operator=(
    const Matrix<T2>& rhs
  ){
    if(nrows != rhs.nrows || ncols != rhs.ncols){
      throw invalid_argument("The rhs matrix must have the same number of rows and columns for this operation.");
    }
    ptr_t nnz = get_nnz();
    for(ptr_t i=0; i < nnz; ++i){
      data[i] = (T) rhs.data[i];
    }
  }

  /**
   * Assignment from same type matrix
   */
  template <typename T>
  Matrix<T>& Matrix<T>::operator=(
      const Matrix<T>& rhs
  ){
    if(nrows != rhs.nrows || ncols != rhs.ncols){
      throw invalid_argument("The rhs matrix must have the same number of rows and columns for this operation.");
    }
    memcpy((void*) data, (void*) rhs.data, sizeof(T) * get_nnz());
  }


  // Setters
  /**
   * @brief Set matrix to have all values equal to val.
   * @param val Value that should be set in all matrix elements.
   */
  template <typename T>
  void Matrix<T>::set(const T val)
  {
    if(val == 0.0){
      memset((void*) data, 0, sizeof(T) * get_nnz());
      return;
    }
    fill_n(data, get_nnz(), val);
  }

  /**
   * @brief Initialize matrix with random values.
   * Values and their locations are allocated using uniform normal distributions.
   * @param min  Minimum value
   * @param max  Maximum value
   * @param fill Percent fill for the matrix.
   * @param seed Seed for the random generator; default is 0, i.e., random seed
   */
  template <typename T>
  void Matrix<T>::set_random(
    const T min,
    const T max,
    const float fill,
    const uint64_t seed
  ){
    ptr_t nnz = get_nnz();
    ptr_t nfill = (ptr_t) (nnz * fill);
    ptr_t nzeros = nnz - nfill;
    if(fill <= 0 || fill > 1.0){
      throw out_of_range("Fill must be in the range (0,1].");
    }
    if(fill < 0.5){
      #pragma omp parallel
      {
        /* distribution for values */
        random_device rd;
        default_random_engine gen(seed ? seed : rd());
        uniform_real_distribution<float> dis((float) min, (float) max);
        auto val = bind(dis, gen);
        /* distribution for indexes */
        idx_t nthreads = omp_get_num_threads();
        idx_t tid = omp_get_thread_num();
        ptr_t zpt = (ptr_t) ((double) nnz / nthreads);
        uniform_int_distribution<ptr_t> pdis(tid * zpt, (tid+1) * zpt);
        auto place = bind(pdis, gen);
        /* set everything to 0 */
        memset((void*) data, 0, sizeof(T) * get_nnz());
        /* fill some of the values randomly */
        #pragma omp for schedule(static)
        for(ptr_t i=0; i < nfill; ++i){
          data[(const ptr_t) place()] = (T) val();
        }
      }
      return;
    }
    #pragma omp parallel
    {
      random_device rd;
      default_random_engine gen(seed ? seed : rd());
      uniform_real_distribution<float> dis((float) min, (float) max);
      auto val = bind(dis, gen);
      #pragma omp for schedule(static, 1024)
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) val();
      }

      /* set some values back to 0 */
      idx_t nthreads = omp_get_num_threads();
      idx_t tid = omp_get_thread_num();
      ptr_t zpt = (ptr_t) ((double) nnz / nthreads);
      uniform_int_distribution<ptr_t> pdis(tid * zpt, (tid+1) * zpt);
      auto place = bind(pdis, gen);
      #pragma omp for schedule(static)
      for(ptr_t i=0; i < nzeros; ++i){
        data[(const ptr_t) place()] = 0;
      }
    }
  }

  // I/O methods
  /**
   * @brief Figure out a type ID for encoding type in binary file.
   */
  template <typename T>
  uint8_t Matrix<T>::get_type_id()
  {
    if(is_same<T, sim_t>::value){ return 0; }
    if(is_same<T, int8_t>::value){ return 1; }
    if(is_same<T, int16_t>::value){ return 2; }
    if(is_same<T, int32_t>::value){ return 3; }
    if(is_same<T, int64_t>::value){ return 4; }
    if(is_same<T, intptr_t>::value){ return 5; }
    if(is_same<T, uint8_t>::value){ return 6; }
    if(is_same<T, uint16_t>::value){ return 7; }
    if(is_same<T, uint32_t>::value){ return 8; }
    if(is_same<T, uint64_t>::value){ return 9; }
    if(is_same<T, uintptr_t>::value){ return 10; }
    if(is_same<T, float>::value){ return 11; }
    if(is_same<T, double>::value){ return 12; }
    if(is_same<T, long double>::value){ return 13; }
    throw invalid_argument("Unsupported type for binary reading or writing.");
  }
  /**
   * @brief Figure out a type ID for encoding type in binary file.
   */
  template <typename T>
  string Matrix<T>::get_type()
  {
    if(is_same<T, sim_t>::value){ return string("sim_t"); }
    if(is_same<T, int8_t>::value){ return string("int8_t"); }
    if(is_same<T, int16_t>::value){ return string("int16_t"); }
    if(is_same<T, int32_t>::value){ return string("int32_t"); }
    if(is_same<T, int64_t>::value){ return string("int64_t"); }
    if(is_same<T, intptr_t>::value){ return string("intptr_t"); }
    if(is_same<T, uint8_t>::value){ return string("uint8_t"); }
    if(is_same<T, uint16_t>::value){ return string("uint16_t"); }
    if(is_same<T, uint32_t>::value){ return string("uint32_t"); }
    if(is_same<T, uint64_t>::value){ return string("uint64_t"); }
    if(is_same<T, uintptr_t>::value){ return string("uintptr_t"); }
    if(is_same<T, float>::value){ return string("float"); }
    if(is_same<T, double>::value){ return string("double"); }
    if(is_same<T, long double>::value){ return string("long double"); }
    throw invalid_argument("Unsupported type for binary reading or writing.");
  }

  /**
   * @brief Read the data array from a bmat file.
   * @param fin Input file pointer
   * @param type_id The type of data stored in the file.
   * @param nnz Number of elements that should be read.
   * @return array of data elements read from the file.
   */
  template <typename T>
  T* Matrix<T>::read_binary_data(
    FILE * fin,
    const uint8_t type_id,
    const uint64_t nnz
  ){
    ptr_t ret;
    T * data = (T *) malloc(sizeof(T) * nnz);
    if(!data){
      fclose(fin);
      throw runtime_error(fmt::format("Could not allocate {} elements.", nnz));
    }
    if(Matrix<T>::get_type_id() == type_id){
      ret = fread(data, sizeof(T), nnz, fin);
      if(ret != nnz){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not read data from matrix input file. Expected {} elements, read only {}.",
          nnz, ret));
      }
      return data;
    }
    if(Matrix<T>::get_type_id() == 0){
      throw runtime_error("Trying to read a similarity matrix into a structure that is not of type Matrix<sim_t>.");
    }

    if(type_id == 1){
      int8_t * d = (int8_t *) malloc(sizeof(int8_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type int8_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 2){
      int16_t * d = (int16_t *) malloc(sizeof(int16_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type int16_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 3){
      int32_t * d = (int32_t *) malloc(sizeof(int32_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type int32_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 4){
      int64_t * d = (int64_t *) malloc(sizeof(int64_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type int64_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 5){
      intptr_t * d = (intptr_t *) malloc(sizeof(intptr_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type intptr_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 6){
      uint8_t * d = (uint8_t *) malloc(sizeof(uint8_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type uint8_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 7){
      uint16_t * d = (uint16_t *) malloc(sizeof(uint16_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type uint16_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 8){
      uint32_t * d = (uint32_t *) malloc(sizeof(uint32_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type uint32_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 9){
      uint64_t * d = (uint64_t *) malloc(sizeof(uint64_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type uint64_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 10){
      uintptr_t * d = (uintptr_t *) malloc(sizeof(uintptr_t) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type uintptr_t.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 11){
      float * d = (float *) malloc(sizeof(float) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type float.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 12){
      double * d = (double *) malloc(sizeof(double) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type double.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    if(type_id == 13){
      long double * d = (long double *) malloc(sizeof(long double) * nnz);
      if(!d){
        fclose(fin);
        free(data);
        throw runtime_error(fmt::format("Could not allocate {} temporary elements of type long double.", nnz));
      }
      for(ptr_t i=0; i < nnz; ++i){
        data[i] = (T) d[i];
      }
      free(d);
      return data;
    }
    throw invalid_argument("Unsupported type for binary reading or writing.");
  }


  /**
   * Check if a line in the file is a header
   * @param  line Line in ASCII file
   * @return      Whether this is a header line
   */
  bool is_header(const string &line)
  {
    for(auto &ch : line){
      if(std::isalpha(static_cast<unsigned char>(ch)) && ch != 'e' && ch != 'E'){
        return true;
      }
    }
    return false;
  }

  /**
   * Find the delimiter in the line string. Delimiter is assumed to be the first
   * character after the first number in the first data line.
   * @param line Line in ASCII file
   */
  const char find_delimiter(const string &line)
  {
    istringstream iss(line);
    long double d;
    iss >> d;
    return iss.peek();
  }

  /**
   * Extract numbers from a line of text with the given delimiter
   * @param line      Line in ASCII file
   * @param lvd       Vector to store numbers in
   * @param delimiter Delimiter character (e.g., comma, semicolon, or tab)
   */
  void get_numbers(const string &line, vector<long double> &lvd, const char delimiter)
  {
    istringstream iss(line);
    long double d;
    char c;
    lvd.clear();
    while(iss >> d){
      lvd.push_back(d);
      do{
        iss.ignore(); // expected delimiter right after number
        while(iss.peek() == ' '){ // allow spaces after commas
          iss.ignore();
        }
        if (iss.peek() == delimiter){ // empty cell, assume 0
          lvd.push_back(0);
        }
      } while(iss.peek() == delimiter);
    }
  }

  /**
   * Read a matrix stored in a binary matrix file.
   * Note that file pointer is not closed after read and should be closed
   * by calling function.
   * @param fin File pointer where matrix is to be read from.
   * @return Matrix that was read.
   */
  template <typename T>
  Matrix<T>* Matrix<T>::read_binary(FILE* fin)
  {
    uint64_t nrows;
    uint64_t ncols;
    ptr_t ret;
    uint8_t type_id;
    ret  = fread(&nrows, sizeof(uint64_t), 1, fin);
    if(ret != 1){
      fclose(fin);
      throw runtime_error("Could not read number of rows from file.");
    }
    ret  = fread(&ncols, sizeof(uint64_t), 1, fin);
    if(ret != 1){
      fclose(fin);
      throw runtime_error("Could not read number of columns from file.");
    }
    ret  = fread(&type_id, sizeof(uint8_t), 1, fin);
    if(ret != 1){
      fclose(fin);
      throw runtime_error("Could not read the data type id from file.");
    }
    uint64_t nnz = nrows * ncols;
    // read data based on the encoded type
    auto data = Matrix<T>::read_binary_data(fin, type_id, nnz);
    return new Matrix<T>(nrows, ncols, data, true);
  }


  /**
   * @brief Read a matrix from a file. File can be ASCII or binary.
   * ASCII file should be a standard CSV, TSV, or space delimited table, with
   * or without a header. See the \ref{header} parameter.
   * Binary file should be one written by this library.
   * File format can be deduced from the file extension, with expected
   * ASCII format extensions being .csv, .tsv, .txt, .mat, and binary format
   * extension being .bmat.
   *
   * @param fpath  File path from which to read
   * @param format File format:
   *               'd' = deduce from file extension
   *               't' = ascii text delimited by space, tab, or comma
   *               'b' = binary
   */
  template <typename T>
  Matrix<T>* Matrix<T>::read(
    const string fpath,
    const char format
  )
  {
    uint64_t nrows;
    uint64_t ncols;
    if(format == 'b' || format == 'B' || std::string_view(fpath).ends_with(".bmat")){
      // Binary matrix
      FILE * fin  = fopen (fpath.c_str(), "rb");
      if(!fin){
        throw invalid_argument(fmt::format("Could not open file {}.", fpath));
      }
      auto mat = Matrix<T>::read_binary(fin);
      fclose(fin);
      return mat;
    }
    // ASCII file
    if(is_same<T, sim_t>::value){
      throw invalid_argument("read not yet implemented for sim_t matrices.");
    }
    ifstream ifile (fpath, ifstream::in);
    string line;
    if(!ifile.is_open()){
      ifile.close();
      throw invalid_argument(fmt::format("Could not open file {}.", fpath));
    }
    // find nrows and ncols
    if(!getline(ifile, line)){
      ifile.close();
      throw invalid_argument(fmt::format("Could not read first line from {}.", fpath));
    }
    bool has_header;
    if(has_header = is_header(line)){ // allow and skip one header file and empty lines
      do {
        getline(ifile, line);
      } while(line.empty() && !ifile.eof());
      if(line.empty() || ifile.eof()){
        ifile.close();
        throw invalid_argument(fmt::format("Could not read first data line from {}.", fpath));
      }
    }
    nrows = 1;
    std::vector<long double> vld;
    auto delimiter = find_delimiter(line);
    get_numbers(line, vld, delimiter);
    ncols = vld.size();
    if (ncols == 0){
      ifile.close();
      throw runtime_error(fmt::format("Invalid input file {}. The file seems to have {} header but 0 data columns.",
        fpath, has_header ? "a" : "no"));
    }
    while(getline(ifile, line)){
      if(!line.empty()){ // deals with windows line endings in linux system
        nrows++;
      }
    }
    // allocate matrix
    auto mat = new Matrix<T>(nrows, ncols);
    // reset file to start
    ifile.clear();
    ifile.seekg(0, ifile.beg);
    if(has_header){
      getline(ifile, line);
    }
    uint64_t i=0;
    while(getline(ifile, line)){
      if(!line.empty()){
        get_numbers(line, vld, delimiter);
        if(vld.size() != ncols){
          delete mat;
          ifile.close();
          throw runtime_error(fmt::format("Invalid input file {}. Data row {} returned {} values. "
            "Expected {}. Delimiter is '{}'. File has {} header.",
            fpath, i+1, vld.size(), ncols, delimiter, has_header ? "a" : "no"));
        }
        for(uint64_t j=0; j < ncols; ++j){
          mat->operator()(i, j) = (T) vld[j];
        }
        assert(i++ < nrows);
      }
    }
    ifile.close();
    return mat;
  }


  /**
   * Write matrix to a binary file.
   * Note that file pointer is not closed after write and should be closed
   * by calling function.
   * @param mat  Matrix to write
   * @param fout File pointer that matrix should be written to.
   */
  template <typename T>
  void Matrix<T>::write_binary(const Matrix<T>& mat, FILE* fout)
  {
    uint8_t type_id = Matrix<T>::get_type_id();
    // Binary matrix file
    uint64_t wnrows = mat.nrows;
    uint64_t wncols = mat.ncols;
    /* write out nrows, ncols, and type id */
    if(fwrite(&wnrows, sizeof(uint64_t), 1, fout) != 1){
      fclose(fout);
      throw runtime_error("Could not write nrows to file.");
    }
    if(fwrite(&wncols, sizeof(uint64_t), 1, fout) != 1){
      fclose(fout);
      throw runtime_error("Could not write ncols to file.");
    }
    if(fwrite(&type_id, sizeof(uint8_t), 1, fout) != 1){
      fclose(fout);
      throw runtime_error("Could not write type id to file.");
    }
    if(fwrite(mat.data, sizeof(T), mat.get_nnz(), fout) != mat.get_nnz()){
      fclose(fout);
      throw runtime_error("Could not write data values to file.");
    }
  }

  /**
   * @brief Write a given matrix to a file.
   * @param mat Matrix to write out.
   * @param fpath  File path from which to read
   * @param format File format:
   *               'd' = deduce from file extension
   *               't' = ascii text delimited by space, tab, or comma
   *               'b' = binary
   * @param delimiter Delimiter for ASCII output
   */
  template <typename T>
  void Matrix<T>::write(
    const Matrix<T>& mat,
    const string fpath,
    const char format,
    const char delimiter
  ){
    uint8_t type_id = Matrix<T>::get_type_id();
    if(format == 'b' || format == 'B' || std::string_view(fpath).ends_with(".bmat")){
      FILE * fout = fopen (fpath.c_str(), "wb");
      if(!fout){
        throw invalid_argument(fmt::format("Could not open file {} for writing.", fpath));
      }
      Matrix<T>::write_binary(mat, fout);
      fclose(fout);
      return;
    }
    // ASCII file
    fstream ofile;
    ofile.open(fpath, ios::out);
    if(!ofile){
      throw invalid_argument(fmt::format("Could not open file {} for writing.", fpath));
    }
    for(idx_t i=0; i < mat.nrows; ++i){
      for(idx_t j=0; j < mat.ncols-1; ++j){
        ofile << mat(i,j) << delimiter;
      }
      ofile << mat(i, mat.ncols-1) << endl;
    }
    ofile.close();
  }

  /**
   * @brief Write a given matrix to a file.
   * @param fpath  File path from which to read
   * @param format File format:
   *               'd' = deduce from file extension
   *               't' = ascii text delimited by space, tab, or comma
   *               'b' = binary
   * @param delimiter Delimiter for ASCII output
   */
  template <typename T>
  void Matrix<T>::write(
    const string fpath,
    const char format,
    const char delimiter
  ) const
  {
    Matrix<T>::write(*this, fpath, format, delimiter);
  }

  /**
   * Print the matrix to a stream
   * @param stream What stream to print to
   */
  template <typename T>
  void Matrix<T>::print(ostream& stream) const {
    stream << "Matrix<" << this->get_info() << endl;
    for(idx_t i=0; i < nrows; ++i){
      for(idx_t j=0; j < ncols; ++j){
        stream << this->operator()(i,j) << " ";
      }
      stream << endl;
    }
  }


  /**
   * Human readable information about this matrix
   * @return string describing the matrix
   */
  template <typename T>
  string Matrix<T>::get_info() const
  {
    return "Matrix<" + this->get_type() + ">(" + to_string(nrows) + ", " + to_string(ncols) + ")";
  }

  // Normalization and Scaling
  /**
   * @brief Normalize matrix, i.e., scale rows by their norm.
   * @param norm Which norm, one of 1,2, or i for (infinity)
   * @return A copy of the matrix with rows normalized by the given norm.
   */
  template <typename T>
  Matrix<T> Matrix<T>::normalized(const char norm)
  {
    T* d = (T*) malloc(sizeof(T) * get_nnz());
    switch(norm){
    case '1':
      #pragma omp parallel for schedule(static, 1024)
      for(ptr_t i=0; i < nrows; ++i){
        long double nrm = 0;
        for(ptr_t j=0; j < ncols; ++j){
          nrm += data[i * ncols + j];
        }
        if(nrm != 0){ /* ignorning empty rows */
          nrm = 1.0/nrm;
          for(ptr_t j=0; j < ncols; ++j){
            d[i * ncols + j] = data[i * ncols + j] * nrm;
          }
        }
      }
      break;
    case '2':
      #pragma omp parallel for schedule(static, 1024)
      for(ptr_t i=0; i < nrows; ++i){
        long double nrm = 0;
        for(ptr_t j=0; j < ncols; ++j){
          nrm += data[i * ncols + j] * data[i * ncols + j];
        }
        if(nrm != 0){ /* ignorning empty rows */
          nrm = 1.0/std::sqrt(nrm);
          for(ptr_t j=0; j < ncols; ++j){
            d[i * ncols + j] = data[i * ncols + j] * nrm;
          }
        }
      }
      break;
    case 'i':
      #pragma omp parallel for schedule(static, 1024)
      for(ptr_t i=0; i < nrows; ++i){
        long double nrm = 0;
        for(ptr_t j=0; j < ncols; ++j){
          if(std::abs(data[i * ncols + j]) > nrm){
            nrm = std::abs(data[i * ncols + j]);
          }
        }
        if(nrm != 0){ /* ignorning empty rows */
          nrm = 1.0/nrm;
          for(ptr_t j=0; j < ncols; ++j){
            d[i * ncols + j] = data[i * ncols + j] * nrm;
          }
        }
      }
      break;
    default:
      throw invalid_argument("Only L1, L2, and L-infinity norms are supported. Try norm=1, 2, or i.");
      break;
    }
    return Matrix<T>(nrows, ncols, d);
  }


  /**
   * @brief Check if matrix is normalized with the given norm
   * Note: function does a spot check of 2 random rows in the matrix.
   * @param norm Which norm, one of 1,2, or i for (infinity)
   * @param error Allowed error in the norm float comparison, defaults to 1e-5.
   * @return Whether matrix is normalized.
   */
  template <typename T>
  bool Matrix<T>::is_normalized(
    const char norm,
    const long double error
  ) const
  {
    idx_t i1 = rand() % nrows;
    idx_t i2 = nrows - 1; // always check the last row
    long double nrm = 0;
    switch(norm){
    case '1':
      for(idx_t j=0; j < ncols; ++j){
        nrm += data[i1 * ncols + j];
      }
      if(std::abs(nrm - 1.0) > error){
        return false;
      }
      nrm = 0.0;
      for(idx_t j=0; j < ncols; ++j){
        nrm += data[i2 * ncols + j];
      }
      if(std::abs(nrm - 1.0) > error){
        return false;
      }
      break;
    case '2':
      for(idx_t j=0; j < ncols; ++j){
        nrm += data[i1 * ncols + j] * data[i1 * ncols + j];
      }
      if(std::abs(std::sqrt(nrm) - 1.0) > error){
        return false;
      }
      nrm = 0.0;
      for(idx_t j=0; j < ncols; ++j){
        nrm += data[i2 * ncols + j] * data[i2 * ncols + j];
      }
      if(std::abs(std::sqrt(nrm) - 1.0) > error){
        return false;
      }
      break;
    case 'i':
      for(idx_t j=0; j < ncols; ++j){
        if(std::abs(data[i1 * ncols + j]) > 1.0 + error){
          return false;
        }
      }
      for(idx_t j=0; j < ncols; ++j){
        if(std::abs(data[i2 * ncols + j]) > 1.0 + error){
          return false;
        }
      }
      break;
    default:
      throw invalid_argument("Only L1, L2, and L-infinity norms are supported. Try norm=1, 2, or i.");
      break;
    }
    return true;
  }

  /**
   * @brief Normalize matrix, i.e., scale rows by their norm.
   * @param norm Which norm, one of 1,2, or i for (infinity)
   * @return Normalizes this matrix and returns a reference to this.
   */
  template <typename T>
  Matrix<T>& Matrix<T>::normalize(const char norm)
  {
    switch(norm){
    case '1':
      #pragma omp parallel for schedule(static, 1024)
      for(idx_t i=0; i < nrows; ++i){
        long double nrm = 0.0;
        for(idx_t j=0; j < ncols; ++j){
          nrm += data[i * ncols + j];
        }
        if(nrm != 0.0){ /* ignorning empty rows */
          nrm = 1.0/nrm;
          for(idx_t j=0; j < ncols; ++j){
            data[i * ncols + j] = (T) (data[i * ncols + j] * nrm);
          }
        }
      }
      break;
    case '2':
      #pragma omp parallel for schedule(static, 1024)
      for(idx_t i=0; i < nrows; ++i){
        long double nrm = 0.0;
        for(idx_t j=0; j < ncols; ++j){
          nrm += data[i * ncols + j] * data[i * ncols + j];
        }
        if(nrm != 0.0){ /* ignorning empty rows */
          nrm = 1.0/std::sqrt(nrm);
          for(idx_t j=0; j < ncols; ++j){
            data[i * ncols + j] = (T) (data[i * ncols + j] * nrm);
          }
        }
      }
      break;
    case 'i':
      #pragma omp parallel for schedule(static, 1024)
      for(idx_t i=0; i < nrows; ++i){
        long double nrm = 0.0;
        for(idx_t j=0; j < ncols; ++j){
          if(std::abs(data[i * ncols + j]) > nrm){
            nrm = std::abs(data[i * ncols + j]);
          }
        }
        if(nrm != 0.0){ /* ignorning empty rows */
          nrm = 1.0/nrm;
          for(idx_t j=0; j < ncols; ++j){
            data[i * ncols + j] = (T) (data[i * ncols + j] * nrm);
          }
        }
      }
      break;
    default:
      throw invalid_argument("Only L1, L2, and L-infinity norms are supported. Try norm=1, 2, or i.");
      break;
    }
    return *this;
  }



  /**
   * Scale the columns of the matrix in `x` by the chosen `mode`. Mode can be one of:
   *   `stdev` or `zscore`: Standardize the columns, i.e., x' = (x-mean)/stdev (unbiased, divide by n-1)
   *   `stdevp`: Standardize the columns, i.e., x' = (x-mean)/stdev (biased, divide by n)
   *   `minmax` or 'min-max': x' = (x-min)/(max-min)
   *   `mean`: x' = (x-mean)/(max-min)
   *   `sqrt`: x' = sqrt(x)
   *   `log`: x' = 1 + log_2(x)
   *   `maxtf`: x' = .5 + .5 * x/max
   * `x` is assumed to have the same number of columns as this matrix.
   *
   * @param n       Number of vectors to scale.
   * @param x       Data array of size n * ncols.
   * @param mode    How to scale.
   * @param factors Optional pointer for return factors (they will be allocated within).
   */
  template <typename T>
  void Matrix<T>::compute_scaling_factors(const string mode, double** factors)
  {
    if(!factors){
      return;
    }
    uint64_t nfactors = 0;
    if(mode == "stdev" || mode == "stdevp" || mode == "zscore" ||
      mode == "minmax" || mode == "min-max"
    ){
      nfactors = 2;
    } else if(mode == "maxtf"){
      nfactors = 1;
    } else if(mode == "mean"){
      nfactors = 3;
    } else if(mode != "log" && mode != "sqrt") {
      throw invalid_argument("Scaling mode must be one of: stdev, stdevp, zscore, minmax, maxtf, mean, log, or sqrt.");
    }
    double * cscale = nullptr;
    if(nfactors){
      cscale = (double*) calloc(sizeof(double), nfactors * ncols);
      if(!cscale){
        throw runtime_error("Could not allocate array for storing scaling factors.");
      }
    }
    if(nrows == 0){
      if(factors){
        *factors = cscale;
      } else {
        free(cscale);
      }
      return;
    }

    // compute the scaling factors
    T v;
    double s, s2;
    if(mode == "stdev" || mode == "stdevp" || mode == "zscore"){
      // find the mean and stdev for each column
      for(idx_t i=0; i < nrows; ++i){
        #pragma omp parallel for if(ncols > 1e+5) private(v) schedule(static, 32)
        for(idx_t j=0; j < ncols; ++j){
          v = data[i * ncols + j];
          cscale[j] += v;  // sum
          cscale[ncols + j] += v * v; // sum of squares
        }
      }
      double b = mode == "stdevp" ? 0 : nrows == 1 ? 1 - 1e-10 : 1; // avoid division by 0
      #pragma omp parallel for if(ncols > 1e+5) private(s, s2) schedule(static, 32)
      for(idx_t j=0; j < ncols; ++j){
        s = cscale[j];
        s2 = cscale[ncols + j];
        cscale[ncols + j] = sqrt((s2 - s * s / (nrows - b)) / (nrows - b)); // stddev
        if(cscale[ncols + j] == 0){
          cscale[ncols + j] = 1e-10; // avoid division by 0
        }
        cscale[j] = s / nrows; // mean
      }
    } else if(mode == "minmax" || mode == "min-max") {
      // find min and max for each column
      memcpy((void*) cscale, (void*) data, sizeof(T) * ncols);  // first row
      memcpy((void*) (cscale + ncols), (void*) data, sizeof(T) * ncols);  // first row
      for(idx_t i=1; i < nrows; ++i){
        #pragma omp parallel for if(ncols > 1e+5) private(v) schedule(static, 32)
        for(idx_t j=0; j < ncols; ++j){
          v = data[i * ncols + j];
          if(v < cscale[j]){
            cscale[j] = v; // min
          }
          if(v > cscale[ncols + j]){
            cscale[ncols + j] = v; // max
          }
        }
      }
      #pragma omp parallel for if(ncols > 1e+5) schedule(static, 32)
      for(idx_t j=0; j < ncols; ++j){
        if(cscale[j] == cscale[ncols + j]){
          cscale[ncols + j] += 1e-10; // avoid division by 0
        }
      }
    } else if(mode == "mean"){
      // find min, max, and mean for each column
      memcpy((void*) cscale, (void*) data, sizeof(T) * ncols);  // first row
      memcpy((void*) (cscale + ncols), (void*) data, sizeof(T) * ncols);  // first row
      memcpy((void*) (cscale + 2 * ncols), (void*) data, sizeof(T) * ncols);  // first row
      for(idx_t i=1; i < nrows; ++i){
        #pragma omp parallel for if(ncols > 1e+5) private(v) schedule(static, 32)
        for(idx_t j=0; j < ncols; ++j){
          v = data[i * ncols + j];
          if(v < cscale[j]){
            cscale[j] = v; // min
          }
          if(v > cscale[ncols + j]){
            cscale[ncols + j] = v; // max
          }
          cscale[2 * ncols + j] += v;
        }
      }
      #pragma omp parallel for if(ncols > 1e+5) schedule(static, 32)
      for(idx_t j=0; j < ncols; ++j){
        if(cscale[j] == cscale[ncols + j]){
          cscale[ncols + j] += 1e-10; // avoid division by 0
        }
      }
    } else if(mode == "maxtf"){
      // find max for each column
      memcpy((void*) cscale, (void*) data, sizeof(T) * ncols);  // first row
      #pragma omp parallel for if(ncols > 1e+5) schedule(static, 32)
      for(idx_t j=0; j < ncols; ++j){
        if(cscale[j] < 0){
          cscale[j] *= -1;  // absolute value
        }
      }
      for(idx_t i=1; i < nrows; ++i){
        #pragma omp parallel for if(ncols > 1e+5) private(v) schedule(static, 32)
        for(idx_t j=0; j < ncols; ++j){
          v = data[i * ncols + j] > 0 ? data[i * ncols + j] : (data[i * ncols + j] * -1);
          if(v > cscale[j]){
            cscale[j] = v; // max absolute value
          }
        }
      }
      #pragma omp parallel for if(ncols > 1e+5) schedule(static, 32)
      for(idx_t j=0; j < ncols; ++j){
        if(cscale[j] == 0){
          cscale[j] = 1e-10; // avoid division by 0
        }
      }
    }
    *factors = cscale;
  }

  /**
   * Scale the columns of the matrix in `x` by the chosen `mode`. Mode can be one of:
   *   `stdev` or `zscore`: Standardize the columns, i.e., x' = (x-mean)/stdev (unbiased, divide by n-1)
   *   `stdevp`: Standardize the columns, i.e., x' = (x-mean)/stdev (biased, divide by n)
   *   `minmax` or 'min-max': x' = (x-min)/(max-min)
   *   `mean`: x' = (x-mean)/(max-min)
   *   `sqrt`: x' = sqrt(x)
   *   `log`: x' = 1 + log_2(x)
   *   `maxtf`: x' = .5 + .5 * x/max
   * `x` is assumed to have the same number of columns as this matrix.
   *
   * @param n       Number of vectors to scale.
   * @param x       Data array of size n * ncols.
   * @param mode    How to scale.
   * @param factors Factors to scale by (should be obtained via `compute_scaling_factors`).
   */
  template <typename T>
  void Matrix<T>::scale(const idx_t n, T* x, const string mode, const double* factors) const
  {
    if(!factors){
      return;
    }
    // apply scaling
    if(mode == "stdev" || mode == "stdevp" || mode == "zscore"){
      #pragma omp parallel for if(n > 1e+5) schedule(static)
      for(idx_t i=0; i < n; ++i){
        for(idx_t j=0; j < ncols; ++j){
          x[i * ncols + j] = (x[i * ncols + j] - factors[j]) / factors[ncols + j];
        }
      }
    } else if(mode == "minmax" || mode == "min-max") {
      #pragma omp parallel for if(n > 1e+5) schedule(static)
      for(idx_t i=0; i < n; ++i){
        for(idx_t j=0; j < ncols; ++j){
          x[i * ncols + j] = (x[i * ncols + j] - factors[j]) / (factors[ncols + j] - factors[j]);
        }
      }
    } else if(mode == "mean"){
      #pragma omp parallel for if(n > 1e+5) schedule(static)
      for(idx_t i=0; i < n; ++i){
        for(idx_t j=0; j < ncols; ++j){
          x[i * ncols + j] = (x[i * ncols + j] - factors[2 * ncols + j]) / (factors[ncols + j] - factors[j]);
        }
      }
    } else if(mode == "maxtf"){
      #pragma omp parallel for if(n > 1e+5) schedule(static)
      for(idx_t i=0; i < n; ++i){
        for(idx_t j=0; j < ncols; ++j){
          x[i * ncols + j] = 0.5 + 0.5 * x[i * ncols + j] / factors[j];
        }
      }
    } else if(mode == "log") {
      double logscale = 1.0 / log(2.0);
      #pragma omp parallel for if(n > 1e+5) schedule(static)
      for(idx_t i=0; i < n; ++i){
        for(idx_t j=0; j < ncols; ++j){
          // scale value but keep sign
          x[i * ncols + j] = 1 + (x[i * ncols + j] > 0 ? log(x[i * ncols + j]) : -log(-x[i * ncols + j])) * logscale;
        }
      }
    } else if(mode == "sqrt") {
      #pragma omp parallel for if(n > 1e+5) schedule(static)
      for(idx_t i=0; i < n; ++i){
        for(idx_t j=0; j < ncols; ++j){
          // scale value but keep sign
          x[i * ncols + j] = x[i * ncols + j] >= 0 ? sqrt(x[i * ncols + j]) : (-sqrt(-x[i * ncols + j]));
        }
      }
    }
  }

  /**
   * Undo the scaling of the columns of the vectors in `x`, given vectors were
   * scaled by `mode`. Mode can be one of:
   *   `stdev` or `zscore`: Standardize the columns, i.e., x' = (x-mean)/stdev (unbiased, divide by n-1)
   *   `stdevp`: Standardize the columns, i.e., x' = (x-mean)/stdev (biased, divide by n)
   *   `minmax` or 'min-max': x' = (x-min)/(max-min)
   *   `mean`: x' = (x-mean)/(max-min)
   *   `sqrt`: x' = sqrt(x)
   *   `log`: x' = 1 + log_2(x)
   *   `maxtf`: x' = .5 + .5 * x/max
   * `x` is assumed to have the same number of columns as this matrix.
   *
   * @param n       Number of vectors to unscale.
   * @param x       Data array of size n * ncols.
   * @param mode    How the vectors were scaled.
   * @param factors Factors used when scaling.
   */
  template <typename T>
  void Matrix<T>::unscale(const idx_t n, T* x, const string mode, const double* factors) const
  {
    throw runtime_error("Not yet implemented.");
  }


  /**
   * Scale the columns of the matrix by the chosen `mode`. Mode can be one of:
   *   `stdev` or `zscore`: Standardize the columns, i.e., x' = (x-mean)/stdev (unbiased, divide by n-1)
   *   `stdevp`: Standardize the columns, i.e., x' = (x-mean)/stdev (biased, divide by n)
   *   `minmax` or 'min-max': x' = (x-min)/(max-min)
   *   `mean`: x' = (x-mean)/(max-min)
   *   `sqrt`: x' = sqrt(x)
   *   `log`: x' = 1 + log_2(x)
   *   `maxtf`: x' = .5 + .5 * x/max
   *
   * @param mode    How to scale.
   * @param factors Optional pointer for return factors (they will be allocated within).
   */
  template <typename T>
  void Matrix<T>::scale(const string mode, double** factors)
  {
    double * cscale;
    compute_scaling_factors(mode, &cscale);
    scale(nrows, data, mode, cscale);
    if(factors){
      *factors = cscale;
    } else {
      free(cscale);
    }
  }

  /**
   * Scale the columns of the matrix by the chosen `mode` using the given factors. Mode can be one of:
   *   `stdev` or `zscore`: Standardize the columns, i.e., x' = (x-mean)/stdev (unbiased, divide by n-1)
   *   `stdevp`: Standardize the columns, i.e., x' = (x-mean)/stdev (biased, divide by n)
   *   `minmax` or 'min-max': x' = (x-min)/(max-min)
   *   `mean`: x' = (x-mean)/(max-min)
   *   `sqrt`: x' = sqrt(x)
   *   `log`: x' = 1 + log_2(x)
   *   `maxtf`: x' = .5 + .5 * x/max
   *
   * @param mode    How to scale.
   * @param factors Factors to scale by (should be obtained via `compute_scaling_factors`).
   */
  template <typename T>
  void Matrix<T>::scale(const string mode, const double* factors)
  {
    scale(nrows, data, mode, factors);
  }


  // Permutation operations
  /**
   * Permute the columns of the matrix given the permutation array perm
   * @param perm Permutation array of size ncols containing new column ids for each column.
   */
  template <typename T>
  void Matrix<T>::permute_columns(const idx_t * perm)
  {
    #pragma omp parallel
    {
      T* row = (T*) malloc(sizeof(T) * ncols);
      if(!row){
        throw runtime_error("Could not allocate temporary array for storing result row.");
      }
      #pragma omp for schedule(static, 20)
      for(idx_t i=0; i < nrows; ++i){
        for(idx_t j=0; j < ncols; ++j){
          row[perm[j]] = data[i*ncols + j];
        }
        memcpy(data + (ncols * i), row, sizeof(T) * ncols);
      }
      free(row);
    }
  }

  // Access operators
  /**
   * @brief Subscript operator
   * @param  i Row ID
   * @param  j Column ID
   * @return Reference to value at MAT[i,j]
   * @thows out_of_range
   */
  template <typename T>
  inline T& Matrix<T>::operator() (
    const idx_t i,
    const idx_t j
  ){
    return data[i * ncols + j];
  }

  /**
   * @brief Subscript operator
   * @param  i Row ID
   * @param  j Column ID
   * @return Value at MAT[i,j]
   * @thows out_of_range
   */
  template <typename T>
  inline T  Matrix<T>::operator() (
    const idx_t i,
    const idx_t j
  ) const
  {
    return data[i * ncols + j];
  }

  // Mathematical operations
  /**
   * @brief Matrix addition
   * @param rhs Matrix to be added
   * @return new matrix consisting of this + rhs
   */
  template <typename T>
  Matrix<T> Matrix<T>::operator+(
    const Matrix<T>& rhs
  ){
    if(nrows != rhs.nrows || ncols != rhs.ncols){
      throw invalid_argument("The rhs matrix must have the same number of rows and columns for this operation.");
    }
    Matrix<T> other = Matrix<T>(this);
    ptr_t nnz = get_nnz();
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      other.data[i] += rhs.data[i];
    }
    return other;
  }

  /**
   * @brief Matrix addition
   * @param rhs Matrix to be added to this one
   * @return reference to this matrix after rhs was added
   */
  template <typename T>
  Matrix<T>& Matrix<T>::operator+=(
    const Matrix<T>& rhs
  ){
    if(nrows != rhs.nrows || ncols != rhs.ncols){
      throw invalid_argument("The rhs matrix must have the same number of rows and columns for this operation.");
    }
    ptr_t nnz = get_nnz();
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      data[i] += rhs.data[i];
    }
    return *this;
  }

  /**
   * @brief Matrix subtraction
   * @param rhs Matrix to be subtracted
   * @return new matrix consisting of this - rhs
   */
  template <typename T>
  Matrix<T> Matrix<T>::operator-(
    const Matrix<T>& rhs
  ){
    if(nrows != rhs.nrows || ncols != rhs.ncols){
      throw invalid_argument("The rhs matrix must have the same number of rows and columns for this operation.");
    }
    Matrix<T> other = Matrix<T>(this);
    ptr_t nnz = get_nnz();
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      other.data[i] -= rhs.data[i];
    }
    return other;
  }

  /**
   * @brief Matrix subtraction
   * @param rhs Matrix to be subtracted from this one
   * @return reference to this matrix after rhs was subtracted
   */
  template <typename T>
  Matrix<T>& Matrix<T>::operator-=(
    const Matrix<T>& rhs
  ){
    if(nrows != rhs.nrows || ncols != rhs.ncols){
      throw invalid_argument("The rhs matrix must have the same number of rows and columns for this operation.");
    }
    ptr_t nnz = get_nnz();
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      data[i] -= rhs.data[i];
    }
    return *this;
  }

  /**
   * @brief Matrix Matrix multiplication
   * @param rhs Matrix to be multiplied with this one
   * @return new matrix consisting of this * rhs
   */
  template <typename T>
  Matrix<T> Matrix<T>::operator*(
    const Matrix<T>& rhs
  ){
    if(ncols != rhs.nrows){
      throw invalid_argument("The lhs ncols must equal rhs nrows for matrix matrix multiplication.");
    }
    Matrix<T> other = Matrix<T>(nrows, other.ncols);
    #pragma omp parallel
    {
      T v;
      const idx_t ncols2 = rhs.ncols;
      #pragma omp for schedule(static, 32)
      for(idx_t i=0; i < nrows; ++i){
        for(idx_t j=0; j < ncols2; ++j){
          // c_{i,j} = \sum_{k=1}^m a_{i,k} b_{k,j}
          v = 0;
          for(idx_t k=0; k < ncols; ++k){
            v += data[i * ncols + k] * rhs.data[k * ncols2 + j];
          }
          other(i,j) = v;
        }
      }
    }
    return other;
  }

  /**
   * @brief Matrix Matrix multiplication of square matrices of same size.
   * @param rhs Matrix to be multiplied with this one, storing the result in this matrix
   * @return reference of this matrix which also stores the result of the multiplication
   */
  template <typename T>
  Matrix<T>& Matrix<T>::operator*=(
    const Matrix<T>& rhs
  ){
    if(nrows != ncols){
      throw invalid_argument("The *= operation is only valid for square matrices.");
    }
    if(ncols != rhs.nrows){
      throw invalid_argument("The lhs dimensions must equal rhs dimensions for the matrix *= operation.");
    }
    #pragma omp parallel
    {
      T* tmp = (T*) malloc(nrows * sizeof(T));
      if(!tmp){
        throw runtime_error("Could not allocate temporary array for storing result row.");
      }
      T v;
      const idx_t ncols2 = rhs.ncols;
      #pragma omp for schedule(static, 32)
      for(idx_t i=0; i < nrows; ++i){
        for(idx_t j=0; j < ncols2; ++j){
          // c_{i,j} = \sum_{k=1}^m a_{i,k} b_{k,j}
          v = 0;
          for(idx_t k=0; k < ncols; ++k){
            v += data[i * ncols + k] * rhs.data[k * ncols2 + j];
          }
          tmp[j] = v;
        }
        for(idx_t j=0; j < ncols2; ++j){
          data[i * ncols + j] = tmp[j];
        }
      }
      free(tmp);
    }
    return *this;
  }

  /**
   * @brief Transpose matrix, i.e., rotate it by 90 degrees.
   * Rows become columns and columns become rows.
   * @return A copy matrix that is the transpose of this one.
   */
  template <typename T>
  Matrix<T>* Matrix<T>::transposed()
  {
    auto other = new Matrix<T>(ncols, nrows);
    #pragma omp parallel for schedule(static, 20)
    for(ptr_t i=0; i < nrows; ++i){
      for(ptr_t j=0; j < ncols; ++j){
        other->data[j * nrows + i] = data[i * ncols + j];
      }
    }
    return other;
  }

  /**
   * @brief Transpose matrix, i.e., rotate it by 90 degrees.
   * Rows become columns and columns become rows.
   * @return Transposes this matrix and returns a reference to this.
   */
  template <typename T>
  Matrix<T>& Matrix<T>::transpose()
  {
    T* d = (T*) malloc(sizeof(T) * get_nnz());
    #pragma omp parallel for schedule(static, 20)
    for(ptr_t i=0; i < nrows; ++i){
      for(ptr_t j=0; j < ncols; ++j){
        d[j * nrows + i] = data[i * ncols + j];
      }
    }
    T* tmp = data;
    data = d;
    free(tmp);
    idx_t n = nrows;
    nrows = ncols;
    ncols = n;
    return *this;
  }

  // Matrix/scalar operations

  /**
   * @brief Matrix scalar addition
   * @param val Value to be added to all elements of matrix
   * @return new matrix consisting of this + val added to all elements
   */
  template <typename T>
  Matrix<T> Matrix<T>::operator+(
    const T val
  ){
    ptr_t nnz = get_nnz();
    Matrix<T> other = Matrix<T>(this);
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      other.data[i] += val;
    }
    return other;
  }

  /**
   * @brief Matrix scalar addition
   * @param val Value to be added to all elements of matrix
   * @return reference to this matrix after val was added
   */
  template <typename T>
  Matrix<T>& Matrix<T>::operator+=(
    const T val
  ){
    ptr_t nnz = get_nnz();
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      data[i] += val;
    }
    return *this;
  }

  /**
   * @brief Matrix scalar subtraction
   * @param val Value to be subtracted from all elements of matrix
   * @return new matrix consisting of elements of this - val
   */
  template <typename T>
  Matrix<T> Matrix<T>::operator-(
    const T val
  ){
    ptr_t nnz = get_nnz();
    Matrix<T> other = Matrix<T>(this);
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      other.data[i] -= val;
    }
    return other;
  }

  /**
   * @brief Matrix scalar subtraction
   * @param val Value to be subtracted from all elements of matrix
   * @return reference to this matrix after val was subtracted
   */
  template <typename T>
  Matrix<T>& Matrix<T>::operator-=(
    const T val
  ){
    ptr_t nnz = get_nnz();
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      data[i] -= val;
    }
    return *this;
  }

  /**
   * @brief Matrix scalar multiplication
   * @param val Value that all elements of matrix should be multiplied by
   * @return new matrix consisting of elements of this multiplied by val
   */
  template <typename T>
  Matrix<T> Matrix<T>::operator*(
    const T val
  ){
    ptr_t nnz = get_nnz();
    Matrix<T> other = Matrix<T>(this);
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      other.data[i] *= val;
    }
    return other;
  }

  /**
   * @brief Matrix scalar multiplication
   * @param val Value that all elements of matrix should be multiplied by
   * @return reference to this matrix after elements were multiplied by val
   */
  template <typename T>
  Matrix<T>& Matrix<T>::operator*=(
    const T val
  ){
    ptr_t nnz = get_nnz();
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      data[i] *= val;
    }
    return *this;
  }

  /**
   * @brief Matrix scalar division
   * @param val Value that all elements of matrix should be devided by
   * @return new matrix consisting of elements of this divided by val
   */
  template <typename T>
  Matrix<T> Matrix<T>::operator/(
    const T val
  ){
    Matrix<T> other = Matrix<T>(this);
    ptr_t nnz = get_nnz();
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      other.data[i] /= val;
    }
    return other;
  }

  /**
   * @brief Matrix scalar division
   * @param val Value that all elements of matrix should be divided by
   * @return reference to this matrix after elements were divided by val
   */
  template <typename T>
  Matrix<T>& Matrix<T>::operator/=(
    const T val
  ){
    ptr_t nnz = get_nnz();
    #pragma omp parallel for schedule(static, 1024)
    for(ptr_t i=0; i < nnz; ++i){
      data[i] /= val;
    }
    return *this;
  }

  /**
   * @brief check if this matrix is equal to rhs
   * @param rhs Comparison matrix
   * @return true if matrices are equal
   */
  template <typename T>
  bool Matrix<T>::operator==(
    const Matrix<T>& rhs
  ) const
  {
    if(nrows != rhs.nrows || ncols != rhs.ncols){
      return false;
    }
    uint8_t t = Matrix<T>::get_type_id();
    ptr_t nnz = get_nnz();
    if(t > 10){
      for(ptr_t i=0; i < nnz; ++i){
        if(data[i] > rhs.data[i] ? data[i] - rhs.data[i] > FASTA_MATRIX_EPSILON : \
          rhs.data[i] - data[i] > FASTA_MATRIX_EPSILON ){
          return false;
        }
      }
      return true;
    }
    for(ptr_t i=0; i < nnz; ++i){
      if(data[i] != rhs.data[i]){
        return false;
      }
    }
    return true;
  }

  template <typename T>
  std::ostream& operator<< (std::ostream &out, const Matrix<T>& mat) {
    for(idx_t i=0; i < mat.nrows; ++i){
      for(idx_t j=0; j < mat.ncols; ++j){
        out << mat(i,j) << " ";
      }
      out << endl;
    }
    return out;
  }

} // namespace fasta
