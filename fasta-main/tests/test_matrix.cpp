/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#include <iostream>
#include <string>
#include <gtest/gtest.h>

#include <fasta/Matrix.h>

using namespace fasta;

TEST(Matrix, binary_io_and_equality) {
  Matrix<float> mat(3, 5);
  mat.set_random();
  mat.write("mat.bmat");
  // std::cout << "Mat: " << std::endl << mat;
  auto mat2 = Matrix<float>::read("mat.bmat");
  // std::cout << "Mat2: " << std::endl << *mat2;
  ASSERT_TRUE(mat2->nrows == 3);
  ASSERT_TRUE(mat2->ncols == 5);
  ASSERT_TRUE(mat == *mat2);
  delete mat2;
}

TEST(Matrix, normalize) {
  Matrix<float> mat(100, 10000);
  mat.set_random();
  mat.normalize();
  for(uint32_t i=0; i < mat.nrows; ++i){
    double norm = 0.0;
    for(uint32_t j=0; j < mat.ncols; ++j){
      norm += mat(i,j) * mat(i,j);
    }
    norm = std::sqrt(norm);
    ASSERT_NEAR(norm, 1.0, 1e-5);
  }
}

TEST(Matrix, is_normalized) {
  Matrix<float> mat(10, 10000);
  for(uint32_t i=0; i < 10; ++i){
    mat.set_random();
    mat.normalize();
    ASSERT_TRUE(mat.is_normalized());
  }
}

TEST(Matrix, transpose) {
  Matrix<int> mat(305, 434);
  mat.set_random();
  auto mat2 = mat;
  mat2.transpose();
  ASSERT_EQ(mat.nrows, mat2.ncols);
  ASSERT_EQ(mat.ncols, mat2.nrows);
  for(idx_t i=0; i < mat.nrows; ++i){
    for(idx_t j=0; j < mat.ncols; ++j){
      ASSERT_EQ(mat(i,j), mat2(j,i));
    }
  }
}

TEST(Matrix, transposed) {
  Matrix<int> mat(305, 434);
  mat.set_random();
  auto mat2 = mat.transposed();
  ASSERT_EQ(mat.nrows, mat2->ncols);
  ASSERT_EQ(mat.ncols, mat2->nrows);
  for(idx_t i=0; i < mat.nrows; ++i){
    for(idx_t j=0; j < mat.ncols; ++j){
      ASSERT_EQ(mat(i,j), mat2->operator()(j,i));
    }
  }
  delete mat2;
}

TEST(Matrix, permute) {
  idx_t nrows = 305;
  idx_t ncols = 434;
  Matrix<int> mat(nrows, ncols);
  // create permutation vector
  idx_t perm[ncols];
  for (idx_t i=0; i < ncols; ++i){
    perm[i] = i;
  }
  std::random_shuffle (perm, perm + ncols);
  // make a copy of mat
  auto mat2 = mat;
  for(idx_t i=0; i < mat.nrows; ++i){
    for(idx_t j=0; j < mat.ncols; ++j){
      ASSERT_EQ(mat(i,j), mat2(i,j));
    }
  }
  //permute mat2
  mat2.permute_columns(perm);
  ASSERT_EQ(mat.nrows, mat2.nrows);
  ASSERT_EQ(mat.ncols, mat2.ncols);
  for(idx_t i=0; i < mat.nrows; ++i){
    for(idx_t j=0; j < mat.ncols; ++j){
      ASSERT_EQ(mat(i,j), mat2(i,perm[i]));
    }
  }
}

TEST(Matrix, add_rows) {
  Matrix<int> mat(305, 434);
  mat.set_random();
  Matrix<int> mat2(23, 434);
  mat.add(mat2);
  ASSERT_EQ(mat.nrows, 328);
  ASSERT_EQ(mat.get_max_nrows(), 328);
  ASSERT_EQ(mat.ncols, 434);
  for(idx_t i=305; i < mat.nrows; ++i){
    for(idx_t j=0; j < mat.ncols; ++j){
      ASSERT_EQ(mat(i,j), 0);
    }
  }
  mat2.add(1, mat.data);
  ASSERT_EQ(mat2.nrows, 24);
  ASSERT_EQ(mat2.get_max_nrows(), 46);
  for(idx_t j=0; j < mat.ncols; ++j){
    ASSERT_EQ(mat(0,j), mat2(23,j));
  }
}

TEST(Matrix, copy_constructor) {
  Matrix<int> mat(305, 434);
  mat.set_random();
  auto mat2 = new Matrix<int>(mat);
  ASSERT_EQ(mat2->nrows, 305);
  ASSERT_EQ(mat2->get_max_nrows(), 305);
  ASSERT_EQ(mat2->ncols, 434);
  for(idx_t i=0; i < mat.nrows; ++i){
    for(idx_t j=0; j < mat.ncols; ++j){
      ASSERT_EQ(mat(i,j), mat2->operator()(i,j));
    }
  }
  delete mat2;
}

TEST(Matrix, scale) {
  Matrix<int> mat(5, 104754);
  mat.set_random();
  double * factors = nullptr;
  mat.scale("minmax", &factors);
  ASSERT_NE(nullptr, factors);
  for(idx_t i=0; i < mat.nrows; ++i){
    for(idx_t j=0; j < mat.ncols; ++j){
      ASSERT_LE(mat(i,j), 1.00000001);
      ASSERT_GE(mat(i,j), -0.00000001);
    }
  }
  free(factors);
}

