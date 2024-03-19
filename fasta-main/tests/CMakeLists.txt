# Copyright (c) David C. Anastasiu & Ghoerghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

cmake_minimum_required(VERSION 3.17...3.21 FATAL_ERROR)

enable_testing()

# file(GLOB FASTA_TEST_SRC CONFIGURE_DEPENDS
#   ${CMAKE_CURRENT_LIST_DIR}/test_*.cpp
# )
SET(FASTA_TEST_SRC
  test_bit_vector.cpp
  test_kv_vector.cpp
  test_kvv_view.cpp
  test_matrix.cpp
)
message(STATUS "Fasta tests: ${FASTA_TEST_SRC}")

add_executable(fasta_test ${FASTA_TEST_SRC})
set_target_properties(fasta_test PROPERTIES
  CXX_STANDARD 20
)
target_link_libraries(fasta_test PRIVATE fmt spdlog fasta)

include(FetchContent)
FetchContent_Declare(googletest
  GIT_REPOSITORY    https://github.com/google/googletest.git
  GIT_TAG           389cb68b87193358358ae87cc56d257fd0d80189
  SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-src"
  BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)
set(BUILD_GMOCK CACHE BOOL OFF)
set(INSTALL_GTEST CACHE BOOL OFF)
FetchContent_MakeAvailable(googletest)

find_package(OpenMP REQUIRED)

target_link_libraries(fasta_test PRIVATE
  OpenMP::OpenMP_CXX
  gtest_main
)

# Defines `gtest_discover_tests()`.
include(GoogleTest)
gtest_discover_tests(fasta_test)
