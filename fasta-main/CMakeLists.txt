# Copyright (c) David C. Anastasiu & Gheorghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

cmake_minimum_required(VERSION 3.17...3.21 FATAL_ERROR)

# Get the version of the library
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(fasta)
fasta_getversion(VERSION)

project(fasta
  VERSION ${VERSION}
  DESCRIPTION "Efficient exact nearest neighor search of dense data using cosine similarity."
  HOMEPAGE_URL "https://github.com/davidanastasiu/fasta"
  LANGUAGES CXX)
include(GNUInstallDirs)

# Include guards
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

#--- Add dependencies ---
include(FetchContent)

# Lightweight C++ command line option parser   https://github.com/jarro2783/cxxopts
FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
  GIT_TAG        v2.2.1
  GIT_SHALLOW    TRUE
)
set(CXXOPTS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(CXXOPTS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(CXXOPTS_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
set(CXXOPTS_ENABLE_WARNINGS OFF CACHE BOOL "" FORCE)
# FetchContent_MakeAvailable(cxxopts)
FetchContent_GetProperties(cxxopts)
if(NOT cxxopts_POPULATED)
    FetchContent_Populate(cxxopts)
endif()
add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR})
include_directories(${cxxopts_SOURCE_DIR}/include)

# Logging
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.9.2
		GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(spdlog)
if (NOT spdlog_POPULATED)
    FetchContent_Populate(spdlog)
    add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR})
endif()

# String formatting
FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG        6.1.2
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(fmt)

#--- Project options ---
option(FASTA_ENABLE_WARNINGS "Add warnings to CMAKE_CXX_FLAGS." ON)
option(FASTA_ENABLE_DEMOS "Build demos." ON)
option(FASTA_ENABLE_PYTHON "Build Python extension." ON)
option(FASTA_BUILD_TESTING "Build tests." ON)

# Set CXX standard and optionally enable warnings
fasta_set_cxx_standard()
if (FASTA_ENABLE_WARNINGS)
    fasta_enable_warnings()
endif()

# Add targets
message(STATUS "Adding target: fasta")
add_subdirectory(fasta)

if(FASTA_ENABLE_PYTHON)
  add_subdirectory(fasta/python)
endif()
if(FASTA_ENABLE_DEMOS)
  message(STATUS "Adding target: demos")
  add_subdirectory(demos)
endif()

# CTest must be included in the top level to enable the `make test` target.
include(CTest)
if(FASTA_BUILD_TESTING)
  message(STATUS "Adding target: tests")
  add_subdirectory(tests)
endif()

message(STATUS "Disabling CXX conversion warning due to cxxopts.")
if (FASTA_ENABLE_WARNINGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-conversion")
endif()
message(STATUS "CMAKE_CXX_FLAGS before exec target: ${CMAKE_CXX_FLAGS}")

message(STATUS "Adding target: exec")
add_subdirectory(exec)

#--- Some status information ---
message(STATUS "OS:                ${CMAKE_SYSTEM_NAME}")
message(STATUS "Compiler:          ${CMAKE_CXX_COMPILER_ID} ${CXX_COMPILER_VERSION}")
message(STATUS "CXX flags release: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CXX flags debug:   ${CMAKE_CXX_FLAGS_DEBUG}")
