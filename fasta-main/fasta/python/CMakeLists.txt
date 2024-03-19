# Copyright (c) David C. Anastasiu & Gheorghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

project(pyfasta
  DESCRIPTION "Python bindings for fasta."
  HOMEPAGE_URL "https://github.com/davidanastasiu/fasta"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(SWIG REQUIRED COMPONENTS python)
include(${SWIG_USE_FILE})

set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)
set(SWIG_SOURCE_FILE_EXTENSIONS swig)

macro(configure_swigfasta source)
  set_source_files_properties(${source} PROPERTIES
    CPLUSPLUS ON
    USE_TARGET_INCLUDE_DIRECTORIES TRUE
  )
  if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND NOT WIN32)
    set_source_files_properties(${source} PROPERTIES
      SWIG_FLAGS -DSWIGWORDSIZE64
    )
  endif()
  if(WIN32)
    set_source_files_properties(${source} PROPERTIES
      SWIG_FLAGS -DSWIGWIN
    )
  endif()
endmacro()

configure_swigfasta(swigfasta.swig)

if(TARGET fasta)
  # Manually add headers as extra dependencies of swigfasta.
  set(SWIG_MODULE_swigfasta_EXTRA_DEPS)
  foreach(h ${FAISS_HEADERS})
    list(APPEND SWIG_MODULE_swigfasta_EXTRA_DEPS "${fasta_SOURCE_DIR}/fasta/${h}")
  endforeach()
else()
  find_package(fasta REQUIRED)
endif()

swig_add_library(swigfasta
  TYPE SHARED
  LANGUAGE python
  SOURCES swigfasta.swig
)

if(NOT WIN32)
  # NOTE: Python does not recognize the dylib extension.
  set_target_properties(swigfasta PROPERTIES SUFFIX .so)
endif()

find_package(OpenMP REQUIRED)

target_link_libraries(swigfasta PRIVATE
  fasta
  Python::Module
  Python::NumPy
  OpenMP::OpenMP_CXX
)

# Hack so that python_callbacks.h can be included as
# `#include <fasta/python/python_callbacks.h>`.
target_include_directories(swigfasta PRIVATE ${PROJECT_SOURCE_DIR}/../..)

find_package(Python REQUIRED
  COMPONENTS Development NumPy
)

add_library(fasta_python_callbacks EXCLUDE_FROM_ALL
  python_callbacks.cpp
)
set_property(TARGET fasta_python_callbacks
  PROPERTY POSITION_INDEPENDENT_CODE ON
)

# Hack so that python_callbacks.h can be included as
# `#include <fasta/python/python_callbacks.h>`.
target_include_directories(fasta_python_callbacks PRIVATE ${PROJECT_SOURCE_DIR}/../..)
target_include_directories(fasta_python_callbacks PRIVATE ${Python_INCLUDE_DIRS})

target_link_libraries(swigfasta PRIVATE fasta_python_callbacks)

configure_file(setup.py setup.py COPYONLY)
configure_file(__init__.py __init__.py COPYONLY)

file(GLOB files "${PROJECT_SOURCE_DIR}/../../pytools/*.py")
file(COPY ${files} DESTINATION pytools/)
