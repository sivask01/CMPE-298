# Copyright (c) David C. Anastasiu & Gheorghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

project(exec
  DESCRIPTION "Command-line executable for fasta."
  HOMEPAGE_URL "https://github.com/davidanastasiu/fasta"
  LANGUAGES CXX)

# Add dependencies

# Set sources
# file(GLOB FASTA_MAIN_SRC CONFIGURE_DEPENDS
#   ${CMAKE_CURRENT_LIST_DIR}/*.cpp
# )
SET(FASTA_MAIN_SRC
  main.cpp
)

# add executable target

add_executable(${PROJECT_NAME} EXCLUDE_FROM_ALL ${FASTA_MAIN_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES
  CXX_STANDARD 20
  OUTPUT_NAME fasta
)
target_link_libraries(${PROJECT_NAME} fmt spdlog fasta)
