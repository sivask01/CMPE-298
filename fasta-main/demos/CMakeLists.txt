# Copyright (c) David C. Anastasiu & Gheorghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

file(GLOB FASTA_DEMOS CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_LIST_DIR}/*.cpp
)

foreach(program ${FASTA_DEMOS})
  get_filename_component(program_name ${program} NAME_WE )
  message(STATUS "  Setting up demo: ${program_name}")
  add_executable(${program_name} EXCLUDE_FROM_ALL ${program})
  target_link_libraries(${program_name} PRIVATE fasta)
endforeach()
