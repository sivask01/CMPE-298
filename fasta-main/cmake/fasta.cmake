# Copyright (c) David C. Anastasiu & Gheorghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.
#
# Helpter functions
#

include_guard()

# Request C++20 without gnu extensions
macro(fasta_set_cxx_standard)
  if (FASTA_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD ${FASTA_CXX_STANDARD})
  else()
    set(CMAKE_CXX_STANDARD 11)
  endif()
  set(CMAKE_CXX_EXTENSIONS OFF)
endmacro()

# Helper function to enable warnings
function(fasta_enable_warnings)
  if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W2")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "[Cc]lang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wextra -Wno-suggest-override -Wno-shadow -Wsign-compare -Wwrite-strings -Wpointer-arith -Winit-self -Wconversion -Wno-sign-conversion")
  endif()

  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} PARENT_SCOPE)
endfunction()

function(fasta_getversion version_arg)
  # Parse the current version from the fasta header
  file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/fasta/version.h" fasta_version_defines
    REGEX "#define FASTA_VERSION_(MAJOR|MINOR|PATCH)")
  foreach(ver ${fasta_version_defines})
    if(ver MATCHES "#define FASTA_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
      set(FASTA_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
  endforeach()
  set(VERSION ${FASTA_VERSION_MAJOR}.${FASTA_VERSION_MINOR}.${FASTA_VERSION_PATCH})

  message(STATUS "fasta version ${VERSION}")

  # Return the information to the caller
  set(${version_arg} ${VERSION} PARENT_SCOPE)
endfunction()
