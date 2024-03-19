# Copyright (c) David C. Anastasiu & Gheorghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# file(GLOB FASTA_HEADERS CONFIGURE_DEPENDS
#   ${CMAKE_CURRENT_SOURCE_DIR}/*.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/impl/*.h
# )
SET(FASTA_HEADERS
  Fasta.h
  GemmIndex.h
  Index.h
  Matrix.h
  types.h
  version.h
  impl/Fasta.h
  impl/Index.h
  impl/Matrix.h
  utils/BitVector.h
  utils/Counters.h
  utils/KVVector.h
  utils/KVVectorView.h
  utils/sort.h
  utils/Timers.h
)

# Export FASTA_HEADERS variable to parent scope.
set(FASTA_HEADERS ${FASTA_HEADERS} PARENT_SCOPE)

# Add the library
add_library(fasta INTERFACE)

# Handle `#include <fasta/foo.h>`.
target_include_directories(fasta SYSTEM INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

find_package(OpenMP REQUIRED)
# set(CMAKE_CXX_FLAGS "-g ${CMAKE_CXX_FLAGS} -Wall -Werror -Wextra -Wno-suggest-override -Wno-shadow -Wsign-compare -Wwrite-strings -Wpointer-arith -Winit-self -Wconversion -Wno-sign-conversion ${COMPILER_SPECIFIC_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
target_link_libraries(fasta INTERFACE OpenMP::OpenMP_CXX fmt spdlog)

### Installation instructions ###

install(TARGETS fasta fmt spdlog
  EXPORT fasta-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
foreach(header ${FASTA_HEADERS})
  get_filename_component(dir ${header} DIRECTORY )
  install(FILES ${header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fasta/${dir}
  )
endforeach()

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/cmake/fasta-config-version.cmake"
  VERSION ${CMAKE_PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_file(${PROJECT_SOURCE_DIR}/cmake/fasta-config.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/fasta-config.cmake
  COPYONLY
)

install(FILES ${PROJECT_BINARY_DIR}/cmake/fasta-config.cmake
  ${PROJECT_BINARY_DIR}/cmake/fasta-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/fasta
)

install(EXPORT fasta-targets
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/fasta
)
