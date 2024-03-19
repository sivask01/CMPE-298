/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#pragma once
#include <cstdint>

#define FASTA_VERSION_MAJOR 0
#define FASTA_VERSION_MINOR 0
#define FASTA_VERSION_PATCH 1

namespace fasta
{
  static constexpr struct {
    uint8_t major, minor, patch;
  } version = {
    FASTA_VERSION_MAJOR,
    FASTA_VERSION_MINOR,
    FASTA_VERSION_PATCH
  };
} // namespace cxxopts
