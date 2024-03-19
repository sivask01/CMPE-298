# Copyright (c) David C. Anastasiu & Gheorghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

import logging

logger = logging.getLogger(__name__)

# we import * so that the symbol X can be accessed as fasta.X
logger.info("Loading fasta.")
from .swigfaiss import *
logger.info("Successfully loaded fasta.")
