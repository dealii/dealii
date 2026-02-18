## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2018 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Configuration for complex value support
#

if(${DEAL_II_WITH_COMPLEX_VALUES})
   set(DEAL_II_EXPAND_COMPLEX_SCALARS
       "std::complex<double>"
       "std::complex<float>"
      )
   set(DEAL_II_EXPAND_COMPLEX_VECTORS
       "Vector<std::complex<double> >"
       "Vector<std::complex<float> >"
      )
   set(DEAL_II_EXPAND_COMPLEX_BLOCK_VECTORS
       "BlockVector<std::complex<double> >"
       "BlockVector<std::complex<float> >"
      )
   set(DEAL_II_EXPAND_COMPLEX_LA_PARALLEL_VECTORS
       "LinearAlgebra::distributed::Vector<std::complex<double> >"
       "LinearAlgebra::distributed::Vector<std::complex<float> >"
      )
   set(DEAL_II_EXPAND_COMPLEX_LA_PARALLEL_BLOCK_VECTORS
       "LinearAlgebra::distributed::BlockVector<std::complex<double> >"
       "LinearAlgebra::distributed::BlockVector<std::complex<float> >"
      )
endif()
