## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------
 
#
# Configuration for complex value support
#

IF(${DEAL_II_WITH_COMPLEX_VALUES})
   SET(DEAL_II_EXPAND_COMPLEX_SCALARS 
       "std::complex<double>"
       "std::complex<float>"
      )
   SET(DEAL_II_EXPAND_COMPLEX_VECTORS 
       "Vector<std::complex<double> >"
       "Vector<std::complex<float> >"
      )
   SET(DEAL_II_EXPAND_COMPLEX_BLOCK_VECTORS
       "BlockVector<std::complex<double> >"
       "BlockVector<std::complex<float> >"
      )
   SET(DEAL_II_EXPAND_COMPLEX_LA_VECTORS
       "LinearAlgebra::Vector<std::complex<double> >"
       "LinearAlgebra::Vector<std::complex<float> >"
      )
   SET(DEAL_II_EXPAND_COMPLEX_LA_PARALLEL_VECTORS
       "LinearAlgebra::distributed::Vector<std::complex<double> >"
       "LinearAlgebra::distributed::Vector<std::complex<float> >"
      )
   SET(DEAL_II_EXPAND_COMPLEX_LA_PARALLEL_BLOCK_VECTORS
       "LinearAlgebra::distributed::BlockVector<std::complex<double> >"
       "LinearAlgebra::distributed::BlockVector<std::complex<float> >"
      )
ENDIF()

