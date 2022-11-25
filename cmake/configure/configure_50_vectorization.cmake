## ---------------------------------------------------------------------
##
## Copyright (C) 2019 - 2020 by the deal.II authors
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
# Configuration for real scalar vectorization
#


set(DEAL_II_EXPAND_REAL_SCALARS_VECTORIZED "VectorizedArray<double,1>" "VectorizedArray<float,1>")
set(DEAL_II_EXPAND_FLOAT_VECTORIZED "VectorizedArray<float,1>")

if(${DEAL_II_VECTORIZATION_WIDTH_IN_BITS} GREATER 0)
   set(DEAL_II_EXPAND_REAL_SCALARS_VECTORIZED
      "${DEAL_II_EXPAND_REAL_SCALARS_VECTORIZED}" "VectorizedArray<double,2>" "VectorizedArray<float,4>")
   set(DEAL_II_EXPAND_FLOAT_VECTORIZED  "${DEAL_II_EXPAND_FLOAT_VECTORIZED}" "VectorizedArray<float,4>")
endif()

if((${DEAL_II_VECTORIZATION_WIDTH_IN_BITS} GREATER 128))
   set(DEAL_II_EXPAND_REAL_SCALARS_VECTORIZED
      "${DEAL_II_EXPAND_REAL_SCALARS_VECTORIZED}" "VectorizedArray<double,4>" "VectorizedArray<float,8>")
   set(DEAL_II_EXPAND_FLOAT_VECTORIZED  "${DEAL_II_EXPAND_FLOAT_VECTORIZED}" "VectorizedArray<float,8>")
endif()

if((${DEAL_II_VECTORIZATION_WIDTH_IN_BITS} GREATER 256))
   set(DEAL_II_EXPAND_REAL_SCALARS_VECTORIZED
      "${DEAL_II_EXPAND_REAL_SCALARS_VECTORIZED}" "VectorizedArray<double,8>" "VectorizedArray<float,16>")
   set(DEAL_II_EXPAND_FLOAT_VECTORIZED  "${DEAL_II_EXPAND_FLOAT_VECTORIZED}" "VectorizedArray<float,16>")
endif()
