// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_SACADO

#  include <deal.II/differentiation/ad/sacado_number_types.h>

DEAL_II_NAMESPACE_OPEN

/*---------------------- Explicit Instantiations ----------------------*/

#  include "sacado_number_types.inst1"
#  ifdef DEAL_II_TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
#    include "sacado_number_types.inst2"
#  endif

DEAL_II_NAMESPACE_CLOSE

#endif
