//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

// Driver file for MGTools::Make_*_sparsity_* routines.

// Call this file after defining PATTERN to the desired sparsity
// pattern type.


template void
MGTools::make_sparsity_pattern<deal_II_dimension> (
  const MGDoFHandler<deal_II_dimension> &,
  PATTERN &,
  const unsigned int);

template void
MGTools::make_flux_sparsity_pattern<deal_II_dimension> (
  const MGDoFHandler<deal_II_dimension> &,
  PATTERN &,
  const unsigned int);

template void
MGTools::make_flux_sparsity_pattern_edge<deal_II_dimension> (
  const MGDoFHandler<deal_II_dimension> &,
  PATTERN &,
  const unsigned int);

#if deal_II_dimension > 1

template void
MGTools::make_flux_sparsity_pattern<deal_II_dimension> (
  const MGDoFHandler<deal_II_dimension> &,
  PATTERN &,
  const unsigned int,
  const Table<2,DoFTools::Coupling>&,
  const Table<2,DoFTools::Coupling>&);

template void
MGTools::make_flux_sparsity_pattern_edge<deal_II_dimension> (
  const MGDoFHandler<deal_II_dimension> &,
  PATTERN &,
  const unsigned int,
  const Table<2,DoFTools::Coupling>&);

#endif

