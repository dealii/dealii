//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <dofs/block_info.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_tools.h>
#include <fe/fe.h>
#include <fe/fe_tools.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
void
BlockInfo::initialize(const DoFHandler<dim, spacedim>& dof)
{
  const FiniteElement<dim, spacedim>& fe = dof.get_fe();
  std::vector<unsigned int> sizes(fe.n_blocks());
  DoFTools::count_dofs_per_block(dof, sizes);
  bi_global.reinit(sizes);
}


template <int dim, int spacedim>
void
BlockInfo::initialize_local(const DoFHandler<dim, spacedim>& dof)
{
  const FiniteElement<dim, spacedim>& fe = dof.get_fe();
  std::vector<unsigned int> sizes(fe.n_blocks());

  base_element.resize(fe.n_blocks());
  
  for (unsigned int i=0;i<base_element.size();++i)
    base_element[i] = fe.block_to_base_index(i).first;
  
  local_renumbering.resize(fe.n_dofs_per_cell());
  FETools::compute_block_renumbering(fe,
				     local_renumbering,
				     sizes, false);
  bi_local.reinit(sizes);
}


template <int dim, int spacedim>
void
BlockInfo::initialize(const MGDoFHandler<dim, spacedim>& dof)
{
  initialize((DoFHandler<dim, spacedim>&) dof);
  
  std::vector<std::vector<unsigned int> > sizes;
  MGTools::count_dofs_per_block(dof, sizes);
  
  levels.resize(sizes.size());
  for (unsigned int i=0;i<sizes.size();++i)
    levels[i].reinit(sizes[i]);
}


template void BlockInfo::initialize(const DoFHandler<deal_II_dimension,deal_II_dimension>&);
template void BlockInfo::initialize(const MGDoFHandler<deal_II_dimension,deal_II_dimension>&);
template void BlockInfo::initialize_local(const DoFHandler<deal_II_dimension,deal_II_dimension>&);


DEAL_II_NAMESPACE_CLOSE
