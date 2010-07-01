//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
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

  base_elements.resize(fe.n_blocks());

  for (unsigned int i=0;i<base_elements.size();++i)
    base_elements[i] = fe.block_to_base_index(i).first;

  local_renumbering.resize(fe.n_dofs_per_cell());
  FETools::compute_block_renumbering(fe,
				     local_renumbering,
				     sizes, false);
  bi_local.reinit(sizes);
}


template <int dim, int spacedim>
void
BlockInfo::initialize(const MGDoFHandler<dim, spacedim>& dof, bool levels_only)
{
  if (!levels_only)
    initialize((DoFHandler<dim, spacedim>&) dof);

  std::vector<std::vector<unsigned int> > sizes (dof.get_tria().n_levels());
  for (unsigned int i=0; i<sizes.size(); ++i)
    sizes[i].resize (dof.get_fe().n_blocks());
  MGTools::count_dofs_per_block(dof, sizes);

  levels.resize(sizes.size());
  for (unsigned int i=0;i<sizes.size();++i)
    levels[i].reinit(sizes[i]);
}


template void BlockInfo::initialize(const DoFHandler<deal_II_dimension,deal_II_dimension>&);
template void BlockInfo::initialize(const MGDoFHandler<deal_II_dimension,deal_II_dimension>&, bool);
template void BlockInfo::initialize_local(const DoFHandler<deal_II_dimension,deal_II_dimension>&);

#if deal_II_dimension==1 || deal_II_dimension==2
template void BlockInfo::initialize(const DoFHandler<deal_II_dimension,deal_II_dimension+1>&);
template void BlockInfo::initialize(const MGDoFHandler<deal_II_dimension,deal_II_dimension+1>&, bool);
template void BlockInfo::initialize_local(const DoFHandler<deal_II_dimension,deal_II_dimension+1>&);
#endif


DEAL_II_NAMESPACE_CLOSE
