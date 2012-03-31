//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2009, 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/meshworker/dof_info.h>
#include <deal.II/base/quadrature_lib.h>

DEAL_II_NAMESPACE_OPEN


namespace MeshWorker
{
  template <int dim, int spacedim, typename number>
  DoFInfo<dim,spacedim,number>::DoFInfo(const BlockInfo& info)
                  : block_info(&info, typeid(*this).name())
  {
    indices_by_block.resize(info.local().size());
    for (unsigned int i=0;i<indices_by_block.size();++i)
      indices_by_block[i].resize(info.local().block_size(i));
  }


  template <int dim, int spacedim, typename number>
  DoFInfo<dim,spacedim,number>::DoFInfo()
  {}


  template <int dim, int spacedim, typename number>
  void
  DoFInfo<dim,spacedim,number>::set_block_indices()
  {
    for (unsigned int i=0;i<indices.size();++i)
      {
        const std::pair<unsigned int, unsigned int>
          bi = block_info->local().global_to_local(this->block_info->renumber(i));
        indices_by_block[bi.first][bi.second] = indices_org[i];
      }
                                     // Remove this after
                                     // changing block codes
    for (unsigned int i=0;i<indices.size();++i)
      indices[this->block_info->renumber(i)] = indices_org[i];
  }


  template <int dim, int spacedim, typename number>
  void
  DoFInfo<dim,spacedim,number>::get_indices(const typename DoFHandler<dim, spacedim>::cell_iterator& c)
  {
    if (!c->has_children())
      {
        indices.resize(c->get_fe().dofs_per_cell);

        if (block_info == 0 || block_info->local().size() == 0)
          c->get_dof_indices(indices);
        else
          {
            indices_org.resize(c->get_fe().dofs_per_cell);
            c->get_dof_indices(indices_org);
            set_block_indices();
          }
      }
    else
      indices.resize(0);

    level_cell = false;
  }


  template <int dim, int spacedim, typename number>
  void
  DoFInfo<dim,spacedim,number>::get_indices(const typename MGDoFHandler<dim, spacedim>::cell_iterator& c)
  {
    indices.resize(c->get_fe().dofs_per_cell);

    if (block_info == 0 || block_info->local().size() == 0)
      c->get_mg_dof_indices(indices);
    else
      {
        indices_org.resize(c->get_fe().dofs_per_cell);
        c->get_mg_dof_indices(indices_org);
        set_block_indices();
      }
    level_cell = true;
  }

}


DEAL_II_NAMESPACE_CLOSE

