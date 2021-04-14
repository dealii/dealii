// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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


#include <deal.II/distributed/cell_weights.h>

#include <deal.II/dofs/dof_accessor.h>


DEAL_II_NAMESPACE_OPEN


namespace parallel
{
  template <int dim, int spacedim>
  CellWeights<dim, spacedim>::CellWeights(
    const dealii::DoFHandler<dim, spacedim> &dof_handler,
    const WeightingFunction &                weighting_function)
  {
    reinit(dof_handler, weighting_function);
  }



  template <int dim, int spacedim>
  CellWeights<dim, spacedim>::~CellWeights()
  {
    connection.disconnect();
  }



  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::reinit(
    const DoFHandler<dim, spacedim> &dof_handler,
    const typename CellWeights<dim, spacedim>::WeightingFunction
      &weighting_function)
  {
    connection.disconnect();
    connection = dof_handler.get_triangulation().signals.cell_weight.connect(
      make_weighting_callback(dof_handler, weighting_function));
  }



  // ---------- static member functions ----------

  // ---------- selection of weighting functions ----------

  template <int dim, int spacedim>
  typename CellWeights<dim, spacedim>::WeightingFunction
  CellWeights<dim, spacedim>::constant_weighting(const unsigned int factor)
  {
    return [factor](const typename DoFHandler<dim, spacedim>::cell_iterator &,
                    const FiniteElement<dim, spacedim> &) -> unsigned int {
      return factor;
    };
  }



  template <int dim, int spacedim>
  typename CellWeights<dim, spacedim>::WeightingFunction
  CellWeights<dim, spacedim>::ndofs_weighting(
    const std::pair<float, float> &coefficients)
  {
    return [coefficients](
             const typename DoFHandler<dim, spacedim>::cell_iterator &,
             const FiniteElement<dim, spacedim> &future_fe) -> unsigned int {
      const float result =
        std::trunc(coefficients.first *
                   std::pow(future_fe.n_dofs_per_cell(), coefficients.second));

      Assert(result >= 0. &&
               result <=
                 static_cast<float>(std::numeric_limits<unsigned int>::max()),
             ExcMessage(
               "Cannot cast determined weight for this cell to unsigned int!"));

      return static_cast<unsigned int>(result);
    };
  }



  template <int dim, int spacedim>
  typename CellWeights<dim, spacedim>::WeightingFunction
  CellWeights<dim, spacedim>::ndofs_weighting(
    const std::vector<std::pair<float, float>> &coefficients)
  {
    return [coefficients](
             const typename DoFHandler<dim, spacedim>::cell_iterator &,
             const FiniteElement<dim, spacedim> &future_fe) -> unsigned int {
      float result = 0;
      for (const auto &pair : coefficients)
        result +=
          pair.first * std::pow(future_fe.n_dofs_per_cell(), pair.second);
      result = std::trunc(result);

      Assert(result >= 0. &&
               result <=
                 static_cast<float>(std::numeric_limits<unsigned int>::max()),
             ExcMessage(
               "Cannot cast determined weight for this cell to unsigned int!"));

      return static_cast<unsigned int>(result);
    };
  }



  // ---------- handling callback functions ----------

  template <int dim, int spacedim>
  std::function<unsigned int(
    const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
    const typename dealii::Triangulation<dim, spacedim>::CellStatus     status)>
  CellWeights<dim, spacedim>::make_weighting_callback(
    const DoFHandler<dim, spacedim> &dof_handler,
    const typename CellWeights<dim, spacedim>::WeightingFunction
      &weighting_function)
  {
    const parallel::TriangulationBase<dim, spacedim> *tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &(dof_handler.get_triangulation()));

    Assert(
      tria != nullptr,
      ExcMessage(
        "parallel::CellWeights requires a parallel::TriangulationBase object."));

    return
      [&dof_handler, tria, weighting_function](
        const typename dealii::Triangulation<dim, spacedim>::cell_iterator
          &                                                             cell,
        const typename dealii::Triangulation<dim, spacedim>::CellStatus status)
        -> unsigned int {
        return CellWeights<dim, spacedim>::weighting_callback(
          cell,
          status,
          std::cref(dof_handler),
          std::cref(*tria),
          weighting_function);
      };
  }



  template <int dim, int spacedim>
  unsigned int
  CellWeights<dim, spacedim>::weighting_callback(
    const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell_,
    const typename dealii::Triangulation<dim, spacedim>::CellStatus     status,
    const DoFHandler<dim, spacedim> &                 dof_handler,
    const parallel::TriangulationBase<dim, spacedim> &triangulation,
    const typename CellWeights<dim, spacedim>::WeightingFunction
      &weighting_function)
  {
    // Check if we are still working with the correct combination of
    // Triangulation and DoFHandler.
    Assert(&triangulation == &(dof_handler.get_triangulation()),
           ExcMessage(
             "Triangulation associated with the DoFHandler has changed!"));
    (void)triangulation;

    // Skip if the DoFHandler has not been initialized yet.
    if (dof_handler.get_fe_collection().size() == 0)
      return 0;

    // Convert cell type from Triangulation to DoFHandler to be able
    // to access the information about the degrees of freedom.
    const typename DoFHandler<dim, spacedim>::cell_iterator cell(*cell_,
                                                                 &dof_handler);

    // Determine which FiniteElement object will be present on this cell after
    // refinement and will thus specify the number of degrees of freedom.
    unsigned int fe_index = numbers::invalid_unsigned_int;
    switch (status)
      {
        case Triangulation<dim, spacedim>::CELL_PERSIST:
        case Triangulation<dim, spacedim>::CELL_REFINE:
        case Triangulation<dim, spacedim>::CELL_INVALID:
          fe_index = cell->future_fe_index();
          break;

        case Triangulation<dim, spacedim>::CELL_COARSEN:
#ifdef DEBUG
          for (const auto &child : cell->child_iterators())
            Assert(child->is_active() && child->coarsen_flag_set(),
                   typename dealii::Triangulation<
                     dim>::ExcInconsistentCoarseningFlags());
#endif

          fe_index = dealii::internal::hp::DoFHandlerImplementation::
            dominated_future_fe_on_children<dim, spacedim>(cell);
          break;

        default:
          Assert(false, ExcInternalError());
          break;
      }

    // Return the cell weight determined by the function of choice.
    return weighting_function(cell, dof_handler.get_fe(fe_index));
  }



  // ---------- deprecated functions ----------

  template <int dim, int spacedim>
  CellWeights<dim, spacedim>::CellWeights(
    const dealii::DoFHandler<dim, spacedim> &dof_handler)
    : dof_handler(&dof_handler)
    , triangulation(
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &(dof_handler.get_triangulation())))
  {
    Assert(
      triangulation != nullptr,
      ExcMessage(
        "parallel::CellWeights requires a parallel::TriangulationBase object."));

    register_constant_weighting();
  }



  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::register_constant_weighting(
    const unsigned int factor)
  {
    connection.disconnect();

    connection = triangulation->signals.cell_weight.connect(
      [this,
       factor](const typename Triangulation<dim, spacedim>::cell_iterator &cell,
               const typename Triangulation<dim, spacedim>::CellStatus status)
        -> unsigned int {
        return this->weighting_callback(cell,
                                        status,
                                        std::cref(*(this->dof_handler)),
                                        std::cref(*(this->triangulation)),
                                        this->constant_weighting(factor));
      });
  }



  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::register_ndofs_weighting(
    const unsigned int factor)
  {
    connection.disconnect();

    connection = triangulation->signals.cell_weight.connect(
      [this,
       factor](const typename Triangulation<dim, spacedim>::cell_iterator &cell,
               const typename Triangulation<dim, spacedim>::CellStatus status)
        -> unsigned int {
        return this->weighting_callback(cell,
                                        status,
                                        std::cref(*(this->dof_handler)),
                                        std::cref(*(this->triangulation)),
                                        this->ndofs_weighting({factor, 1}));
      });
  }



  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::register_ndofs_squared_weighting(
    const unsigned int factor)
  {
    connection.disconnect();

    connection = triangulation->signals.cell_weight.connect(
      [this,
       factor](const typename Triangulation<dim, spacedim>::cell_iterator &cell,
               const typename Triangulation<dim, spacedim>::CellStatus status)
        -> unsigned int {
        return this->weighting_callback(cell,
                                        status,
                                        std::cref(*(this->dof_handler)),
                                        std::cref(*(this->triangulation)),
                                        this->ndofs_weighting({factor, 2}));
      });
  }



  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::register_custom_weighting(
    const std::function<
      unsigned int(const FiniteElement<dim, spacedim> &,
                   const typename DoFHandler<dim, spacedim>::cell_iterator &)>
      custom_function)
  {
    connection.disconnect();

    const std::function<
      unsigned int(const typename DoFHandler<dim, spacedim>::cell_iterator &,
                   const FiniteElement<dim, spacedim> &)>
      converted_function =
        [&custom_function](
          const typename DoFHandler<dim, spacedim>::cell_iterator &cell,
          const FiniteElement<dim, spacedim> &future_fe) -> unsigned int {
      return custom_function(future_fe, cell);
    };

    connection = triangulation->signals.cell_weight.connect(
      [this, converted_function](
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const typename Triangulation<dim, spacedim>::CellStatus     status)
        -> unsigned int {
        return this->weighting_callback(cell,
                                        status,
                                        std::cref(*(this->dof_handler)),
                                        std::cref(*(this->triangulation)),
                                        converted_function);
      });
  }
} // namespace parallel


// explicit instantiations
#include "cell_weights.inst"

DEAL_II_NAMESPACE_CLOSE
