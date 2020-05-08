// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#ifndef dealii_base_parsed_convergence_table_h
#define dealii_base_parsed_convergence_table_h

#include <deal.II/base/config.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools_integrate_difference.h>

DEAL_II_NAMESPACE_OPEN

/**
 * @brief The ParsedConvergenceTable class
 *
 * This class simplifies the construction of convergence tables, reading the
 * options for the generation of the table from a parameter file. It provides a
 * series of methods that can be used to compute the error given a reference
 * exact solution, or the difference between two numerical solutions, or any
 * other custom computation of the error, given via std::function objects.
 *
 * An example usage of this class is given by
 * @code
 * ParsedConvergenceTable table;
 *
 * ParameterHandler prm;
 * table.add_parameters(prm);
 *
 * for (unsigned int i = 0; i < n_cycles; ++i)
 *   {
 *     ... // do some computations
 *     table.error_from_exact(dof_handler, solution, exact_solution);
 *   }
 * table.output_table(std::cout);
 * @endcode
 *
 * The above code constructs a ParsedConvergenceTable that works for
 * scalar problems, and will produce an error table with `H1_norm`, `L2_norm`,
 * and `Linfty_norm` norms of the error.
 *
 * Whenever a call to the methods error_from_exact() or difference() is made,
 * the instance of this class inspects its parameters, computes all norms
 * specified by the parameter given at construction time, possibly modified
 * via a parameter file, computes all extra column entries specified using the
 * method add_extra_column(), and writes one row of the convergence table.
 *
 * Once you have finished with the computations, a call to output_table() will
 * generate a formatted convergence table on the provided stream, and to the
 * file (if any) specified in the parameter file.
 *
 * With a small modification, the same code can be used to estimate the errors
 * of mixed or multi-physics problems, e.g.:
 * @code
 * using namespace VectorTools;
 * ParsedConvergenceTable table({"u,u,p"},{{H1_norm, L2_norm}, {L2_norm}});
 *
 * ParameterHandler prm;
 * table.add_parameters(prm);
 *
 * for (unsigned int i = 0; i < n_cycles; ++i)
 *   {
 *     ... // do some computations
 *     table.error_from_exact(dof_handler, solution, exact_solution);
 *   }
 * table.output_table(std::cout);
 * @endcode
 *
 * The above code assumes that you are solving a Stokes problem with three
 * components. Two components for the vector velocity field `u`, and one
 * component for the pressure field `p`, and will produce an error
 * table with `H1` and `L2` norm of the error in the velocity field (first two,
 * components) and `L2` error in the pressure field.
 *
 * You may also call `table.output_table()` without arguments, to write the
 * table only to the file specified in the parameter file.
 *
 * By calling the method add_parameters() passing a ParameterHandler object,
 * the following options will be defined in the given ParameterHandler object
 * (in the current level of the ParameterHandler object, i.e., whatever level
 * you have entered with the ParamterHandler::enter_subsection() method),
 * and can be modified at run time through a parameter file:
 * @code
 * set Enable computation of the errors = true
 * set Error file name                  =
 * set Error precision                  = 3
 * set Exponent for p-norms             = 2
 * set Extra columns                    = dofs, cells
 * set List of error norms to compute   = Linfty_norm, L2_norm, H1_norm
 * set Rate key                         = dofs
 * set Rate mode                        = reduction_rate_log2
 * @endcode
 *
 * When using this class, please cite
 * @code{.bib}
 * @article{SartoriGiulianiBardelloni-2018-a,
 *  Author = {Sartori, Alberto and Giuliani, Nicola and
 *            Bardelloni, Mauro and Heltai, Luca},
 *  Journal = {SoftwareX},
 *  Pages = {318--327},
 *  Title = {{deal2lkit: A toolkit library for high performance
 *            programming in deal.II}},
 *  Doi = {10.1016/j.softx.2018.09.004},
 *  Volume = {7},
 *  Year = {2018}}
 * @endcode
 *
 * @author Luca Heltai, 2019.
 */
class ParsedConvergenceTable
{
public:
  /**
   * Minimal constructor for ParsedConvergenceTable objects.
   *
   * The number of components must match the number of components of the
   * finite element space that is used to compute the errors. If
   * a component name is repeated, than it is interpreted as a vector field,
   * and the errors of the repeated components are grouped together.
   *
   * The size of the vector @p list_of_error_norms must match the number of
   * unique component names, and may contain zero or more comma separated
   * identifiers for the norm to compute for each component (see the
   * documentation of VectorTools::NormType for the available options).
   *
   * For example, the following constructor
   * @code
   * using namespace VectorTools;
   * ParsedConvergenceTable table({"u", "v", "v"},
   *                              {{Linfty_norm}, {L2_norm, H1_norm}});
   * @endcode
   * would produce (if the parameter file is left untouched) a table similar to
   * @code
   * cells dofs u_Linfty_norm    v_L2_norm      v_H1_norm
   * 4     9    1.183e-01 -    5.156e-02 -    2.615e-01 -
   * 16    25   3.291e-02 2.50 1.333e-02 2.65 1.272e-01 1.41
   * 64    81   8.449e-03 2.31 3.360e-03 2.34 6.313e-02 1.19
   * 256   289  2.126e-03 2.17 8.418e-04 2.18 3.150e-02 1.09
   * 1024  1089 5.325e-04 2.09 2.106e-04 2.09 1.574e-02 1.05
   * @endcode
   *
   * See the other constructor for a documentation of all the parameters you can
   * change.
   *
   * @param component_names Specify the names of the components;
   * @param list_of_error_norms Specify what error norms to compute for each
   * unique component name.
   */
  ParsedConvergenceTable(
    const std::vector<std::string> &                    component_names = {"u"},
    const std::vector<std::set<VectorTools::NormType>> &list_of_error_norms = {
      {VectorTools::H1_norm, VectorTools::L2_norm, VectorTools::Linfty_norm}});

  /**
   * Full constructor for ParsedConvergenceTable.
   *
   * @param component_names Names of the components. Repeated consecutive
   * names are interpreted as components of a vector valued field;
   * @param list_of_error_norms Specify what error norms to compute for each
   * unique component name;
   * @param exponent The exponent to use in p-norms;
   * @param extra_columns Extra columns to add. These may be "cells" or "dofs";
   * @param rate_key Specify the extra column by which we will compute the
   * error rates. This key can either be one of "cells" or "dofs", or, if you
   * add extra columns to the table via the method add_extra_column(), it may
   * be one of the extra columns you added;
   * @param rate_mode Specify the rate mode to use when computing error rates.
   * This maybe either "reduction_rate", "reduction_rate_log2", or "none". See
   * the documentation of ConvergenceTable::RateMode for an explanation of how
   * each of this mode behaves;
   * @param error_file_name Name of error output file (with extension txt,
   * gpl, tex, or org). If different from the empty string, than
   * output_table() also writes in this file in the format deduced from its
   * extension;
   * @param precision How many digits to use when writing the error;
   * @param compute_error Control whether the filling of the table is enabled
   * or not. This flag may be used to disable at run time any error computation;
   *
   * The parameters you specify with this constructor can be written to a
   * ParameterHandler object by calling the add_parameters() method. Once you
   * call the add_parameters() method, the following options will be defined in
   * the given ParameterHandler object, and the parameters of the instance of
   * this class will follow the modification you make to the ParameterHandler
   * object at run time:
   * @code
   * # Listing of Parameters
   * # ---------------------
   * # When set to false, no computations are performed.
   * set Enable computation of the errors = true
   *
   * # Set this to a filename with extension .txt, .gpl, .org, or .tex to enable
   * # writing the convergence table to a file.
   * set Error file name                  =
   *
   * # Number of digits to use when printing the error.
   * set Error precision                  = 3
   *
   * # Extra columns to add to the table. Available options are dofs and cells.
   * set Extra columns                    = dofs, cells
   *
   * # The exponent to use when computing p-norms.
   * set Exponent for p-norms             = 2
   *
   * # Each component is separated by a semicolon and each norm by a comma. See
   * # the documentation of VectorTools::NormType for a list of implemented
   * # norms. If you want to skip a component, leave its entry empty.
   * set List of error norms to compute   = Linfty_norm, L2_norm, H1_norm
   *
   * # Key to use when computing convergence rates. If this is set to a
   * # column that is not present, or to none, then no error rates are computed.
   * set Rate key                         = dofs
   *
   * # What type of error rate to compute. Available options are
   * # reduction_rate_log2, reduction_rate, and none.
   * set Rate mode                        = reduction_rate_log2
   * @endcode
   */
  ParsedConvergenceTable(
    const std::vector<std::string> &                    component_names,
    const std::vector<std::set<VectorTools::NormType>> &list_of_error_norms,
    const double                                        exponent,
    const std::set<std::string> &                       extra_columns,
    const std::string &                                 rate_key,
    const std::string &                                 rate_mode,
    const std::string &                                 error_file_name,
    const unsigned int                                  precision,
    const bool                                          compute_error);

  /**
   * Attach all the parameters in this class to entries of the parameter
   * handler @p prm. Whenever the content of @p prm changes, the parameters
   * of this class will be updated.
   */
  void
  add_parameters(ParameterHandler &prm);

  /**
   * Add a row to the error table, containing the error between @p solution and
   * the @p exact function, in the norm(s) specified in the parameter file.
   *
   * If you specify a @p weight function during this call, then this is used
   * to compute weighted errors. The weight function can be either a scalar
   * function (which will be used for all components), or a vector function.
   * When it is a vector function, an assertion is triggered if the number of
   * components does not coincide with the number of components of the
   * underlying finite element space.
   */
  template <typename DoFHandlerType, typename VectorType>
  void
  error_from_exact(
    const DoFHandlerType &                           vspace,
    const VectorType &                               solution,
    const Function<DoFHandlerType::space_dimension> &exact,
    const Function<DoFHandlerType::space_dimension> *weight = nullptr);

  /**
   * Same as above, with a different mapping.
   */
  template <typename DoFHandlerType, typename VectorType>
  void
  error_from_exact(
    const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
      &                                              mapping,
    const DoFHandlerType &                           vspace,
    const VectorType &                               solution,
    const Function<DoFHandlerType::space_dimension> &exact,
    const Function<DoFHandlerType::space_dimension> *weight = nullptr);

  /**
   * Add an additional column (with name @p column_name) to the table, by invoking
   * the function @p custom_function, when calling error_from_exact() or
   * difference().
   *
   * You can call this method as many times as you want. If @p column_name was
   * already used in a previous call, then calling this method with the same
   * name will overwrite whatever function you had previously specified. If
   * you use a lambda function for this call, make sure that the variables used
   * internally in the lambda function remain valid until the call to
   * error_from_exact() or difference().
   *
   * Make sure you add all extra columns before the first call to
   * error_from_exact() or difference(). Adding additional columns to the
   * convergence table after you already started filling the table will trigger
   * an exception.
   *
   * This method may be used, for example, to compute the error w.r.t. to
   * time step increments in time, for example:
   * @code
   * using namespace VectorTools;
   * ParsedConvergenceTable table({"u"}, {{L2_norm}});
   *
   * double dt = .5;
   * auto dt_function = [&]() {
   *        return dt;
   * };
   *
   * table.add_extra_column("dt", dt_function, false);
   *
   * for (unsigned int i = 0; i < n_cycles; ++i)
   *   {
   *     // ... compute solution at the current dt
   *
   *     table.error_from_exact(dof_handler, solution, exact_solution);
   *     dt /= 2.0;
   *   }
   * table.output_table(std::cout);
   * @endcode
   * will produce a table similar to
   * @code
   *    dt        u_L2_norm
   * 5.000e-1    5.156e-02 -
   * 2.500e-2    1.333e-02 2.65
   * 1.250e-2    3.360e-03 2.34
   * 6.250e-3    8.418e-04 2.18
   * @endcode
   * provided that you use the following parameter file (only non default
   * entries are shown here):
   * @code
   * set Extra columns                  =
   * set List of error norms to compute = L2_norm
   * set Rate key                       = dt
   * @endcode
   *
   *
   * @param column_name Name of the column to add;
   * @param custom_function Function that will be called to fill the given
   * entry. You need to make sure that the scope of this function is valid
   * up to the call to error_from_exact() or difference();
   * @param compute_rate If set to true, then this column will be included in
   * the list of columns for which error rates are computed. You may want to set
   * this to false if you want to compute error rates with respect to this
   * column. In this case, you should also specify @p column_name as the rate
   * key in the parameter file.
   */
  void
  add_extra_column(const std::string &            column_name,
                   const std::function<double()> &custom_function,
                   const bool                     compute_rate = true);

  /**
   * Difference between two solutions in the same vector space.
   */
  template <typename DoFHandlerType, typename VectorType>
  void
  difference(const DoFHandlerType &,
             const VectorType &,
             const VectorType &,
             const Function<DoFHandlerType::space_dimension> *weight = nullptr);

  /**
   * Same as above, with a non default mapping.
   */
  template <typename DoFHandlerType, typename VectorType>
  void
  difference(const Mapping<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &mapping,
             const DoFHandlerType &,
             const VectorType &,
             const VectorType &,
             const Function<DoFHandlerType::space_dimension> *weight = nullptr);

  /**
   * Write the error table to the @p out stream (in text format), and
   * (possibly) to the file stream specified in the parameters (with the format
   * deduced from the file name extension).
   */
  void
  output_table(std::ostream &out);

  /**
   * Write the error table to the file stream specified in the parameters.
   *
   * If the "Error file name" option in the parameter file is set to the empty
   * string, no output is written.
   */
  void
  output_table();

private:
  /**
   * Add rates to the output table.
   */
  void
  prepare_table_for_output();

  /**
   * Names of the solution components.
   */
  const std::vector<std::string> component_names;

  /**
   * Same as above, but containing repeated component names only once.
   */
  const std::vector<std::string> unique_component_names;

  /**
   * Masks for each unique component name.
   */
  const std::vector<ComponentMask> unique_component_masks;

  /**
   * Additional methods to call when adding rows to the table.
   */
  std::map<std::string, std::pair<std::function<double()>, bool>>
    extra_column_functions;

  /**
   * Type of error to compute per components.
   */
  std::vector<std::set<VectorTools::NormType>> norms_per_unique_component;

  /**
   * Exponent to use in p-norm types.
   */
  double exponent;

  /**
   * The actual table
   */
  ConvergenceTable table;

  /**
   * Extra columns to add to the table.
   */
  std::set<std::string> extra_columns;

  /**
   * The name of column with respect to which we compute convergence rates.
   */
  std::string rate_key;

  /**
   * Reduction rate mode. See ConvergenceTable::RateMode for a documentation.
   */
  std::string rate_mode;

  /**
   * The precision used to output the table.
   */
  unsigned int precision;

  /**
   * Filename to use when writing to file.
   */
  std::string error_file_name;

  /**
   * Compute the error. If this is false, all methods that perform the
   * computation of the error are disabled and don't do anything.
   */
  bool compute_error;
};



#ifndef DOXYGEN
// ============================================================
// Template functions
// ============================================================
template <typename DoFHandlerType, typename VectorType>
void
ParsedConvergenceTable::difference(
  const DoFHandlerType &                           dh,
  const VectorType &                               solution1,
  const VectorType &                               solution2,
  const Function<DoFHandlerType::space_dimension> *weight)
{
  AssertThrow(solution1.size() == solution2.size(),
              ExcDimensionMismatch(solution1.size(), solution2.size()));
  VectorType solution(solution1);
  solution -= solution2;
  error_from_exact(dh,
                   solution,
                   ConstantFunction<DoFHandlerType::space_dimension>(
                     0, component_names.size()),
                   weight);
}



template <typename DoFHandlerType, typename VectorType>
void
ParsedConvergenceTable::difference(
  const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
    &                                              mapping,
  const DoFHandlerType &                           dh,
  const VectorType &                               solution1,
  const VectorType &                               solution2,
  const Function<DoFHandlerType::space_dimension> *weight)
{
  AssertThrow(solution1.size() == solution2.size(),
              ExcDimensionMismatch(solution1.size(), solution2.size()));
  VectorType solution(solution1);
  solution -= solution2;
  error_from_exact(mapping,
                   dh,
                   solution,
                   ConstantFunction<DoFHandlerType::space_dimension>(
                     0, component_names.size()),
                   weight);
}



template <typename DoFHandlerType, typename VectorType>
void
ParsedConvergenceTable::error_from_exact(
  const DoFHandlerType &                           dh,
  const VectorType &                               solution,
  const Function<DoFHandlerType::space_dimension> &exact,
  const Function<DoFHandlerType::space_dimension> *weight)
{
  error_from_exact(StaticMappingQ1<DoFHandlerType::dimension,
                                   DoFHandlerType::space_dimension>::mapping,
                   dh,
                   solution,
                   exact,
                   weight);
}



template <typename DoFHandlerType, typename VectorType>
void
ParsedConvergenceTable::error_from_exact(
  const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
    &                                              mapping,
  const DoFHandlerType &                           dh,
  const VectorType &                               solution,
  const Function<DoFHandlerType::space_dimension> &exact,
  const Function<DoFHandlerType::space_dimension> *weight)
{
  const int  dim          = DoFHandlerType::dimension;
  const int  spacedim     = DoFHandlerType::space_dimension;
  const auto n_components = component_names.size();

  if (compute_error)
    {
      AssertDimension(exact.n_components, n_components);
      AssertDimension(dh.get_fe().n_components(), n_components);

      const types::global_cell_index n_active_cells =
        dh.get_triangulation().n_global_active_cells();
      const unsigned int n_dofs = dh.n_dofs();

      for (const auto &col : extra_columns)
        if (col == "cells")
          {
            table.add_value("cells", n_active_cells);
            table.set_tex_caption("cells", "\\# cells");
            table.set_tex_format("cells", "r");
          }
        else if (col == "dofs")
          {
            table.add_value("dofs", n_dofs);
            table.set_tex_caption("dofs", "\\# dofs");
            table.set_tex_format("dofs", "r");
          }

      // A vector of zero std::functions with n_components components
      const std::vector<std::function<double(const Point<spacedim> &)>>
        zero_components(n_components,
                        [](const Point<spacedim> &) { return 0.0; });

      // The default weight function, with n_components components
      std::vector<std::function<double(const Point<spacedim> &)>>
        weight_components(n_components,
                          [](const Point<spacedim> &) { return 1.0; });

      if (weight != nullptr)
        {
          if (weight->n_components == 1)
            {
              for (auto &f : weight_components)
                f = [&](const Point<spacedim> &p) { return weight->value(p); };
            }
          else
            {
              AssertDimension(weight->n_components, n_components);
              for (unsigned int i = 0; i < n_components; ++i)
                weight_components[i] = [&](const Point<spacedim> &p) {
                  return weight->value(p, i);
                };
            }
        }

      for (unsigned int i = 0; i < norms_per_unique_component.size(); ++i)
        {
          std::map<VectorTools::NormType, double> errors;

          const auto &norms = norms_per_unique_component[i];
          const auto &mask  = unique_component_masks[i];

          // Simple case first
          if (norms.empty())
            continue;

          auto components_expr = zero_components;
          for (unsigned int i = 0; i < n_components; ++i)
            if (mask[i] == true)
              components_expr[i] = weight_components[i];

          FunctionFromFunctionObjects<spacedim> select_component(
            components_expr);

          Vector<float> difference_per_cell(
            dh.get_triangulation().n_global_active_cells());

          QGauss<dim> q_gauss((dh.get_fe().degree + 1) * 2);

          for (const auto &norm : norms)
            {
              difference_per_cell = 0;
              VectorTools::integrate_difference(mapping,
                                                dh,
                                                solution,
                                                exact,
                                                difference_per_cell,
                                                q_gauss,
                                                norm,
                                                &select_component,
                                                exponent);

              errors[norm] = VectorTools::compute_global_error(
                dh.get_triangulation(), difference_per_cell, norm, exponent);

              std::string name = unique_component_names[i] + "_" +
                                 Patterns::Tools::to_string(norm);
              std::string latex_name = "$\\| " + unique_component_names[i] +
                                       " - " + unique_component_names[i] +
                                       "_h \\|_{" +
                                       Patterns::Tools::to_string(norm) + "}$";

              table.add_value(name, errors[norm]);
              table.set_precision(name, precision);
              table.set_scientific(name, true);
              table.set_tex_caption(name, latex_name);
            }
        }

      for (const auto &extra_col : extra_column_functions)
        {
          const double custom_error = extra_col.second.first();

          std::string name = extra_col.first;
          table.add_value(name, custom_error);
          table.set_precision(name, precision);
          table.set_scientific(name, true);
        }
    }
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
