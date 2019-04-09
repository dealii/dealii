//-----------------------------------------------------------
//
//    Copyright (C) 2019 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal2lkit library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal2lkit distribution.
//
//-----------------------------------------------------------

#ifndef dealii_base_parsed_convergence_table_h
#define dealii_base_parsed_convergence_table_h

#include <deal.II/base/config.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_convergence_table_flags.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

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
 * for(unsigned int i=0; i<n_cycles; ++i) {
 *   ... // do some computations
 *   table.error_from_exact(dof_handler, solution, exact_solution);
 * }
 * table.output_table(std::cout);
 * @endcode
 *
 * The above code constructs a ParsedConvergenceTable that works for
 * scalar problems, and will produce an error table with `H1`, `L2`, and
 * `Linfty` norms of the error.
 *
 * Whenever a call to the methods error_from_exact() or difference() is made,
 * the instance of this class inspects its parameters, and computes all norms
 * specified by the parameter given at construction time, possibly modified
 * via a ParameterFile.
 *
 * With a small modification, the same code can be used to estimate the errors
 * of mixed or multi-physics problems, e.g.:
 * @code
 * using namespace ParsedConvergenceTableFlags;
 * ParsedConvergenceTable table({"u,u,p"},{{H1, L2}, {L2}});
 *
 * ParameterHandler prm;
 * table.add_parameters(prm);
 *
 * for(unsigned int i=0; i<n_cycles; ++i) {
 *   ... // do some computations
 *   table.error_from_exact(dof_handler, solution, exact_solution);
 * }
 * table.output_table(std::cout);
 * @endcode
 *
 * The above code assumes that you are solving a Stokes problem with three
 * components. Two components for the vector velocity field `u`, and one
 * component for the pressure field `p`, and will produce an error
 * table with `H1` and `L2` norm of the error in the velocity field (first two,
 * components) and `L2` error in the pressure field.
 *
 * By calling the method add_parameters() passing a ParameterHandler object,
 * the following options will be defined in the given ParameterHandler object,
 * and can be modified at run time through a parameter file:
 * @code
 * set Enable computation of the errors = true
 * set Enable output to streams         = true
 * set Error file name                  =
 * set Error precision                  = 3
 * set Extra columns                    = dofs, cells
 * set List of error norms to compute   = Linfty, L2, H1
 * set Rate key                         = dofs
 * @endcode
 *
 * When using this class, please cite
 * @code{.bib}
 * @article{SartoriGiulianiBardelloni-2018-a,
 * 	Author = {Sartori, Alberto and Giuliani, Nicola and
 *            Bardelloni, Mauro and Heltai, Luca},
 * 	Journal = {SoftwareX},
 * 	Pages = {318--327},
 * 	Title = {{deal2lkit: A toolkit library for high performance
 *            programming in deal.II}},
 *  Doi = {10.1016/j.softx.2018.09.004},
 * 	Volume = {7},
 * 	Year = {2018}}
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
   * strings that identify the norm types to compute for each component.
   *
   * For example, the following constructor
   * @code
   * using namespace ParsedConvergenceTableFlags;
   * ParsedConvergenceTable table({"u", "v", "v"}, {{Linfty}, {L2, H1}});
   * @endcode
   *
   * Will produce (if the parameter file is left untouched) the following table:
   * @code
   * cells dofs    u_Linfty         v_L2           v_H1
   * 4     9    1.183e-01 -    5.156e-02 -    2.615e-01 -
   * 16    25   3.291e-02 2.50 1.333e-02 2.65 1.272e-01 1.41
   * 64    81   8.449e-03 2.31 3.360e-03 2.34 6.313e-02 1.19
   * 256   289  2.126e-03 2.17 8.418e-04 2.18 3.150e-02 1.09
   * 1024  1089 5.325e-04 2.09 2.106e-04 2.09 1.574e-02 1.05
   * @endcode
   *
   * Se the documentation of NormType for a list of
   * available norms you can compute using this class, and see the other
   * constructor for a documentation of all the parameters you can change.
   *
   * @param component_names Specify the names of the components
   * @param list_of_error_norms Specify what error norms to compute for each
   * unique component name
   */
  ParsedConvergenceTable(
    const std::vector<std::string> &component_names = {"u"},
    const std::vector<std::set<ParsedConvergenceTableFlags::NormType>>
      &list_of_error_norms = {{ParsedConvergenceTableFlags::H1,
                               ParsedConvergenceTableFlags::L2,
                               ParsedConvergenceTableFlags::Linfty}});

  /**
   * Full constructor for ParsedConvergenceTable.
   *
   * @param component_names Names of the components. Repeated consecutive
   * names are interpreted as components of a vector valued field;
   * @param list_of_error_norms Specify what error norms to compute for each
   * unique component name;
   * @param extra_columns Extra columns to add;
   * @param rate_key Specify the extra column by which we will compute the
   * error rates;
   * @param error_file_name Name of error output file (with extension txt,
   * gpl, tex, or org). If different from the empty string, than
   * output_table() also writes in this file in the format deduced from its
   * extension;
   * @param precision How many digits to use when writing the error;
   * @param compute_error If set to `false`, then no error computation is
   * performed, and no output is produced;
   * @param output_error If set to `false`, then the call to output_table()
   * will produce only an output file (if @p error_file_name is not empty),
   * but no output will be written on the given stream.
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
   * # When set to false, printing of the convergence table to the stream
   * # specified as input to output_table() is disabled.
   * set Enable output to streams         = true
   *
   * # Set this to a filename with extension .txt, .gpl, .org, or .tex to enable
   * # writing the convergence table to a file.
   * set Error file name                  =
   *
   * # Number of digits to use when printing the error.
   * set Error precision                  = 3
   *
   * # Extra columns to add to the table. Availabl options are dofs, cells, and
   * # dt.
   * set Extra columns                    = dofs, cells
   *
   * # Each component is separated by a semicolon, and each norm by a comma.
   * # Implemented norms are Linfty, L2, W1infty, H1, and custom. If you want to
   * # skip a component, use none.
   * set List of error norms to compute   = Linfty, L2, H1
   *
   * # Key to use when computing convergence rates. If this is set to a
   * # column that is not present, or to none, then no error rates are computed.
   * set Rate key                         = dofs
   * @endcode
   */
  ParsedConvergenceTable(
    const std::vector<std::string> &component_names,
    const std::vector<std::set<ParsedConvergenceTableFlags::NormType>>
      &list_of_error_norms,
    const std::set<ParsedConvergenceTableFlags::ExtraColumns> &extra_columns,
    const ParsedConvergenceTableFlags::ExtraColumns &          rate_key,
    const std::string &                                        error_file_name,
    const unsigned int &                                       precision,
    const bool &                                               compute_error,
    const bool &                                               output_error);

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
   * If you specify a @p weight expression during this call, then this is used
   * to compute weighted errors. If the weight is different from "1.0", then
   * you have to make sure that muparser is installed and enabled, otherwise an
   * exception will be thrown.
   */
  template <typename DoFHandlerType, typename VectorType>
  void
  error_from_exact(const DoFHandlerType &                           vspace,
                   const VectorType &                               solution,
                   const Function<DoFHandlerType::space_dimension> &exact,
                   double                                           dt = 0.,
                   const std::string &weight                           = "1.0");

  /**
   * Same as above, with a different mapping.
   */
  template <typename DoFHandlerType, typename VectorType>
  void
  error_from_exact(const Mapping<DoFHandlerType::dimension,
                                 DoFHandlerType::space_dimension> & mapping,
                   const DoFHandlerType &                           vspace,
                   const VectorType &                               solution,
                   const Function<DoFHandlerType::space_dimension> &exact,
                   double                                           dt = 0.,
                   const std::string &weight                           = "1.0");


  /**
   * Call the given custom function to compute a custom error for any component
   * for which you specify the ConvergenceTableFlags::custom norm. This function
   * may be called several times to add different columns with different names,
   * provided you use an appropriate function every time, and
   * that you update @p column_name everytime.
   *
   * Should you wish to do so, then make sure the parameter
   * @p add_table_extras is set to false for each call, except one, so
   * that you only add extra columns informations once.
   *
   * An example usage of this class with a custom error function is the
   * following:
   * @code
   * using namespace ParsedConvergenceTableFlags;
   * ParsedConvergenceTable table({"u"}, {{"custom"}});
   *
   * for(unsigned int i=0; i<n_cycles; ++i) {
   *   // ... refine and distribute dofs
   *   auto cycle_function = [&](const unsigned int component) {
   *      return i;
   *   };
   *   table.custom_error(cycle_function, dof_handler, "cycle", true);
   * }
   * table.output_table(std::cout);
   * @endcode
   *
   * Will produce (if the parameter file is left untouched) the following table:
   * @code
   * cells dofs  cycle
   * 4     9     0
   * 16    25    1
   * 64    81    2
   * 256   289   3
   * 1024  1089  4
   * @endcode
   *
   * Should you wish to add another column, then the above code should be
   * replaced by:
   * @code
   * using namespace ParsedConvergenceTableFlags;
   * ParsedConvergenceTable table({"u"}, {"custom"});
   *
   * for(unsigned int i=0; i<n_cycles; ++i) {
   *   // ... refine and distribute dofs
   *   auto cycle_function = [&](const unsigned int) {
   *      return i;
   *   };
   *   auto cycle_square_function = [&](const unsigned int) {
   *      return i*i;
   *   };
   *   table.custom_error(cycle_function, dof_handler, "cycle", false);
   *   table.custom_error(cycle_square_function, dof_handler,
   *                      "cycle_square", true);
   * }
   * table.output_table(std::cout);
   * @endcode
   *
   * This will produce the table:
   * @code
   * cells dofs  cycle cycle_square
   * 4     9     0     0
   * 16    25    1     1
   * 64    81    2     4
   * 256   289   3     9
   * 1024  1089  4     16
   * @endcode
   */
  template <typename DoFHandlerType>
  void
  custom_error(const std::function<double(const unsigned int component)>
                 &                   custom_error_function,
               const DoFHandlerType &dh,
               const std::string &   column_name      = "custom",
               const bool            add_table_extras = false,
               const double          dt               = 0.);

  /**
   * Difference between two solutions in the same vector space.
   */
  template <typename DoFHandlerType, typename VectorType>
  void
  difference(const DoFHandlerType &,
             const VectorType &,
             const VectorType &,
             double dt = 0.);

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
             double dt = 0.);

  /**
   * Write the error table to the @p out stream, and (possibly) to the file
   * stream specified in the parameters.
   */
  void
  output_table(std::ostream &out = std::cout);

private:
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
   * Type of error to compute per components.
   */
  std::vector<std::set<ParsedConvergenceTableFlags::NormType>>
    norms_per_unique_component;

  /**
   * The actual table
   */
  ConvergenceTable table;
  /**
   * Extra columns to add to the table.
   */
  std::set<ParsedConvergenceTableFlags::ExtraColumns> extra_columns;

  /**
   * Wether or not to calculate the rates according to the given keys.
   */
  ParsedConvergenceTableFlags::ExtraColumns rate_key;

  /**
   * The precision used to output the table.
   */
  unsigned int precision;

  /**
   * Filename to use when writing to file.
   */
  std::string error_file_name;

  /**
   * Compute the error. If this is false, all functions regarding
   * errors are disabled and don't do anything.
   */
  bool compute_error;

  /**
   * Output the error file also on screen.
   */
  bool output_error;
};



#ifndef DOXYGEN
// ============================================================
// Template instantiations
// ============================================================
template <typename DoFHandlerType, typename VectorType>
void
ParsedConvergenceTable::difference(const DoFHandlerType &dh,
                                   const VectorType &    solution1,
                                   const VectorType &    solution2,
                                   double                dt)
{
  AssertThrow(solution1.size() == solution2.size(),
              ExcDimensionMismatch(solution1.size(), solution2.size()));
  VectorType solution(solution1);
  solution -= solution2;
  error_from_exact(dh,
                   solution,
                   ConstantFunction<DoFHandlerType::space_dimension>(
                     0, component_names.size()),
                   dt);
}



template <typename DoFHandlerType, typename VectorType>
void
ParsedConvergenceTable::difference(
  const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
    &                   mapping,
  const DoFHandlerType &dh,
  const VectorType &    solution1,
  const VectorType &    solution2,
  double                dt)
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
                   dt);
}



template <typename DoFHandlerType, typename VectorType>
void
ParsedConvergenceTable::error_from_exact(
  const DoFHandlerType &                           dh,
  const VectorType &                               solution,
  const Function<DoFHandlerType::space_dimension> &exact,
  double                                           dt,
  const std::string &                              weight)
{
  error_from_exact(StaticMappingQ1<DoFHandlerType::dimension,
                                   DoFHandlerType::space_dimension>::mapping,
                   dh,
                   solution,
                   exact,
                   dt,
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
  double                                           dt,
  const std::string &                              weight)
{
  const int dim      = DoFHandlerType::dimension;
  const int spacedim = DoFHandlerType::space_dimension;
  if (compute_error)
    {
      AssertDimension(exact.n_components, component_names.size());

      const unsigned int n_active_cells =
        dh.get_triangulation().n_global_active_cells();
      const unsigned int n_dofs = dh.n_dofs();

      for (const auto &col : extra_columns)
        if (col == ParsedConvergenceTableFlags::cells)
          {
            table.add_value("cells", n_active_cells);
            table.set_tex_caption("cells", "\\# cells");
            table.set_tex_format("cells", "r");
          }
        else if (col == ParsedConvergenceTableFlags::dofs)
          {
            table.add_value("dofs", n_dofs);
            table.set_tex_caption("dofs", "\\# dofs");
            table.set_tex_format("dofs", "r");
          }
        else if (col == ParsedConvergenceTableFlags::dt)
          {
            table.add_value("dt", dt);
            table.set_tex_caption("dt", "\\Delta t");
            table.set_tex_format("dt", "r");
          }


      std::map<ParsedConvergenceTableFlags::NormType, unsigned int>
        norm_to_index;
      norm_to_index[ParsedConvergenceTableFlags::L2]      = 0;
      norm_to_index[ParsedConvergenceTableFlags::H1]      = 1;
      norm_to_index[ParsedConvergenceTableFlags::Linfty]  = 2;
      norm_to_index[ParsedConvergenceTableFlags::W1infty] = 3;
      norm_to_index[ParsedConvergenceTableFlags::custom]  = 4;

      const auto n_components = component_names.size();
      for (unsigned int i = 0; i < unique_component_names.size(); ++i)
        {
          std::vector<double> error(5, 0.0);
          const auto &        norms = norms_per_unique_component[i];
          const auto &        mask  = unique_component_masks[i];

#  ifndef DEAL_II_WITH_MUPARSER
          Assert(weight == "1.0",
                 ExcMessage("If you want to use a weight which is "
                            "different from 1.0, you have to "
                            "install deal.II with muparser enabled."));
          ComponentSelectFunction<spacedim> select_component(
            {mask.first_selected_component(),
             mask.first_selected_component() + mask.n_selected_components()},
            n_components);
#  else
          std::string comp_sel;
          {
            // Select all required components
            std::vector<std::string> expr(n_components, "0");

            for (unsigned int i = 0; i < n_components; ++i)
              if (mask[i] == true)
                expr[i] = weight;

            comp_sel = expr[0];
            for (unsigned int i = 1; i < expr.size(); ++i)
              comp_sel += "; " + expr[i];
          }
          FunctionParser<spacedim> select_component(comp_sel);
#  endif

          Vector<float> difference_per_cell(
            dh.get_triangulation().n_global_active_cells());

          QGauss<dim> q_gauss((dh.get_fe().degree + 1) * 2);

          for (const auto &norm : norms)
            {
              if (norm == ParsedConvergenceTableFlags::L2)
                {
                  VectorTools::integrate_difference(mapping,
                                                    dh,
                                                    solution,
                                                    exact,
                                                    difference_per_cell,
                                                    q_gauss,
                                                    VectorTools::L2_norm,
                                                    &select_component);
                }

              error[norm_to_index[ParsedConvergenceTableFlags::L2]] =
                VectorTools::compute_global_error(dh.get_triangulation(),
                                                  difference_per_cell,
                                                  VectorTools::L2_norm);
              difference_per_cell = 0;

              if (norm == ParsedConvergenceTableFlags::H1)
                {
                  VectorTools::integrate_difference(mapping,
                                                    dh,
                                                    solution,
                                                    exact,
                                                    difference_per_cell,
                                                    q_gauss,
                                                    VectorTools::H1_norm,
                                                    &select_component);
                }
              error[norm_to_index[ParsedConvergenceTableFlags::H1]] =
                VectorTools::compute_global_error(dh.get_triangulation(),
                                                  difference_per_cell,
                                                  VectorTools::H1_norm);
              difference_per_cell = 0;

              if (norm == ParsedConvergenceTableFlags::W1infty)
                {
                  VectorTools::integrate_difference(mapping,
                                                    dh,
                                                    solution,
                                                    exact,
                                                    difference_per_cell,
                                                    q_gauss,
                                                    VectorTools::W1infty_norm,
                                                    &select_component);
                }

              error[norm_to_index[ParsedConvergenceTableFlags::W1infty]] =
                VectorTools::compute_global_error(
                  dh.get_triangulation(),
                  difference_per_cell,
                  VectorTools::W1infty_seminorm);
              difference_per_cell = 0;

              if (norm == ParsedConvergenceTableFlags::Linfty)
                {
                  VectorTools::integrate_difference(mapping,
                                                    dh,
                                                    solution,
                                                    exact,
                                                    difference_per_cell,
                                                    q_gauss,
                                                    VectorTools::Linfty_norm,
                                                    &select_component);
                }

              error[norm_to_index[ParsedConvergenceTableFlags::Linfty]] =
                VectorTools::compute_global_error(dh.get_triangulation(),
                                                  difference_per_cell,
                                                  VectorTools::Linfty_norm);
              difference_per_cell = 0;

              // Now add what we just computed.
              if (norm != ParsedConvergenceTableFlags::none)
                {
                  std::string name = unique_component_names[i] + "_" +
                                     Patterns::Tools::to_string(norm);
                  std::string latex_name =
                    "$\\| " + unique_component_names[i] + " - " +
                    unique_component_names[i] + "_h \\|_{" +
                    Patterns::Tools::to_string(norm) + "}$";

                  table.add_value(name, error[norm_to_index[norm]]);
                  table.set_precision(name, precision);
                  table.set_scientific(name, true);
                  table.set_tex_caption(name, latex_name);
                }
            }
        }
    }
}


template <typename DoFHandlerType>
void
ParsedConvergenceTable::custom_error(
  const std::function<double(const unsigned int component)>
    &                   custom_error_function,
  const DoFHandlerType &dh,
  const std::string &   column_name,
  const bool            add_table_extras,
  const double          dt)
{
  if (compute_error)
    {
      const unsigned int  n_components = norms_per_unique_component.size();
      std::vector<double> c_error(norms_per_unique_component.size());
      const unsigned int  n_active_cells =
        dh.get_triangulation().n_global_active_cells();
      const unsigned int n_dofs = dh.n_dofs();
      if (add_table_extras)
        {
          for (const auto &col : extra_columns)
            {
              if (col == ParsedConvergenceTableFlags::cells)
                {
                  table.add_value("cells", n_active_cells);
                  table.set_tex_caption("cells", "\\# cells");
                  table.set_tex_format("cells", "r");
                }
              if (col == ParsedConvergenceTableFlags::dofs)
                {
                  table.add_value("dofs", n_dofs);
                  table.set_tex_caption("dofs", "\\# dofs");
                  table.set_tex_format("dofs", "r");
                }
              if (col == ParsedConvergenceTableFlags::dt)
                {
                  table.add_value("dt", dt);
                  table.set_tex_caption("dt", "\\Delta t");
                  table.set_tex_format("dt", "r");
                }
            }
        }

      for (unsigned int j = 0; j < n_components; ++j)
        {
          const auto &norms = norms_per_unique_component[j];

          for (const auto &norm : norms)
            if (norm == ParsedConvergenceTableFlags::custom)
              {
                const double custom_error = custom_error_function(j);

                std::string name =
                  unique_component_names[j] + "_" + column_name;
                std::string latex_name = "$\\| " + unique_component_names[j] +
                                         " - " + unique_component_names[j] +
                                         "_h \\|_{\\text{" + column_name +
                                         "} $";

                table.add_value(name, custom_error);
                table.set_precision(name, precision);
                table.set_scientific(name, true);
                table.set_tex_caption(name, latex_name);
              }
        }
    }
}
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
