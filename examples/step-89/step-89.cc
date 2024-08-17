/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2023 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 *
 * Authors: Johannes Heinz, TU Wien, 2023
 *          Maximilian Bergbauer, TUM, 2023
 *          Marco Feder, SISSA, 2023
 *          Peter Munch, University of Augsburg/Uppsala University, 2023
 */

// @sect3{Include files}
//
// The program starts with including all the relevant header files.
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/non_matching/mapping_info.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
// The following header file provides the class FERemoteEvaluation, which allows
// to access values and/or gradients at remote triangulations similar to
// FEEvaluation.
#include <deal.II/matrix_free/fe_remote_evaluation.h>

// We pack everything that is specific for this program into a namespace
// of its own.
namespace Step89
{
  using namespace dealii;

  // @sect3{Initial conditions for vibrating membrane}
  //
  // In the following, let us first define a function that provides
  // the initial condition for the vibrating membrane test case.  It
  // implements both the initial pressure (component 0) and velocity
  // (components 1 to 1 + dim).
  //
  // There is also a function that computes the duration of one
  // oscillation.
  template <int dim>
  class InitialConditionVibratingMembrane : public Function<dim>
  {
  public:
    InitialConditionVibratingMembrane(const double modes);

    double value(const Point<dim> &p, const unsigned int comp) const final;

    double get_period_duration(const double speed_of_sound) const;

  private:
    const double M;
  };

  template <int dim>
  InitialConditionVibratingMembrane<dim>::InitialConditionVibratingMembrane(
    const double modes)
    : Function<dim>(dim + 1, 0.0)
    , M(modes)
  {
    static_assert(dim == 2, "Only implemented for dim==2");
  }

  template <int dim>
  double
  InitialConditionVibratingMembrane<dim>::value(const Point<dim>  &p,
                                                const unsigned int comp) const
  {
    if (comp == 0)
      return std::sin(M * numbers::PI * p[0]) *
             std::sin(M * numbers::PI * p[1]);

    return 0.0;
  }

  template <int dim>
  double InitialConditionVibratingMembrane<dim>::get_period_duration(
    const double speed_of_sound) const
  {
    return 2.0 / (M * std::sqrt(dim) * speed_of_sound);
  }

  // @sect3{Gauss pulse}
  //
  // Next up is a function that provides the values of a pressure
  // Gauss pulse.  As with the previous function, it implements both
  // the initial pressure (component 0) and velocity (components 1 to
  // 1 + dim).
  template <int dim>
  class GaussPulse : public Function<dim>
  {
  public:
    GaussPulse(const double shift_x, const double shift_y);

    double value(const Point<dim> &p, const unsigned int comp) const final;

  private:
    const double shift_x;
    const double shift_y;
  };

  template <int dim>
  GaussPulse<dim>::GaussPulse(const double shift_x, const double shift_y)
    : Function<dim>(dim + 1, 0.0)
    , shift_x(shift_x)
    , shift_y(shift_y)
  {
    static_assert(dim == 2, "Only implemented for dim==2");
  }

  template <int dim>
  double GaussPulse<dim>::value(const Point<dim>  &p,
                                const unsigned int comp) const
  {
    if (comp == 0)
      return std::exp(-1000.0 * ((std::pow(p[0] - shift_x, 2)) +
                                 (std::pow(p[1] - shift_y, 2))));

    return 0.0;
  }

  // @sect3{Helper functions}
  //
  // The following namespace contains free helper functions that are used in the
  // tutorial.
  namespace HelperFunctions
  {
    // Helper function to check if a boundary ID is related to a non-matching
    // face. A @c std::set that contains all non-matching boundary IDs is
    // handed over in addition to the face ID under question.
    bool is_non_matching_face(
      const std::set<types::boundary_id> &non_matching_face_ids,
      const types::boundary_id            face_id)
    {
      return non_matching_face_ids.find(face_id) != non_matching_face_ids.end();
    }

    // Helper function to set the initial conditions for the vibrating membrane
    // test case.
    template <int dim, typename Number, typename VectorType>
    void set_initial_condition(MatrixFree<dim, Number> matrix_free,
                               const Function<dim>    &initial_solution,
                               VectorType             &dst)
    {
      VectorTools::interpolate(*matrix_free.get_mapping_info().mapping,
                               matrix_free.get_dof_handler(),
                               initial_solution,
                               dst);
    }

    // Helper function to compute the time step size according to the CFL
    // condition.
    double
    compute_dt_cfl(const double hmin, const unsigned int degree, const double c)
    {
      return hmin / (std::pow(degree, 1.5) * c);
    }

    // Helper function that writes vtu output.
    template <typename VectorType, int dim>
    void write_vtu(const VectorType      &solution,
                   const DoFHandler<dim> &dof_handler,
                   const Mapping<dim>    &mapping,
                   const unsigned int     degree,
                   const std::string     &name_prefix)
    {
      DataOut<dim>          data_out;
      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;
      data_out.set_flags(flags);

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation(
          dim + 1, DataComponentInterpretation::component_is_part_of_vector);
      std::vector<std::string> names(dim + 1, "U");

      interpretation[0] = DataComponentInterpretation::component_is_scalar;
      names[0]          = "P";

      data_out.add_data_vector(dof_handler, solution, names, interpretation);

      data_out.build_patches(mapping, degree, DataOut<dim>::curved_inner_cells);
      data_out.write_vtu_in_parallel(name_prefix + ".vtu",
                                     dof_handler.get_mpi_communicator());
    }
  } // namespace HelperFunctions

  // @sect3{Material parameter description}
  //
  // The following class stores the information if the fluid is
  // homogeneous as well as the material properties at every cell.
  // This class helps access the correct values without accessing a
  // large vector of materials in the homogeneous case.
  //
  // The class is provided a map from material ids to material
  // properties (given as a pair of values for the speed of sound and
  // the density).  If the map has only one entry, the material is
  // homogeneous -- using the same values everywhere -- and we can
  // remember that fact in the `homogeneous` member variable and use it
  // to optimize some code paths below. If the material is not
  // homogeneous, we will fill a vector with the correct materials for
  // each batch of cells; this information can then be access via
  // FEEvaluationData::read_cell_data().
  //
  // As is usual when working with the MatrixFree framework, we will
  // not only access material parameters at a single site, but for
  // whole "batches" of cells. As a consequence, the functions below
  // returned VectorizedArray objects for a batch at a time.
  template <typename Number>
  class CellwiseMaterialData
  {
  public:
    template <int dim>
    CellwiseMaterialData(
      const MatrixFree<dim, Number, VectorizedArray<Number>> &matrix_free,
      const std::map<types::material_id, std::pair<double, double>>
        &material_id_map)
      : homogeneous(material_id_map.size() == 1)
    {
      Assert(material_id_map.size() > 0,
             ExcMessage("No materials given to CellwiseMaterialData."));

      if (homogeneous)
        {
          speed_of_sound_homogeneous = material_id_map.begin()->second.first;
          density_homogeneous        = material_id_map.begin()->second.second;
        }
      else
        {
          const auto n_cell_batches =
            matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();

          speed_of_sound.resize(n_cell_batches);
          density.resize(n_cell_batches);

          for (unsigned int cell = 0; cell < n_cell_batches; ++cell)
            {
              speed_of_sound[cell] = 1.;
              density[cell]        = 1.;
              for (unsigned int v = 0;
                   v < matrix_free.n_active_entries_per_cell_batch(cell);
                   ++v)
                {
                  const auto material_id =
                    matrix_free.get_cell_iterator(cell, v)->material_id();

                  speed_of_sound[cell][v] =
                    material_id_map.at(material_id).first;
                  density[cell][v] = material_id_map.at(material_id).second;
                }
            }
        }
    }

    bool is_homogeneous() const
    {
      return homogeneous;
    }

    const AlignedVector<VectorizedArray<Number>> &get_speed_of_sound() const
    {
      Assert(!homogeneous, ExcMessage("Use get_homogeneous_speed_of_sound()."));
      return speed_of_sound;
    }

    const AlignedVector<VectorizedArray<Number>> &get_density() const
    {
      Assert(!homogeneous, ExcMessage("Use get_homogeneous_density()."));
      return density;
    }

    VectorizedArray<Number> get_homogeneous_speed_of_sound() const
    {
      Assert(homogeneous, ExcMessage("Use get_speed_of_sound()."));
      return speed_of_sound_homogeneous;
    }

    VectorizedArray<Number> get_homogeneous_density() const
    {
      Assert(homogeneous, ExcMessage("Use get_density()."));
      return density_homogeneous;
    }

  private:
    const bool homogeneous;

    /* Materials in the inhomogeneous case. */
    AlignedVector<VectorizedArray<Number>> speed_of_sound;
    AlignedVector<VectorizedArray<Number>> density;

    /* Materials in the homogeneous case. */
    VectorizedArray<Number> speed_of_sound_homogeneous;
    VectorizedArray<Number> density_homogeneous;
  };

  // To be able to access the material data in every cell in a thread
  // safe way, the following class @c MaterialEvaluation is
  // used. Similar to @c FEEvaluation, functions below will create
  // their own instances of this class; thus, there can be no race
  // conditions even if these functions run multiple times in
  // parallel. For inhomogeneous materials, a @c reinit_cell() or @c
  // reinit_face() function is used to set the correct material at the
  // current cell batch. In the homogeneous case the @c _reinit()
  // functions don't have to reset the materials.
  template <int dim, typename Number>
  class MaterialEvaluation
  {
  public:
    MaterialEvaluation(
      const MatrixFree<dim, Number, VectorizedArray<Number>> &matrix_free,
      const CellwiseMaterialData<Number>                     &material_data)
      : phi(matrix_free)
      , phi_face(matrix_free, true)
      , material_data(material_data)
    {
      if (material_data.is_homogeneous())
        {
          speed_of_sound = material_data.get_homogeneous_speed_of_sound();
          density        = material_data.get_homogeneous_density();
        }
    }

    bool is_homogeneous() const
    {
      return material_data.is_homogeneous();
    }

    // The following two functions update the data for the current
    // cell or face, given a cell batch index. If the material is
    // homogeneous, there is nothing to do.  Otherwise, we reinit the
    // FEEvaluation object and store the data for the current object.
    void reinit_cell(const unsigned int cell)
    {
      if (!material_data.is_homogeneous())
        {
          phi.reinit(cell);
          speed_of_sound =
            phi.read_cell_data(material_data.get_speed_of_sound());
          density = phi.read_cell_data(material_data.get_density());
        }
    }

    void reinit_face(const unsigned int face)
    {
      if (!material_data.is_homogeneous())
        {
          phi_face.reinit(face);
          speed_of_sound =
            phi_face.read_cell_data(material_data.get_speed_of_sound());
          density = phi_face.read_cell_data(material_data.get_density());
        }
    }

    // The following functions then return the speed of sound and
    // density for the current cell batch.
    VectorizedArray<Number> get_speed_of_sound() const
    {
      return speed_of_sound;
    }

    VectorizedArray<Number> get_density() const
    {
      return density;
    }

  private:
    /* Members needed for the inhomogeneous case. */
    FEEvaluation<dim, -1, 0, 1, Number>     phi;
    FEFaceEvaluation<dim, -1, 0, 1, Number> phi_face;

    /* Material defined at every cell. */
    const CellwiseMaterialData<Number> &material_data;

    /* Materials at current cell. */
    VectorizedArray<Number> speed_of_sound;
    VectorizedArray<Number> density;
  };


  // @sect3{Boundary conditions}
  //
  // To be able to use the same kernel, for all face integrals we define
  // a class that returns the needed values at boundaries. In this tutorial
  // homogeneous pressure Dirichlet boundary conditions are applied via
  // the mirror principle, i.e. $p_h^+=-p_h^- + 2g$ with $g=0$.
  template <int dim, typename Number>
  class BCEvaluationP
  {
  public:
    BCEvaluationP(const FEFaceEvaluation<dim, -1, 0, 1, Number> &pressure_m)
      : pressure_m(pressure_m)
    {}

    typename FEFaceEvaluation<dim, -1, 0, 1, Number>::value_type
    get_value(const unsigned int q) const
    {
      return -pressure_m.get_value(q);
    }

  private:
    const FEFaceEvaluation<dim, -1, 0, 1, Number> &pressure_m;
  };

  // We don't have to apply boundary conditions for the velocity, i.e.
  // $\mathbf{u}_h^+=\mathbf{u}_h^-$.
  template <int dim, typename Number>
  class BCEvaluationU
  {
  public:
    BCEvaluationU(const FEFaceEvaluation<dim, -1, 0, dim, Number> &velocity_m)
      : velocity_m(velocity_m)
    {}

    typename FEFaceEvaluation<dim, -1, 0, dim, Number>::value_type
    get_value(const unsigned int q) const
    {
      return velocity_m.get_value(q);
    }

  private:
    const FEFaceEvaluation<dim, -1, 0, dim, Number> &velocity_m;
  };

  // @sect3{Acoustic operator}
  //
  // The following class then defines the acoustic operator. The class is
  // heavily based on matrix-free methods. For a better understanding in
  // matrix-free methods please refer to step-67.
  //
  // At the top we define a flag that decides whether we want to use
  // mortaring. If the remote evaluators are set up with a
  // VectorizedArray we are using point-to-point interpolation;
  // otherwise we make use of Nitsche-type mortaring. The decision is
  // made using `std::is_floating_point_v`, which is a variable that
  // is true if the template argument is a floating point type (such
  // as `float` or `double`) and false otherwise (in particular if the
  // template argument is a vectorized array type).
  template <int dim, typename Number, typename remote_value_type>
  class AcousticOperator
  {
    static constexpr bool use_mortaring =
      std::is_floating_point_v<remote_value_type>;

  public:
    // In case of Nitsche-type mortaring, `NonMatching::MappingInfo` has to
    // be provided in the constructor.
    AcousticOperator(
      const MatrixFree<dim, Number>                &matrix_free,
      std::shared_ptr<CellwiseMaterialData<Number>> material_data,
      const std::set<types::boundary_id>           &remote_face_ids,
      std::shared_ptr<FERemoteEvaluation<dim, 1, remote_value_type>>
        pressure_r_eval,
      std::shared_ptr<FERemoteEvaluation<dim, dim, remote_value_type>>
        velocity_r_eval,
      std::shared_ptr<FERemoteEvaluation<dim, 1, remote_value_type>> c_r_eval,
      std::shared_ptr<FERemoteEvaluation<dim, 1, remote_value_type>> rho_r_eval,
      std::shared_ptr<NonMatching::MappingInfo<dim, dim, Number>>    nm_info =
        nullptr)
      : matrix_free(matrix_free)
      , material_data(material_data)
      , remote_face_ids(remote_face_ids)
      , pressure_r_eval(pressure_r_eval)
      , velocity_r_eval(velocity_r_eval)
      , c_r_eval(c_r_eval)
      , rho_r_eval(rho_r_eval)
      , nm_mapping_info(nm_info)
    {
      if (use_mortaring)
        Assert(nm_info,
               ExcMessage(
                 "In case of Nitsche-type mortaring NonMatching::MappingInfo \
                  has to be provided."));
    }

    // The following function then evaluates the acoustic operator.
    // It first updates the precomputed values in corresponding the
    // FERemoteEvaluation objects. The material parameters do not change and
    // thus, we do not have to update precomputed values in @c c_r_eval and @c
    // rho_r_eval.
    //
    // It then either performs a matrix-free loop with Nitsche-type
    // mortaring at non-matching faces (if `use_mortaring` is true) or
    // with point-to-point interpolation at non-matching faces (in the
    // `else` branch). The difference is only in the third argument to
    // the loop function, denoting the function object that is
    // executed at boundary faces.
    template <typename VectorType>
    void evaluate(VectorType &dst, const VectorType &src) const
    {
      pressure_r_eval->gather_evaluate(src, EvaluationFlags::values);
      velocity_r_eval->gather_evaluate(src, EvaluationFlags::values);

      if constexpr (use_mortaring)
        {
          matrix_free.loop(
            &AcousticOperator::local_apply_cell<VectorType>,
            &AcousticOperator::local_apply_face<VectorType>,
            &AcousticOperator::local_apply_boundary_face_mortaring<VectorType>,
            this,
            dst,
            src,
            true,
            MatrixFree<dim, Number>::DataAccessOnFaces::values,
            MatrixFree<dim, Number>::DataAccessOnFaces::values);
        }
      else
        {
          matrix_free.loop(
            &AcousticOperator::local_apply_cell<VectorType>,
            &AcousticOperator::local_apply_face<VectorType>,
            &AcousticOperator::local_apply_boundary_face_point_to_point<
              VectorType>,
            this,
            dst,
            src,
            true,
            MatrixFree<dim, Number>::DataAccessOnFaces::values,
            MatrixFree<dim, Number>::DataAccessOnFaces::values);
        }
    }

    // In the `private` section of the class, we then define the
    // functions that evaluate volume, interior face, and boundary
    // face integrals. The concrete terms these functions evaluate are
    // stated in the documentation at the top of this tutorial
    // program. Each of these functions has its own object of type
    // `MaterialEvaluation` that provides access to the material at
    // each cell or face.
  private:
    template <typename VectorType>
    void local_apply_cell(
      const MatrixFree<dim, Number>               &matrix_free,
      VectorType                                  &dst,
      const VectorType                            &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
    {
      FEEvaluation<dim, -1, 0, 1, Number>   pressure(matrix_free, 0, 0, 0);
      FEEvaluation<dim, -1, 0, dim, Number> velocity(matrix_free, 0, 0, 1);

      MaterialEvaluation material(matrix_free, *material_data);

      for (unsigned int cell = cell_range.first; cell < cell_range.second;
           ++cell)
        {
          velocity.reinit(cell);
          pressure.reinit(cell);

          pressure.gather_evaluate(src, EvaluationFlags::gradients);
          velocity.gather_evaluate(src, EvaluationFlags::gradients);

          // Get the materials on the corresponding cell. Since we
          // introduced @c MaterialEvaluation we can write the code
          // independent of whether the material is homogeneous or
          // inhomogeneous.
          material.reinit_cell(cell);
          const auto c   = material.get_speed_of_sound();
          const auto rho = material.get_density();

          for (const unsigned int q : pressure.quadrature_point_indices())
            {
              pressure.submit_value(rho * c * c * velocity.get_divergence(q),
                                    q);
              velocity.submit_value(1.0 / rho * pressure.get_gradient(q), q);
            }

          pressure.integrate_scatter(EvaluationFlags::values, dst);
          velocity.integrate_scatter(EvaluationFlags::values, dst);
        }
    }

    // This next function evaluates the fluxes at faces between cells with the
    // same material. If boundary faces are under consideration fluxes into
    // neighboring faces do not have to be considered which is enforced via
    // `weight_neighbor=false`. For non-matching faces the fluxes into
    // neighboring faces are not considered as well. This is because we iterate
    // over each side of the non-matching face separately (similar to a cell
    // centric loop).
    //
    // In this and following functions, we also introduce the factors
    // $\tau$ and $\gamma$ that are computed from $\rho$ and $c$ along
    // interfaces and that appear in the bilinear forms. Their
    // definitions are provided in the introduction.
    template <bool weight_neighbor,
              typename InternalFaceIntegratorPressure,
              typename InternalFaceIntegratorVelocity,
              typename ExternalFaceIntegratorPressure,
              typename ExternalFaceIntegratorVelocity>
    inline DEAL_II_ALWAYS_INLINE void evaluate_face_kernel(
      InternalFaceIntegratorPressure                           &pressure_m,
      InternalFaceIntegratorVelocity                           &velocity_m,
      ExternalFaceIntegratorPressure                           &pressure_p,
      ExternalFaceIntegratorVelocity                           &velocity_p,
      const typename InternalFaceIntegratorPressure::value_type c,
      const typename InternalFaceIntegratorPressure::value_type rho) const
    {
      const auto tau   = 0.5 * rho * c;
      const auto gamma = 0.5 / (rho * c);

      for (const unsigned int q : pressure_m.quadrature_point_indices())
        {
          const auto n  = pressure_m.normal_vector(q);
          const auto pm = pressure_m.get_value(q);
          const auto um = velocity_m.get_value(q);

          const auto pp = pressure_p.get_value(q);
          const auto up = velocity_p.get_value(q);

          // Compute homogeneous local Lax-Friedrichs fluxes and submit the
          // corresponding values to the integrators.
          const auto momentum_flux =
            0.5 * (pm + pp) + 0.5 * tau * (um - up) * n;
          velocity_m.submit_value(1.0 / rho * (momentum_flux - pm) * n, q);
          if constexpr (weight_neighbor)
            velocity_p.submit_value(1.0 / rho * (momentum_flux - pp) * (-n), q);

          const auto mass_flux = 0.5 * (um + up) + 0.5 * gamma * (pm - pp) * n;
          pressure_m.submit_value(rho * c * c * (mass_flux - um) * n, q);
          if constexpr (weight_neighbor)
            pressure_p.submit_value(rho * c * c * (mass_flux - up) * (-n), q);
        }
    }

    // This function evaluates the fluxes at faces between cells with different
    // materials. This can only happen over non-matching interfaces. Therefore,
    // it is implicitly known that `weight_neighbor=false` and we can omit the
    // parameter.
    template <typename InternalFaceIntegratorPressure,
              typename InternalFaceIntegratorVelocity,
              typename ExternalFaceIntegratorPressure,
              typename ExternalFaceIntegratorVelocity,
              typename MaterialIntegrator>
    void evaluate_face_kernel_inhomogeneous(
      InternalFaceIntegratorPressure                           &pressure_m,
      InternalFaceIntegratorVelocity                           &velocity_m,
      const ExternalFaceIntegratorPressure                     &pressure_p,
      const ExternalFaceIntegratorVelocity                     &velocity_p,
      const typename InternalFaceIntegratorPressure::value_type c,
      const typename InternalFaceIntegratorPressure::value_type rho,
      const MaterialIntegrator                                 &c_r,
      const MaterialIntegrator                                 &rho_r) const
    {
      const auto tau_m   = 0.5 * rho * c;
      const auto gamma_m = 0.5 / (rho * c);

      for (const unsigned int q : pressure_m.quadrature_point_indices())
        {
          const auto c_p           = c_r.get_value(q);
          const auto rho_p         = rho_r.get_value(q);
          const auto tau_p         = 0.5 * rho_p * c_p;
          const auto gamma_p       = 0.5 / (rho_p * c_p);
          const auto tau_sum_inv   = 1.0 / (tau_m + tau_p);
          const auto gamma_sum_inv = 1.0 / (gamma_m + gamma_p);

          const auto n  = pressure_m.normal_vector(q);
          const auto pm = pressure_m.get_value(q);
          const auto um = velocity_m.get_value(q);

          const auto pp = pressure_p.get_value(q);
          const auto up = velocity_p.get_value(q);


          // Compute inhomogeneous fluxes and submit the corresponding values
          // to the integrators.
          const auto momentum_flux =
            pm - tau_m * tau_sum_inv * (pm - pp) +
            tau_m * tau_p * tau_sum_inv * (um - up) * n;
          velocity_m.submit_value(1.0 / rho * (momentum_flux - pm) * n, q);


          const auto mass_flux =
            um - gamma_m * gamma_sum_inv * (um - up) +
            gamma_m * gamma_p * gamma_sum_inv * (pm - pp) * n;

          pressure_m.submit_value(rho * c * c * (mass_flux - um) * n, q);
        }
    }

    // This function evaluates the inner face integrals.
    template <typename VectorType>
    void local_apply_face(
      const MatrixFree<dim, Number>               &matrix_free,
      VectorType                                  &dst,
      const VectorType                            &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
    {
      FEFaceEvaluation<dim, -1, 0, 1, Number> pressure_m(
        matrix_free, true, 0, 0, 0);
      FEFaceEvaluation<dim, -1, 0, 1, Number> pressure_p(
        matrix_free, false, 0, 0, 0);
      FEFaceEvaluation<dim, -1, 0, dim, Number> velocity_m(
        matrix_free, true, 0, 0, 1);
      FEFaceEvaluation<dim, -1, 0, dim, Number> velocity_p(
        matrix_free, false, 0, 0, 1);

      MaterialEvaluation material(matrix_free, *material_data);

      for (unsigned int face = face_range.first; face < face_range.second;
           ++face)
        {
          velocity_m.reinit(face);
          velocity_p.reinit(face);

          pressure_m.reinit(face);
          pressure_p.reinit(face);

          pressure_m.gather_evaluate(src, EvaluationFlags::values);
          pressure_p.gather_evaluate(src, EvaluationFlags::values);

          velocity_m.gather_evaluate(src, EvaluationFlags::values);
          velocity_p.gather_evaluate(src, EvaluationFlags::values);

          material.reinit_face(face);
          evaluate_face_kernel<true>(pressure_m,
                                     velocity_m,
                                     pressure_p,
                                     velocity_p,
                                     material.get_speed_of_sound(),
                                     material.get_density());

          pressure_m.integrate_scatter(EvaluationFlags::values, dst);
          pressure_p.integrate_scatter(EvaluationFlags::values, dst);
          velocity_m.integrate_scatter(EvaluationFlags::values, dst);
          velocity_p.integrate_scatter(EvaluationFlags::values, dst);
        }
    }


    // @sect4{Matrix-free boundary function for point-to-point interpolation}
    //
    // The following function then evaluates the boundary face integrals and the
    // non-matching face integrals using point-to-point interpolation.
    template <typename VectorType>
    void local_apply_boundary_face_point_to_point(
      const MatrixFree<dim, Number>               &matrix_free,
      VectorType                                  &dst,
      const VectorType                            &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
    {
      FEFaceEvaluation<dim, -1, 0, 1, Number> pressure_m(
        matrix_free, true, 0, 0, 0);
      FEFaceEvaluation<dim, -1, 0, dim, Number> velocity_m(
        matrix_free, true, 0, 0, 1);

      BCEvaluationP pressure_bc(pressure_m);
      BCEvaluationU velocity_bc(velocity_m);

      MaterialEvaluation material(matrix_free, *material_data);

      // Remote evaluators.
      auto pressure_r = pressure_r_eval->get_data_accessor();
      auto velocity_r = velocity_r_eval->get_data_accessor();
      auto c_r        = c_r_eval->get_data_accessor();
      auto rho_r      = rho_r_eval->get_data_accessor();

      for (unsigned int face = face_range.first; face < face_range.second;
           ++face)
        {
          velocity_m.reinit(face);
          pressure_m.reinit(face);

          pressure_m.gather_evaluate(src, EvaluationFlags::values);
          velocity_m.gather_evaluate(src, EvaluationFlags::values);

          if (HelperFunctions::is_non_matching_face(
                remote_face_ids, matrix_free.get_boundary_id(face)))
            {
              // If @c face is non-matching we have to query values via the
              // FERemoteEvaluaton objects. This is done by passing the
              // corresponding FERemoteEvaluaton objects to the function that
              // evaluates the kernel. As mentioned above, each side of the
              // non-matching interface is traversed separately and we do not
              // have to consider the neighbor in the kernel. Note, that the
              // values in the FERemoteEvaluaton objects are already updated at
              // this point.

              // For point-to-point interpolation we simply use the
              // corresponding FERemoteEvaluaton objects in combination with the
              // standard FEFaceEvaluation objects.
              velocity_r.reinit(face);
              pressure_r.reinit(face);

              material.reinit_face(face);

              // If we are considering a homogeneous material, do not use the
              // inhomogeneous fluxes. While it would be possible
              // to use the inhomogeneous fluxes they are more expensive to
              // compute.
              if (material.is_homogeneous())
                {
                  evaluate_face_kernel<false>(pressure_m,
                                              velocity_m,
                                              pressure_r,
                                              velocity_r,
                                              material.get_speed_of_sound(),
                                              material.get_density());
                }
              else
                {
                  c_r.reinit(face);
                  rho_r.reinit(face);
                  evaluate_face_kernel_inhomogeneous(
                    pressure_m,
                    velocity_m,
                    pressure_r,
                    velocity_r,
                    material.get_speed_of_sound(),
                    material.get_density(),
                    c_r,
                    rho_r);
                }
            }
          else
            {
              // If @c face is a standard boundary face, evaluate the integral
              // as usual in the matrix free context. To be able to use the same
              // kernel as for inner faces we pass the boundary condition
              // objects to the function that evaluates the kernel. As detailed
              // above `weight_neighbor=false`.
              material.reinit_face(face);
              evaluate_face_kernel<false>(pressure_m,
                                          velocity_m,
                                          pressure_bc,
                                          velocity_bc,
                                          material.get_speed_of_sound(),
                                          material.get_density());
            }

          pressure_m.integrate_scatter(EvaluationFlags::values, dst);
          velocity_m.integrate_scatter(EvaluationFlags::values, dst);
        }
    }

    // @sect4{Matrix-free boundary function for Nitsche-type mortaring}
    //
    // This function evaluates the boundary face integrals and the
    // non-matching face integrals using Nitsche-type mortaring.
    template <typename VectorType>
    void local_apply_boundary_face_mortaring(
      const MatrixFree<dim, Number>               &matrix_free,
      VectorType                                  &dst,
      const VectorType                            &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
    {
      FEFaceEvaluation<dim, -1, 0, 1, Number> pressure_m(
        matrix_free, true, 0, 0, 0);
      FEFaceEvaluation<dim, -1, 0, dim, Number> velocity_m(
        matrix_free, true, 0, 0, 1);

      // For Nitsche-type mortaring we are evaluating the integrals over
      // intersections. This is why, quadrature points are arbitrarily
      // distributed on every face. Thus, we can not make use of face batches
      // and FEFaceEvaluation but have to consider each face individually and
      // make use of @c FEFacePointEvaluation to evaluate the integrals in the
      // arbitrarily distributed quadrature points.
      // Since the setup of FEFacePointEvaluation is more expensive than that of
      // FEEvaluation we do the setup only once. For this we are using the
      // helper function @c get_thread_safe_fe_face_point_evaluation_object().
      FEFacePointEvaluation<1, dim, dim, Number> &pressure_m_mortar =
        get_thread_safe_fe_face_point_evaluation_object<1>(
          thread_local_pressure_m_mortar, 0);
      FEFacePointEvaluation<dim, dim, dim, Number> &velocity_m_mortar =
        get_thread_safe_fe_face_point_evaluation_object<dim>(
          thread_local_velocity_m_mortar, 1);

      BCEvaluationP pressure_bc(pressure_m);
      BCEvaluationU velocity_bc(velocity_m);

      MaterialEvaluation material(matrix_free, *material_data);

      auto pressure_r_mortar = pressure_r_eval->get_data_accessor();
      auto velocity_r_mortar = velocity_r_eval->get_data_accessor();
      auto c_r               = c_r_eval->get_data_accessor();
      auto rho_r             = rho_r_eval->get_data_accessor();

      for (unsigned int face = face_range.first; face < face_range.second;
           ++face)
        {
          if (HelperFunctions::is_non_matching_face(
                remote_face_ids, matrix_free.get_boundary_id(face)))
            {
              material.reinit_face(face);

              // First fetch the DoF values with standard FEFaceEvaluation
              // objects.
              pressure_m.reinit(face);
              velocity_m.reinit(face);

              pressure_m.read_dof_values(src);
              velocity_m.read_dof_values(src);

              // Project the internally stored values into the face DoFs
              // of the current face.
              pressure_m.project_to_face(EvaluationFlags::values);
              velocity_m.project_to_face(EvaluationFlags::values);

              // For mortaring, we have to consider every face from the face
              // batches separately and have to use the FEFacePointEvaluation
              // objects to be able to evaluate the integrals with the
              // arbitrarily distributed quadrature points.
              for (unsigned int v = 0;
                   v < matrix_free.n_active_entries_per_face_batch(face);
                   ++v)
                {
                  constexpr unsigned int n_lanes =
                    VectorizedArray<Number>::size();
                  velocity_m_mortar.reinit(face * n_lanes + v);
                  pressure_m_mortar.reinit(face * n_lanes + v);

                  // Evaluate using FEFacePointEvaluation. As buffer,
                  // simply use the internal buffers from the
                  // FEFaceEvaluation objects.
                  velocity_m_mortar.evaluate_in_face(
                    &velocity_m.get_scratch_data().begin()[0][v],
                    EvaluationFlags::values);

                  pressure_m_mortar.evaluate_in_face(
                    &pressure_m.get_scratch_data().begin()[0][v],
                    EvaluationFlags::values);

                  velocity_r_mortar.reinit(face * n_lanes + v);
                  pressure_r_mortar.reinit(face * n_lanes + v);

                  // As above, if we are considering a homogeneous
                  // material, do not use the inhomogeneous
                  // fluxes. Since we are operating on face @c v we
                  // call @c material.get_density()[v].
                  if (material.is_homogeneous())
                    {
                      evaluate_face_kernel<false>(
                        pressure_m_mortar,
                        velocity_m_mortar,
                        pressure_r_mortar,
                        velocity_r_mortar,
                        material.get_speed_of_sound()[v],
                        material.get_density()[v]);
                    }
                  else
                    {
                      c_r.reinit(face * n_lanes + v);
                      rho_r.reinit(face * n_lanes + v);

                      evaluate_face_kernel_inhomogeneous(
                        pressure_m_mortar,
                        velocity_m_mortar,
                        pressure_r_mortar,
                        velocity_r_mortar,
                        material.get_speed_of_sound()[v],
                        material.get_density()[v],
                        c_r,
                        rho_r);
                    }

                  // Integrate using FEFacePointEvaluation. As buffer,
                  // simply use the internal buffers from the
                  // FEFaceEvaluation objects.
                  velocity_m_mortar.integrate_in_face(
                    &velocity_m.get_scratch_data().begin()[0][v],
                    EvaluationFlags::values);

                  pressure_m_mortar.integrate_in_face(
                    &pressure_m.get_scratch_data().begin()[0][v],
                    EvaluationFlags::values);
                }

              // Collect the contributions from the face DoFs to
              // the internal cell DoFs to be able to use the
              // member function @c distribute_local_to_global().
              pressure_m.collect_from_face(EvaluationFlags::values);
              velocity_m.collect_from_face(EvaluationFlags::values);

              pressure_m.distribute_local_to_global(dst);
              velocity_m.distribute_local_to_global(dst);
            }
          else
            {
              velocity_m.reinit(face);
              pressure_m.reinit(face);

              pressure_m.gather_evaluate(src, EvaluationFlags::values);
              velocity_m.gather_evaluate(src, EvaluationFlags::values);

              material.reinit_face(face);
              evaluate_face_kernel<false>(pressure_m,
                                          velocity_m,
                                          pressure_bc,
                                          velocity_bc,
                                          material.get_speed_of_sound(),
                                          material.get_density());

              pressure_m.integrate_scatter(EvaluationFlags::values, dst);
              velocity_m.integrate_scatter(EvaluationFlags::values, dst);
            }
        }
    }

    const MatrixFree<dim, Number> &matrix_free;

    const std::shared_ptr<CellwiseMaterialData<Number>> material_data;

    const std::set<types::boundary_id> remote_face_ids;

    // FERemoteEvaluation objects are stored as shared pointers. This way,
    // they can also be used for other operators without precomputing the values
    // multiple times.
    const std::shared_ptr<FERemoteEvaluation<dim, 1, remote_value_type>>
      pressure_r_eval;
    const std::shared_ptr<FERemoteEvaluation<dim, dim, remote_value_type>>
      velocity_r_eval;

    const std::shared_ptr<FERemoteEvaluation<dim, 1, remote_value_type>>
      c_r_eval;
    const std::shared_ptr<FERemoteEvaluation<dim, 1, remote_value_type>>
      rho_r_eval;

    const std::shared_ptr<NonMatching::MappingInfo<dim, dim, Number>>
      nm_mapping_info;

    // We store FEFacePointEvaluation objects as members in a thread local
    // way, since its creation is more expensive compared to FEEvaluation
    // objects.
    mutable Threads::ThreadLocalStorage<
      std::unique_ptr<FEFacePointEvaluation<1, dim, dim, Number>>>
      thread_local_pressure_m_mortar;

    mutable Threads::ThreadLocalStorage<
      std::unique_ptr<FEFacePointEvaluation<dim, dim, dim, Number>>>
      thread_local_velocity_m_mortar;

    // Helper function to create and get FEFacePointEvaluation objects in a
    // thread safe way. On each thread, FEFacePointEvaluation is created if it
    // has not been created by now. After that, simply return the object
    // corresponding to the thread under consideration.
    template <int n_components>
    FEFacePointEvaluation<n_components, dim, dim, Number> &
    get_thread_safe_fe_face_point_evaluation_object(
      Threads::ThreadLocalStorage<
        std::unique_ptr<FEFacePointEvaluation<n_components, dim, dim, Number>>>
                  &fe_face_point_eval_thread_local,
      unsigned int fist_selected_comp) const
    {
      if (fe_face_point_eval_thread_local.get() == nullptr)
        {
          fe_face_point_eval_thread_local = std::make_unique<
            FEFacePointEvaluation<n_components, dim, dim, Number>>(
            *nm_mapping_info,
            matrix_free.get_dof_handler().get_fe(),
            true,
            fist_selected_comp);
        }
      return *fe_face_point_eval_thread_local.get();
    }
  };

  // @sect3{Inverse mass operator}
  //
  // For the time stepping methods below, we need the inverse mass
  // operator. We apply it via a loop over all (batches of) cells as
  // always, where the contribution of each cell is computed in a
  // matrix-free way:
  template <int dim, typename Number>
  class InverseMassOperator
  {
  public:
    InverseMassOperator(const MatrixFree<dim, Number> &matrix_free)
      : matrix_free(matrix_free)
    {}

    template <typename VectorType>
    void apply(VectorType &dst, const VectorType &src) const
    {
      dst.zero_out_ghost_values();
      matrix_free.cell_loop(&InverseMassOperator::local_apply_cell<VectorType>,
                            this,
                            dst,
                            src);
    }

  private:
    template <typename VectorType>
    void local_apply_cell(
      const MatrixFree<dim, Number>               &mf,
      VectorType                                  &dst,
      const VectorType                            &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
    {
      FEEvaluation<dim, -1, 0, dim + 1, Number> phi(mf);
      MatrixFreeOperators::CellwiseInverseMassMatrix<dim, -1, dim + 1, Number>
        minv(phi);

      for (unsigned int cell = cell_range.first; cell < cell_range.second;
           ++cell)
        {
          phi.reinit(cell);
          phi.read_dof_values(src);
          minv.apply(phi.begin_dof_values(), phi.begin_dof_values());
          phi.set_dof_values(dst);
        }
    }

    const MatrixFree<dim, Number> &matrix_free;
  };

  // @sect3{Runge-Kutta time-stepping}
  //
  // This class implements a Runge-Kutta scheme of order 2.
  template <int dim, typename Number, typename remote_value_type>
  class RungeKutta2
  {
    using VectorType = LinearAlgebra::distributed::Vector<Number>;

  public:
    RungeKutta2(
      const std::shared_ptr<InverseMassOperator<dim, Number>>
        inverse_mass_operator,
      const std::shared_ptr<AcousticOperator<dim, Number, remote_value_type>>
        acoustic_operator)
      : inverse_mass_operator(inverse_mass_operator)
      , acoustic_operator(acoustic_operator)
    {}

    // The `run()` function of this class contains the time loop. It
    // starts by initializing some member variables (such as
    // short-cuts to objects that describe the finite element, its
    // properties, and the mapping; an by initializing vectors). It
    // then computes the time step size via minimum element edge
    // length.  Assuming non-distorted elements, we can compute the
    // edge length as the distance between two vertices.  From this,
    // we can then obtain the time step size via the CFL condition.
    void run(const MatrixFree<dim, Number> &matrix_free,
             const double                   cr,
             const double                   end_time,
             const double                   speed_of_sound,
             const Function<dim>           &initial_condition,
             const std::string             &vtk_prefix)
    {
      const auto &dof_handler = matrix_free.get_dof_handler();
      const auto &mapping     = *matrix_free.get_mapping_info().mapping;
      const auto  degree      = dof_handler.get_fe().degree;

      VectorType solution;
      matrix_free.initialize_dof_vector(solution);
      VectorType solution_temp;
      matrix_free.initialize_dof_vector(solution_temp);

      HelperFunctions::set_initial_condition(matrix_free,
                                             initial_condition,
                                             solution);

      double h_local_min = std::numeric_limits<double>::max();
      for (const auto &cell : dof_handler.active_cell_iterators())
        h_local_min =
          std::min(h_local_min,
                   (cell->vertex(1) - cell->vertex(0)).norm_square());
      h_local_min = std::sqrt(h_local_min);
      const double h_min =
        Utilities::MPI::min(h_local_min, dof_handler.get_mpi_communicator());

      const double dt =
        cr * HelperFunctions::compute_dt_cfl(h_min, degree, speed_of_sound);

      // We can then perform the time integration loop:
      double       time     = 0.0;
      unsigned int timestep = 0;
      while (time < end_time)
        {
          HelperFunctions::write_vtu(solution,
                                     matrix_free.get_dof_handler(),
                                     mapping,
                                     degree,
                                     "step_89-" + vtk_prefix +
                                       std::to_string(timestep));

          std::swap(solution, solution_temp);
          time += dt;
          ++timestep;
          perform_time_step(dt, solution, solution_temp);
        }
    }

    // The main work of this class is done by a `private` member
    // function that performs one Runge-Kutta 2 time step. Recall that
    // this method requires two sub-steps ("stages") computing
    // intermediate values `k1` and `k2`, after which the intermediate
    // values are summed with weights to obtain the new solution at
    // the end of the time step. The RK2 method allows for the
    // elimination of the intermediate vector `k2` by utilizing the
    // `dst` vector as temporary storage.
  private:
    void
    perform_time_step(const double dt, VectorType &dst, const VectorType &src)
    {
      VectorType k1 = src;

      /* First stage. */
      evaluate_stage(k1, src);

      /* Second stage. */
      k1.sadd(0.5 * dt, 1.0, src);
      evaluate_stage(dst, k1);

      /* Summing things into the output vector. */
      dst.sadd(dt, 1.0, src);
    }

    // Evaluating a single Runge-Kutta stage is a straightforward step
    // that really only requires applying the operator once, and then
    // applying the inverse of the mass matrix.
    void evaluate_stage(VectorType &dst, const VectorType &src)
    {
      acoustic_operator->evaluate(dst, src);
      dst *= -1.0;
      inverse_mass_operator->apply(dst, dst);
    }

    const std::shared_ptr<InverseMassOperator<dim, Number>>
      inverse_mass_operator;
    const std::shared_ptr<AcousticOperator<dim, Number, remote_value_type>>
      acoustic_operator;
  };


  // @sect3{Construction of non-matching triangulations}
  //
  // Let us now make our way to the higher-level functions of this program.
  //
  // The first of these functions creates a two dimensional square
  // triangulation that spans from $(0,0)$ to $(1,1)$. It consists of
  // two sub-domains.  The left sub-domain spans from $(0,0)$ to
  // $(0.525,1)$. The right sub-domain spans from $(0.525,0)$ to
  // $(1,1)$. The left sub-domain has elements that are three times
  // smaller compared to the ones for the right sub-domain.
  //
  // At non-matching interfaces, we need to provide different boundary
  // IDs for the cells that make up the two parts (because, while they
  // may be physically adjacent, they are not logically neighbors
  // given that the faces of cells on both sides do not match, and the
  // Triangulation class will therefore treat the interface between
  // the two parts as a "boundary"). These boundary IDs have to differ
  // because later on RemotePointEvaluation has to search for remote
  // points for each face, that are defined in the same mesh (since we
  // merge the mesh) but not on the same side of the non-matching
  // interface. As a consequence, we declare at the top symbolic names
  // for these boundary indicators, and ensure that we return a set
  // with these values to the caller for later use.
  //
  // The actual mesh is then constructed by first constructing the
  // left and right parts separately (setting material ids to zero and
  // one, respectively), and using the appropriate boundary ids for
  // all parts of the mesh. We then use
  // GridGenerator::merge_triangulations() to combine them into one
  // (non-matching) mesh. We have to pay attention that should the two
  // sub-triangulations have vertices at the same locations, that they
  // are not merged (connecting the two triangulations logically)
  // since we want the interface to be an actual boundary. We achieve
  // this by providing a tolerance of zero for the merge, see the
  // documentation of the function
  // GridGenerator::merge_triangulations().
  template <int dim>
  void build_non_matching_triangulation(
    Triangulation<dim>           &tria,
    std::set<types::boundary_id> &non_matching_faces,
    const unsigned int            refinements)
  {
    const double length = 1.0;

    const types::boundary_id non_matching_id_left  = 98;
    const types::boundary_id non_matching_id_right = 99;

    non_matching_faces = {non_matching_id_left, non_matching_id_right};

    /* Construct left part of mesh. */
    Triangulation<dim> tria_left;
    const unsigned int subdiv_left = 11;
    GridGenerator::subdivided_hyper_rectangle(tria_left,
                                              {subdiv_left, 2 * subdiv_left},
                                              {0.0, 0.0},
                                              {0.525 * length, length});

    for (const auto &cell : tria_left.active_cell_iterators())
      cell->set_material_id(0);
    for (const auto &face : tria_left.active_face_iterators())
      if (face->at_boundary())
        {
          face->set_boundary_id(0);
          if (face->center()[0] > 0.525 * length - 1e-6)
            face->set_boundary_id(non_matching_id_left);
        }

    /* Construct right part of mesh. */
    Triangulation<dim> tria_right;
    const unsigned int subdiv_right = 4;
    GridGenerator::subdivided_hyper_rectangle(tria_right,
                                              {subdiv_right, 2 * subdiv_right},
                                              {0.525 * length, 0.0},
                                              {length, length});

    for (const auto &cell : tria_right.active_cell_iterators())
      cell->set_material_id(1);
    for (const auto &face : tria_right.active_face_iterators())
      if (face->at_boundary())
        {
          face->set_boundary_id(0);
          if (face->center()[0] < 0.525 * length + 1e-6)
            face->set_boundary_id(non_matching_id_right);
        }

    /* Merge triangulations. */
    GridGenerator::merge_triangulations(tria_left,
                                        tria_right,
                                        tria,
                                        /*tolerance*/ 0.,
                                        /*copy_manifold_ids*/ false,
                                        /*copy_boundary_ids*/ true);

    /* Refine the result. */
    tria.refine_global(refinements);
  }

  // @sect3{Set up and running of the two schemes}

  // @sect4{Setup and running of the point-to-point interpolation scheme}
  //
  // We are now at the two functions that run the overall schemes (the
  // point-to-point and the mortaring schemes).  The first of these
  // functions fills a FERemoteEvaluationCommunicator object that is
  // needed for point-to-point interpolation. Additionally, the
  // corresponding remote evaluators are set up using this remote
  // communicator. Eventually, the operators are handed to the time
  // integrator that runs the simulation.
  template <int dim, typename Number>
  void run_with_point_to_point_interpolation(
    const MatrixFree<dim, Number>      &matrix_free,
    const std::set<types::boundary_id> &non_matching_faces,
    const std::map<types::material_id, std::pair<double, double>> &materials,
    const double                                                   end_time,
    const Function<dim> &initial_condition,
    const std::string   &vtk_prefix)
  {
    const auto &dof_handler = matrix_free.get_dof_handler();
    const auto &tria        = dof_handler.get_triangulation();

    // Communication objects know about the communication pattern. That is,
    // they know about the cells and quadrature points that have to be
    // evaluated at remote faces. This information is given via
    // RemotePointEvaluation. Additionally, the communication objects
    // have to be able to match the quadrature points of the remote
    // points (that provide exterior information) to the quadrature points
    // defined at the interior cell. In case of point-to-point interpolation
    // a vector of pairs with face batch Ids and the number of faces in the
    // batch is needed. @c FERemoteCommunicationObjectEntityBatches
    // is a container to store this information.
    //
    // The information is filled outside of the actual class since in some cases
    // the information is available from some heuristic and
    // it is possible to skip some expensive operations. This is for example
    // the case for sliding rotating interfaces with equally spaced elements on
    // both sides of the non-matching interface @cite duerrwaechter2021an.
    //
    // For the standard case of point to point-to-point interpolation without
    // any heuristic we make use of the utility function
    // Utilities::compute_remote_communicator_faces_point_to_point_interpolation().
    // Please refer to the documentation of this function to see how to manually
    // set up the remote communicator from outside.
    //
    // Among the inputs for the remote communicator we need a list of
    // boundary ids for the non-matching faces, along with a function
    // object for each boundary id that returns a vector of true/false
    // flags in which exactly the vertices of cells of the
    // triangulation are marked that have a face at the boundary id in
    // question.
    std::vector<
      std::pair<types::boundary_id, std::function<std::vector<bool>()>>>
      non_matching_faces_marked_vertices;
    for (const auto &nm_face : non_matching_faces)
      {
        auto marked_vertices = [&nm_face, &tria]() -> std::vector<bool> {
          std::vector<bool> mask(tria.n_vertices(), true);

          for (const auto &cell : tria.active_cell_iterators())
            for (auto const &f : cell->face_indices())
              if (cell->face(f)->at_boundary() &&
                  cell->face(f)->boundary_id() == nm_face)
                for (const auto v : cell->vertex_indices())
                  mask[cell->vertex_index(v)] = false;

          return mask;
        };

        non_matching_faces_marked_vertices.emplace_back(
          std::make_pair(nm_face, marked_vertices));
      }

    auto remote_communicator =
      Utilities::compute_remote_communicator_faces_point_to_point_interpolation(
        matrix_free, non_matching_faces_marked_vertices);

    // We are using point-to-point interpolation and can therefore
    // easily access the corresponding data at face batches. This
    // is why we use a @c VectorizedArray as @c remote_value_type:
    using remote_value_type = VectorizedArray<Number>;

    // We then set up FERemoteEvaluation objects that access the
    // pressure and velocity at remote faces, along with an object to
    // describe cell-wise material data.
    const auto pressure_r =
      std::make_shared<FERemoteEvaluation<dim, 1, remote_value_type>>(
        remote_communicator, dof_handler, /*first_selected_component*/ 0);

    const auto velocity_r =
      std::make_shared<FERemoteEvaluation<dim, dim, remote_value_type>>(
        remote_communicator, dof_handler, /*first_selected_component*/ 1);

    const auto material_data =
      std::make_shared<CellwiseMaterialData<Number>>(matrix_free, materials);

    // If we have an inhomogeneous problem, we have to set up the
    // material handler that accesses the materials at remote faces.
    const auto c_r =
      std::make_shared<FERemoteEvaluation<dim, 1, remote_value_type>>(
        remote_communicator,
        matrix_free.get_dof_handler().get_triangulation(),
        /*first_selected_component*/ 0);
    const auto rho_r =
      std::make_shared<FERemoteEvaluation<dim, 1, remote_value_type>>(
        remote_communicator,
        matrix_free.get_dof_handler().get_triangulation(),
        /*first_selected_component*/ 0);

    // If the domain is not homogeneous, i.e., if material parameters
    // change from cell to cell, we initialize and fill DoF vectors
    // that contain the material properties.  Materials do not change
    // during the simulation, therefore there is no need to ever
    // compute the values after the first @c gather_evaluate() (at the
    // end of the following block) again.
    if (!material_data->is_homogeneous())
      {
        Vector<Number> c(
          matrix_free.get_dof_handler().get_triangulation().n_active_cells());
        Vector<Number> rho(
          matrix_free.get_dof_handler().get_triangulation().n_active_cells());

        for (const auto &cell : matrix_free.get_dof_handler()
                                  .get_triangulation()
                                  .active_cell_iterators())
          {
            c[cell->active_cell_index()] =
              materials.at(cell->material_id()).first;
            rho[cell->active_cell_index()] =
              materials.at(cell->material_id()).second;
          }

        c_r->gather_evaluate(c, EvaluationFlags::values);
        rho_r->gather_evaluate(rho, EvaluationFlags::values);
      }


    // Next, we set up the inverse mass operator and the acoustic
    // operator. Using `remote_value_type=VectorizedArray<Number>`
    // makes the operator use point-to-point interpolation. These two
    // objects are then used to create a `RungeKutta2` object to
    // perform the time integration.
    //
    // We also compute the maximum speed of sound, needed for the
    // computation of the time-step size, and then run the time integrator.
    // For the examples considered here, we found a limiting Courant number of
    // $\mathrm{Cr}\approx 0.36$ to maintain stability. To ensure, the
    // error of the temporal discretization is small, we use a considerably
    // smaller Courant number of $0.2$.
    const auto inverse_mass_operator =
      std::make_shared<InverseMassOperator<dim, Number>>(matrix_free);

    const auto acoustic_operator =
      std::make_shared<AcousticOperator<dim, Number, remote_value_type>>(
        matrix_free,
        material_data,
        non_matching_faces,
        pressure_r,
        velocity_r,
        c_r,
        rho_r);

    RungeKutta2<dim, Number, remote_value_type> time_integrator(
      inverse_mass_operator, acoustic_operator);

    double speed_of_sound_max = 0.0;
    for (const auto &mat : materials)
      speed_of_sound_max = std::max(speed_of_sound_max, mat.second.first);


    time_integrator.run(matrix_free,
                        /*Cr*/ 0.2,
                        end_time,
                        speed_of_sound_max,
                        initial_condition,
                        vtk_prefix);
  }

  // @sect4{Setup and running of the Nitsche-type mortaring scheme}
  //
  // The alternative to the previous function is to use the mortaring
  // scheme -- implemented in the following function.  This function
  // can only be run when deal.II is configured using CGAL (but the
  // previous function can still be used without CGAL), and so errors
  // out if CGAL is not available.
  template <int dim, typename Number>
  void run_with_nitsche_type_mortaring(
    const MatrixFree<dim, Number>      &matrix_free,
    const std::set<types::boundary_id> &non_matching_faces,
    const std::map<types::material_id, std::pair<double, double>> &materials,
    const double                                                   end_time,
    const Function<dim> &initial_condition,
    const std::string   &vtk_prefix)
  {
#ifndef DEAL_II_WITH_CGAL
    (void)matrix_free;
    (void)non_matching_faces;
    (void)materials;
    (void)end_time;
    (void)initial_condition;
    (void)vtk_prefix;

    ConditionalOStream pcout(
      std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

    pcout << "In this function, mortars are computed using CGAL. "
             "Configure deal.II with DEAL_II_WITH_CGAL to run this function.\n";

    return;
#else

    const auto &dof_handler       = matrix_free.get_dof_handler();
    const auto &tria              = dof_handler.get_triangulation();
    const auto &mapping           = *matrix_free.get_mapping_info().mapping;
    const auto  n_quadrature_pnts = matrix_free.get_quadrature().size();

    std::vector<
      std::pair<types::boundary_id, std::function<std::vector<bool>()>>>
      non_matching_faces_marked_vertices;

    for (const auto &nm_face : non_matching_faces)
      {
        auto marked_vertices = [&]() {
          std::vector<bool> mask(tria.n_vertices(), true);

          for (const auto &cell : tria.active_cell_iterators())
            for (auto const &f : cell->face_indices())
              if (cell->face(f)->at_boundary() &&
                  cell->face(f)->boundary_id() == nm_face)
                for (const auto v : cell->vertex_indices())
                  mask[cell->vertex_index(v)] = false;

          return mask;
        };

        non_matching_faces_marked_vertices.emplace_back(
          std::make_pair(nm_face, marked_vertices));
      }

    // The only parts in this function that are functionally different
    // from the previous one follow here.
    //
    // First, quadrature points are arbitrarily distributed on each non-matching
    // face. Therefore, we have to make use of FEFacePointEvaluation.
    // FEFacePointEvaluation needs NonMatching::MappingInfo to work at the
    // correct quadrature points that are in sync with used FERemoteEvaluation
    // object. Using
    // `compute_remote_communicator_faces_nitsche_type_mortaring()` to reinit
    // NonMatching::MappingInfo ensures this. In the case of mortaring, we have
    // to use the weights provided by the quadrature rules that are used to set
    // up NonMatching::MappingInfo. Therefore we set the flag @c
    // use_global_weights.
    typename NonMatching::MappingInfo<dim, dim, Number>::AdditionalData
      additional_data;
    additional_data.use_global_weights = true;

    // Set up NonMatching::MappingInfo with needed update flags and
    // @c additional_data.
    auto nm_mapping_info =
      std::make_shared<NonMatching::MappingInfo<dim, dim, Number>>(
        mapping,
        update_values | update_JxW_values | update_normal_vectors |
          update_quadrature_points,
        additional_data);

    auto remote_communicator =
      Utilities::compute_remote_communicator_faces_nitsche_type_mortaring(
        matrix_free,
        non_matching_faces_marked_vertices,
        n_quadrature_pnts,
        0,
        nm_mapping_info.get());

    // Second, since quadrature points are arbitrarily distributed we
    // have to consider each face in a batch separately and can not
    // make use of @c VecorizedArray.
    using remote_value_type = Number;

    // The rest of the code is then identical to what we had in the
    // previous function (though it functions differently because of
    // the difference in `remote_value_type`).
    const auto pressure_r =
      std::make_shared<FERemoteEvaluation<dim, 1, remote_value_type>>(
        remote_communicator, dof_handler, /*first_selected_component*/ 0);

    const auto velocity_r =
      std::make_shared<FERemoteEvaluation<dim, dim, remote_value_type>>(
        remote_communicator, dof_handler, /*first_selected_component*/ 1);

    const auto material_data =
      std::make_shared<CellwiseMaterialData<Number>>(matrix_free, materials);

    const auto c_r =
      std::make_shared<FERemoteEvaluation<dim, 1, remote_value_type>>(
        remote_communicator,
        matrix_free.get_dof_handler().get_triangulation(),
        /*first_selected_component*/ 0);
    const auto rho_r =
      std::make_shared<FERemoteEvaluation<dim, 1, remote_value_type>>(
        remote_communicator,
        matrix_free.get_dof_handler().get_triangulation(),
        /*first_selected_component*/ 0);

    if (!material_data->is_homogeneous())
      {
        Vector<Number> c(
          matrix_free.get_dof_handler().get_triangulation().n_active_cells());
        Vector<Number> rho(
          matrix_free.get_dof_handler().get_triangulation().n_active_cells());

        for (const auto &cell : matrix_free.get_dof_handler()
                                  .get_triangulation()
                                  .active_cell_iterators())
          {
            c[cell->active_cell_index()] =
              materials.at(cell->material_id()).first;
            rho[cell->active_cell_index()] =
              materials.at(cell->material_id()).second;
          }

        c_r->gather_evaluate(c, EvaluationFlags::values);
        rho_r->gather_evaluate(rho, EvaluationFlags::values);
      }

    const auto inverse_mass_operator =
      std::make_shared<InverseMassOperator<dim, Number>>(matrix_free);

    const auto acoustic_operator =
      std::make_shared<AcousticOperator<dim, Number, remote_value_type>>(
        matrix_free,
        material_data,
        non_matching_faces,
        pressure_r,
        velocity_r,
        c_r,
        rho_r,
        nm_mapping_info);

    RungeKutta2<dim, Number, remote_value_type> time_integrator(
      inverse_mass_operator, acoustic_operator);

    double speed_of_sound_max = 0.0;
    for (const auto &mat : materials)
      speed_of_sound_max = std::max(speed_of_sound_max, mat.second.first);

    time_integrator.run(matrix_free,
                        /*Cr*/ 0.2,
                        end_time,
                        speed_of_sound_max,
                        initial_condition,
                        vtk_prefix);
#endif
  }
} // namespace Step89


// @sect3{main()}
//
// Finally, the `main()` function executes the different versions of handling
// non-matching interfaces.
//
// Similar to step-87, the minimum requirement of this tutorial is MPI.
// The parallel::distributed::Triangulation class is used if deal.II is
// configured with p4est. Otherwise parallel::shared::Triangulation is used.
int main(int argc, char *argv[])
{
  using namespace dealii;
  constexpr int dim = 2;
  using Number      = double;

  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  std::cout.precision(5);
  ConditionalOStream pcout(std::cout,
                           (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                            0));

  const unsigned int refinements = 1;
  const unsigned int degree      = 3;

#ifdef DEAL_II_WITH_P4EST
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
#else
  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD);
#endif

  pcout << "Create non-matching grid..." << std::endl;

  std::set<types::boundary_id> non_matching_faces;
  Step89::build_non_matching_triangulation(tria,
                                           non_matching_faces,
                                           refinements);

  pcout << " - Refinement level: " << refinements << std::endl;
  pcout << " - Number of cells: " << tria.n_cells() << std::endl;

  pcout << "Create DoFHandler..." << std::endl;
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(FE_DGQ<dim>(degree) ^ (dim + 1)));
  pcout << " - Number of DoFs: " << dof_handler.n_dofs() << std::endl;

  AffineConstraints<Number> constraints;
  constraints.close();

  pcout << "Set up MatrixFree..." << std::endl;
  typename MatrixFree<dim, Number>::AdditionalData data;
  data.mapping_update_flags             = update_gradients | update_values;
  data.mapping_update_flags_inner_faces = update_values;
  data.mapping_update_flags_boundary_faces =
    update_quadrature_points | update_values;

  MatrixFree<dim, Number> matrix_free;
  matrix_free.reinit(
    MappingQ1<dim>(), dof_handler, constraints, QGauss<dim>(degree + 1), data);


  // @sect4{Run vibrating membrane test case} Homogeneous pressure}
  // Dirichlet boundary conditions are applied for
  // simplicity. Therefore, modes can not be chosen arbitrarily.
  pcout << "Run vibrating membrane test case..." << std::endl;
  const double                                            modes = 10.0;
  std::map<types::material_id, std::pair<double, double>> homogeneous_material;
  homogeneous_material[numbers::invalid_material_id] = std::make_pair(1.0, 1.0);
  const auto initial_solution_membrane =
    Step89::InitialConditionVibratingMembrane<dim>(modes);

  /* Run vibrating membrane test case using point-to-point interpolation: */
  pcout << " - Point-to-point interpolation: " << std::endl;
  Step89::run_with_point_to_point_interpolation(
    matrix_free,
    non_matching_faces,
    homogeneous_material,
    8.0 * initial_solution_membrane.get_period_duration(
            homogeneous_material.begin()->second.first),
    initial_solution_membrane,
    "vm-p2p");

  /* Run vibrating membrane test case using Nitsche-type mortaring: */
  pcout << " - Nitsche-type mortaring: " << std::endl;
  Step89::run_with_nitsche_type_mortaring(
    matrix_free,
    non_matching_faces,
    homogeneous_material,
    8.0 * initial_solution_membrane.get_period_duration(
            homogeneous_material.begin()->second.first),
    initial_solution_membrane,
    "vm-nitsche");

  // @sect4{Run test case with inhomogeneous material}
  // Run simple test case with inhomogeneous material and Nitsche-type
  // mortaring:
  pcout << "Run test case with inhomogeneous material..." << std::endl;
  std::map<types::material_id, std::pair<double, double>>
    inhomogeneous_material;
  inhomogeneous_material[0] = std::make_pair(1.0, 1.0);
  inhomogeneous_material[1] = std::make_pair(3.0, 1.0);
  Step89::run_with_nitsche_type_mortaring(matrix_free,
                                          non_matching_faces,
                                          inhomogeneous_material,
                                          /*runtime*/ 0.3,
                                          Step89::GaussPulse<dim>(0.3, 0.5),
                                          "inhomogeneous");


  return 0;
}
