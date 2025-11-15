#include <deal.II/base/function.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_hermite.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_interface_values.h>

#include <deal.II/non_matching/fe_immersed_values.h>
#include <deal.II/non_matching/fe_values.h>
#include <deal.II/non_matching/mesh_classifier.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <new>

#define DIM 2
#define PLANE 0
#define BESSEL_ZERO 2.40483

std::string dtos(const double d)
{
    std::ostringstream s;
    s << d;
    return s.str();
}

namespace CutHermiteWaveCircle 
{
    using namespace dealii;

    template <int dim>
    class WaveSystemValues
    {
    public:
        WaveSystemValues() :
            frequency(1.),
            wavespeed(1.)
                {
                    wavenumber = Point<dim>();
                    wavenumber(0) = 1.0;
                    for (int i = 1; i < dim; ++i)
                        wavenumber(i) = 0.;
                };
        WaveSystemValues(
            const Point<dim>& wavenumber, 
            const double frequency) : 
                wavenumber(wavenumber),
                frequency(frequency)
                {
                    wavespeed = frequency / wavenumber.norm();
                };

        Point<dim> wavenumber;
        double frequency;
        double wavespeed;
    }; //WaveSystemValues

    class LevelSet : public Function<DIM>
    {
    public:
        LevelSet(const double radius = 1.) :
            radius(radius) {}

        virtual double
        value (const Point<DIM> &p, const unsigned int component) const override
        {
            (void) component;
            double temp = 0;
            for (unsigned int i = 0; i < DIM; ++i)
                temp += p(i) * p(i);
            return temp - radius * radius;
        }

        const double radius;
    };

    class Solution : public Function<DIM>
    {
    public:
        Solution(
            const WaveSystemValues<DIM>* vals = nullptr,
            const double start_time = 0.) :
                systemvals(vals),
                current_time(start_time) 
                {
                    Assert(
                        systemvals != nullptr, 
                        ExcMessage("ERROR: Solution function needs a non-null pointer to system values.")
                        );

                    double lower = 2, upper = 3, middle = 0.;
                    double result_upper = 0., result_lower = 0., result_middle = 0.;

                    double TOL = 10e-15;
                    int n = 0, NN = 50000;

                    result_lower = std::cyl_bessel_j(0,lower);
                    result_upper = std::cyl_bessel_j(0, upper);

                    while ((result_lower > TOL) && n < NN)
                    {
                        ++n;
                        middle = 0.5 * (lower + upper);
                        result_middle = std::cyl_bessel_j(0, middle);

                        if (result_middle < 0)
                        {
                            upper = middle;
                            result_upper = result_middle;
                        }
                        else
                        {
                            lower = middle;
                            result_lower = result_middle;
                        }
                    }

                    this->bessel_zero = lower;
                }

        double
        get_time() const
        {
            return this->current_time;
        }

        std::string
        print_time() const
        {
            char print_t[10];
            std::sprintf(print_t, "%.3f", this->current_time);

            for (unsigned int i = 0; i < 10; ++i)
                if (print_t[i] == '.')
                    print_t[i] = '_';

            return std::string(print_t);
        }

        double
        set_current_time(const double new_time)
        {
            return current_time = new_time;
        }

        double
        increment_time(const double delta_time)
        {
            return current_time += delta_time;
        }

        virtual double
        value(
            const Point<DIM> &p,
            const unsigned int component) const override
            {
                (void) component;
#if PLANE 
                double temp = 0.;
                for (unsigned int i = 0; i < DIM; ++i)
                    temp += this->systemvals->wavenumber(i) * p(i);
                temp -= this->systemvals->frequency * this->current_time;
                temp *= M_PI;
//                return std::sin(temp);
                return std::exp( - temp * temp);
#else
                double radius = p.norm();
                double temp = std::cyl_bessel_j(0, radius * this->bessel_zero);
                return temp * std::cos(this->current_time * this->bessel_zero);
#endif
            }

        double
        time_derivative_value(
            const Point<DIM> &p) const
            {
#if PLANE
                double temp = 0.;
                for (unsigned int i = 0; i < DIM; ++i)
                    temp += this->systemvals->wavenumber(i) * p(i);
                temp -= this->systemvals->frequency * this->current_time;
                temp *= M_PI;
                double dtemp_dt = - this->systemvals->frequency * M_PI;
//                return dtemp_dt * std::cos(temp);
                return dtemp_dt * -2 * temp * std::exp( - temp * temp);
#else
                double radius = p.norm();
                double temp = std::cyl_bessel_j(0, radius * this->bessel_zero);
                return - temp * this->bessel_zero * std::sin(this->current_time * this->bessel_zero);
#endif
            }

    private:
        const WaveSystemValues<DIM>* systemvals;
        double current_time;
        double bessel_zero;
    }; //Solution

    class WaveSolver
    {
        public:
        WaveSolver(
            const WaveSystemValues<DIM>& vals = WaveSystemValues<DIM>(),
            const unsigned int fe_deg = 1) :
            wave_values(vals),
            fe_degree(fe_deg),
            levelset_mapping(),
            levelset_fe(levelset_degree),
            levelset_quadrature(levelset_degree + 1),
            levelset_function(1.),
            levelset_dof_handler(grid),
            solution_face_quadrature(fe_degree + 1),
            solution_dof_handler(grid),
            mesh_classifier(levelset_dof_handler, levelset),
            analytical_solution(&wave_values)
            {};
        void run(
            const double start_time = 0.,
            const double end_time = 5.,
            const unsigned int refinement = 0, 
            const bool plot_solution = false);

        private:
        void make_grid(const unsigned int refine_level = 0);
        void setup_levelset();
        void distribute_dofs();
        void initialise_quadratures();

        void initialise_matrices();
        void assembly_preprocessing();
        void assemble_system();
        void prepare_system_matrices();

        void make_boundary_term();
        void make_initial_conditions();
        double compute_l2_error();

        void make_initial_timestep();
        void make_timestep();
        void solve_in_time(
            const double end_time,
            std::vector<double> &l2_error_timeseries,
            const bool plot_sol = false,
            const std::string &run_name = ""
        );
        void print_errors(
            const std::vector<double>& errs,
            const unsigned int degree,
            const bool to_file = true) const;
        void plot_solution(const std::string &run_name);

        void clean_memory();

        const WaveSystemValues<DIM> wave_values;

        const unsigned int levelset_degree = 2;
        unsigned int fe_degree;

        Triangulation<DIM> grid;
        double dx;
        double dt;

        MappingCartesian<DIM> levelset_mapping;
        FE_Q<DIM> levelset_fe;
        QGauss<DIM> levelset_quadrature;
        LevelSet levelset_function;
        DoFHandler<DIM> levelset_dof_handler;
        Vector<double> levelset;

        hp::MappingCollection<DIM> solution_mapping_collection;
        hp::FECollection<DIM> solution_fe_collection;
        hp::QCollection<DIM> solution_quadrature;
        QGauss<DIM-1> solution_face_quadrature;
        hp::QCollection<1> solution_quadrature_1d;
        DoFHandler<DIM> solution_dof_handler;
        Vector<double> solution;
        Vector<double> initial_tdev_and_tstepping_temp_vector;
        Vector<double> temp_2;
        Vector<double> previous_timestep;
        Vector<double> next_timestep;

        NonMatching::MeshClassifier<DIM> mesh_classifier;
        NonMatching::RegionUpdateFlags update_flags;
        NonMatching::FEValues<DIM>* nonmatching_fe_values;

        SparsityPattern sparsity_pattern;
        SparseMatrix<double> mass_matrix;
        SparseMatrix<double> stiffness_matrix;
        SparseMatrix<double> stabilisation_term;
        SparseMatrix<double> tstep_noninvert;
        SparseDirectUMFPACK mass_inv;
        Vector<double> boundaries;
        double nitsche_weight_dx;

        Solution analytical_solution;

        DataOut<DIM> data_out;
    }; //WaveSolver

    void 
    WaveSolver::make_grid(const unsigned int refine_level)
    {
        std::cout << "Creating background mesh" << std::endl;

        double width = 1.21 * this->levelset_function.radius;
        GridGenerator::hyper_cube(this->grid, -width, width);
        this->grid.refine_global(refine_level + 3);

        this->solution_mapping_collection.push_back(MappingCartesian<DIM>());
    }

    void 
    WaveSolver::setup_levelset()
    {
        std::cout << "Setting up discrete level set function" << std::endl;

        this->levelset_dof_handler.distribute_dofs(this->levelset_fe);
        this->levelset.reinit(this->levelset_dof_handler.n_dofs());

        VectorTools::interpolate(
            this->levelset_mapping,
            this->levelset_dof_handler,
            this->levelset_function,
            this->levelset);

        this->mesh_classifier.reclassify();
    }

    enum
    ActiveFEIndex
    {
        active = 0,
        nothing = 1
    };

    void 
    WaveSolver::distribute_dofs()
    {
        std::cout << "Distributing degrees of freedom" << std::endl;

        this->solution_fe_collection.push_back(FE_Hermite<DIM>(this->fe_degree));
        this->solution_fe_collection.push_back(FE_Nothing<DIM>());

        unsigned int cell_count = 0, active_cell_count = 0;

        this->dx = this->grid.begin_active()->minimum_vertex_distance();

        for (const auto &cell : this->solution_dof_handler.active_cell_iterators())
        {
            ++cell_count;
            const NonMatching::LocationToLevelSet cell_location =
                this->mesh_classifier.location_to_level_set(cell);

            if (cell_location == NonMatching::LocationToLevelSet::outside)
                cell->set_active_fe_index(ActiveFEIndex::nothing);
            else
            {
                cell->set_active_fe_index(ActiveFEIndex::active);
                ++active_cell_count;
            }
        }

        std::cout << "Total number of elements: " << cell_count << std::endl;
        std::cout << "Total active elements: " << active_cell_count << std::endl;
        this->solution_dof_handler.distribute_dofs(this->solution_fe_collection);
    }

    void
    WaveSolver::initialise_quadratures()
    {
        this->solution_quadrature.push_back(QGauss<DIM>(this->fe_degree + 1));
        this->solution_face_quadrature = QGauss<DIM-1>(this->fe_degree + 1);
        this->solution_quadrature_1d.push_back(QGauss<1>(this->fe_degree + 1));
    }

    bool 
    face_has_penalties(
            const NonMatching::MeshClassifier<DIM> &mesh_classifier,
            const Triangulation<DIM>::active_cell_iterator &cell,
            const unsigned int face_index)
    {
        if (cell->at_boundary(face_index))
            return false;

        const NonMatching::LocationToLevelSet 
        this_location = mesh_classifier.location_to_level_set(cell),
        that_location = mesh_classifier.location_to_level_set(cell->neighbor(face_index));

        Assert(
            (this_location != NonMatching::LocationToLevelSet::unassigned ||
                that_location != NonMatching::LocationToLevelSet::unassigned),
            ExcMessage("Error: There is a mesh cell with an unassigned location relative"
                        " to the level-set. Check that all cells are correctly assigned."));

        if (this_location == NonMatching::LocationToLevelSet::inside &&
            that_location == NonMatching::LocationToLevelSet::inside)
            return false;
        
        if (this_location == NonMatching::LocationToLevelSet::outside ||
            that_location == NonMatching::LocationToLevelSet::outside)
            return false;

        return true;
    }

    void 
    WaveSolver::initialise_matrices()
    {
        std::cout << "Initialising matrix structures" << std::endl;

        const auto
        face_has_penalty = [&](const auto &cell, const unsigned int face)
        {
            return face_has_penalties(this->mesh_classifier, cell, face);
        };

        const unsigned int n_dofs_total = this->solution_dof_handler.n_dofs();
        std::cout << "Total DoFs in system: " << n_dofs_total << std::endl;
        DynamicSparsityPattern matrix_shape(n_dofs_total, n_dofs_total);
        this->solution.reinit(n_dofs_total);
        this->initial_tdev_and_tstepping_temp_vector.reinit(n_dofs_total);
        this->temp_2.reinit(n_dofs_total);
        this->previous_timestep.reinit(n_dofs_total);
        this->next_timestep.reinit(n_dofs_total);
        this->boundaries.reinit(n_dofs_total);

        const AffineConstraints<double> constraints;
        const bool                      keep_constrained_dofs = true;

        const unsigned int n_components = this->solution_fe_collection.n_components();

        Table<2, DoFTools::Coupling> cell_coupling(n_components, n_components);
        Table<2, DoFTools::Coupling> face_coupling(n_components, n_components);
        cell_coupling[0][0] = DoFTools::always;
        face_coupling[0][0] = DoFTools::always;

        DoFTools::make_flux_sparsity_pattern(
            this->solution_dof_handler,
            matrix_shape,
            constraints,
            keep_constrained_dofs,
            cell_coupling,
            face_coupling,
            numbers::invalid_subdomain_id,
            face_has_penalty);

        this->sparsity_pattern.copy_from(matrix_shape);

        this->mass_matrix.reinit(sparsity_pattern);
        this->stiffness_matrix.reinit(sparsity_pattern);
        this->stabilisation_term.reinit(sparsity_pattern);
        this->tstep_noninvert.reinit(sparsity_pattern);
    }

    void
    obtain_normalised_face_derivatives(
        const Mapping<DIM>& mapping,
        const DoFHandler<DIM>& dof_handler,
        const Point<DIM>& cell_side_lengths,
        Table<3, double>& normalised_face_derivatives)
    {
        Assert(
            (dynamic_cast<const MappingCartesian<DIM>*>(&mapping) != nullptr),
            ExcMessage("ERROR: Current version of code only works with rectangular cells (MappingCartesian)")
        );

        // Extract information on polynomial basis from DoFHandler object
        Assert(
            (dynamic_cast<const FE_Poly<DIM>*>(&(dof_handler.get_fe())) != nullptr),
            ExcMessage("ERROR: Current version of code requires an FE class with FE_Poly as a base class.")
        );
        const FE_Poly<DIM>& fe_object = dynamic_cast<const FE_Poly<DIM>&>(dof_handler.get_fe());
        const unsigned int poly_degree = fe_object.get_degree();

        Assert(
            (dynamic_cast<const TensorProductPolynomials<DIM>*>(&(fe_object.get_poly_space())) != nullptr),
            ExcMessage("ERROR: Current version of code requires an underlying polynomial space constructed from tensor product polynomials.")
        );
        const std::vector<Polynomials::Polynomial<double>> poly_vector =
            dynamic_cast<const TensorProductPolynomials<DIM>&>(fe_object.get_poly_space()).get_underlying_polynomials();
        const unsigned int n_basis_polys = static_cast<unsigned int>(poly_vector.size());
        AssertDimension(n_basis_polys, poly_degree + 1);

        // Check whether basis is Hermite (ie has special rescaling) or not
        unsigned int herm_regularity = 0;
        const FE_Hermite<DIM>* fe_herm_ptr = dynamic_cast<const FE_Hermite<DIM>*>(&fe_object);
        if (fe_herm_ptr != nullptr)
            herm_regularity = fe_herm_ptr->get_regularity;

        // Set up arrays for storing derivative information
        std::vector<double> derivative_values_at_0(n_basis_polys);
        std::vector<double> derivative_values_at_1(n_basis_polys);

        normalised_face_derivatives.resize(n_basis_polys, n_basis_polys, 2 * DIM);

        // Calculate derivatives of the 1D polynomials running normal to each face
        Point<DIM> h_factors, h_factors_out;
        for (unsigned int d = 0; d < DIM; ++d)
            h_factors(d) = h_factors_out(d) = 1.;
        unsigned int poly_count = 0;

        for (const auto &p : poly_vector)
        {
            h_factors = h_factors_out;

            p.value(0., derivative_values_at_0);
            p.value(1., derivative_values_at_1);

            for (unsigned int i = 0; i <= poly_degree; ++i)
                for (unsigned int d = 0; d < DIM; ++d)
                {
                    normalised_face_derivatives[poly_count][i][2*d] = h_factors(d) * derivative_values_at_0[i];
                    normalised_face_derivatives[poly_count][i][2*d+1] = h_factors(d) * derivative_values_at_1[i];

                    h_factors(d) /= cell_side_lengths(d);
                }
            
            ++poly_count;

            for (unsigned int d = 0; d < DIM; ++d)
                h_factors_out(d) =
                    ((herm_regularity == 0) || (poly_count == herm_regularity + 1)) ?
                    1. :
                    h_factors_out(d) * cell_side_length(d);
        }
    }

    void
    create_stabilisation_matrix(
        const Mapping<DIM>& mapping,
        const DoFHandler<DIM>& dof_handler,
        const Quadrature<DIM-1>& face_quadrature,
        const NonMatching::MeshClassifier<DIM>& mesh_classifier,
        const std::vector<double>& weights,
        const Point<DIM>& cell_side_lengths,
        SparseMatrix<double>& stabilisation_matrix)
    {
        Assert(
            (dynamic_cast<const MappingCartesian<DIM> *>(&mapping) != nullptr),
            ExcMessage("Warning: current code for calculating ghost penalty term only valid "
                       "for MappingCartesian. The provided mapping is not compatible.")
                       );

        const unsigned int poly_degree = dynamic_cast<const FE_Poly<DIM>&>(dof_handler.get_fe())
                                            .get_degree();
        AssertDimension(weights.size(), poly_degree + 1);

        // Obtain a look-up table of normalised derivative values normal to each face
        Table<3, double> poly_basis_derivatives;
        obtain_normalised_face_derivatives(mapping, dof_handler, cell_side_lengths, poly_basis_derivatives);
        poly_basis_derivatives.size(0);

        // Create arrays to store DoF information, pre-calculating where possible
        const unsigned int n_dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

        std::vector<unsigned int> cell_to_interface(n_dofs_per_cell),
                                neighbor_to_interface(n_dofs_per_cell);
        std::vector<unsigned int> interface_dof_locations;

        const std::vector<unsigned int> fe_lexicographic_to_hierarchic =
                dynamic_cast<const FE_Hermite<DIM>&>(dof_handler.get_fe())
                    .get_lexicographic_to_hierarchic_numbering();

        Table<2, unsigned int> dofs_of_nonzero_values_on_faces =
            dynamic_cast<const FE_Hermite<DIM> &>(dof_handler.get_fe())
                .get_dofs_corresponding_to_outward_normal_derivatives(0);

        // Iterate over faces 
        FEInterfaceValues<DIM> fe_stabilisation_int_vals(
            mapping,
            dof_handler.get_fe(),
            face_quadrature,
            update_values | update_JxW_values | update_quadrature_points | update_normal_vectors | update_gradients);

        for (auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::active))
            for (unsigned int f : cell->face_indices())
                if (face_has_penalties(mesh_classifier, cell, f))
                {
                    const unsigned int invalid_subface =
                        numbers::invalid_unsigned_int;

                    fe_stabilisation_int_vals.reinit(
                        cell,
                        f,
                        invalid_subface,
                        cell->neighbor(f),
                        cell->neighbor_of_neighbor(f),
                        invalid_subface);

                    const unsigned int n_interface_dofs = fe_stabilisation_int_vals.n_current_interface_dofs();
                    interface_dof_locations.resize(n_interface_dofs);

                    // Calculate mapping from cell and neighbour to interface
                    std::fill(cell_to_interface.begin(), cell_to_interface.end(), numbers::invalid_unsigned_int);
                    std::fill(neighbor_to_interface.begin(), neighbor_to_interface.end(), numbers::invalid_unsigned_int);
                    std::fill(interface_dof_locations.begin(), interface_dof_locations.end(), 0);

                    for (unsigned int i = 0; i < n_interface_dofs; ++i)
                    {
                        std::array<unsigned int, 2> temp_map =
                            fe_stabilisation_int_vals.interface_dof_to_dof_indices(i);

                        bool on_cell = (temp_map[0] != numbers::invalid_unsigned_int),
                            on_neighbor = (temp_map[1] != numbers::invalid_unsigned_int);

                        if (on_cell)
                            cell_to_interface[temp_map[0]] = i;
                        if (on_neighbor)
                            neighbor_to_interface[temp_map[1]] = i;

                        interface_dof_locations[i] = (on_cell ? 1 : 0) + (on_neighbor ? 2 : 0);
                    }

                    Assert(
                        (std::find(
                            cell_to_interface.begin(), 
                            cell_to_interface.end(), 
                            numbers::invalid_unsigned_int) == cell_to_interface.end()),
                        ExcMessage("Error: failed to generate complete cell to interface mapping.")
                        );
                    Assert(
                        (std::find(
                            neighbor_to_interface.begin(), 
                            neighbor_to_interface.end(), 
                            numbers::invalid_unsigned_int) == neighbor_to_interface.end()),
                        ExcMessage("Error: failed to generate complete neighbor to interface mapping.")
                        );
                    Assert(
                        (std::find(
                            interface_dof_locations.begin(), 
                            interface_dof_locations.end(), 
                            0) == interface_dof_locations.end()),
                        ExcMessage("Error: failed to find locations of all DoFs on interface.")
                        );

                    // For FE_Hermite, the value of the polynomial component normal to the interface of a basis
                    // function with non-zero values on the interface will be 1 at the interface. Therefore
                    // for a Cartesian grid, the value of a normal derivative to the interface can
                    // be found by multiplying the correct 1D derivative with the values of the non-zero shape
                    // functions on that face.
                    // Note: this is effectively a sleight of hand needed to get normal derivatives of 4th order
                    // or higher, which are necessary to stailise p>=4 polynomial bases in cutFEM
                    const unsigned int normal_derivative_axis = static_cast<unsigned int>(f / 2);
                    const unsigned int cell_direction = f - 2 * derivative_axis;
                    const unsigned int neighbor_direction = 1 - cell_direction;

                    const unsigned int no_consecutive_nonzero_dofs = Utilities::pow(poly_degree + 1, derivative_axis);

                    // Iterate over quadrature points
                    FullMatrix<double> local_stabilization(
                        n_interface_dofs,
                        n_interface_dofs);

                    for (
                        unsigned int q = 0;
                        q < fe_stabilisation_int_vals.n_quadrature_points;
                        ++q)
                    {
                        // Pre-calculate normal derivatives at the given quadrature point
                        Table<3, double> normal_derivatives_on_selected_faces(n_dofs_per_cell, poly_degree + 1, 2);

                        // Loop over indices in lexicographic ordering, and convert to hierarchic when needed
                        for (unsigned int current_cell_dof = 0; current_cell_dof < n_dofs_per_cell; ++current_cell_dof)
                        {
                            const unsigned int hierarchic_index = fe_lexicographic_to_hierarchic[current_cell_dof];

                            // Calculate index of differentiated polynomial (same at both faces)
                                // The poly_index will increment each time a batch of consecutive non-zero DoFs
                                // is completed, and reset once it goes through all polynomials in the basis
                            const unsigned int internal_batch_index = current_cell_dof % no_consecutive_nonzero_dofs;
                            const unsigned int intermediate_p_index = (current_cell_dof - internal_batch_index) / no_consecutive_nonzero_dofs;
                            const unsigned int poly_index = intermediate_p_index % (poly_degree + 1);

                            // Calculate indices of nonzero value DoFs (separated by offset at each face)
                            const unsigned int external_batch_index = (intermediate_p_index - poly_index) / (poly_degree + 1);
                            const unsigned int hierarchic_shape_index = 
                                nonzero_face_dofs[2 * derivative_axis][internal_batch_index + no_consecutive_nonzero_dofs * external_batch_index];

                            // Obtain value of shape function at given quadrature point
                                // Both nonzero DoFs will take the same value, since the perpendicular
                                // component will be 1, and the parallel component will use the same polys
                                // => only need to calculate this value once on one face
                            const double current_shape_value = fe_stabilisation_int_vals.shape_value(
                                true, 
                                cell_to_interface[hierarchic_shape_index], 
                                q);

                            for (unsigned int i = 0; i <= poly_degree; ++i)
                            {
                                normal_derivatives_on_selected_faces[hierarchic_index][i][0] = 
                                    current_shape_value * poly_basis_derivatives[poly_index][i][2*normal_derivative_axis];
                                normal_derivatives_on_selected_faces[hierarchic_index][i][1] =
                                    current_shape_value * poly_basis_derivatives[poly_index][i][2*normal_derivative_axis+1];
                            }
                        }

                        // Obtain jump values of derivatives across interface
                        Table<2, double> normal_derivative_jump_terms(n_interface_dofs, poly_degree + 1);

                        for (unsigned int i = 0; i < n_interface_dofs; ++i)
                            for (unsigned int p = 0; p <= poly_degree; ++p)
                            {
                                normal_derivative_jump_terms[i][p] = 0;
                                if (interface_dof_locations[i] % 2 == 1)
                                {
                                    // non-zero on cell
                                    unsigned int cell_location = 0;
                                    for (const unsigned int j : cell_to_interface)
                                    {
                                        if (j == i) break;
                                        ++cell_location;
                                    }

                                    double normal_value =
                                        normal_derivatives_on_selected_faces[cell_location][p][cell_direction];
                                    normal_derivative_jump_terms[i][p] += (cell_direction == 0) ? normal_value : - normal_value;
                                }
                                if (interface_dof_locations[i] > 1)
                                {
                                    // non-zero on neighbor
                                    unsigned int neighbor_location = 0;
                                    for (const unsigned int j : neighbor_to_interface)
                                    {
                                        if (j == i) break;
                                        ++neighbor_location;
                                    }

                                    double normal_value =
                                        normal_derivatives_on_selected_faces[neighbor_location][p][neighbor_direction];
                                    normal_derivative_jump_terms[i][p] += (neighbor_direction == 0) ? normal_value : - normal_value;
                                }
                            }

                        for (unsigned int i = 0; i < n_interface_dofs; ++i)
                            for (unsigned int j = 0; j < n_interface_dofs; ++j)
                            {
                                double len_factor = cell_side_length;

                                for (unsigned int p = 0; p <= poly_degree; ++p)
                                {
                                    local_stabilization(i, j) += 
                                        len_factor * weights[p] *
                                        normal_derivative_jump_terms[i][p] *
                                        normal_derivative_jump_terms[j][p] *
                                        fe_stabilisation_int_vals.JxW(q);
                                    len_factor *= cell_side_length * cell_side_length;
                                }
                            }
                    }

                    const std::vector<types::global_dof_index>
                    local_interface_dof_indices =
                        fe_stabilisation_int_vals.get_interface_dof_indices();

                    stabilisation_matrix.add(
                        local_interface_dof_indices, 
                        local_stabilization);
                }
    }

    static inline double
    get_nitsche_parameter(const unsigned int fe_degree)
    {
        return static_cast<double>(5 * fe_degree * (fe_degree + 1));
    }

    static inline double
    get_cfl(const double wavespeed, const unsigned int fe_degree)
    {
        return 0.75 / (std::sqrt(3 * DIM) * (fe_degree + 1) * wavespeed);
    }

    void
    WaveSolver::assembly_preprocessing()
    {
        this->update_flags.inside =
            update_values | update_gradients | update_JxW_values | update_quadrature_points;
        this->update_flags.surface = this->update_flags.inside | update_normal_vectors;

        this->nonmatching_fe_values =
            new NonMatching::FEValues<DIM>(
                this->solution_mapping_collection,
                this->solution_fe_collection,
                this->solution_quadrature,
                this->solution_quadrature_1d,
                this->update_flags,
                this->mesh_classifier,
                this->levelset_dof_handler,
                this->levelset);

        this->nitsche_weight_dx = get_nitsche_parameter(this->fe_degree) / this->dx;
        this->dt = get_cfl(this->wave_values.wavespeed, this->fe_degree) * this->dx;
    }

    void 
    WaveSolver::assemble_system()
    {
        std::cout << "Assembling matrices" << std::endl;

        const unsigned int n_dofs_per_cell = this->solution_fe_collection[0].dofs_per_cell;

        std::cout << "DoFs per cell: " << n_dofs_per_cell << std::endl;
        std::cout << "DPO vector for FE object: ";
        std::cout << this->solution_fe_collection[0].n_dofs_per_vertex() << " "
                  << this->solution_fe_collection[0].n_dofs_per_line() << " "
                  << this->solution_fe_collection[0].n_dofs_per_quad();
        std::cout << std::endl;

        FullMatrix<double> local_mass(n_dofs_per_cell, n_dofs_per_cell);
        FullMatrix<double> local_stiffness(n_dofs_per_cell, n_dofs_per_cell);
        Vector<double> local_initial_conds(n_dofs_per_cell);
        Vector<double> local_initial_tdevs(n_dofs_per_cell);
        std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

        for (const auto &cell :
         this->solution_dof_handler.active_cell_iterators() |
           IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::active))
        {
            local_mass = 0;
            local_stiffness = 0;
            local_initial_conds = 0;
            local_initial_tdevs = 0;

            this->nonmatching_fe_values->reinit(cell);

            const /*std_cxx17*/std::optional<FEValues<DIM>>& inside_fe_values =
                this->nonmatching_fe_values->get_inside_fe_values();
            if (inside_fe_values)
                for (const unsigned int q : inside_fe_values->quadrature_point_indices())
                {
                    const Point<DIM> & current_point = inside_fe_values->quadrature_point(q);

                    for (const unsigned int i : inside_fe_values->dof_indices())
                    {
                        for (const unsigned int j : inside_fe_values->dof_indices())
                        {
                            local_mass(i,j) += inside_fe_values->JxW(q) *
                                inside_fe_values->shape_value(i,q) *
                                inside_fe_values->shape_value(j,q);
                            local_stiffness(i,j) += inside_fe_values->JxW(q) *
                                inside_fe_values->shape_grad(i,q) *
                                inside_fe_values->shape_grad(j,q);
                        }

                        local_initial_conds(i) += inside_fe_values->JxW(q) *
                            inside_fe_values->shape_value(i,q) *
                            this->analytical_solution.value(current_point, 0);
                        local_initial_tdevs(i) += inside_fe_values->JxW(q) *
                            inside_fe_values->shape_value(i,q) *
                            this->analytical_solution.time_derivative_value(current_point);
                    }
                }

            const /*std_cxx17*/std::optional<NonMatching::FEImmersedSurfaceValues<DIM>> &
                surface_fe_values = this->nonmatching_fe_values->get_surface_fe_values();
            if (surface_fe_values)
                for (const unsigned int q : surface_fe_values->quadrature_point_indices())
                {
                    const Tensor<1, DIM>& normal = surface_fe_values->normal_vector(q);

                    for (const unsigned int i : surface_fe_values->dof_indices())
                        for (const unsigned int j : surface_fe_values->dof_indices())
                        {
                            local_stiffness(i,j) += surface_fe_values->JxW(q) *
                                ( this->nitsche_weight_dx * surface_fe_values->shape_value(i,q) *
                                    surface_fe_values->shape_value(j,q) -
                                 normal * surface_fe_values->shape_grad(i,q) *
                                    surface_fe_values->shape_value(j,q) -
                                 normal * surface_fe_values->shape_grad(j,q) *
                                    surface_fe_values->shape_value(i,q) );
                        }
                }

            cell->get_dof_indices(local_dof_indices);
            this->mass_matrix.add(local_dof_indices, local_mass);
            this->stiffness_matrix.add(local_dof_indices, local_stiffness);
            this->solution.add(local_dof_indices, local_initial_conds);
            this->initial_tdev_and_tstepping_temp_vector.add(local_dof_indices, local_initial_tdevs);
        }

        const double stabfactor = 0.0625;
        std::vector<double> stab_weights(fe_degree + 1, 1.0);
        for (unsigned int i = 0; i < fe_degree; ++i)
            stab_weights[i+1] = stabfactor * stab_weights[i];

        create_stabilisation_matrix(
            this->solution_mapping_collection[0],
            this->solution_dof_handler,
            this->solution_face_quadrature,
            this->mesh_classifier,
            stab_weights,
            this->dx,
            this->stabilisation_term);
    }

    void
    WaveSolver::prepare_system_matrices()
    {
        const double mass_ghost_parameter = 1.;
        const double stiffness_ghost_parameter = 0.25 / (this->dx * this->dx);

        std::cout << "Mass matrix Linfty norm: " << this->mass_matrix.linfty_norm() << std::endl;
        std::cout << "Stiffness matrix Linfty norm: " << this->stiffness_matrix.linfty_norm() << std::endl;
        std::cout << "Stabilisation Linfty norm: " << this->stabilisation_term.linfty_norm() << std::endl;

        this->mass_matrix.add( mass_ghost_parameter, this->stabilisation_term );
        this->stiffness_matrix.add( stiffness_ghost_parameter, this->stabilisation_term);

        std::cout << "Stabilised mass Linfty norm: "<< this->mass_matrix.linfty_norm() << std::endl;
        std::cout << "Stabilised stiffness Linfty norm: " << this->stiffness_matrix.linfty_norm() << std::endl;

        this->tstep_noninvert.copy_from(this->mass_matrix);
        this->tstep_noninvert.add( -0.5 * this->dt * this->dt, this->stiffness_matrix);

        this->mass_inv.initialize(this->mass_matrix);
        std::cout << "m=" << this->mass_inv.m() << " n=" << this->mass_inv.n() << std::endl;
    }

    void
    WaveSolver::make_boundary_term()
    {
        this->boundaries = 0;

        const unsigned int n_dofs_per_cell = this->solution_fe_collection[0].dofs_per_cell;
        Vector<double> local_boundary_contribution(n_dofs_per_cell);
        std::vector<types::global_dof_index> local_to_global(n_dofs_per_cell);

        for (const auto &cell : this->solution_dof_handler.active_cell_iterators() |
                                    IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::active))
        {
            local_boundary_contribution = 0;
            this->nonmatching_fe_values->reinit(cell);
            cell->get_dof_indices(local_to_global);

            const /*std_cxx17*/std::optional<NonMatching::FEImmersedSurfaceValues<DIM>>
                &boundary_fe_values = nonmatching_fe_values->get_surface_fe_values();

            if (boundary_fe_values)
                for (const unsigned int q : boundary_fe_values->quadrature_point_indices())
                {
                    const Tensor<1, DIM> &normal = boundary_fe_values->normal_vector(q);
                    const double boundary_val = this->analytical_solution.value(
                        boundary_fe_values->quadrature_point(q), 0) *
                        boundary_fe_values->JxW(q);

                    for (const unsigned int i : boundary_fe_values->dof_indices())
                    {
                        local_boundary_contribution(i) +=
                            boundary_val * (this->nitsche_weight_dx * boundary_fe_values->shape_value(i,q) -
                                            normal * boundary_fe_values->shape_grad(i,q));
                    }
                }

            this->boundaries.add(local_to_global, local_boundary_contribution);
        }
    }

    void 
    WaveSolver::make_initial_conditions()
    {
        this->next_timestep = this->solution;
        this->solution = 0;
        this->mass_inv.vmult(this->solution, this->next_timestep);

        this->next_timestep = 0;
        
        // initial_tdev_and_tstepping_vector contains M u_t, which we
        // use for an initial timestep so does not need to be changed

        this->make_boundary_term();
    }

    double 
    WaveSolver::compute_l2_error()
    {
        std::vector<types::global_dof_index> local_to_global(this->solution_fe_collection[0].dofs_per_cell);
        double error_total_sq = 0., local_error = 0.;

        for (const auto &cell : this->solution_dof_handler.active_cell_iterators() |
                                    IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::active))
        {
            this->nonmatching_fe_values->reinit(cell);
            cell->get_dof_indices(local_to_global);

            const /*std_cxx17*/std::optional<FEValues<DIM>>& inside_fe_values =
                this->nonmatching_fe_values->get_inside_fe_values();

            if (inside_fe_values)
                for (const unsigned int q : inside_fe_values->quadrature_point_indices())
                {
                    local_error = this->analytical_solution.value(inside_fe_values->quadrature_point(q), 0);

                    for (const unsigned int i : inside_fe_values->dof_indices())
                        local_error -= this->solution(local_to_global[i]) 
                            * inside_fe_values->shape_value(i,q);
                    
                    error_total_sq += local_error * local_error * inside_fe_values->JxW(q);
                }
        }

        return std::sqrt(error_total_sq);
    }

    void 
    WaveSolver::make_initial_timestep()
    {
        // Start by calculating boundary contributions, which appear as a forcing term
        this->next_timestep = this->boundaries;

        this->analytical_solution.increment_time(this->dt);
        this->make_boundary_term();

        //this->next_timestep += this->boundaries;
        //this->next_timestep *= 0.5;

        this->next_timestep *= 0.5 * this->dt * this->dt;

        // Now add the initial time derivative
        this->initial_tdev_and_tstepping_temp_vector *= this->dt;
        this->next_timestep += this->initial_tdev_and_tstepping_temp_vector;
        this->initial_tdev_and_tstepping_temp_vector = 0;

        // Add the initial values and contributions from the Laplacian
        this->tstep_noninvert.vmult(this->initial_tdev_and_tstepping_temp_vector, this->solution);
        this->next_timestep += this->initial_tdev_and_tstepping_temp_vector;

#if PLANE
#else
        this->stiffness_matrix.vmult(this->initial_tdev_and_tstepping_temp_vector, this->solution);
        this->initial_tdev_and_tstepping_temp_vector *= - 0.5 * this->dt * this->dt;

        for (unsigned int i = 2; i < this->fe_degree; i += 2)
        {
            this->mass_inv.vmult(this->temp_2, this->initial_tdev_and_tstepping_temp_vector);
            this->stiffness_matrix.vmult(this->initial_tdev_and_tstepping_temp_vector, this->temp_2);
            this->initial_tdev_and_tstepping_temp_vector *= - this->dt * this->dt / ((i + 1) * (i + 2));
            this->next_timestep += this->initial_tdev_and_tstepping_temp_vector;
        }
#endif

        this->previous_timestep = this->solution;
        this->solution = 0;
        this->mass_inv.vmult(this->solution, this->next_timestep);
    }

    void 
    WaveSolver::make_timestep()
    {
        this->next_timestep = this->boundaries;

        this->analytical_solution.increment_time(this->dt);
        this->make_boundary_term();

        //this->next_timestep += this->boundaries;
        //this->next_timestep *= 0.5;

        this->next_timestep *= this->dt * this->dt;

        this->initial_tdev_and_tstepping_temp_vector = 0;
        this->tstep_noninvert.vmult(this->initial_tdev_and_tstepping_temp_vector, this->solution);
        this->initial_tdev_and_tstepping_temp_vector *= 2.;

        this->next_timestep += this->initial_tdev_and_tstepping_temp_vector;
        this->initial_tdev_and_tstepping_temp_vector = 0;   //This was added for clarity, not because it could be arbitrarily increasing order of accuracy

#if PLANE
#else
        this->stiffness_matrix.vmult(this->initial_tdev_and_tstepping_temp_vector, this->solution);
        this->initial_tdev_and_tstepping_temp_vector *= - this->dt * this->dt;

        for (unsigned int i = 2; i < this->fe_degree; i += 2)
        {
            this->mass_inv.vmult(this->temp_2, this->initial_tdev_and_tstepping_temp_vector);
            this->stiffness_matrix.vmult(this->initial_tdev_and_tstepping_temp_vector, this->temp_2);
            this->initial_tdev_and_tstepping_temp_vector *= - this->dt * this->dt / ((i + 1) * (i + 2));
            this->next_timestep += this->initial_tdev_and_tstepping_temp_vector;
        }
#endif

        this->initial_tdev_and_tstepping_temp_vector = this->solution;
        this->solution = 0;
        this->mass_inv.vmult(this->solution, this->next_timestep);
        this->solution -= this->previous_timestep;

        this->previous_timestep = this->initial_tdev_and_tstepping_temp_vector;
        this->initial_tdev_and_tstepping_temp_vector = 0;
    }

    void WaveSolver::solve_in_time(
        const double end_time,
        std::vector<double> &l2_error_timeseries,
        const bool plot_sol,
        const std::string &run_name)
    {
        unsigned int n_timesteps = std::ceil((end_time - this->analytical_solution.get_time()) / this->dt);
        l2_error_timeseries.reserve(n_timesteps + 1);

        // initial conditions
        make_initial_conditions();
        l2_error_timeseries.push_back(compute_l2_error());
        if (plot_sol)
            this->data_out.add_data_vector(
                this->solution_dof_handler,
                this->solution, 
                run_name + std::string("_t_") + this->analytical_solution.print_time()
                );

        // initial timestep
        make_initial_timestep();
        l2_error_timeseries.push_back(compute_l2_error());
        if (plot_sol)
            this->data_out.add_data_vector(
                this->solution_dof_handler,
                this->solution, 
                run_name + std::string("_t_") + this->analytical_solution.print_time()
                );

        // timestepping solver
        for (unsigned int i = 1; i < n_timesteps; ++i)
        {
            make_timestep();
            l2_error_timeseries.push_back(compute_l2_error());
            if (plot_sol)
                this->data_out.add_data_vector(
                    this->solution_dof_handler,
                    this->solution, 
                    run_name + std::string("_t_") + this->analytical_solution.print_time()
                    );
        }
    }

    void WaveSolver::print_errors(
            const std::vector<double>& errs,
            const unsigned int degree,
            const bool to_file) const
            {
                auto print_run = [&](std::ostream& output_stream) -> std::ostream&
                {
                    for (const double l2_err : errs)
                        output_stream << l2_err << std::endl;
                    return output_stream;
                };
                if (to_file)
                {
                    std::ostringstream filename_errors;
                    filename_errors << std::setprecision(2);
                    filename_errors << "error_series_degree_" << degree;
                    filename_errors << "_gridsize_" << this->dx;
                    filename_errors << ".txt";
                    std::ofstream error_series(filename_errors.str());
                    print_run(error_series);
                    error_series.close();
                }
                else
                    print_run(std::cout);
            }

    void 
    WaveSolver::plot_solution(const std::string &run_name)
    {
        this->data_out.set_cell_selection(
            [this](const typename Triangulation<DIM>::cell_iterator &cell) 
                {
                    return cell->is_active() &&
                        this->mesh_classifier.location_to_level_set(cell) !=
                        NonMatching::LocationToLevelSet::outside;
                }
            );
        
        this->data_out.build_patches(
            this->solution_mapping_collection[0],
            this->fe_degree,
            DataOut<DIM>::CurvedCellRegion::curved_inner_cells);

        std::ofstream solution_plots(
            std::string("wave_equation_") + run_name + std::string(".vtu"));
        this->data_out.write_vtu(solution_plots);
        solution_plots.close();
    }

    void
    WaveSolver::clean_memory()
    {
        this->grid.clear();

        this->solution_mapping_collection = hp::MappingCollection<DIM>();
        this->solution_fe_collection = hp::FECollection<DIM>();
        this->solution_quadrature = hp::QCollection<DIM>();
        this->solution_quadrature_1d = hp::QCollection<1>();
        
        delete this->nonmatching_fe_values;
    }

    void
    WaveSolver::run(
        const double start_time,
        const double end_time,
        const unsigned int refinement, 
        const bool plot_sol)
    {
        Assert(
            start_time < end_time,
            ExcMessage("ERROR: Simulation must have a finite duration."));
        Assert(
            (this->fe_degree % 2 == 1), 
            ExcMessage("ERROR: Provided finite element degree for Hermite "
                        "interpolation polynomials should be an odd number."));

        this->analytical_solution.set_current_time(start_time);

        make_grid(refinement);
        setup_levelset();

        distribute_dofs();
        initialise_quadratures();
        initialise_matrices();
        assembly_preprocessing();
        assemble_system();
        prepare_system_matrices();

        std::vector<double> l2_errors;

        std::ostringstream name_frag;
        name_frag << "solution_p_" << this->fe_degree << "_h_" << dx;
        std::string name_prefix = name_frag.str();
        for (auto & s : name_prefix)
            if (s == '.')
                s = '_';

        solve_in_time(end_time, l2_errors, plot_sol, name_prefix);
        
        print_errors(l2_errors, this->fe_degree, true);
        if (plot_sol)
            plot_solution(name_prefix);

        clean_memory();
    }
} // CutHermiteWaveCircle



int main() 
{
    CutHermiteWaveCircle::WaveSystemValues<DIM> system_parameters;

    for (unsigned int j = 5; j < 6; j += 2)
    {
        CutHermiteWaveCircle::WaveSolver system(system_parameters, j);

        for (unsigned int i = 0; i < 4; ++i)
            system.run(0., 3., i, (i == 1));
    }

    return 0;
}