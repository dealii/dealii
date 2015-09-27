/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2003 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Guido Kanschat, University of Heidelberg, 2003
 *         Baerbel Janssen, University of Heidelberg, 2010
 *         Wolfgang Bangerth, Texas A&M University, 2010
 */


// parallel geometric multigrid. work in progress!

#define DEBUG_RENUMBERING

// As discussed in the introduction, most of
// this program is copied almost verbatim
// from step-6, which itself is only a slight
// modification of step-5. Consequently, a
// significant part of this program is not
// new if you've read all the material up to
// step-6, and we won't comment on that part
// of the functionality that is
// unchanged. Rather, we will focus on those
// aspects of the program that have to do
// with the multigrid functionality which
// forms the new aspect of this tutorial
// program.

// @sect3{Include files}

// Again, the first few include files
// are already known, so we won't
// comment on them:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

// These, now, are the include necessary for
// the multilevel methods. The first two
// declare classes that allow us to enumerate
// degrees of freedom not only on the finest
// mesh level, but also on intermediate
// levels (that's what the MGDoFHandler class
// does) as well as allow to access this
// information (iterators and accessors over
// these cells).
//
// The rest of the include files deals with
// the mechanics of multigrid as a linear
// operator (solver or preconditioner).
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

// This is C++:
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// The last step is as in all
// previous programs:
namespace Step50
{
  using namespace dealii;





  // send to proc 0 and compute&print checksum
  void check(const ConstraintMatrix & cm)
  {
      unsigned int myid=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

      std::string rank0data;
      {
          const IndexSet & is = cm.get_local_lines();
          unsigned int num_lines = cm.n_constraints();
          std::ostringstream oss;
          {
          boost::archive::text_oarchive oa(oss);
          oa << num_lines;
          for (unsigned int n=0;n<is.n_elements();++n)
            {
              types::global_dof_index idx = is.nth_index_in_set(n);
              if (cm.is_constrained(idx))
              {
                  const std::vector<std::pair<types::global_dof_index,double> > * v
                          = cm.get_constraint_entries(idx);
                  oa << idx;
                  unsigned int s = v->size();
                  oa << s;
                  for (unsigned int i=0;i<v->size();++i)
                      oa << (*v)[i].first << (*v)[i].second;
                  double inhom = cm.get_inhomogeneity(idx);
                  oa << inhom;
                  //std::cout << "check() idx: " << idx << " v: " << v->size() << " inhom: " << inhom << std::endl;
              }
            }
          }
//          std::cout << "check(): myid=" << myid
//                    << " data length: " << oss.str().size()
//                    << " lines: " << num_lines
//                    << std::endl;
          std::string data = oss.str();
          char * ptr = (char*)data.c_str();
          if (myid!=0)
            MPI_Send(ptr,data.size()+1,MPI_CHAR, 0,7001,MPI_COMM_WORLD);
          else
              rank0data = data;
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if (myid==0)
      {
          ConstraintMatrix all;
          unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);
          std::vector<char> receive;
          for (unsigned int p=0;p<numprocs;++p)
          {
              MPI_Status status;
              status.MPI_SOURCE=0;
              std::string instr;
              if (p!=0)
              {
                  int len;
              MPI_Probe(p, 7001, MPI_COMM_WORLD, &status);
              MPI_Get_count(&status, MPI_BYTE, &len);
              receive.resize(len);
              char * ptr = &receive[0];
              MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                  MPI_COMM_WORLD, &status);
              instr = ptr;
              }
              else
              {
                  instr = rank0data;
//                  std::cout << " my length is of course " << instr.length()<< std::endl;
              }
              std::istringstream iss(instr);
              boost::archive::text_iarchive ia(iss);

              unsigned int num_lines;
              ia >> num_lines;
              std::cout << "check(): got "
                        << " lines: " << num_lines
                        << " from " << status.MPI_SOURCE
                        << std::endl;

              for (unsigned int line=0;line<num_lines;++line)
              {
                  types::global_dof_index idx;
                  std::vector<std::pair<types::global_dof_index,double> > d;
                  double inhom;
                  unsigned int ll;

                  ia >> idx >> ll;
                  for (unsigned int i=0;i<ll;++i)
                  {
                      types::global_dof_index idx;
                      double val;
                      ia >> idx >> val;
                      d.push_back(std::make_pair(idx, val));
                  }
                  ia >> inhom;
//                  std::cout << "check() idx: " << idx << " d: " << d.size() << " inhom: " << inhom << std::endl;
                  if (!all.is_constrained(idx))
                  {
                      all.add_line(idx);
                      all.add_entries(idx,d);
                      all.set_inhomogeneity(idx, inhom);
                  }
                  else
                  {
                      const std::vector<std::pair<types::global_dof_index,double> > * v
                              = all.get_constraint_entries(idx);
                      if (v->size() != d.size())
                          std::cout << "check(): DIFFERENT size " << idx << std::endl;
                      else
                      {
                          for (unsigned int i=0;i<d.size();++i)
                          {
                              if (d[i]!=(*v)[i])
                                std::cout << "check(): DIFFERENT " << idx
                                          << " entry " << i
                                          << " a: " << d[i].first
                                          << " b: " << (*v)[i].first
                                             << std::endl;

                          }

                      }




                  }
              }
          }

          all.close();
          all.print(std::cout);
      }

      MPI_Barrier(MPI_COMM_WORLD);


  }




  // @sect3{The <code>LaplaceProblem</code> class template}

  // This main class is basically the same
  // class as in step-6. As far as member
  // functions is concerned, the only addition
  // is the <code>assemble_multigrid</code>
  // function that assembles the matrices that
  // correspond to the discrete operators on
  // intermediate levels:
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem (const unsigned int deg);
    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void assemble_multigrid ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    parallel::distributed::Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>    mg_dof_handler;

    typedef TrilinosWrappers::SparseMatrix matrix_t;
    typedef TrilinosWrappers::MPI::Vector vector_t;

    matrix_t system_matrix;

    IndexSet locally_relevant_set;

    // We need an additional object for the
    // hanging nodes constraints. They are
    // handed to the transfer object in the
    // multigrid. Since we call a condense
    // inside the multigrid these constraints
    // are not allowed to be inhomogeneous so
    // we store them in different ConstraintMatrix
    // objects.
    ConstraintMatrix     hanging_node_constraints;
    ConstraintMatrix     constraints;

    vector_t       solution;
    vector_t       system_rhs;

    const unsigned int degree;

    // The following four objects are the
    // only additional member variables,
    // compared to step-6. They first three
    // represent the
    // operators that act on individual
    // levels of the multilevel hierarchy,
    // rather than on the finest mesh as do
    // the objects above while the last object
    // stores information about the boundary
    // indices on each level and information
    // about indices lying on a refinement
    // edge between two different refinement
    // levels.
    //
    // To facilitate having objects on each
    // level of a multilevel hierarchy,
    // deal.II has the MGLevelObject class
    // template that provides storage for
    // objects on each level. What we need
    // here are matrices on each level, which
    // implies that we also need sparsity
    // patterns on each level. As outlined in
    // the @ref mg_paper, the operators
    // (matrices) that we need are actually
    // twofold: one on the interior of each
    // level, and one at the interface
    // between each level and that part of
    // the domain where the mesh is
    // coarser. In fact, we will need the
    // latter in two versions: for the
    // direction from coarse to fine mesh and
    // from fine to coarse. Fortunately,
    // however, we here have a self-adjoint
    // problem for which one of these is the
    // transpose of the other, and so we only
    // have to build one; we choose the one
    // from coarse to fine.
    MGLevelObject<matrix_t> mg_matrices;
    MGLevelObject<matrix_t> mg_interface_matrices;
    MGConstrainedDoFs                    mg_constrained_dofs;
  };



  // @sect3{Nonconstant coefficients}

  // The implementation of nonconstant
  // coefficients is copied verbatim
  // from step-5 and step-6:

  template <int dim>
  class Coefficient : public Function<dim>
  {
  public:
    Coefficient () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;
  };



  template <int dim>
  double Coefficient<dim>::value (const Point<dim> &p,
                                  const unsigned int) const
  {
    if (p.square() < 0.5*0.5)
      return 20;
    else
      return 1;
  }



  template <int dim>
  void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
                                     std::vector<double>            &values,
                                     const unsigned int              component) const
  {
    const unsigned int n_points = points.size();

    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));

    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Coefficient<dim>::value (points[i]);
  }


  // @sect3{The <code>LaplaceProblem</code> class implementation}

  // @sect4{LaplaceProblem::LaplaceProblem}

  // The constructor is left mostly
  // unchanged. We take the polynomial degree
  // of the finite elements to be used as a
  // constructor argument and store it in a
  // member variable.
  //
  // By convention, all adaptively refined
  // triangulations in deal.II never change by
  // more than one level across a face between
  // cells. For our multigrid algorithms,
  // however, we need a slightly stricter
  // guarantee, namely that the mesh also does
  // not change by more than refinement level
  // across vertices that might connect two
  // cells. In other words, we must prevent the
  // following situation:
  //
  // @image html limit_level_difference_at_vertices.png ""
  //
  // This is achieved by passing the
  // Triangulation::limit_level_difference_at_vertices
  // flag to the constructor of the
  // triangulation class.
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem (const unsigned int degree)
    :
    triangulation (MPI_COMM_WORLD,Triangulation<dim>::
                   limit_level_difference_at_vertices,
                   parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy),
    fe (degree),
    mg_dof_handler (triangulation),
    degree(degree)
  {
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)!=0)
      deallog.depth_console(0);
  }

  std::string id_to_string(const CellId &id)
  {
    std::ostringstream ss;
    ss << id;
    return ss.str();
  }


  // @sect4{LaplaceProblem::setup_system}

  // The following function extends what the
  // corresponding one in step-6 did. The top
  // part, apart from the additional output,
  // does the same:
  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    mg_dof_handler.distribute_dofs (fe);
    mg_dof_handler.distribute_mg_dofs (fe);


#ifdef     DEBUG_RENUMBERING
    static unsigned int n=0;
    ++n;
    unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);
    if (numprocs==1) // renumber DoFs
      {
        std::map<std::string,std::vector<types::global_dof_index> > dofmap;
        std::map<std::string,std::vector<types::global_dof_index> > mgdofmap;

        {
          std::ifstream f(("ordering." + Utilities::int_to_string(n)).c_str());

          while (!f.eof())
            {
              CellId id;
              f >> id;
              if (f.eof())
                break;
              std::vector<types::global_dof_index> &d = dofmap[id_to_string(id)];
              d.reserve(fe.dofs_per_cell);
              for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                {
                  unsigned int temp;
                  f >> temp;
                  d.push_back(temp);
                }
            }
        }
        {

          std::ifstream f(("mgordering." + Utilities::int_to_string(n)).c_str());

          while (!f.eof())
            {
              CellId id;
              f >> id;
              if (f.eof())
                break;
              std::vector<types::global_dof_index> &d = mgdofmap[id_to_string(id)];
              d.reserve(fe.dofs_per_cell);
              for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                {
                  unsigned int temp;
                  f >> temp;
                  d.push_back(temp);
                }
            }
        }
        std::ofstream f(("temp." + Utilities::int_to_string(n)).c_str(), std::ofstream::out);

        for (typename DoFHandler<dim>::active_cell_iterator cell = mg_dof_handler.begin_active();
             cell != mg_dof_handler.end(); ++cell)
          {
            if (!cell->is_locally_owned())
              continue;

            std::vector<types::global_dof_index>   &renumbered = dofmap[id_to_string(cell->id())];
            if (renumbered.size()!=fe.dofs_per_cell)
              {
                std::cout << " ERROR " << cell->id() << " not found!" << std::endl;

              }
            cell->set_dof_indices(renumbered);
            cell->update_cell_dof_indices_cache();

            std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
            f << cell->id() << ' ';
            cell->get_dof_indices(local_dof_indices);
            for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
              f << local_dof_indices[i] << ' ';
            f << std::endl;
          }
        for (typename DoFHandler<dim>::level_cell_iterator cell = mg_dof_handler.begin();
             cell != mg_dof_handler.end(); ++cell)
          {
            if (!cell->is_locally_owned_on_level())
              continue;

            std::vector<types::global_dof_index>   &renumbered = mgdofmap[id_to_string(cell->id())];
            cell->set_mg_dof_indices(renumbered);
            cell->update_cell_dof_indices_cache();
          }

      }
       if (numprocs>1)
         {
           unsigned int myid=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
           unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

           for (unsigned int i=0; i<numprocs; ++i)
             {
               MPI_Barrier(MPI_COMM_WORLD);
               if (myid!=i)
                 continue;
                 {
                   std::ofstream f(("ordering." + Utilities::int_to_string(n)).c_str(), (myid>0)?std::ofstream::app:std::ofstream::out);
                   std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);
                   for (typename DoFHandler<dim>::active_cell_iterator cell = mg_dof_handler.begin_active();
                        cell != mg_dof_handler.end(); ++cell)
                     {
                       if (!cell->is_locally_owned())
                         continue;

                       f << cell->id() << ' ';
                       cell->get_dof_indices(local_dof_indices);
                       for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                         f << local_dof_indices[i] << ' ';
                       f << std::endl;
                     }
                   f.close();
                 }
                 {
                   std::ofstream f(("mgordering." + Utilities::int_to_string(n)).c_str(), (myid>0)?std::ofstream::app:std::ofstream::out);
                   std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);
                   for (typename DoFHandler<dim>::level_cell_iterator cell = mg_dof_handler.begin();
                       cell != mg_dof_handler.end(); ++cell)
                     {
                       if (!cell->is_locally_owned_on_level())
                         continue;

                       f << cell->id() << ' ';
                       cell->get_mg_dof_indices(local_dof_indices);
                       for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                         f << local_dof_indices[i] << ' ';
                       f << std::endl;
                     }
                   f.close();
                 }
             }
         }

#endif



    DoFTools::extract_locally_relevant_dofs (mg_dof_handler,
                                             locally_relevant_set);


    //solution.reinit (mg_dof_handler.n_dofs());
    //system_rhs.reinit (mg_dof_handler.n_dofs());
    solution.reinit(mg_dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
    system_rhs.reinit(mg_dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);

    // But it starts to be a wee bit different
    // here, although this still doesn't have
    // anything to do with multigrid
    // methods. step-6 took care of boundary
    // values and hanging nodes in a separate
    // step after assembling the global matrix
    // from local contributions. This works,
    // but the same can be done in a slightly
    // simpler way if we already take care of
    // these constraints at the time of copying
    // local contributions into the global
    // matrix. To this end, we here do not just
    // compute the constraints do to hanging
    // nodes, but also due to zero boundary
    // conditions. We will
    // use this set of constraints later on to
    // help us copy local contributions
    // correctly into the global linear system
    // right away, without the need for a later
    // clean-up stage:
    constraints.reinit (locally_relevant_set);
    hanging_node_constraints.reinit (locally_relevant_set);
    DoFTools::make_hanging_node_constraints (mg_dof_handler, hanging_node_constraints);
    DoFTools::make_hanging_node_constraints (mg_dof_handler, constraints);

    typename FunctionMap<dim>::type      dirichlet_boundary;
    ConstantFunction<dim>                    homogeneous_dirichlet_bc (1.0);
    dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
    VectorTools::interpolate_boundary_values (mg_dof_handler,
                                              dirichlet_boundary,
                                              constraints);
    constraints.close ();
    hanging_node_constraints.close ();

    //check(constraints);
    //check(hanging_node_constraints);


    CompressedSimpleSparsityPattern csp(mg_dof_handler.n_dofs(), mg_dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (mg_dof_handler, csp, constraints);
    system_matrix.reinit (mg_dof_handler.locally_owned_dofs(), csp, MPI_COMM_WORLD, true);

    // The multigrid constraints have to be
    // initialized. They need to know about
    // the boundary values as well, so we
    // pass the <code>dirichlet_boundary</code>
    // here as well.
    mg_constrained_dofs.clear();
    mg_constrained_dofs.initialize(mg_dof_handler, dirichlet_boundary);


    // Now for the things that concern the
    // multigrid data structures. First, we
    // resize the multilevel objects to hold
    // matrices and sparsity patterns for every
    // level. The coarse level is zero (this is
    // mandatory right now but may change in a
    // future revision). Note that these
    // functions take a complete, inclusive
    // range here (not a starting index and
    // size), so the finest level is
    // <code>n_levels-1</code>.  We first have
    // to resize the container holding the
    // SparseMatrix classes, since they have to
    // release their SparsityPattern before the
    // can be destroyed upon resizing.
    const unsigned int n_levels = triangulation.n_global_levels();

    mg_interface_matrices.resize(0, n_levels-1);
    mg_interface_matrices.clear ();
    mg_matrices.resize(0, n_levels-1);
    mg_matrices.clear ();

    // Now, we have to provide a matrix on each
    // level. To this end, we first use the
    // MGTools::make_sparsity_pattern function
    // to first generate a preliminary
    // compressed sparsity pattern on each
    // level (see the @ref Sparsity module for
    // more information on this topic) and then
    // copy it over to the one we really
    // want. The next step is to initialize
    // both kinds of level matrices with these
    // sparsity patterns.
    //
    // It may be worth pointing out that the
    // interface matrices only have entries for
    // degrees of freedom that sit at or next
    // to the interface between coarser and
    // finer levels of the mesh. They are
    // therefore even sparser than the matrices
    // on the individual levels of our
    // multigrid hierarchy. If we were more
    // concerned about memory usage (and
    // possibly the speed with which we can
    // multiply with these matrices), we should
    // use separate and different sparsity
    // patterns for these two kinds of
    // matrices.
    for (unsigned int level=0; level<n_levels; ++level)
      {
        CompressedSparsityPattern csp;
        csp.reinit(mg_dof_handler.n_dofs(level),
                   mg_dof_handler.n_dofs(level));
        MGTools::make_sparsity_pattern(mg_dof_handler, csp, level);

        mg_matrices[level].reinit(mg_dof_handler.locally_owned_mg_dofs(level),
                                  mg_dof_handler.locally_owned_mg_dofs(level),
                                  csp,
                                  MPI_COMM_WORLD, true);

        mg_interface_matrices[level].reinit(mg_dof_handler.locally_owned_mg_dofs(level),
                                            mg_dof_handler.locally_owned_mg_dofs(level),
                                            csp,
                                            MPI_COMM_WORLD, true);
      }
  }


  // @sect4{LaplaceProblem::assemble_system}

  // The following function assembles the
  // linear system on the finest level of the
  // mesh. It is almost exactly the same as in
  // step-6, with the exception that we don't
  // eliminate hanging nodes and boundary
  // values after assembling, but while copying
  // local contributions into the global
  // matrix. This is not only simpler but also
  // more efficient for large problems.
  //
  // This latter trick is something that only
  // found its way into deal.II over time and
  // wasn't used in the initial version of this
  // tutorial program. There is, however, a
  // discussion of this function in the
  // introduction of step-27.
  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    const QGauss<dim>  quadrature_formula(degree+1);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points  |  update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const Coefficient<dim> coefficient;
    std::vector<double>    coefficient_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = mg_dof_handler.begin_active(),
    endc = mg_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs = 0;

          fe_values.reinit (cell);

          coefficient.value_list (fe_values.get_quadrature_points(),
                                  coefficient_values);

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  cell_matrix(i,j) += (coefficient_values[q_point] *
                                       fe_values.shape_grad(i,q_point) *
                                       fe_values.shape_grad(j,q_point) *
                                       fe_values.JxW(q_point));

                cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                                10.0 *
                                fe_values.JxW(q_point));
              }

          cell->get_dof_indices (local_dof_indices);
          constraints.distribute_local_to_global (cell_matrix, cell_rhs,
                                                  local_dof_indices,
                                                  system_matrix, system_rhs);
        }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }


  // @sect4{LaplaceProblem::assemble_multigrid}

  // The next function is the one that builds
  // the linear operators (matrices) that
  // define the multigrid method on each level
  // of the mesh. The integration core is the
  // same as above, but the loop below will go
  // over all existing cells instead of just
  // the active ones, and the results must be
  // entered into the correct matrix. Note also
  // that since we only do multilevel
  // preconditioning, no right-hand side needs
  // to be assembled here.
  //
  // Before we go there, however, we have to
  // take care of a significant amount of book
  // keeping:
  template <int dim>
  void LaplaceProblem<dim>::assemble_multigrid ()
  {
    QGauss<dim>  quadrature_formula(1+degree);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const Coefficient<dim> coefficient;
    std::vector<double>    coefficient_values (n_q_points);

    // Next a few things that are specific to building the multigrid
    // data structures (since we only need them in the current
    // function, rather than also elsewhere, we build them here
    // instead of the <code>setup_system</code> function). Some of the
    // following may be a bit obscure if you're not familiar with the
    // algorithm actually implemented in deal.II to support multilevel
    // algorithms on adaptive meshes; if some of the things below seem
    // strange, take a look at the @ref mg_paper.
    //
    // Our first job is to identify those degrees of freedom on each level
    // that are located on interfaces between adaptively refined levels, and
    // those that lie on the interface but also on the exterior boundary of
    // the domain. The <code>MGConstraints</code> already computed the
    // information for us when we called initialize in

    // <code>setup_system()</code>.
    // of type IndexSet on each level (get_refinement_edge_indices(),
    // get_refinement_edge_boundary_indices()).


    //    std::vector<std::vector<bool> > interface_dofs
//      = mg_constrained_dofs.get_refinement_edge_indices ();
//    std::vector<std::vector<bool> > boundary_interface_dofs
//      = mg_constrained_dofs.get_refinement_edge_boundary_indices ();


    // now communicate mg_constraint_dofs
    if (false)
    {
      for (unsigned int l=0;l<triangulation.n_global_levels(); ++l)
        {
          std::cout << "level " << l << std::endl;

          std::map<unsigned int, IndexSet> to_send;


          // determine dofs we need to get information for
          {
            IndexSet dofset;
            DoFTools::extract_locally_relevant_mg_dofs (mg_dof_handler,
                dofset, l);
            dofset.subtract_set(mg_dof_handler.locally_owned_mg_dofs(l));

            unsigned int dest_cpu = 0;
            while (dofset.n_elements()>0)
              {
                types::global_dof_index first_idx = dofset.nth_index_in_set(0);

                while (!mg_dof_handler.locally_owned_mg_dofs_per_processor(l)[dest_cpu]
                                                                              .is_element(first_idx))
                  ++dest_cpu;
                Assert(mg_dof_handler.locally_owned_mg_dofs_per_processor(l)[dest_cpu]
                                                                             .is_element(first_idx),
                                                                             ExcInternalError());
                to_send[dest_cpu] = dofset & mg_dof_handler.locally_owned_mg_dofs_per_processor(l)[dest_cpu];
                dofset.subtract_set(mg_dof_handler.locally_owned_mg_dofs_per_processor(l)[dest_cpu]);
              }
          }


          // create list of processor to send to and get a list of things to
          // receive
          std::vector<unsigned int> destinations;
          for (std::map<unsigned int, IndexSet>::iterator it = to_send.begin();
              it != to_send.end(); ++it)
            destinations.push_back(it->first);

          std::vector<unsigned int> to_receive =
              Utilities::MPI::compute_point_to_point_communication_pattern(MPI_COMM_WORLD,
                                                            destinations);

          std::cout << "I want to send " << destinations.size()
              << " and receive " << to_receive.size() << std::endl;

          // send messages
          std::vector<MPI_Request> requests (to_send.size() + to_receive.size());
          unsigned int req_index = 0;
          std::vector<std::string> buffers(to_send.size());
          {
            unsigned int index = 0;
            for (std::map<unsigned int, IndexSet>::iterator it = to_send.begin();
                it!=to_send.end();++it, ++index)
              {
                std::ostringstream oss;
                it->second.block_write(oss);
                buffers[index] = oss.str();

                char *ptr = const_cast<char*>(buffers[index].c_str());
                MPI_Isend(ptr, buffers[index].length(),
                                        MPI_BYTE, it->first,
                                        1000+l, MPI_COMM_WORLD, &requests[req_index]);
                ++req_index;

                std::cout << "send to " << it->first << std::endl;
              }
          }

          // receive
          unsigned int n_to_recieve = to_receive.size();
          std::vector<char> receive;
          std::vector<std::string> buffers2(n_to_recieve);
          unsigned int index = 0;
          while (n_to_recieve>0)
            {
              MPI_Status status;
              int len;
              MPI_Probe(MPI_ANY_SOURCE, 1000+l, MPI_COMM_WORLD, &status);
              MPI_Get_count(&status, MPI_BYTE, &len);
              receive.resize(len);

              char *ptr = &receive[0];
              MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                  MPI_COMM_WORLD, &status);

              IndexSet idxset;
              std::istringstream iss(std::string(ptr, len));
              idxset.block_read(iss);

              std::cout << "got idxset from " << status.MPI_SOURCE << std::endl;

              // compute and send reply
              std::ostringstream oss;
              IndexSet reply_ref_edge = mg_constrained_dofs.refinement_edge_indices[l] & idxset;
	      IndexSet reply_bdry(idxset.size()); 
	      {//= boundary_indices[l]
		
		for (unsigned int n=0;n<idxset.n_elements();++n)
		  {
		    types::global_dof_index idx = idxset.nth_index_in_set(n);
		    if (mg_constrained_dofs.is_boundary_index(l, idx))
		      reply_bdry.add_index(idx);
		  }
		// test
		/*		for (unsigned int n=0;n<;++n)
		  {
		if (mg_constrained_dofs.is_boundary_index(l, n))
		      reply_bdry.add_index(n);
		      }*/
		
	      }

	       // TH: I don't think we need this
	       //	       IndexSet reply_ref_edge_bdry = mg_constrained_dofs.refinement_edge_boundary_indices[l] & idxset;
	       reply_ref_edge.block_write(oss);
	       //reply_ref_edge_bdry.block_write(oss);
	       reply_bdry.block_write(oss);
	       
              buffers2[index] = oss.str();
              char *ptr2 = const_cast<char*>(buffers2[index].c_str());
              MPI_Isend(ptr2, buffers2[index].length(),
                  MPI_BYTE, status.MPI_SOURCE,
                  2000+l, MPI_COMM_WORLD, &requests[req_index]);
              ++req_index;

              ++index;
              --n_to_recieve;
            }

          // receive answers
          {
            std::vector<char> receive;
            for (unsigned int n=0;n<destinations.size();++n)
              {
                MPI_Status status;
                int len;
                MPI_Probe(MPI_ANY_SOURCE, 2000+l, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_BYTE, &len);
                receive.resize(len);

                char *ptr = &receive[0];
                MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                    MPI_COMM_WORLD, &status);

                std::istringstream iss(std::string(ptr,len));
                IndexSet reply_ref_edge;
                //IndexSet reply_ref_edge_bdry;
                reply_ref_edge.block_read(iss);
		IndexSet reply_bdry;
		reply_bdry.block_read(iss);
                //reply_ref_edge_bdry.block_read(iss);
                unsigned int c1 = mg_constrained_dofs.refinement_edge_indices[l].n_elements();
                //unsigned int c2 = mg_constrained_dofs.refinement_edge_boundary_indices[l].n_elements();

                mg_constrained_dofs.refinement_edge_indices[l].add_indices(reply_ref_edge);
                //mg_constrained_dofs.refinement_edge_boundary_indices[l].add_indices(reply_ref_edge_bdry);

		unsigned int c3 = mg_constrained_dofs.boundary_indices[l].size();
		for (unsigned int n=0;n<reply_bdry.n_elements();++n)
		  {
		    types::global_dof_index idx = reply_bdry.nth_index_in_set(n);
		    Assert(mg_constrained_dofs.boundary_indices[l].is_element(idx), ExcInternalError());
		    
		      //		    mg_constrained_dofs.boundary_indices[l].insert(idx);
		  }
		
		/*
                std::cout << "L=" << l
			  << " new " << mg_constrained_dofs.refinement_edge_indices[l].n_elements()-c1
		  //<< " and " << mg_constrained_dofs.refinement_edge_boundary_indices[l].n_elements()-c2
			  << " and " << mg_constrained_dofs.boundary_indices[l].size()-c3
			  << std::endl;*/
		if (mg_constrained_dofs.boundary_indices[l].size()>c3)
		  std::cout << "HEY: boundary indices transfer did something!" << std::endl;
		
              }

          }

          // finish all requests:
          if (requests.size() > 0)
            MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);

          MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    // The indices just identified will later be used to decide where
    // the assembled value has to be added into on each level.  On the
    // other hand, we also have to impose zero boundary conditions on
    // the external boundary of each level. But this the
    // <code>MGConstraints</code> knows it. So we simply ask for them
    // by calling <code>get_boundary_indices ()</code>.  The third
    // step is to construct constraints on all those degrees of
    // freedom: their value should be zero after each application of
    // the level operators. To this end, we construct ConstraintMatrix
    // objects for each level, and add to each of these constraints
    // for each degree of freedom. Due to the way the ConstraintMatrix
    // stores its data, the function to add a constraint on a single
    // degree of freedom and force it to be zero is called
    // Constraintmatrix::add_line(); doing so for several degrees of
    // freedom at once can be done using
    // Constraintmatrix::add_lines():
    std::vector<ConstraintMatrix> boundary_constraints (triangulation.n_global_levels());
    std::vector<ConstraintMatrix> boundary_interface_constraints (triangulation.n_global_levels());
    for (unsigned int level=0; level<triangulation.n_global_levels(); ++level)
      {
        // TODO: here we get missing entries!
        IndexSet dofset;
        DoFTools::extract_locally_relevant_mg_dofs (mg_dof_handler, dofset, level);
        boundary_constraints[level].reinit(dofset);
        boundary_constraints[level].add_lines (mg_constrained_dofs.get_refinement_edge_indices(level));
        boundary_constraints[level].add_lines (mg_constrained_dofs.get_boundary_indices(level));


        boundary_constraints[level].close ();
        //check(boundary_constraints[level]);


        //boundary_interface_constraints[level]
        //.add_lines (mg_constrained_dofs.get_refinement_edge_boundary_indices (level));
        boundary_interface_constraints[level].close ();
        //check(boundary_interface_constraints[level]);
      }

    // Now that we're done with most of our preliminaries, let's start
    // the integration loop. It looks mostly like the loop in
    // <code>assemble_system</code>, with two exceptions: (i) we don't
    // need a right hand side, and more significantly (ii) we don't
    // just loop over all active cells, but in fact all cells, active
    // or not. Consequently, the correct iterator to use is
    // MGDoFHandler::cell_iterator rather than
    // MGDoFHandler::active_cell_iterator. Let's go about it:
    typename DoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
                                            endc = mg_dof_handler.end();

    for (; cell!=endc; ++cell)
      if (cell->level_subdomain_id()==triangulation.locally_owned_subdomain())
        {
          cell_matrix = 0;
          fe_values.reinit (cell);

          coefficient.value_list (fe_values.get_quadrature_points(),
                                  coefficient_values);

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += (coefficient_values[q_point] *
                                     fe_values.shape_grad(i,q_point) *
                                     fe_values.shape_grad(j,q_point) *
                                     fe_values.JxW(q_point));

          // The rest of the assembly is again slightly
          // different. This starts with a gotcha that is easily
          // forgotten: The indices of global degrees of freedom we
          // want here are the ones for current level, not for the
          // global matrix. We therefore need the function
          // MGDoFAccessorLLget_mg_dof_indices, not
          // MGDoFAccessor::get_dof_indices as used in the assembly of
          // the global system:
          cell->get_mg_dof_indices (local_dof_indices);

          // Next, we need to copy local contributions into the level
          // objects. We can do this in the same way as in the global
          // assembly, using a constraint object that takes care of
          // constrained degrees (which here are only boundary nodes,
          // as the individual levels have no hanging node
          // constraints). Note that the
          // <code>boundary_constraints</code> object makes sure that
          // the level matrices contains no contributions from degrees
          // of freedom at the interface between cells of different
          // refinement level.
          boundary_constraints[cell->level()]
          .distribute_local_to_global (cell_matrix,
                                       local_dof_indices,
                                       mg_matrices[cell->level()]);

          // The next step is again slightly more obscure (but
          // explained in the @ref mg_paper): We need the remainder of
          // the operator that we just copied into the
          // <code>mg_matrices</code> object, namely the part on the
          // interface between cells at the current level and cells
          // one level coarser. This matrix exists in two directions:
          // for interior DoFs (index $i$) of the current level to
          // those sitting on the interface (index $j$), and the other
          // way around. Of course, since we have a symmetric
          // operator, one of these matrices is the transpose of the
          // other.
          //
          // The way we assemble these matrices is as follows: since
          // the are formed from parts of the local contributions, we
          // first delete all those parts of the local contributions
          // that we are not interested in, namely all those elements
          // of the local matrix for which not $i$ is an interface DoF
          // and $j$ is not. The result is one of the two matrices
          // that we are interested in, and we then copy it into the
          // <code>mg_interface_matrices</code> object. The
          // <code>boundary_interface_constraints</code> object at the
          // same time makes sure that we delete contributions from
          // all degrees of freedom that are not only on the interface
          // but also on the external boundary of the domain.
          //
          // The last part to remember is how to get the other
          // matrix. Since it is only the transpose, we will later (in
          // the <code>solve()</code> function) be able to just pass
          // the transpose matrix where necessary.

          const IndexSet &interface_dofs_on_level
            = mg_constrained_dofs.get_refinement_edge_indices(cell->level());


          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
	      /** old HEAD:
		  if (!mg_constrained_dofs.at_refinement_edge(cell->level(),local_dof_indices[i])
		  || mg_constrained_dofs.at_refinement_edge(cell->level(),local_dof_indices[j])) */
              if ( !(interface_dofs_on_level.is_element(local_dof_indices[i])==true &&
                     interface_dofs_on_level.is_element(local_dof_indices[j])==false))
                cell_matrix(i,j) = 0;

          boundary_interface_constraints[cell->level()]
          .distribute_local_to_global (cell_matrix,
                                       local_dof_indices,
                                       mg_interface_matrices[cell->level()]);
        }

    for (unsigned int i=0; i<triangulation.n_global_levels(); ++i)
      {
        mg_matrices[i].compress(VectorOperation::add);
        mg_interface_matrices[i].compress(VectorOperation::add);
      }
  }



  // @sect4{LaplaceProblem::solve}

  // This is the other function that is significantly different in
  // support of the multigrid solver (or, in fact, the preconditioner
  // for which we use the multigrid method).
  //
  // Let us start out by setting up two of the components of
  // multilevel methods: transfer operators between levels, and a
  // solver on the coarsest level. In finite element methods, the
  // transfer operators are derived from the finite element function
  // spaces involved and can often be computed in a generic way
  // independent of the problem under consideration. In that case, we
  // can use the MGTransferPrebuilt class that, given the constraints
  // on the global level and an MGDoFHandler object computes the
  // matrices corresponding to these transfer operators.
  //
  // The second part of the following lines deals with the coarse grid
  // solver. Since our coarse grid is very coarse indeed, we decide
  // for a direct solver (a Householder decomposition of the coarsest
  // level matrix), even if its implementation is not particularly
  // sophisticated. If our coarse mesh had many more cells than the
  // five we have here, something better suited would obviously be
  // necessary here.
  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {

    // Create the object that deals with the transfer between
    // different refinement levels. We need to pass it the hanging
    // node constraints.
    MGTransferPrebuilt<vector_t> mg_transfer(hanging_node_constraints, mg_constrained_dofs);
    // Now the prolongation matrix has to be built.  This matrix needs
    // to take the boundary values on each level into account and
    // needs to know about the indices at the refinement edges. The
    // <code>MGConstraints</code> knows about that so pass it as an
    // argument.
    mg_transfer.build_matrices(mg_dof_handler);


    matrix_t &coarse_matrix = mg_matrices[0];
    //coarse_matrix.copy_from (mg_matrices[0]);
    //MGCoarseGridHouseholder<double,vector_t> coarse_grid_solver;
    //coarse_grid_solver.initialize (coarse_matrix);

    SolverControl coarse_solver_control (1000, 1e-10, false, false);
    SolverGMRES<vector_t> coarse_solver(coarse_solver_control);
    PreconditionIdentity id;
    MGCoarseGridLACIteration<SolverGMRES<vector_t>,vector_t> coarse_grid_solver(coarse_solver,
        coarse_matrix,
        id);

    // The next component of a multilevel solver or preconditioner is
    // that we need a smoother on each level. A common choice for this
    // is to use the application of a relaxation method (such as the
    // SOR, Jacobi or Richardson method). The MGSmootherPrecondition
    // class provides support for this kind of smoother. Here, we opt
    // for the application of a single SOR iteration. To this end, we
    // define an appropriate <code>typedef</code> and then setup a
    // smoother object.
    //
    // The last step is to initialize the smoother object with our
    // level matrices and to set some smoothing parameters.  The
    // <code>initialize()</code> function can optionally take
    // additional arguments that will be passed to the smoother object
    // on each level. In the current case for the SOR smoother, this
    // could, for example, include a relaxation parameter. However, we
    // here leave these at their default values. The call to
    // <code>set_steps()</code> indicates that we will use two pre-
    // and two post-smoothing steps on each level; to use a variable
    // number of smoother steps on different levels, more options can
    // be set in the constructor call to the <code>mg_smoother</code>
    // object.
    //
    // The last step results from the fact that
    // we use the SOR method as a smoother -
    // which is not symmetric - but we use the
    // conjugate gradient iteration (which
    // requires a symmetric preconditioner)
    // below, we need to let the multilevel
    // preconditioner make sure that we get a
    // symmetric operator even for nonsymmetric
    // smoothers:
    typedef TrilinosWrappers::PreconditionJacobi Smoother;
    MGSmootherPrecondition<matrix_t, Smoother, vector_t> mg_smoother;
    mg_smoother.initialize(mg_matrices, Smoother::AdditionalData(0.5));
    mg_smoother.set_steps(2);
    //mg_smoother.set_symmetric(false);

    // The next preparatory step is that we
    // must wrap our level and interface
    // matrices in an object having the
    // required multiplication functions. We
    // will create two objects for the
    // interface objects going from coarse to
    // fine and the other way around; the
    // multigrid algorithm will later use the
    // transpose operator for the latter
    // operation, allowing us to initialize
    // both up and down versions of the
    // operator with the matrices we already
    // built:
    mg::Matrix<vector_t> mg_matrix(mg_matrices);
    mg::Matrix<vector_t> mg_interface_up(mg_interface_matrices);
    mg::Matrix<vector_t> mg_interface_down(mg_interface_matrices);

    // Now, we are ready to set up the
    // V-cycle operator and the
    // multilevel preconditioner.
    Multigrid<vector_t > mg(mg_dof_handler,
                            mg_matrix,
                            coarse_grid_solver,
                            mg_transfer,
                            mg_smoother,
                            mg_smoother);
    //mg.set_debug(6);
    mg.set_edge_matrices(mg_interface_down, mg_interface_up);




    PreconditionMG<dim, vector_t, MGTransferPrebuilt<vector_t> >
    preconditioner(mg_dof_handler, mg, mg_transfer);

    if (false)//true)
      {
	static int round=0;
	++round;
	std::cout << " ********* start " << round << std::endl;
	solution = 1.0;
	
	TrilinosWrappers::MPI::Vector temp = solution;
	temp = 0.0;
	preconditioner.vmult(temp, solution);
	  
	std::cout << " ********* end " << round << " " << temp.l2_norm() << " min=" << temp.minimal_value() << std::endl;
      }
    
    
    // With all this together, we can finally
    // get about solving the linear system in
    // the usual way:
    SolverControl solver_control (500, 1e-8*system_rhs.l2_norm(), false);
    SolverGMRES<vector_t>    cg (solver_control);

    //solution = 0;

    if (false)
      {
        // code to optionally compare to Trilinos ML
        TrilinosWrappers::PreconditionAMG prec;

        TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
        //    Amg_data.constant_modes = constant_modes;
        Amg_data.elliptic = true;
        Amg_data.higher_order_elements = true;
        Amg_data.smoother_sweeps = 2;
        Amg_data.aggregation_threshold = 0.02;
        // Amg_data.symmetric = true;

        prec.initialize (system_matrix,
                         Amg_data);
        cg.solve (system_matrix, solution, system_rhs,
                  prec);
      }
    else
      {
        cg.solve (system_matrix, solution, system_rhs,
                  preconditioner);
      }

    constraints.distribute (solution);
  }



  // @sect4{Postprocessing}

  // The following two functions postprocess a
  // solution once it is computed. In
  // particular, the first one refines the mesh
  // at the beginning of each cycle while the
  // second one outputs results at the end of
  // each such cycle. The functions are almost
  // unchanged from those in step-6, with the
  // exception of two minor differences: The
  // KellyErrorEstimator::estimate function
  // wants an argument of type DoFHandler, not
  // MGDoFHandler, and so we have to cast from
  // derived to base class; and we generate
  // output in VTK format, to use the more
  // modern visualization programs available
  // today compared to those that were
  // available when step-6 was written.
  template <int dim>
  void LaplaceProblem<dim>::refine_grid ()
  {
#ifdef     DEBUG_RENUMBERING
    static int n = 0;
    n++;

    if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)==1)
      {
        std::ifstream f(("refine." + Utilities::int_to_string(n)).c_str());
        std::set<CellId> to_refine;
        while(!f.eof())
          {
            CellId temp;
            f >> temp;
            to_refine.insert(temp);
          }

        typename DoFHandler<dim>::active_cell_iterator
        cell = mg_dof_handler.begin_active(),
        endc = mg_dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned()
              && to_refine.find(cell->id()) != to_refine.end())
            {
              cell->set_refine_flag();
                //std::cout << cell->id() << std::endl;
              }
    triangulation.execute_coarsening_and_refinement ();
    return;
      }
#endif
    
      {

    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    TrilinosWrappers::MPI::Vector temp_solution;
    temp_solution.reinit(locally_relevant_set, MPI_COMM_WORLD);
    temp_solution = solution;

    KellyErrorEstimator<dim>::estimate (static_cast<DoFHandler<dim>&>(mg_dof_handler),
                                        QGauss<dim-1>(degree+1),
                                        typename FunctionMap<dim>::type(),
                                        temp_solution,
                                        estimated_error_per_cell);

    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_fraction (triangulation,
                                       estimated_error_per_cell,
                                       0.3, 0.0);

    
    triangulation.prepare_coarsening_and_refinement ();


#ifdef     DEBUG_RENUMBERING
    if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)>1)
    {
      // write refine.n

      std::ofstream f;
      int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      if (myid==0)
	f.open(("refine." + Utilities::int_to_string(n)).c_str());

      {
	// local contribution
	std::ostringstream oss;
	typename DoFHandler<dim>::active_cell_iterator
	  cell = mg_dof_handler.begin_active(),
	  endc = mg_dof_handler.end();
	for (; cell!=endc; ++cell)
	  {
	    if (cell->is_locally_owned() && cell->refine_flag_set())
	      oss << cell->id() << std::endl;
	    
	    cell->clear_coarsen_flag();
	    
	  }
	
	
	if (myid!=0)
	  {
	    std::string s= oss.str();
	    
	    MPI_Send((void*)s.c_str(), s.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	  }
	else
	  f << oss.str();
      }

      if (myid==0)
	for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)-1;++i)
	{
	  MPI_Status status;
	  int len;
	  MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
	  MPI_Get_count(&status, MPI_BYTE, &len);
	  std::vector<char> receive;
	  receive.resize(len+1);

	  char *ptr = &receive[0];
	  MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
		   MPI_COMM_WORLD, &status);
	  f << ptr;
	}
    }
#endif
      }
      
    triangulation.execute_coarsening_and_refinement ();
  }



  template <int dim>
  void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
  {
    DataOut<dim> data_out;

    TrilinosWrappers::MPI::Vector temp_solution;
    temp_solution.reinit(locally_relevant_set, MPI_COMM_WORLD);
    temp_solution = solution;


    TrilinosWrappers::MPI::Vector temp = solution;
    system_matrix.residual(temp,solution,system_rhs);
    TrilinosWrappers::MPI::Vector res_ghosted = temp_solution;
    res_ghosted = temp;

    data_out.attach_dof_handler (mg_dof_handler);
    data_out.add_data_vector (temp_solution, "solution");
    data_out.add_data_vector (res_ghosted, "res");
    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");

    data_out.build_patches (0);

    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, 5) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4) +
                                  ".vtu");
    std::ofstream output (filename.c_str());
    data_out.write_vtu (output);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); ++i)
          filenames.push_back (std::string("solution-") +
                               Utilities::int_to_string (cycle, 5) +
                               "." +
                               Utilities::int_to_string(i, 4) +
                               ".vtu");
        const std::string
        pvtu_master_filename = ("solution-" +
                                Utilities::int_to_string (cycle, 5) +
                                ".pvtu");
        std::ofstream pvtu_master (pvtu_master_filename.c_str());
        data_out.write_pvtu_record (pvtu_master, filenames);

        const std::string
        visit_master_filename = ("solution-" +
                                 Utilities::int_to_string (cycle, 5) +
                                 ".visit");
        std::ofstream visit_master (visit_master_filename.c_str());
        data_out.write_visit_record (visit_master, filenames);

	std::cout << "wrote " << pvtu_master_filename << std::endl;
	
      }
  }


  // @sect4{LaplaceProblem::run}

  // Like several of the functions above, this
  // is almost exactly a copy of of the
  // corresponding function in step-6. The only
  // difference is the call to
  // <code>assemble_multigrid</code> that takes
  // care of forming the matrices on every
  // level that we need in the multigrid
  // method.
  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<20; ++cycle)
      {
        deallog << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube (triangulation);

            triangulation.refine_global (4);
          }
        else
          refine_grid ();

        deallog << "   Number of active cells:       "
                << triangulation.n_global_active_cells()
                << std::endl;

        setup_system ();

        deallog << "   Number of degrees of freedom: "
                << mg_dof_handler.n_dofs()
                << " (by level: ";
        for (unsigned int level=0; level<triangulation.n_global_levels(); ++level)
          deallog << mg_dof_handler.n_dofs(level)
                  << (level == triangulation.n_global_levels()-1
                      ? ")" : ", ");
        deallog << std::endl;


        assemble_system ();
        assemble_multigrid ();

        unsigned int myid=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
        unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);


	if (false)
	  {
	    
        static unsigned int n=0;
        ++n;

        MPI_Barrier(MPI_COMM_WORLD);
        for (unsigned int lvl = 0;lvl<triangulation.n_global_levels();++lvl)
        for (unsigned int i=0; i<numprocs; ++i)
          {
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid!=i)
              continue;

            std::ofstream f(("matrix."+Utilities::int_to_string(n)+"."+Utilities::int_to_string(lvl)).c_str(),
                (myid>0)?std::ofstream::app:std::ofstream::out);
            mg_matrices[lvl].print(f);
            mg_interface_matrices[lvl].print(f);
          }
        MPI_Barrier(MPI_COMM_WORLD);
	  }
	
        deallog << "rhs: " << std::setprecision(15) << system_rhs.l2_norm() << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        solve ();
        output_results (cycle);

	MPI_Barrier(MPI_COMM_WORLD);
	TrilinosWrappers::MPI::Vector temp = solution;
	system_matrix.residual(temp,solution,system_rhs);
	constraints.set_zero(temp);

	for (unsigned int i=0;i<solution.size();++i)
	  if (solution.in_local_range(i) && solution[i]<0.8)
	    std::cout << "WRONG solution at " << i
		      << " is " << solution[i] << std::endl;
	
	deallog << "residual " << temp.l2_norm() << std::endl;
	if (temp.l2_norm()>1e-5)
	  {
	    for (unsigned int i=0;i<temp.size();++i)
	      if (temp.in_local_range(i) && temp[i]>0.06)
                std::cout << "RESIDUAL at " << i
                          << " is " << temp[i] << std::endl;
	    
	  }

      }
  }
}


// @sect3{The main() function}
//
// This is again the same function as
// in step-6:
int main (int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  try
    {
      using namespace dealii;
      using namespace Step50;

      LaplaceProblem<2> laplace_problem(1/*degree*/);
      laplace_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      throw;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      throw;
    }

  return 0;
}
