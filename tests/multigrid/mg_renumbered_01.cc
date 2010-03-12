#define baerbel_mg_test 1

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/utilities.h>

#include <lac/constraint_matrix.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/sparse_ilu.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/precondition.h>

#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

#include <multigrid/multigrid.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_transfer_component.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_coarse.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_matrix.h>

#include <fstream>
#include <sstream>

using namespace dealii;


template <int dim, typename number, int spacedim>
void
reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
	       MGLevelObject<dealii::Vector<number> > &v)
{
  for (unsigned int level=v.get_minlevel();
       level<=v.get_maxlevel();++level)
    {
      unsigned int n = mg_dof.n_dofs (level);
      v[level].reinit(n);
    }
}

template <int dim>
void initialize (const MGDoFHandler<dim> &dof,
    Vector<double> &u)
{
  unsigned int counter=0;
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices(dofs_per_cell);
  for (typename MGDoFHandler<dim>::active_cell_iterator
      cell = dof.begin_active();
      cell != dof.end(); ++cell)
  {
    cell->get_dof_indices(dof_indices);
    for(unsigned int i=0; i<dofs_per_cell; ++i)
      u(dof_indices[i]) = ++counter;
  }
}



template <int dim>
void initialize (const MGDoFHandler<dim> &dof,
    MGLevelObject<Vector<double> > &u)
{
  unsigned int counter=0;
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices(dofs_per_cell);
  typename MGDoFHandler<dim>::cell_iterator
      cell = dof.begin(0);
    cell->get_mg_dof_indices(dof_indices);
    for(unsigned int i=0; i<dofs_per_cell; ++i)
      u[0](dof_indices[i]) = ++counter;
}


template <int dim>
void print_diff (const MGDoFHandler<dim> &dof_1, const MGDoFHandler<dim> &dof_2,
    const Vector<double> &u, const Vector<double> &v)
{
  Vector<double> diff;
  diff.reinit (u);
  const unsigned int dofs_per_cell = dof_1.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices_1(dofs_per_cell);
  std::vector<unsigned int> dof_indices_2(dofs_per_cell);
  for (typename MGDoFHandler<dim>::active_cell_iterator
      cell_1 = dof_1.begin_active(), cell_2 = dof_2.begin_active();
      cell_1 != dof_1.end(); ++cell_1, ++cell_2)
  {
    cell_1->get_dof_indices(dof_indices_1);
    cell_2->get_dof_indices(dof_indices_2);
    for(unsigned int i=0; i<dofs_per_cell; ++i)
      diff(dof_indices_1[i]) = u(dof_indices_1[i]) - v(dof_indices_2[i]);
  }
  std::cout << std::endl;
  std::cout << "diff " << diff.l2_norm() << std::endl;
}

template <int dim>
void print(const MGDoFHandler<dim> &dof, std::vector<std::vector<bool> > &interface_dofs)
{
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices(dofs_per_cell);
  for(unsigned int l=0; l<dof.get_tria().n_levels(); ++l)
  {
    deallog << std::endl;
    deallog << "Level " << l << std::endl;
    for (typename MGDoFHandler<dim>::cell_iterator
        cell = dof.begin(l);
        cell != dof.end(l); ++cell)
    {
      cell->get_mg_dof_indices(dof_indices);
      for(unsigned int i=0; i<dofs_per_cell; ++i)
        deallog << ' ' << interface_dofs[l][dof_indices[i]];
    }
  }
}


template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem (const unsigned int deg);
    void run ();

  private:
    void setup_system ();
    void assemble_multigrid (const MGDoFHandler<dim> &mg_dof, 
        MGLevelObject<SparseMatrix<double> > &matrices,
        std::vector<std::set<unsigned int> > 
        &boundary_indices);
    void test ();
    void test_renumber ();
    void test_interface_dofs ();
    void refine_local ();

    Triangulation<dim>   triangulation;
    FESystem<dim>            fe;
    MGDoFHandler<dim>    mg_dof_handler;
    MGDoFHandler<dim>    mg_dof_handler_renumbered;

    MGLevelObject<SparsityPattern> mg_sparsity_renumbered;
    MGLevelObject<SparseMatrix<double> > mg_matrices_renumbered;
    MGLevelObject<SparsityPattern> mg_sparsity;
    MGLevelObject<SparseMatrix<double> > mg_matrices;

    const unsigned int degree;

    std::vector<std::set<unsigned int> > 
    boundary_indices;

    std::vector<std::set<unsigned int> > 
    boundary_indices_renumbered;
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int deg) :
  triangulation (Triangulation<dim>::limit_level_difference_at_vertices),
  fe (FE_Q<dim> (deg),2, FE_Q<dim> (deg),2),
  mg_dof_handler (triangulation),
  mg_dof_handler_renumbered (triangulation),
  degree(deg)
{}


  template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  mg_dof_handler.distribute_dofs (fe);
  mg_dof_handler_renumbered.distribute_dofs (fe);

  std::vector<unsigned int> block_component (2*dim,0);
  for(unsigned int c=dim; c<2*dim; ++c)
    block_component[c] = 1;

  const unsigned int nlevels = triangulation.n_levels();

  DoFHandler<dim> &dof = mg_dof_handler_renumbered;
  DoFRenumbering::component_wise (dof, block_component);
  //DoFRenumbering::Cuthill_McKee (dof);
  for (unsigned int level=0;level<nlevels;++level)
  {
    DoFRenumbering::component_wise (mg_dof_handler_renumbered, level, block_component);
    //DoFRenumbering::Cuthill_McKee (mg_dof_handler_renumbered, level);
  }

  deallog << "Number of degrees of freedom: "
    << mg_dof_handler.n_dofs();

  for (unsigned int l=0;l<triangulation.n_levels();++l)
    deallog << "   " << 'L' << l << ": "
      << mg_dof_handler.n_dofs(l);
  deallog  << std::endl;

  boundary_indices.resize(triangulation.n_levels());
  boundary_indices_renumbered.resize(triangulation.n_levels());

  mg_matrices.resize(0, nlevels-1);
  mg_matrices.clear ();
  mg_matrices_renumbered.resize(0, nlevels-1);
  mg_matrices_renumbered.clear ();
  mg_sparsity.resize(0, nlevels-1);
  mg_sparsity_renumbered.resize(0, nlevels-1);

  for (unsigned int level=0;level<nlevels;++level)
  {
    mg_sparsity[level].reinit (mg_dof_handler.n_dofs(level),
        mg_dof_handler.n_dofs(level),
        mg_dof_handler.max_couplings_between_dofs());
    mg_sparsity_renumbered[level].reinit (mg_dof_handler_renumbered.n_dofs(level),
        mg_dof_handler_renumbered.n_dofs(level),
        mg_dof_handler_renumbered.max_couplings_between_dofs());
    MGTools::make_sparsity_pattern (mg_dof_handler, mg_sparsity[level], level);
    MGTools::make_sparsity_pattern (mg_dof_handler_renumbered, 
        mg_sparsity_renumbered[level], level);
    mg_sparsity[level].compress();
    mg_sparsity_renumbered[level].compress();
    mg_matrices[level].reinit(mg_sparsity[level]);
    mg_matrices_renumbered[level].reinit(mg_sparsity_renumbered[level]);
  }
}

  template <int dim>
void LaplaceProblem<dim>::assemble_multigrid (const MGDoFHandler<dim> &mgdof,
    MGLevelObject<SparseMatrix<double> > &mgmatrices,
    std::vector<std::set<unsigned int> > &boundaryindices) 
{  
  QGauss<dim>  quadrature_formula(1+degree);

  FEValues<dim> fe_values (fe, quadrature_formula, 
      update_values   | update_gradients |
      update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   dofs_per_face   = fe.dofs_per_face;
  const unsigned int   n_q_points      = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  std::vector<unsigned int> face_dof_indices (dofs_per_face);


  std::vector<std::vector<bool> > interface_dofs;
  std::vector<std::vector<bool> > boundary_interface_dofs;
  for (unsigned int level = 0; level<triangulation.n_levels(); ++level)
  {
    std::vector<bool> tmp (mgdof.n_dofs(level));
    interface_dofs.push_back (tmp);
    boundary_interface_dofs.push_back (tmp);
  }

  MGTools::extract_inner_interface_dofs(mgdof, 
      interface_dofs, boundary_interface_dofs);

  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    dirichlet_bc(fe.n_components());
  dirichlet_boundary[0] =             &dirichlet_bc;

  MGTools::make_boundary_list (mgdof, dirichlet_boundary,
			       boundaryindices);

  std::vector<ConstraintMatrix> boundary_constraints (triangulation.n_levels());
  for (unsigned int level=0; level<triangulation.n_levels(); ++level)
    {
      for (unsigned int i=0; i<boundary_interface_dofs[level].size(); ++i)
        if(!(boundary_interface_dofs[level][i] && 
            boundaryindices[level].find(i) != boundaryindices[level].end()))
          boundary_interface_dofs[level][i] = false;

      boundary_constraints[level].add_lines (interface_dofs[level]);
      boundary_constraints[level].add_lines (boundaryindices[level]);
      boundary_constraints[level].close ();
    }

  typename MGDoFHandler<dim>::cell_iterator cell = mgdof.begin(),
           endc = mgdof.end();

  for (; cell!=endc; ++cell)
  {
    const unsigned int level = cell->level();
    cell_matrix = 0;

    fe_values.reinit (cell);
    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          for (unsigned int d=0; d<fe.n_components(); ++d)
            cell_matrix(i,j) += (fe_values.shape_grad_component(j,q_point,d)
              * fe_values.shape_grad_component(i,q_point,d))
              * fe_values.JxW(q_point);
        }
      }

//    cell->get_mg_dof_indices (local_dof_indices);
//      boundary_constraints[level]
//	.distribute_local_to_global (cell_matrix,
//				     local_dof_indices,
//				     mgmatrices[level]);

  }
}

  template <int dim>
void LaplaceProblem<dim>::test_interface_dofs () 
{
  std::vector<std::vector<bool> > interface_dofs;
  std::vector<std::vector<bool> > boundary_interface_dofs;
  for (unsigned int level = 0; level<triangulation.n_levels(); ++level)
  {
    std::vector<bool> tmp (mg_dof_handler.n_dofs(level));
    interface_dofs.push_back (tmp);
    boundary_interface_dofs.push_back (tmp);
  }

  MGTools::extract_inner_interface_dofs(mg_dof_handler_renumbered, 
      interface_dofs, boundary_interface_dofs);
  deallog << "1. Test" << std::endl;
  print(mg_dof_handler_renumbered, interface_dofs);

  deallog << std::endl;

  MGTools::extract_inner_interface_dofs(mg_dof_handler, 
      interface_dofs, boundary_interface_dofs);

  deallog << "2. Test" << std::endl;
  print(mg_dof_handler, interface_dofs);
}

  template <int dim>
void LaplaceProblem<dim>::test_renumber () 
{
  std::vector<std::vector<bool> > v(triangulation.n_levels());
  for(unsigned int l=0; l<triangulation.n_levels(); ++l)
  {
    v[l].resize(mg_dof_handler_renumbered.n_dofs(l));
    std::fill (v[l].begin(), v[l].end(), false);
  }

  const unsigned int dofs_per_cell = 
    mg_dof_handler_renumbered.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices(dofs_per_cell);
  for (typename MGDoFHandler<dim>::cell_iterator
      cell = mg_dof_handler_renumbered.begin();
      cell != mg_dof_handler_renumbered.end(); ++cell)
  {
    cell->get_mg_dof_indices(dof_indices);
    for(unsigned int i=0; i<dofs_per_cell; ++i)
      v[cell->level()][dof_indices[i]] = true;
  }
  for(unsigned int l=0; l<triangulation.n_levels(); ++l)
  {
    deallog << "Level " << l << std::endl;
    for(unsigned int dof=0; dof<mg_dof_handler_renumbered.n_dofs(l);
        ++dof)
      if(!v[l][dof])
        deallog << "FALSE " << dof << std::endl;
  }
}


  template <int dim>
void LaplaceProblem<dim>::test () 
{
  assemble_multigrid(mg_dof_handler, mg_matrices, 
      boundary_indices);
  assemble_multigrid(mg_dof_handler_renumbered, 
      mg_matrices_renumbered, boundary_indices_renumbered);
  
  const unsigned int min_l = mg_matrices.get_minlevel();
  const unsigned int max_l = mg_matrices.get_maxlevel();
  for(unsigned int l=min_l; l<max_l; ++l)
  {
    mg_matrices[l] = IdentityMatrix(mg_dof_handler.n_dofs(l));
    mg_matrices_renumbered[l] = IdentityMatrix(mg_dof_handler.n_dofs(l));
  }

  GrowingVectorMemory<>   vector_memory;
  GrowingVectorMemory<>   vector_memory_renumbered;

  MGTransferPrebuilt<Vector<double> > mg_transfer;
  mg_transfer.build_matrices(mg_dof_handler, boundary_indices);
  MGTransferPrebuilt<Vector<double> > mg_transfer_renumbered;
  mg_transfer_renumbered.build_matrices(mg_dof_handler_renumbered, boundary_indices_renumbered);
  FullMatrix<double> coarse_matrix;
  coarse_matrix.copy_from (mg_matrices[0]);
  MGCoarseGridHouseholder<double, Vector<double> > mg_coarse;
  mg_coarse.initialize(coarse_matrix);

  FullMatrix<double> coarse_matrix_renumered;
  coarse_matrix_renumered.copy_from (mg_matrices_renumbered[0]);
  MGCoarseGridHouseholder<double, Vector<double> > mg_coarse_renumbered;
  mg_coarse_renumbered.initialize(coarse_matrix_renumered);

  typedef PreconditionIdentity RELAXATION;
  MGSmootherPrecondition<SparseMatrix<double>, RELAXATION, Vector<double> >
    mg_smoother(vector_memory);

  MGSmootherPrecondition<SparseMatrix<double>, RELAXATION, Vector<double> >
    mg_smoother_renumbered(vector_memory_renumbered);

  RELAXATION::AdditionalData smoother_data;
  mg_smoother.initialize(mg_matrices, smoother_data);
  mg_smoother_renumbered.initialize(mg_matrices_renumbered, smoother_data);

  mg_smoother.set_steps(1);
  mg_smoother_renumbered.set_steps(1);

  MGMatrix<SparseMatrix<double>, Vector<double> >
    mg_matrix(&mg_matrices);

  MGMatrix<SparseMatrix<double>, Vector<double> >
    mg_matrix_renumbered(&mg_matrices_renumbered);

  Multigrid<Vector<double> > mg(mg_dof_handler,
      mg_matrix,
      mg_coarse,
      mg_transfer,
      mg_smoother,
      mg_smoother);

  Multigrid<Vector<double> > mg_renumbered(mg_dof_handler_renumbered,
      mg_matrix_renumbered,
      mg_coarse_renumbered,
      mg_transfer_renumbered,
      mg_smoother_renumbered,
      mg_smoother_renumbered);

  PreconditionMG<dim, Vector<double>,
    MGTransferPrebuilt<Vector<double> > >
      preconditioner(mg_dof_handler, mg, mg_transfer);

  PreconditionMG<dim, Vector<double>,
    MGTransferPrebuilt<Vector<double> > >
      preconditioner_renumbered(mg_dof_handler_renumbered, 
          mg_renumbered, mg_transfer_renumbered);

  Vector<double> test, dst, dst_renumbered;
  test.reinit(mg_dof_handler.n_dofs());
  dst.reinit(test);
  dst_renumbered.reinit(test);

  initialize(mg_dof_handler, test);
  std::cout << "1. Test " << test.l2_norm() << std::endl;
  preconditioner.vmult(dst, test);
  std::cout << "1. Test " << dst.l2_norm() << std::endl;

  initialize(mg_dof_handler_renumbered, test);
  std::cout << "2. Test " << test.l2_norm() << std::endl;
  preconditioner_renumbered.vmult(dst_renumbered, test);
  std::cout << "2. Test " << dst_renumbered.l2_norm() << std::endl;

  print_diff (mg_dof_handler_renumbered, mg_dof_handler, dst_renumbered, dst);
}




  template <int dim>
void LaplaceProblem<dim>::refine_local ()
{
  bool cell_refined = false;
  for (typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active();
      cell != triangulation.end(); ++cell)
  {
      for (unsigned int vertex=0;
          vertex < GeometryInfo<dim>::vertices_per_cell;
          ++vertex)
      {
        const Point<dim> p = cell->vertex(vertex);
        const Point<dim> origin = (dim == 2 ?
                                    Point<dim>(0,0) :
                                    Point<dim>(0,0,0));
        const double dist = p.distance(origin);
        if(dist<0.25/M_PI)
        {
          cell->set_refine_flag ();
          cell_refined = true;
          break;
        }
      }
  }
  //Wenn nichts verfeinert wurde bisher, global verfeinern!
  if(!cell_refined)
    for (typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active();
        cell != triangulation.end(); ++cell)
      cell->set_refine_flag();


  triangulation.execute_coarsening_and_refinement ();
}

  template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<7; ++cycle)
  {
    deallog << "Cycle " << cycle << std::endl;

    if (cycle == 0)
    {
      GridGenerator::hyper_cube(triangulation, -1, 1);
      triangulation.refine_global (1);
    }
    //triangulation.refine_global (1);
    refine_local ();
    setup_system ();
    test();
  };
}

int main () 
{
  LaplaceProblem<2> laplace_problem_2d(1);
  laplace_problem_2d.run ();

  return 0;
}
