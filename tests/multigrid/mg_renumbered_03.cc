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
#include <numerics/mesh_worker.h>
#include <numerics/mesh_worker_loop.h>
#include <numerics/mesh_worker_output.h>

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
#include <base/std_cxx1x/bind.h>

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
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices(dofs_per_cell);
  for (typename MGDoFHandler<dim>::active_cell_iterator
      cell = dof.begin_active();
      cell != dof.end(); ++cell)
  {
    cell->get_dof_indices(dof_indices);
    for(unsigned int i=0; i<dofs_per_cell; ++i)
    {
      const unsigned int comp = dof.get_fe().system_to_component_index(i).first;
      u(dof_indices[i]) = comp+1;
    }
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
void diff (Vector<double> &diff, const MGDoFHandler<dim> &dof,
    const Vector<double> &v, const unsigned int level)
{
  diff.reinit (v);
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> mg_dof_indices(dofs_per_cell);
  const unsigned int n_comp = dof.get_fe().n_components();
  for (typename MGDoFHandler<dim>::cell_iterator
      cell= dof.begin(level);
      cell != dof.end(level); ++cell)
  {
    cell->get_mg_dof_indices(mg_dof_indices);
    for(unsigned int i=0; i<dofs_per_cell/n_comp; ++i)
    {
      const unsigned int ni = n_comp* i;
      const unsigned int i0 = mg_dof_indices[ni];
      const unsigned int i1 = mg_dof_indices[ni+1];
      const unsigned int i2 = mg_dof_indices[ni+2];
      diff(i0) = 2*v(i0) - v(i1);
      diff(i1) = 3*v(i1) - 2*v(i2);
      diff(i2) = v(i2) - 3*v(i0);
    }
  }
}

template <int dim>
void diff (Vector<double> &diff, const MGDoFHandler<dim> &dof_1, const MGDoFHandler<dim> &dof_2,
    const Vector<double> &u, const Vector<double> &v, const unsigned int level)
{
  diff.reinit (u);
  const unsigned int dofs_per_cell = dof_1.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices_1(dofs_per_cell);
  std::vector<unsigned int> dof_indices_2(dofs_per_cell);
  for (typename MGDoFHandler<dim>::cell_iterator
      cell_1 = dof_1.begin(level), cell_2 = dof_2.begin(level);
      cell_1 != dof_1.end(level); ++cell_1, ++cell_2)
  {
    cell_1->get_mg_dof_indices(dof_indices_1);
    cell_2->get_mg_dof_indices(dof_indices_2);
    for(unsigned int i=0; i<dofs_per_cell; ++i)
      diff(dof_indices_1[i]) = u(dof_indices_1[i]) - v(dof_indices_2[i]);
  }
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
class OutputCreator : public Subscriptor
{
  public:
    void cell(MeshWorker::DoFInfo<dim>& dinfo,
		     MeshWorker::IntegrationInfo<dim>& info);

};

template <int dim>
void OutputCreator<dim>::cell(
  MeshWorker::DoFInfo<dim>& dinfo,
  MeshWorker::IntegrationInfo<dim>& info)
{
  const FEValuesBase<dim>& fe = info.fe_values();
  const std::vector<std::vector<double> >& uh = info.values[0];

  const unsigned int square_root = std::pow(fe.n_quadrature_points, 1./dim)+.5;
  for (unsigned int k1=0; k1<square_root; ++k1)
  {
    for (unsigned int k2=0; k2<square_root; ++k2)
    {
      for(unsigned int d=0; d<dim; ++d)
        dinfo.quadrature_value(k1*square_root+k2, d) = fe.quadrature_point(k1*square_root+k2)[d];
      for(unsigned int i=0; i<uh.size(); ++i)
        dinfo.quadrature_value(k1*square_root+k2, i+dim) = uh[i][k1*square_root+k2];
    }
  }
}


template <int dim>
class LaplaceProblem
{
  public:
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;

    LaplaceProblem (const unsigned int deg);
    void run ();

  private:
    void setup_system ();
    void test ();
    void output_gpl(const MGDoFHandler<dim> &dof,
        MGLevelObject<Vector<double> > &v);
    void refine_local ();

    Triangulation<dim>   triangulation;
    const MappingQ1<dim>      mapping;
    FESystem<dim>            fe;
    MGDoFHandler<dim>    mg_dof_handler;
    MGDoFHandler<dim>    mg_dof_handler_renumbered;

    const unsigned int degree;
    std::vector<std::set<unsigned int> >
      boundary_indices, boundary_indices_renumbered;

};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int deg) :
  triangulation (Triangulation<dim>::limit_level_difference_at_vertices),
  fe (FE_Q<dim> (deg),3),
  mg_dof_handler (triangulation),
  mg_dof_handler_renumbered (triangulation),
  degree(deg)
{}


  template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  mg_dof_handler.distribute_dofs (fe);
  mg_dof_handler_renumbered.distribute_dofs (fe);

  const unsigned int nlevels = triangulation.n_levels();

  DoFHandler<dim> &dof = mg_dof_handler_renumbered;
  DoFRenumbering::component_wise (dof);
  for (unsigned int level=0;level<nlevels;++level)
    DoFRenumbering::component_wise (mg_dof_handler_renumbered, level);

  boundary_indices.resize(nlevels);
  boundary_indices_renumbered.resize(nlevels);

  for(unsigned int l=0; l<nlevels; ++l)
  {
    boundary_indices_renumbered[l].clear();
    boundary_indices[l].clear();
  }

  deallog << "Number of degrees of freedom: "
    << mg_dof_handler.n_dofs();

  for (unsigned int l=0;l<triangulation.n_levels();++l)
    deallog << "   " << 'L' << l << ": "
      << mg_dof_handler.n_dofs(l);
  deallog  << std::endl;

}

template <int dim>
void
LaplaceProblem<dim>::output_gpl(const MGDoFHandler<dim> &dof,
    MGLevelObject<Vector<double> > &v)
{
  MeshWorker::IntegrationInfoBox<dim> info_box;
  const unsigned int n_gauss_points = dof.get_fe().tensor_degree();
  QTrapez<1> trapez;
  QIterated<dim> quadrature(trapez, n_gauss_points);
  info_box.cell_quadrature = quadrature;
  NamedData<MGLevelObject<Vector<double> >* > data;
  data.add(&v, "mg_vector");
  info_box.cell_selector.add("mg_vector");
  info_box.initialize_update_flags();
  UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
  info_box.add_update_flags(update_flags, true, true, true, true);
  info_box.initialize(fe, mapping, data);  
  
  MeshWorker::DoFInfo<dim> dof_info(dof);

  MeshWorker::Assembler::GnuplotPatch assembler;
  assembler.initialize(dim, quadrature.size(), dim+dof.get_fe().n_components());

  for(unsigned int l=0; l<triangulation.n_levels(); ++l)
  {
    OutputCreator<dim> matrix_integrator;
    MeshWorker::integration_loop<dim, dim> (
      dof.begin(l), dof.end(l),
      dof_info, info_box,
      std_cxx1x::bind(&OutputCreator<dim>::cell, &matrix_integrator, _1, _2),
      0,
      0,
      assembler);
  }
}




  template <int dim>
void LaplaceProblem<dim>::test ()
{
  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    dirichlet_bc(fe.n_components());
  dirichlet_boundary[0] =             &dirichlet_bc;

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(mg_dof_handler, dirichlet_boundary);
  MGConstrainedDoFs mg_constrained_dofs_renumbered;
  mg_constrained_dofs_renumbered.initialize(mg_dof_handler_renumbered, dirichlet_boundary);

  ConstraintMatrix constraints;
  constraints.close();

  MGTransferPrebuilt<Vector<double> > mg_transfer (constraints, mg_constrained_dofs);
  mg_transfer.build_matrices(mg_dof_handler);
  MGTransferPrebuilt<Vector<double> > mg_transfer_renumbered (constraints, mg_constrained_dofs_renumbered);
  mg_transfer_renumbered.build_matrices(mg_dof_handler_renumbered);

  Vector<double> test;
  test.reinit(mg_dof_handler.n_dofs());

  MGLevelObject<Vector<double> > v(0, triangulation.n_levels()-1);
  MGLevelObject<Vector<double> > u(0, triangulation.n_levels()-1);
  MGLevelObject<Vector<double> > d(0, triangulation.n_levels()-1);

  initialize(mg_dof_handler, test);
  mg_transfer.copy_to_mg(mg_dof_handler, v, test);

  initialize(mg_dof_handler_renumbered, test);
  mg_transfer_renumbered.copy_to_mg(mg_dof_handler_renumbered, u, test);
  for(unsigned int l=0; l<triangulation.n_levels(); ++l)
  {
    diff(d[l], mg_dof_handler_renumbered, u[l],l);
    deallog << l << " " << u[l].l2_norm() << '\t' << v[l].l2_norm() << '\t'
      << d[l].l2_norm()<< std::endl;
    for(unsigned int i=0; i<d[l].size(); ++i)
      if(d[l](i)!=0)
        deallog << i << " " << d[l](i) << std::endl;
  }
  output_gpl(mg_dof_handler_renumbered, d);
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
  for (unsigned int cycle=0; cycle<6; ++cycle)
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
  std::ofstream logfile("mg_renumbered_03/output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  LaplaceProblem<2> laplace_problem_2d(1);
  laplace_problem_2d.run ();
}
