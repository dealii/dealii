#include <deal.II/dofs/agglomeration_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/poly_utils.h>
#include <deal.II/grid/reference_cell.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_interpolate.h>

#include <algorithm>

using namespace dealii;

// Check that the R-tree based agglomeration works also in 3D.


enum class GridType
{
  grid_generator,
  unstructured
};

enum class PartitionerType
{
  metis,
  rtree
};


template <int dim>
class Poisson
{
private:
  void
  make_agglomerated_grid();
  void
  setup_agglomeration();

  Triangulation<dim>                         tria;
  MappingQ<dim>                              mapping;
  std::unique_ptr<AgglomerationHandler<dim>> ah;
  std::unique_ptr<GridTools::Cache<dim>>     cached_tria;

public:
  Poisson(const GridType        &grid_type        = GridType::grid_generator,
          const PartitionerType &partitioner_type = PartitionerType::rtree,
          const unsigned int     extraction_level = 2);
  void
  run();

  double          penalty_constant = 10.;
  GridType        grid_type;
  PartitionerType partitioner_type;
  unsigned int    extraction_level;
};



template <int dim>
Poisson<dim>::Poisson(const GridType        &grid_type,
                      const PartitionerType &partitioner_type,
                      const unsigned int     extraction_level)
  : mapping(1)
  , grid_type(grid_type)
  , partitioner_type(partitioner_type)
  , extraction_level(extraction_level)
{}

template <int dim>
void
Poisson<dim>::make_agglomerated_grid()
{
  GridIn<dim> grid_in;
  if (grid_type == GridType::unstructured)
    {
      // This test does not read an external grid
      DEAL_II_NOT_IMPLEMENTED();
    }
  else
    {
      GridGenerator::hyper_ball(tria);
      tria.refine_global(3); // 6
    }
  std::cout << "Size of tria: " << tria.n_active_cells() << std::endl;
  cached_tria = std::make_unique<GridTools::Cache<dim>>(tria, mapping);
  ah          = std::make_unique<AgglomerationHandler<dim>>(*cached_tria);


  if (partitioner_type == PartitionerType::metis)
    {
      // Not done in this test
      DEAL_II_NOT_IMPLEMENTED();
    }
  else if (partitioner_type == PartitionerType::rtree)
    {
      // Partition with Rtree

      namespace bgi                                   = boost::geometry::index;
      static constexpr unsigned int max_elem_per_node = 8; // 2
      std::vector<std::pair<BoundingBox<dim>,
                            typename Triangulation<dim>::active_cell_iterator>>
                   boxes(tria.n_active_cells());
      unsigned int i = 0;
      for (const auto &cell : tria.active_cell_iterators())
        boxes[i++] = std::make_pair(mapping.get_bounding_box(cell), cell);

      auto tree = pack_rtree<bgi::rstar<max_elem_per_node>>(boxes);
      std::cout << "Total number of levels in the tree: " << n_levels(tree)
                << std::endl;

      CellsAgglomerator<dim, decltype(tree)> agglomerator{tree,
                                                          extraction_level};
      const auto vec_agglomerates = agglomerator.extract_agglomerates();
      // Flag elements for agglomeration
      for (const auto &agglo : vec_agglomerates)
        {
          ah->define_agglomerate(agglo);
        }


      std::cout << "N subdomains = " << ah->n_agglomerates() << std::endl;
    }
}

template <int dim>
void
Poisson<dim>::run()
{
  make_agglomerated_grid();
}

int
main()
{
  {
    Poisson<3> poisson_problem{GridType::grid_generator,
                               PartitionerType::rtree,
                               2};
    poisson_problem.run();
  }



  return 0;
}
