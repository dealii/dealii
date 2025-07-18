// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_large_count.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>

#include <deal.II/grid/grid_tools.h>

#include <fstream>
#include <memory>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
namespace GridGenerator
{
  template <int dim, int spacedim>
  void
  hyper_cube(Triangulation<dim, spacedim> &tria,
             const double                  left,
             const double                  right,
             const bool                    colorize);
} // namespace GridGenerator

namespace parallel
{
  namespace fullydistributed
  {
    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    Triangulation<dim, spacedim>::Triangulation(const MPI_Comm mpi_communicator)
      : parallel::DistributedTriangulationBase<dim, spacedim>(mpi_communicator)
      , settings(TriangulationDescription::Settings::default_setting)
      , partitioner([](dealii::Triangulation<dim, spacedim> &tria,
                       const unsigned int                    n_partitions) {
        GridTools::partition_triangulation_zorder(n_partitions, tria);
      })
      , currently_processing_create_triangulation_for_internal_usage(false)
      , currently_processing_prepare_coarsening_and_refinement_for_internal_usage(
          false)
    {}



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::create_triangulation(
      const TriangulationDescription::Description<dim, spacedim>
        &construction_data)
    {
      // check if the communicator of this parallel triangulation has been used
      // to construct the TriangulationDescription::Description
      Assert(construction_data.comm == this->mpi_communicator,
             ExcMessage("MPI communicators do not match!"));

      // store internally the settings
      settings = construction_data.settings;

      // set the smoothing properties
      if (settings &
          TriangulationDescription::Settings::construct_multigrid_hierarchy)
        this->set_mesh_smoothing(
          static_cast<
            typename dealii::Triangulation<dim, spacedim>::MeshSmoothing>(
            dealii::Triangulation<dim, spacedim>::none |
            Triangulation<dim, spacedim>::limit_level_difference_at_vertices));
      else
        this->set_mesh_smoothing(
          static_cast<
            typename dealii::Triangulation<dim, spacedim>::MeshSmoothing>(
            dealii::Triangulation<dim, spacedim>::none));

      this->set_mesh_smoothing(construction_data.smoothing);

      // clear internal data structures
      this->coarse_cell_id_to_coarse_cell_index_vector.clear();
      this->coarse_cell_index_to_coarse_cell_id_vector.clear();

      // check if no locally relevant coarse-grid cells have been provided
      if (construction_data.coarse_cell_vertices.empty())
        {
          // 1) create a dummy hypercube
          currently_processing_create_triangulation_for_internal_usage = true;
          GridGenerator::hyper_cube(*this, 0, 1, false);
          currently_processing_create_triangulation_for_internal_usage = false;

          // 2) mark cell as artificial
          auto cell = this->begin();
          cell->set_subdomain_id(dealii::numbers::artificial_subdomain_id);
          cell->set_level_subdomain_id(
            dealii::numbers::artificial_subdomain_id);

          // 3) set up dummy mapping between locally relevant coarse-grid cells
          //    and global cells
          this->coarse_cell_id_to_coarse_cell_index_vector.emplace_back(
            numbers::invalid_coarse_cell_id, 0);
          this->coarse_cell_index_to_coarse_cell_id_vector.emplace_back(
            numbers::invalid_coarse_cell_id);
        }
      else
        {
          // 1) store `coarse-cell index to coarse-cell id`-mapping
          this->coarse_cell_index_to_coarse_cell_id_vector =
            construction_data.coarse_cell_index_to_coarse_cell_id;

          // 2) set up `coarse-cell id to coarse-cell index`-mapping
          std::map<types::coarse_cell_id, unsigned int>
            coarse_cell_id_to_coarse_cell_index_vector;
          for (unsigned int i = 0;
               i < construction_data.coarse_cell_index_to_coarse_cell_id.size();
               ++i)
            coarse_cell_id_to_coarse_cell_index_vector
              [construction_data.coarse_cell_index_to_coarse_cell_id[i]] = i;

          for (auto i : coarse_cell_id_to_coarse_cell_index_vector)
            this->coarse_cell_id_to_coarse_cell_index_vector.emplace_back(i);

          // create locally-relevant
          currently_processing_prepare_coarsening_and_refinement_for_internal_usage =
            true;
          currently_processing_create_triangulation_for_internal_usage = true;
          dealii::Triangulation<dim, spacedim>::create_triangulation(
            construction_data);
          currently_processing_prepare_coarsening_and_refinement_for_internal_usage =
            false;
          currently_processing_create_triangulation_for_internal_usage = false;

          // create a copy of cell_infos such that we can sort them
          auto cell_infos = construction_data.cell_infos;

          // sort cell_infos on each level separately (as done in
          // dealii::Triangulation::create_triangulation())
          for (auto &cell_info : cell_infos)
            std::sort(cell_info.begin(),
                      cell_info.end(),
                      [&](const TriangulationDescription::CellData<dim> &a,
                          const TriangulationDescription::CellData<dim> &b) {
                        const CellId a_id(a.id);
                        const CellId b_id(b.id);

                        const auto a_coarse_cell_index =
                          this->coarse_cell_id_to_coarse_cell_index(
                            a_id.get_coarse_cell_id());
                        const auto b_coarse_cell_index =
                          this->coarse_cell_id_to_coarse_cell_index(
                            b_id.get_coarse_cell_id());

                        // according to their coarse-cell index and if that is
                        // same according to their cell id (the result is that
                        // cells on each level are sorted according to their
                        // index on that level - what we need in the following
                        // operations)
                        if (a_coarse_cell_index != b_coarse_cell_index)
                          return a_coarse_cell_index < b_coarse_cell_index;
                        else
                          return a_id < b_id;
                      });

          // 4a) set all cells artificial (and set the actual
          //     (level_)subdomain_ids in the next step)
          for (const auto &cell : this->cell_iterators())
            {
              if (cell->is_active())
                cell->set_subdomain_id(
                  dealii::numbers::artificial_subdomain_id);

              cell->set_level_subdomain_id(
                dealii::numbers::artificial_subdomain_id);
            }

          // 4b) set actual (level_)subdomain_ids
          for (unsigned int level = 0;
               level < cell_infos.size() && !cell_infos[level].empty();
               ++level)
            {
              auto cell      = this->begin(level);
              auto cell_info = cell_infos[level].begin();
              for (; cell_info != cell_infos[level].end(); ++cell_info)
                {
                  // find cell that has the correct cell
                  while (cell_info->id != cell->id().template to_binary<dim>())
                    ++cell;

                  // subdomain id
                  if (cell->is_active())
                    cell->set_subdomain_id(cell_info->subdomain_id);

                  // level subdomain id
                  if (settings & TriangulationDescription::Settings::
                                   construct_multigrid_hierarchy)
                    cell->set_level_subdomain_id(cell_info->level_subdomain_id);
                }
            }
        }

      this->update_number_cache();
      this->update_cell_relations();
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::create_triangulation(
      const std::vector<Point<spacedim>>       &vertices,
      const std::vector<dealii::CellData<dim>> &cells,
      const SubCellData                        &subcelldata)
    {
      Assert(
        currently_processing_create_triangulation_for_internal_usage,
        ExcMessage(
          "You have called the overload of\n"
          "\n"
          "    parallel::fullydistributed::Triangulation::"
          "create_triangulation()\n"
          "\n"
          "which takes 3 arguments. This function is not yet implemented for "
          "this class. If you have not called this function directly, it "
          "might have been called via a function from the GridGenerator or "
          "GridIn namespace. To set up a fully-distributed Triangulation with "
          "these utility functions, please start by using the same process to "
          "set up a serial Triangulation, parallel::shared::Triangulation, or "
          "a parallel::distributed::Triangulation. Once that is complete use "
          "the copy_triangulation() member function to finish setting up the "
          "original fully distributed Triangulation. Alternatively, you can "
          "use TriangulationDescription::Utilities::"
          "create_description_from_triangulation() or "
          "create_description_from_triangulation_in_groups() to create the "
          "description of the local partition, and pass that description to "
          "parallel::fullydistributed::Triangulation::create_triangulation()."));

      dealii::Triangulation<dim, spacedim>::create_triangulation(vertices,
                                                                 cells,
                                                                 subcelldata);
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::copy_triangulation(
      const dealii::Triangulation<dim, spacedim> &other_tria)
    {
      // pointer to the triangulation for which the construction data
      // should be created (normally it is the input triangulation but
      // in the case of a serial triangulation we create a copy which should
      // be used)
      const dealii::Triangulation<dim, spacedim> *other_tria_ptr = &other_tria;

      // temporary serial triangulation (since the input triangulation is const
      // and we might modify its subdomain_ids and level_subdomain_ids during
      // partitioning)
      dealii::Triangulation<dim, spacedim> serial_tria;

      // check if other triangulation is not a parallel one, which needs to be
      // partitioned
      if (dynamic_cast<const dealii::parallel::TriangulationBase<dim, spacedim>
                         *>(&other_tria) == nullptr)
        {
          // actually copy the serial triangulation
          serial_tria.copy_triangulation(other_tria);

          // partition triangulation
          this->partitioner(serial_tria,
                            dealii::Utilities::MPI::n_mpi_processes(
                              this->mpi_communicator));

          // partition multigrid levels
          if (this->is_multilevel_hierarchy_constructed())
            GridTools::partition_multigrid_levels(serial_tria);

          // use the new serial triangulation to create the construction data
          other_tria_ptr = &serial_tria;
        }

      // create construction data
      const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(*other_tria_ptr,
                                              this->mpi_communicator,
                                              this->settings);

      // finally create triangulation
      this->create_triangulation(construction_data);
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::set_partitioner(
      const std::function<void(dealii::Triangulation<dim, spacedim> &,
                               const unsigned int)> &partitioner,
      const TriangulationDescription::Settings      &settings)
    {
      this->partitioner = partitioner;
      this->settings    = settings;
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::set_partitioner(
      const RepartitioningPolicyTools::Base<dim, spacedim> &partitioner,
      const TriangulationDescription::Settings             &settings)
    {
      this->partitioner_distributed = &partitioner;
      this->settings                = settings;
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::repartition()
    {
      // signal that repartitioning has started
      this->signals.pre_distributed_repartition();

      // create construction_data with the help of the partitioner
      const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(
          *this,
          this->partitioner_distributed->partition(*this),
          this->settings);

      // clear old content
      this->clear();
      this->coarse_cell_id_to_coarse_cell_index_vector.clear();
      this->coarse_cell_index_to_coarse_cell_id_vector.clear();

      // use construction_data to set up new triangulation
      this->create_triangulation(construction_data);

      // signal that repartitioning has completed
      this->signals.post_distributed_repartition();
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::execute_coarsening_and_refinement()
    {
      DEAL_II_NOT_IMPLEMENTED();
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    bool Triangulation<dim, spacedim>::prepare_coarsening_and_refinement()
    {
      Assert(
        currently_processing_prepare_coarsening_and_refinement_for_internal_usage,
        ExcMessage("No coarsening and refinement is supported!"));

      return dealii::Triangulation<dim, spacedim>::
        prepare_coarsening_and_refinement();
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    std::size_t Triangulation<dim, spacedim>::memory_consumption() const
    {
      const std::size_t mem =
        this->dealii::parallel::TriangulationBase<dim, spacedim>::
          memory_consumption() +
        MemoryConsumption::memory_consumption(
          coarse_cell_id_to_coarse_cell_index_vector) +
        MemoryConsumption::memory_consumption(
          coarse_cell_index_to_coarse_cell_id_vector);
      return mem;
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    bool Triangulation<dim, spacedim>::is_multilevel_hierarchy_constructed()
      const
    {
      return (
        settings &
        TriangulationDescription::Settings::construct_multigrid_hierarchy);
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    unsigned int Triangulation<dim, spacedim>::
      coarse_cell_id_to_coarse_cell_index(
        const types::coarse_cell_id coarse_cell_id) const
    {
      const auto coarse_cell_index = std::lower_bound(
        coarse_cell_id_to_coarse_cell_index_vector.begin(),
        coarse_cell_id_to_coarse_cell_index_vector.end(),
        coarse_cell_id,
        [](const std::pair<types::coarse_cell_id, unsigned int> &pair,
           const types::coarse_cell_id &val) { return pair.first < val; });
      if (coarse_cell_index !=
          coarse_cell_id_to_coarse_cell_index_vector.cend())
        return coarse_cell_index->second;
      else
        return numbers::invalid_unsigned_int; // cell could no be found
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    types::coarse_cell_id
      Triangulation<dim, spacedim>::coarse_cell_index_to_coarse_cell_id(
        const unsigned int coarse_cell_index) const
    {
      AssertIndexRange(coarse_cell_index,
                       coarse_cell_index_to_coarse_cell_id_vector.size());

      const auto coarse_cell_id =
        coarse_cell_index_to_coarse_cell_id_vector[coarse_cell_index];
      AssertThrow(coarse_cell_id != numbers::invalid_coarse_cell_id,
                  ExcMessage("You are trying to access a dummy cell!"));
      return coarse_cell_id;
    }


    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::update_cell_relations()
    {
      // Reorganize memory for local_cell_relations.
      this->local_cell_relations.clear();
      this->local_cell_relations.reserve(this->n_locally_owned_active_cells());

      for (const auto &cell : this->active_cell_iterators())
        if (cell->is_locally_owned())
          this->local_cell_relations.emplace_back(
            cell, ::dealii::CellStatus::cell_will_persist);
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::save(const std::string &filename) const
    {
#ifdef DEAL_II_WITH_MPI

      Assert(
        this->cell_attached_data.n_attached_deserialize == 0,
        ExcMessage(
          "Not all SolutionTransfer objects have been deserialized after the last call to load()."));
      Assert(this->n_cells() > 0,
             ExcMessage("Can not save() an empty Triangulation."));

      const int myrank =
        Utilities::MPI::this_mpi_process(this->mpi_communicator);
      const int mpisize =
        Utilities::MPI::n_mpi_processes(this->mpi_communicator);

      // Compute global offset for each rank.
      unsigned int n_locally_owned_cells = this->n_locally_owned_active_cells();

      unsigned int global_first_cell = 0;

      int ierr = MPI_Exscan(&n_locally_owned_cells,
                            &global_first_cell,
                            1,
                            MPI_UNSIGNED,
                            MPI_SUM,
                            this->mpi_communicator);
      AssertThrowMPI(ierr);

      global_first_cell *= sizeof(unsigned int);


      if (myrank == 0)
        {
          std::string   fname = std::string(filename) + ".info";
          std::ofstream f(fname);
          f << "version nproc n_attached_fixed_size_objs n_attached_variable_size_objs n_global_active_cells"
            << std::endl
            << ::dealii::internal::CellAttachedDataSerializer<dim, spacedim>::
                 version_number
            << " " << Utilities::MPI::n_mpi_processes(this->mpi_communicator)
            << " " << this->cell_attached_data.pack_callbacks_fixed.size()
            << " " << this->cell_attached_data.pack_callbacks_variable.size()
            << " " << this->n_global_active_cells() << std::endl;
        }

      // Save cell attached data.
      this->save_attached_data(global_first_cell,
                               this->n_global_active_cells(),
                               filename);

      // Save triangulation description.
      {
        MPI_Info info;
        int      ierr = MPI_Info_create(&info);
        AssertThrowMPI(ierr);

        const std::string fname_tria = filename + "_triangulation.data";

        // Open file.
        MPI_File fh;
        ierr = MPI_File_open(this->mpi_communicator,
                             fname_tria.c_str(),
                             MPI_MODE_CREATE | MPI_MODE_WRONLY,
                             info,
                             &fh);
        AssertThrowMPI(ierr);

        ierr = MPI_File_set_size(fh, 0); // delete the file contents
        AssertThrowMPI(ierr);
        // this barrier is necessary, because otherwise others might already
        // write while one core is still setting the size to zero.
        ierr = MPI_Barrier(this->mpi_communicator);
        AssertThrowMPI(ierr);
        ierr = MPI_Info_free(&info);
        AssertThrowMPI(ierr);
        // ------------------

        // Create construction data.
        const auto construction_data = TriangulationDescription::Utilities::
          create_description_from_triangulation(*this,
                                                this->mpi_communicator,
                                                this->settings);

        // Pack.
        std::vector<char> buffer;
        dealii::Utilities::pack(construction_data, buffer, false);

        // Write offsets to file.
        const std::uint64_t buffer_size = buffer.size();

        std::uint64_t offset = 0;

        ierr = MPI_Exscan(
          &buffer_size,
          &offset,
          1,
          Utilities::MPI::mpi_type_id_for_type<decltype(buffer_size)>,
          MPI_SUM,
          this->mpi_communicator);
        AssertThrowMPI(ierr);

        // Write offsets to file.
        ierr = MPI_File_write_at(
          fh,
          myrank * sizeof(std::uint64_t),
          &buffer_size,
          1,
          Utilities::MPI::mpi_type_id_for_type<decltype(buffer_size)>,
          MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        // global position in file
        const std::uint64_t global_position =
          mpisize * sizeof(std::uint64_t) + offset;

        // Write buffers to file.
        ierr = dealii::Utilities::MPI::LargeCount::File_write_at_c(
          fh,
          global_position,
          buffer.data(),
          buffer.size(), // local buffer
          MPI_CHAR,
          MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        ierr = MPI_File_close(&fh);
        AssertThrowMPI(ierr);
      }
#else
      (void)filename;

      AssertThrow(false, ExcNeedsMPI());
#endif
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::load(const std::string &filename)
    {
#ifdef DEAL_II_WITH_MPI
      Assert(this->n_cells() == 0,
             ExcMessage("load() only works if the Triangulation is empty!"));


      unsigned int version, numcpus, attached_count_fixed,
        attached_count_variable, n_global_active_cells;
      {
        std::string   fname = std::string(filename) + ".info";
        std::ifstream f(fname);
        AssertThrow(f.fail() == false, ExcIO());
        std::string firstline;
        getline(f, firstline);
        f >> version >> numcpus >> attached_count_fixed >>
          attached_count_variable >> n_global_active_cells;
      }

      const auto expected_version = ::dealii::internal::
        CellAttachedDataSerializer<dim, spacedim>::version_number;

      AssertThrow(version == expected_version,
                  ExcMessage("Incompatible version found in .info file."));

      // Load description and construct the triangulation.
      {
        const int myrank =
          Utilities::MPI::this_mpi_process(this->mpi_communicator);
        const int mpisize =
          Utilities::MPI::n_mpi_processes(this->mpi_communicator);

        AssertDimension(numcpus, mpisize);

        // Open file.
        MPI_Info info;
        int      ierr = MPI_Info_create(&info);
        AssertThrowMPI(ierr);

        const std::string fname_tria = filename + "_triangulation.data";

        MPI_File fh;
        ierr = MPI_File_open(this->mpi_communicator,
                             fname_tria.c_str(),
                             MPI_MODE_RDONLY,
                             info,
                             &fh);
        AssertThrowMPI(ierr);

        ierr = MPI_Info_free(&info);
        AssertThrowMPI(ierr);

        // Read offsets from file.
        std::uint64_t buffer_size;

        ierr = MPI_File_read_at(
          fh,
          myrank * sizeof(std::uint64_t),
          &buffer_size,
          1,
          Utilities::MPI::mpi_type_id_for_type<decltype(buffer_size)>,
          MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        std::uint64_t offset = 0;

        ierr = MPI_Exscan(
          &buffer_size,
          &offset,
          1,
          Utilities::MPI::mpi_type_id_for_type<decltype(buffer_size)>,
          MPI_SUM,
          this->mpi_communicator);
        AssertThrowMPI(ierr);

        // global position in file
        const std::uint64_t global_position =
          mpisize * sizeof(std::uint64_t) + offset;

        // Read buffers from file.
        std::vector<char> buffer(buffer_size);
        ierr = dealii::Utilities::MPI::LargeCount::File_read_at_c(
          fh,
          global_position,
          buffer.data(),
          buffer.size(), // local buffer
          MPI_CHAR,
          MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        ierr = MPI_File_close(&fh);
        AssertThrowMPI(ierr);

        auto construction_data = dealii::Utilities::template unpack<
          TriangulationDescription::Description<dim, spacedim>>(buffer, false);

        // WARNING: serialization cannot handle the MPI communicator
        // which is the reason why we have to set it here explicitly
        construction_data.comm = this->mpi_communicator;

        this->create_triangulation(construction_data);
      }

      // Compute global offset for each rank.
      unsigned int n_locally_owned_cells = this->n_locally_owned_active_cells();

      unsigned int global_first_cell = 0;

      int ierr = MPI_Exscan(&n_locally_owned_cells,
                            &global_first_cell,
                            1,
                            MPI_UNSIGNED,
                            MPI_SUM,
                            this->mpi_communicator);
      AssertThrowMPI(ierr);

      global_first_cell *= sizeof(unsigned int);

      Assert(this->n_global_active_cells() == n_global_active_cells,
             ExcMessage("Number of global active cells differ!"));

      // clear all of the callback data, as explained in the documentation of
      // register_data_attach()
      this->cell_attached_data.n_attached_data_sets = 0;
      this->cell_attached_data.n_attached_deserialize =
        attached_count_fixed + attached_count_variable;

      // Load attached cell data, if any was stored.
      this->load_attached_data(global_first_cell,
                               this->n_global_active_cells(),
                               this->n_locally_owned_active_cells(),
                               filename,
                               attached_count_fixed,
                               attached_count_variable);

      this->update_cell_relations();
      this->update_periodic_face_map();
      this->update_number_cache();
#else
      (void)filename;

      AssertThrow(false, ExcNeedsMPI());
#endif
    }



    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    void Triangulation<dim, spacedim>::update_number_cache()
    {
      dealii::parallel::TriangulationBase<dim, spacedim>::update_number_cache();

      // additionally update the number of global coarse cells
      types::coarse_cell_id number_of_global_coarse_cells = 0;

      for (const auto &cell : this->active_cell_iterators())
        if (!cell->is_artificial())
          number_of_global_coarse_cells =
            std::max(number_of_global_coarse_cells,
                     cell->id().get_coarse_cell_id());

      number_of_global_coarse_cells =
        Utilities::MPI::max(number_of_global_coarse_cells,
                            this->mpi_communicator) +
        1;

      this->number_cache.number_of_global_coarse_cells =
        number_of_global_coarse_cells;
    }


  } // namespace fullydistributed
} // namespace parallel



/*-------------- Explicit Instantiations -------------------------------*/
#include "distributed/fully_distributed_tria.inst"


DEAL_II_NAMESPACE_CLOSE
