// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_distributed_grid_tools_h
#define dealii_distributed_grid_tools_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/distributed/tria_base.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/optional.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filter/gzip.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{

  /**
   * This namespace defines parallel algorithms that operate on meshes.
   */
  namespace GridTools
  {
    /**
     * Exchange arbitrary data of type @p DataType provided by the function
     * objects from locally owned cells to ghost cells on other processors.
     *
     * After this call, you typically will have received data from @p unpack on
     * every ghost cell as it was given by @p pack on the owning processor.
     * Whether you do or do not receive information to @p unpack on a given
     * ghost cell depends on whether the @p pack function decided that
     * something needs to be sent. It does so using the boost::optional
     * mechanism: if the boost::optional return object of the @p pack
     * function is empty, then this implies that no data has to be sent for
     * the locally owned cell it was called on. In that case, @p unpack will
     * also not be called on the ghost cell that corresponds to it on the
     * receiving side. On the other hand, if the boost::optional object is
     * not empty, then the data stored within it will be sent to the received
     * and the @p unpack function called with it.
     *
     * @tparam DataType The type of the data to be communicated. It is assumed
     *   to be serializable by boost::serialization. In many cases, this
     *   data type can not be deduced by the compiler, e.g., if you provide
     *   lambda functions for the second and third argument
     *   to this function. In this case, you have to explicitly specify
     *   the @p DataType as a template argument to the function call.
     * @tparam MeshType The type of @p mesh.
     *
     * @param mesh A variable of a type that satisfies the requirements of the
     * @ref ConceptMeshType "MeshType concept".
     * @param pack The function that will be called on each locally owned cell
     *   that is a ghost cell somewhere else. As mentioned above, the function
     *   may return a regular data object of type @p DataType to indicate
     *   that data should be sent, or an empty
     *   <code>boost::optional@<DataType@></code> to indicate that nothing has
     *   to be sent for this cell.
     * @param unpack The function that will be called for each ghost cell
     *   for which data was sent, i.e., for which the @p pack function
     *   on the sending side returned a non-empty boost::optional object.
     *   The @p unpack function is then called with the data sent by the
     *   processor that owns that cell.
     *
     *
     * <h4> An example </h4>
     *
     * Here is an example that shows how this function is to be used
     * in a concrete context. It is taken from the code that makes
     * sure that the @p active_fe_index (a single unsigned integer) is
     * transported from locally owned cells where one can set it in
     * hp::DoFHandler objects, to the corresponding ghost cells on
     * other processors to ensure that one can query the right value
     * also on those processors:
     * @code
     *    auto pack
     *    = [] (const typename dealii::hp::DoFHandler<dim,spacedim>::active_cell_iterator &cell) -> unsigned int
     *    {
     *      return cell->active_fe_index();
     *    };
     *
     *    auto unpack
     *      = [] (const typename dealii::hp::DoFHandler<dim,spacedim>::active_cell_iterator &cell,
     *            const unsigned int                                                        &active_fe_index) -> void
     *    {
     *      cell->set_active_fe_index(active_fe_index);
     *    };
     *
     *    parallel::GridTools::exchange_cell_data_to_ghosts<unsigned int, dealii::hp::DoFHandler<dim,spacedim>>
     *        (dof_handler, pack, unpack);
     * @endcode
     *
     * You will notice that the @p pack lambda function returns an `unsigned int`,
     * not a `boost::optional<unsigned int>`. The former converts automatically
     * to the latter, implying that data will always be transported to the
     * other processor.
     *
     * (In reality, the @p unpack function needs to be a bit more
     * complicated because it is not allowed to call
     * DoFAccessor::set_active_fe_index() on ghost cells. Rather, the
     * @p unpack function directly accesses internal data structures. But
     * you get the idea -- the code could, just as well, have exchanged
     * material ids, user indices, boundary indictors, or any kind of other
     * data with similar calls as the ones above.)
     */
    template <typename DataType, typename MeshType>
    void
    exchange_cell_data_to_ghosts (const MeshType &mesh,
                                  std::function<boost::optional<DataType> (const typename MeshType::active_cell_iterator &)> pack,
                                  std::function<void (const typename MeshType::active_cell_iterator &, const DataType &)> unpack);
  }


}

#ifndef DOXYGEN

namespace parallel
{
  namespace GridTools
  {
    namespace internal
    {
      /**
       * A structure that allows the transfer of data of type T from one processor
       * to another. It corresponds to a packed buffer that stores a list of
       * cells and an array of type T.
       *
       * The vector @p data is the same size as @p cell_ids.
       */
      template <int dim, typename T>
      struct CellDataTransferBuffer
      {
        std::vector<CellId> cell_ids;
        std::vector<T> data;

        /**
         * Write the data of this object to a stream for the purpose of
         * serialization.
         */
        template <class Archive>
        void save (Archive &ar,
                   const unsigned int /*version*/) const
        {
          Assert(cell_ids.size() == data.size(), ExcInternalError());
          // archive the cellids in an efficient binary format
          const size_t n_cells = cell_ids.size();
          ar &n_cells;
          for (auto &it : cell_ids)
            {
              CellId::binary_type binary_cell_id = it.template to_binary<dim>();
              ar &binary_cell_id;
            }

          ar &data;
        }

        /**
         * Read the data of this object from a stream for the purpose of
         * serialization. Throw away the previous content.
         */
        template <class Archive>
        void load (Archive &ar,
                   const unsigned int /*version*/)
        {
          size_t n_cells;
          ar &n_cells;
          cell_ids.clear();
          cell_ids.reserve(n_cells);
          for (unsigned int c=0; c<n_cells; ++c)
            {
              CellId::binary_type value;
              ar &value;
              cell_ids.emplace_back(std::move(value));
            }
          ar &data;
        }

        BOOST_SERIALIZATION_SPLIT_MEMBER()


        /**
         * Pack the data that corresponds to this object into a buffer in
         * the form of a vector of chars and return it.
         */
        std::vector<char>
        pack_data () const
        {
          // set up a buffer and then use it as the target of a compressing
          // stream into which we serialize the current object
          std::vector<char> buffer;
          {
            boost::iostreams::filtering_ostream out;
            out.push(boost::iostreams::gzip_compressor
                     (boost::iostreams::gzip_params
                      (boost::iostreams::gzip::best_compression)));
            out.push(boost::iostreams::back_inserter(buffer));

            boost::archive::binary_oarchive archive(out);

            archive << *this;
            out.flush();
          }

          return buffer;
        }


        /**
         * Given a buffer in the form of an array of chars, unpack it and
         * restore the current object to the state that it was when
         * it was packed into said buffer by the pack_data() function.
         */
        void unpack_data (const std::vector<char> &buffer)
        {
          std::string decompressed_buffer;

          // first decompress the buffer
          {
            boost::iostreams::filtering_ostream decompressing_stream;
            decompressing_stream.push(boost::iostreams::gzip_decompressor());
            decompressing_stream.push(boost::iostreams::back_inserter(decompressed_buffer));

            decompressing_stream.write (&buffer[0], buffer.size());
          }

          // then restore the object from the buffer
          std::istringstream in(decompressed_buffer);
          boost::archive::binary_iarchive archive(in);

          archive >> *this;
        }
      };

    }


    template <typename DataType, typename MeshType>
    void
    exchange_cell_data_to_ghosts (const MeshType &mesh,
                                  std::function<boost::optional<DataType> (const typename MeshType::active_cell_iterator &)> pack,
                                  std::function<void (const typename MeshType::active_cell_iterator &, const DataType &)> unpack)
    {
#ifndef DEAL_II_WITH_MPI
      (void)mesh;
      (void)pack;
      (void)unpack;
      Assert(false, ExcMessage("parallel::GridTools::exchange_cell_data_to_ghosts() requires MPI."));
#else
      constexpr int dim = MeshType::dimension;
      constexpr int spacedim = MeshType::space_dimension;
      auto tria =
        static_cast<const parallel::Triangulation<dim, spacedim>*>(&mesh.get_triangulation());
      Assert (tria != nullptr,
              ExcMessage("The function exchange_cell_data_to_ghosts() only works with parallel triangulations."));

      // map neighbor_id -> data_buffer where we accumulate the data to send
      typedef std::map<dealii::types::subdomain_id, internal::CellDataTransferBuffer<dim, DataType> >
      DestinationToBufferMap;
      DestinationToBufferMap destination_to_data_buffer_map;

      std::map<unsigned int, std::set<dealii::types::subdomain_id> >
      vertices_with_ghost_neighbors = tria->compute_vertices_with_ghost_neighbors();

      for (auto cell : tria->active_cell_iterators())
        if (cell->is_locally_owned())
          {
            std::set<dealii::types::subdomain_id> send_to;
            for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
              {
                const std::map<unsigned int, std::set<dealii::types::subdomain_id> >::const_iterator
                neighbor_subdomains_of_vertex
                  = vertices_with_ghost_neighbors.find (cell->vertex_index(v));

                if (neighbor_subdomains_of_vertex ==
                    vertices_with_ghost_neighbors.end())
                  continue;

                Assert(neighbor_subdomains_of_vertex->second.size()!=0,
                       ExcInternalError());

                send_to.insert(neighbor_subdomains_of_vertex->second.begin(),
                               neighbor_subdomains_of_vertex->second.end());
              }

            if (send_to.size() > 0)
              {
                // this cell's data needs to be sent to someone
                typename MeshType::active_cell_iterator
                mesh_it (tria, cell->level(), cell->index(), &mesh);

                const boost::optional<DataType> data = pack(mesh_it);

                if (data)
                  {
                    const CellId cellid = cell->id();

                    for (auto it : send_to)
                      {
                        const dealii::types::subdomain_id subdomain = it;

                        // find the data buffer for proc "subdomain" if it exists
                        // or create an empty one otherwise
                        typename DestinationToBufferMap::iterator p
                          = destination_to_data_buffer_map.insert (std::make_pair(subdomain,
                                                                                  internal::CellDataTransferBuffer<dim, DataType>()))
                            .first;

                        p->second.cell_ids.emplace_back(cellid);
                        p->second.data.emplace_back(data.get());
                      }
                  }
              }
          }


      // 2. send our messages
      std::set<dealii::types::subdomain_id> ghost_owners = tria->ghost_owners();
      const unsigned int n_ghost_owners = ghost_owners.size();
      std::vector<std::vector<char> > sendbuffers (n_ghost_owners);
      std::vector<MPI_Request> requests (n_ghost_owners);

      unsigned int idx=0;
      for (auto it = ghost_owners.begin();
           it!=ghost_owners.end();
           ++it, ++idx)
        {
          internal::CellDataTransferBuffer<dim, DataType> &data = destination_to_data_buffer_map[*it];

          // pack all the data into the buffer for this recipient and send it.
          // keep data around till we can make sure that the packet has been
          // received
          sendbuffers[idx] = data.pack_data ();
          const int ierr = MPI_Isend(sendbuffers[idx].data(), sendbuffers[idx].size(),
                                     MPI_BYTE, *it,
                                     786, tria->get_communicator(), &requests[idx]);
          AssertThrowMPI(ierr);
        }

      // 3. receive messages
      std::vector<char> receive;
      for (unsigned int idx=0; idx<n_ghost_owners; ++idx)
        {
          MPI_Status status;
          int len;
          int ierr = MPI_Probe(MPI_ANY_SOURCE, 786, tria->get_communicator(), &status);
          AssertThrowMPI(ierr);
          ierr = MPI_Get_count(&status, MPI_BYTE, &len);
          AssertThrowMPI(ierr);

          receive.resize(len);

          char *ptr = &receive[0];
          ierr = MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                          tria->get_communicator(), &status);
          AssertThrowMPI(ierr);

          internal::CellDataTransferBuffer<dim, DataType> cellinfo;
          cellinfo.unpack_data(receive);

          DataType *data = cellinfo.data.data();
          for (unsigned int c=0; c<cellinfo.cell_ids.size(); ++c, ++data)
            {
              const typename Triangulation<dim,spacedim>::cell_iterator
              tria_cell = cellinfo.cell_ids[c].to_cell(*tria);

              const typename MeshType::active_cell_iterator
              cell (tria, tria_cell->level(), tria_cell->index(), &mesh);

              unpack(cell, *data);
            }
        }

      // make sure that all communication is finished
      // when we leave this function.
      if (requests.size())
        {
          const int ierr = MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }
#endif // DEAL_II_WITH_MPI
    }

  }
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
