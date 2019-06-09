
#include <RFAStDFT/block_csr_matrix.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector_operations_internal.h>
#include <numeric>
#include <set>

namespace RealFAStDFT
{

  void DoFInfo::clear()
  {
    row_starts.clear();
    dof_indices.clear();
    block_indices.clear();
  }



  void
  DoFInfo::initialize(const std::vector<unsigned int> &my_rows,
                      const std::shared_ptr<const BlockIndices> &row_blocks)
  {
    Assert(my_rows.size() > 0,
           ExcMessage("selected rows should be non-empty."));

    clear();

    using TUP = std::tuple<unsigned int, unsigned int, unsigned int>;
    std::vector<TUP> cell_dof_indices_local(my_rows.size());
    dof_indices.reserve(my_rows.size());

    // resizes and assign the default element
    row_starts.reserve(2);
    row_starts.push_back({0, 0});

    for (unsigned int i = 0; i < my_rows.size(); ++i)
      {
        const auto pair = row_blocks->global_to_local(my_rows[i]);

        TUP &tup = cell_dof_indices_local[i];
        std::get<0>(tup) = pair.first;
        std::get<1>(tup) = pair.second;
        std::get<2>(tup) = i;
      }

    process_cell_indices(cell_dof_indices_local);
    AssertDimension(dof_indices.size(), my_rows.size());
  }

  template <int dim, typename NumberType>
  void DoFInfo::initialize(
    const DoFHandler<dim> &dof_handler,
    const std::shared_ptr<const dealii::Utilities::MPI::Partitioner> &partitioner,
    const std::shared_ptr<const MatrixFree<dim,NumberType>> &matrix_free,
    const std::shared_ptr<const BlockIndices> &row_blocks,
    const std::vector<unsigned int> &cell_index_permutation)
  {
    const unsigned int n_subcells =
      VectorizedArray<NumberType>::n_array_elements;
    const unsigned int n_cells = matrix_free->n_macro_cells();

    const auto &fe = dof_handler.get_fe();
    AssertDimension(cell_index_permutation.size(), fe.dofs_per_cell);

    const parallel::Triangulation<dim> *tria =
      dynamic_cast<const parallel::Triangulation<dim> *>(
        &dof_handler.get_triangulation());

    Assert(tria != NULL, ExcNotImplemented());

    std::vector<types::global_dof_index> cell_dof_indices(fe.dofs_per_cell);

    // we will translate each local dof into:
    // 1. index within the cell
    // 2. block
    // 3. row within the block
    using TUP = std::tuple<unsigned int, unsigned int, unsigned int>;
    std::vector<TUP> cell_dof_indices_local(fe.dofs_per_cell);

    clear();
    dof_indices.reserve(tria->n_locally_owned_active_cells() *
                        fe.dofs_per_cell);

    // resizes and assign the default element
    row_starts.reserve(n_cells * n_subcells + 1);
    row_starts.push_back({0, 0});

    // loop through cells in the order consistent with MatrixFree class
    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        const auto n_filled = matrix_free->n_components_filled(cell);
        for (unsigned int subcell = 0; subcell < n_filled; ++subcell)
          {
            const auto cell_it =
              matrix_free->get_cell_iterator(cell, subcell);
            // get DoFs on this cell.
            cell_it->get_dof_indices(cell_dof_indices);
            for (unsigned int i = 0; i < cell_dof_indices.size(); ++i)
              {
                const auto dof =
                  partitioner->global_to_local(cell_dof_indices[i]);
                const auto pair = row_blocks->global_to_local(dof);

                TUP &tup = cell_dof_indices_local[i];
                std::get<0>(tup) = pair.first;
                std::get<1>(tup) = pair.second;
                std::get<2>(tup) = cell_index_permutation[i];
              }

            process_cell_indices(cell_dof_indices_local);
          }
        // fill-in non-filled subcells with something that shall never be used.
        // need to be careful here. After the loop above is done, we have correct
        // row_starts in slot one-after the filled. We need to copy it at the
        // beginning of the next SIMD cell chunk and fill in the rest on this chunk
        // with something meaningless
        if (n_filled < n_subcells)
          {
            // copy the last valid entry
            const auto last_valid = row_starts.back();
            for (unsigned int subcell = n_filled + 1; subcell < n_subcells; ++subcell)
              row_starts.push_back({numbers::invalid_unsigned_int,
                                    numbers::invalid_unsigned_int});

            row_starts.push_back(last_valid);
          }
      }

    // check the final size of row_starts
    AssertDimension(row_starts.size(), n_cells * n_subcells + 1);
  }



  void DoFInfo::process_cell_indices(
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
      &cell_dof_indices_local)
  {
    Assert(cell_dof_indices_local.size() > 0, ExcInternalError());
    // sort in lexicographic order (block, row, index)
    std::sort(cell_dof_indices_local.begin(), cell_dof_indices_local.end());

#ifdef DEBUG
    // in debug mode count the total number of elements we collect
    // from the vector of tuples
    unsigned int total_n_elements = 0;
#endif

    // put into our data structures and setup row_starts
    unsigned int current_block = std::get<0>(cell_dof_indices_local[0]);
    unsigned int n_elements = 0;
    unsigned int n_blocks = 0;
    for (const auto &tup : cell_dof_indices_local)
      {
        dof_indices.push_back({std::get<1>(tup), std::get<2>(tup)});
        const auto b = std::get<0>(tup);
        // flash stored data we see a new block:
        if (b != current_block)
          {
#ifdef DEBUG
            total_n_elements += n_elements;
#endif
            block_indices.push_back({current_block, n_elements});
            ++n_blocks;
            n_elements = 0;
            current_block = b;
          }

        ++n_elements;
      }

      // flash the last block:
#ifdef DEBUG
    total_n_elements += n_elements;
#endif
    block_indices.push_back({current_block, n_elements});
    ++n_blocks;

    AssertDimension(total_n_elements, cell_dof_indices_local.size());

    row_starts.push_back(
      {row_starts.back().first + cell_dof_indices_local.size(),
       row_starts.back().second + n_blocks});
  }

  const SparsityPatternStandard::size_type
    SparsityPatternStandard::invalid_entry = numbers::invalid_size_type;



  SparsityPatternStandard::SparsityPatternStandard() : SparsityPatternBase()
  {
    reinit(0, 0, 0);
  }



  void SparsityPatternStandard::reinit(
    const size_type m,
    const size_type n,
    const ArrayView<const unsigned int> &row_lengths)
  {
    AssertDimension(row_lengths.size(), m);

    rows = m;
    cols = n;

    // delete empty matrices
    if ((m == 0) || (n == 0))
      {
        rowstart.reset();
        colnums.reset();

        max_vec_len = max_dim = rows = cols = 0;
        // if dimension is zero: ignore max_per_row
        max_row_length = 0;
        compressed = false;

        return;
      }

    // find out how many entries we need in the @p{colnums} array. if this
    // number is larger than @p{max_vec_len}, then we will need to reallocate
    // memory
    //
    // note that the number of elements per row is bounded by the number of
    // columns
    //
    std::size_t vec_len = 0;
    for (size_type i = 0; i < m; ++i)
      vec_len += std::min(static_cast<size_type>(row_lengths[i]), n);

    // sometimes, no entries are requested in the matrix (this most often
    // happens when blocks in a block matrix are simply zero). in that case,
    // allocate exactly one element, to have a valid pointer to some memory
    if (vec_len == 0)
      {
        vec_len = 1;
        max_vec_len = vec_len;
        colnums = std_cxx14::make_unique<size_type[]>(max_vec_len);
      }

    max_row_length = (row_lengths.size() == 0
                        ? 0
                        : std::min(static_cast<size_type>(*std::max_element(
                                     row_lengths.begin(), row_lengths.end())),
                                   n));

    // allocate memory for the rowstart values, if necessary. even though we
    // re-set the pointers again immediately after deleting their old content,
    // set them to zero in between because the allocation might fail, in which
    // case we get an exception and the destructor of this object will be
    // called
    // -- where we look at the non-nullness of the (now invalid) pointer again
    // and try to delete the memory a second time.
    if (rows > max_dim)
      {
        max_dim = rows;
        rowstart = std_cxx14::make_unique<std::size_t[]>(max_dim + 1);
      }

    // allocate memory for the column numbers if necessary
    if (vec_len > max_vec_len)
      {
        max_vec_len = vec_len;
        colnums = std_cxx14::make_unique<size_type[]>(max_vec_len);
      }

    // set the rowstart array
    rowstart[0] = 0;
    for (size_type i = 1; i <= rows; ++i)
      rowstart[i] = rowstart[i - 1] +
                    std::min(static_cast<size_type>(row_lengths[i - 1]), n);
    Assert((rowstart[rows] == vec_len) ||
             ((vec_len == 1) && (rowstart[rows] == 0)),
           ExcInternalError());

    // preset the column numbers by a value indicating it is not in use
    std::fill_n(colnums.get(), vec_len, invalid_entry);

    compressed = false;
  }

  SparsityPatternStandard::size_type SparsityPatternStandard::
  operator()(const size_type i, const size_type j) const
  {
    Assert((rowstart != nullptr) && (colnums != nullptr), ExcEmptyObject());
    Assert(i < rows, ExcIndexRange(i, 0, rows));
    Assert(j < cols, ExcIndexRange(j, 0, cols));
    Assert(compressed, ExcNotCompressed());

    // let's see whether there is something in this line
    if (rowstart[i] == rowstart[i + 1])
      return invalid_entry;

    // all other entries are sorted, so we can use a binary search algorithm
    //
    // note that the entries are only sorted upon compression, so this would
    // fail for non-compressed sparsity patterns; however, that is why the
    // Assertion is at the top of this function, so it may not be called for
    // noncompressed structures.
    const size_type *const p =
      dealii::Utilities::lower_bound<const size_type *>(
        &colnums[rowstart[i]], &colnums[rowstart[i + 1]], j);
    if ((p != &colnums[rowstart[i + 1]]) && (*p == j))
      return (p - colnums.get());
    else
      return invalid_entry;
  }

  void SparsityPatternStandard::copy_from(const DynamicSparsityPattern &dsp)
  {
    std::vector<unsigned int> row_lengths(dsp.n_rows());
    for (size_type i = 0; i < dsp.n_rows(); ++i)
      row_lengths[i] = dsp.row_length(i);

    reinit(dsp.n_rows(), dsp.n_cols(), row_lengths);

    if (n_rows() != 0 && n_cols() != 0)
      for (size_type row = 0; row < dsp.n_rows(); ++row)
        {
          size_type *cols = &colnums[rowstart[row]];
          const unsigned int row_length = dsp.row_length(row);
          for (unsigned int index = 0; index < row_length; ++index)
            {
              const size_type col = dsp.column_number(row, index);
              *cols++ = col;
            }
        }

    compressed = true;
  }


  //---------------------------- BCSR ---------------------------------

  template <typename NumberType>
  std::size_t BlockCSRMatrix<NumberType>::memory_consumption() const
  {
    std::size_t memory =
      sp.memory_consumption() + values.memory_consumption() +
      MemoryConsumption::memory_consumption(data_start) +
      diagonal.memory_consumption() +
      MemoryConsumption::memory_consumption(diagonal_data_start) +
      diagonal.memory_consumption() +
      MemoryConsumption::memory_consumption(ghost_targets) +
      MemoryConsumption::memory_consumption(import_targets) +
      MemoryConsumption::memory_consumption(import_indices) +
      MemoryConsumption::memory_consumption(ghost_data);

    if (row_blocks.use_count() > 0)
      memory += row_blocks->memory_consumption() / row_blocks.use_count();

    if (col_blocks.use_count() > 0)
      memory += col_blocks->memory_consumption() / col_blocks.use_count();

    if (partitioner.use_count() > 0)
      memory += partitioner->memory_consumption() / partitioner.use_count();

    if (n_import_data > 0)
      memory += sizeof(NumberType) * static_cast<std::size_t>(n_import_data);

    return memory;
  }

  template <typename NumberType>
  void BlockCSRMatrix<NumberType>::print(std::ostream &s,
                                         const unsigned int w,
                                         const unsigned int p) const
  {
    // save the state of out stream
    const std::streamsize old_precision = s.precision(p);
    const std::streamsize old_width = s.width(w);

    // to make the vector write out all the information in order, use as
    // many barriers as there are processors and start writing when it's our
    // turn
#ifdef DEAL_II_WITH_MPI
    if (partitioner->n_mpi_processes() > 1)
      for (unsigned int i = 0; i < partitioner->this_mpi_process(); i++)
        {
          const int ierr = MPI_Barrier(partitioner->get_mpi_communicator());
          AssertThrowMPI(ierr);
        }
#endif

    s << "Process    #" << partitioner->this_mpi_process() << std::endl
      << "Rows:       " << this->m() << std::endl
      << "Columns:    " << this->n() << std::endl
      << "Local range: [" << owned_row_start << ", "
      << owned_row_start + n_owned_rows << ")" << std::endl;

    auto print_row = [&](const size_type r) {
      const auto &M = row_blocks->block_size(r);
      std::stringstream rows_stream;
      if (r < n_owned_row_blocks)
        {
          const auto start = owned_row_start + row_blocks->block_start(r);
          rows_stream << " [" << start << ", " << start + M - 1 << "]";
        }
      else
        {
          const auto &data = ghost_data[r - n_owned_row_blocks];
          const auto &start = data.first;
          for (const auto &range : data.second)
            rows_stream << " [" << start + range.first << ", "
                        << start + range.second - 1 << "]";
        }

      for (auto it = this->begin_local(r); it != this->end_local(r); ++it)
        {
          const auto &c = it->column();
          const auto &N = col_blocks->block_size(c);
          s << "Block (" << r << "," << c << ") " << M << "x" << N << " "
            << rows_stream.str() << " x [" << col_blocks->block_start(c) << ", "
            << col_blocks->block_start(c) + N - 1 << "]" << std::endl;
          for (unsigned int i = 0; i < M; ++i)
            {
              for (unsigned int j = 0; j < N; ++j)
                {
                  s.width(w);
                  s.precision(p);
                  s << *(it->data() +
                         BlockCSRMatrix<NumberType>::local_index(i, j, M, N));
                }
              s << std::endl;
            }
        }
    };

    for (size_type r = 0; r < n_owned_row_blocks; ++r)
      print_row(r);

    if (matrix_is_ghosted)
      {
        s << "Ghost blocks:" << std::endl;
        for (size_type r = n_owned_row_blocks; r < sp.n_rows(); ++r)
          print_row(r);
      }

#ifdef DEAL_II_WITH_MPI
      if (partitioner->n_mpi_processes() > 1)
        {
          int ierr = MPI_Barrier(partitioner->get_mpi_communicator());
          AssertThrowMPI(ierr);

          for (unsigned int i = partitioner->this_mpi_process() + 1;
               i < partitioner->n_mpi_processes();
               i++)
            {
              ierr = MPI_Barrier(partitioner->get_mpi_communicator());
              AssertThrowMPI(ierr);
            }
        }
#endif

    // reset output format
    s.precision(old_precision);
    s.width(old_width);
  }


  namespace internal
  {
    struct BlockData
    {
      unsigned int size;
      types::global_dof_index global_block_start;
      std::vector<std::pair<unsigned int, unsigned int>> ranges;
      std::vector<types::global_dof_index> columns;

      template <class Archive>
      void serialize(Archive &ar, const unsigned int /*version*/)
      {
        ar &size;
        ar &global_block_start;
        ar &ranges;
        ar &columns;
      }
    };
  } // namespace internal



  template <typename NumberType>
  void BlockCSRMatrix<NumberType>::setup_ghosts(
    const DynamicSparsityPattern &locally_owned_sparsity)
  {
    // now we need to setup blocks and sparsity for ghosts row blocks or rows.
    // go through the MPI::Partitioner::import_targets() and
    // MPI::Partitioner::import_indices() to figure out where and what shall
    // we send:
    const std::vector<std::pair<unsigned int, unsigned int>>
      &import_indices_part = partitioner->import_indices();
    const std::vector<std::pair<unsigned int, unsigned int>>
      &import_targets_part = partitioner->import_targets();

    // before modifying row blocks according to sizes on ther processors, record
    // the data we need
    n_owned_rows =
      block_partitioner ? row_blocks->total_size() : partitioner->local_size();

    if (block_partitioner)
      {
        // if we are using block partitioner, we need to calculate our starting
        // point based on local sizes in all cores.
        const std::vector<unsigned int> local_sizes =
          dealii::Utilities::MPI::all_gather(
            partitioner->get_mpi_communicator(), n_owned_rows);

        // now set my start of the local range according to the rank
        owned_row_start = std::accumulate(
          local_sizes.begin(),
          local_sizes.begin() + partitioner->this_mpi_process(),
          types::global_dof_index(0));
      }
    else
      {
        owned_row_start = partitioner->local_range().first;

        // rework import_indices into our own data structure

        // will do split of [b,e) recursively, just in case there are more
        // than two regions inside.
        // last bool parameter forces addition of new range regardless of
        // the state of this->import_indices . We have to force this
        // at the root of recursion when targets switch
        std::function<void(
          const unsigned int beg, const unsigned int end, const bool force_new)>
          add_import;

        add_import = [&](const unsigned int beg,
                         const unsigned int end,
                         const bool force_new) -> void {
          const auto block_index = this->row_blocks->global_to_local(beg);
          const auto &b = block_index.first;
          const auto &start = this->row_blocks->block_start(b);
          Assert(beg >= start, ExcInternalError());
          const auto block_end = start + this->row_blocks->block_size(b);

          const auto range_end =
            end <= block_end ? end - start : block_end - start;

          if (this->import_indices.empty() ||
              this->import_indices.back().first != b || force_new)
            {
              this->import_indices.push_back({b, {{beg - start, range_end}}});
            }
          else
            {
              AssertDimension(this->import_indices.back().first, b);
              this->import_indices.back().second.push_back(
                {beg - start, range_end});
            }

          if (end > block_end)
            add_import(block_end, end, false);
        };

        auto target = import_targets_part.begin();
        unsigned int total_count = 0;
        bool force_new = false;
        for (auto ind : import_indices_part)
          {
            add_import(ind.first, ind.second, force_new);

            total_count += (ind.second - ind.first);
            if (total_count == target->second)
              {
                total_count = 0;
                ++target;
                force_new = true;
              }
            else
              {
                force_new = false;
              }
          }

#ifdef DEBUG
        // sanity check:
        // we should end up having the same total number of elements
        unsigned int part_sum = 0;
        for (auto &ind : import_targets_part)
          part_sum += ind.second;

        unsigned int my_sum = 0;
        for (auto &pair : import_indices)
          for (auto &range : pair.second)
            my_sum += range.second - range.first;

        AssertDimension(my_sum, part_sum);
#endif
      }

    // we will go through the import indices (that are ghosts elsewhere)
    // and gather data from our local sparsity...
    // we need to have 3 things:
    // 1. global row that defines start of ghost block
    // 2. rows in this block relative to start:  [b1,e1), [b2,e2),...   => block
    // size
    // 3. non-zero columns
    // This is exactly what internal::BlockData is for.

    // That's what we shall send to each target.
    std::map<unsigned int, std::vector<internal::BlockData>> objects_to_send;

    auto target = import_targets_part.begin();
    unsigned int total_count = 0;
    // while traversing the data, also setup our own import_targets
    unsigned int this_import_targets = 0;

    if (block_partitioner)
      // in case of block partitioner we can directly use import_indices_part
      {
        for (const auto &ind : import_indices_part)
          {
            for (unsigned int br = ind.first; br < ind.second; ++br)
              {
                const auto &block_size = row_blocks->block_size(br);
                internal::BlockData data;
                data.global_block_start =
                  owned_row_start + row_blocks->block_start(br);
                // FIXME: can keep ranges empty, already have all info
                data.ranges.push_back({0, block_size});
                data.size = block_size;

                for (auto r = locally_owned_sparsity.begin(br);
                     r != locally_owned_sparsity.end(br);
                     ++r)
                  {
                    this_import_targets +=
                      block_size * internal::padded_size<NumberType>(
                                     col_blocks->block_size(r->column()));
                    data.columns.push_back(r->column());
                  }

                objects_to_send[target->first].push_back(data);
              }

            // increase the counter to know when we have all we need to send
            total_count += (ind.second - ind.first);

            if (total_count == target->second)
              {
                // push back collected import size and reset counter
                import_targets.push_back(this_import_targets);
                this_import_targets = 0;
                // need to switch target and move on
                total_count = 0;
                // go to the next target
                ++target;
              }
          } // loop over import_indices_part
      }
    else
      // otherwise use our reworked import_indices
      {
        for (const auto &ind : this->import_indices)
          {
            const unsigned int &br = ind.first;
            internal::BlockData data;
            data.ranges = ind.second;
            data.global_block_start = owned_row_start +
                                      row_blocks->block_start(br);
            unsigned int block_size = 0;
            for (const auto &el : data.ranges)
              block_size += el.second - el.first;

            data.size = block_size;

            for (auto r = locally_owned_sparsity.begin(br);
                 r != locally_owned_sparsity.end(br);
                 ++r)
              {
                this_import_targets +=
                  block_size * internal::padded_size<NumberType>(
                                 col_blocks->block_size(r->column()));
                data.columns.push_back(r->column());
              }

            objects_to_send[target->first].push_back(data);

            // increase the counter to know when we have all we need to send
            total_count += block_size;

            if (total_count == target->second)
              {
                // push back collected import size and reset counter
                import_targets.push_back(this_import_targets);
                this_import_targets = 0;
                // need to switch target and move on
                total_count = 0;
                // go to the next target
                ++target;
              }
          } // loop over import_indices
      }

    // once we through all import indices, target should point to the end
    Assert(target == import_targets_part.end(), ExcInternalError());

    // own import_targets should be consistent with partitioner
    AssertDimension(import_targets.size(), import_targets_part.size());

    // now exchange the data between MPI cores:
    const auto objects_received = dealii::Utilities::MPI::some_to_some(
      partitioner->get_mpi_communicator(), objects_to_send);

    AssertDimension(objects_received.size(),
                    partitioner->ghost_targets().size());

    // now we know sparsity on ghost rows or row blocks
    // add it to DynamicSparsityPattern and finally copy to our local
    // SparsityPattern object

    // we shall know how many ghost rows we will add to initialize
    // DynamicSparsityPattern. for block partitioner that's exactly how many
    // ghost blocks we have. the row partitioner case is more delicate.
    // generally we have complete freedom to block ghost rows however we want as
    // they are only relevant for matrix-free operators. We can do a block of
    // size 1 for each ghost row or make a union for ghost rows coming from a
    // single process. in order to avoid sending too much data and simplify
    // writing into auxiliary arrays during compress() and update_ghosts() we
    // will block according to the blocking on the owning process. this is
    // already reflected in the way we store import_indices.
    unsigned int extra_rows = partitioner->n_ghost_indices();
    if (!block_partitioner)
      {
        extra_rows = 0;
        for (auto &obj : objects_received)
          extra_rows += obj.second.size();
      }

    DynamicSparsityPattern dsp(locally_owned_sparsity.n_rows() + extra_rows,
                               locally_owned_sparsity.n_cols(),
                               locally_owned_sparsity.row_index_set());
    for (DynamicSparsityPattern::size_type r = 0;
         r < locally_owned_sparsity.n_rows();
         ++r)
      for (auto it = locally_owned_sparsity.begin(r);
           it != locally_owned_sparsity.end(r);
           ++it)
        dsp.add(r, it->column());

    // now extend with ghosts
    // note that the ghosts are arranged in the increasing MPI rank, i.e.
    // all ghost elements from proc 0 (if any), then all from proc 1 (if any),
    // etc... so the way we traverse received objects is consistent with ghost
    // indexing
    unsigned int ghost_block_number = 0;
    for (auto obj : objects_received)
      {
        unsigned int this_target_size = 0;
        for (auto &data : obj.second)
          {
            row_blocks->push_back(data.size);
            dsp.add_entries(n_owned_row_blocks + ghost_block_number,
                            data.columns.begin(),
                            data.columns.end(),
                            true);

            for (const auto &c : data.columns)
              this_target_size += data.size * internal::padded_size<NumberType>(
                                                col_blocks->block_size(c));

            ghost_data.push_back({data.global_block_start, data.ranges});

            ++ghost_block_number;
          }

        ghost_targets.push_back(this_target_size);
      }

    AssertDimension(ghost_block_number, extra_rows);

    sp.copy_from(dsp);

    AssertDimension(sp.n_rows(), row_blocks->size());
    // make sure we collected sizes of ghost correctly,
    // size of import_targets is checked above
    AssertDimension(ghost_targets.size(), partitioner->ghost_targets().size());

    n_import_data =
      std::accumulate(import_targets.begin(), import_targets.end(), 0);
  }



  template <typename NumberType>
  void BlockCSRMatrix<NumberType>::reinit(
    const DynamicSparsityPattern &dsp,
    const std::shared_ptr<const BlockIndices> &row_blocks_,
    const std::shared_ptr<const BlockIndices> &col_blocks_,
    const std::shared_ptr<const dealii::Utilities::MPI::Partitioner> row_partitioner_,
    const bool symmetric_)
  {
    // FIXME: we need to check that the sparsity pattern was symmetrized
    symmetric = symmetric_;
    row_blocks = std::make_shared<BlockIndices>(*row_blocks_);
    col_blocks = col_blocks_;
    partitioner = row_partitioner_;

    thread_loop_partitioner =
        std::make_shared<::dealii::parallel::internal::TBBPartitioner>();

    Assert((partitioner->local_size() == row_blocks->total_size()) ||
             (partitioner->local_size() == row_blocks->size()),
           ExcMessage("Local size of partitioner <" +
                      std::to_string(partitioner->local_size()) +
                      "> should be consistent with either "
                      "local number of rows <" +
                      std::to_string(row_blocks->total_size()) +
                      "> or local number of row blocks <" +
                      std::to_string(row_blocks->size()) + ">."));

    block_partitioner = row_blocks->size() == partitioner->local_size();

    n_owned_row_blocks = row_blocks->size();
    n_row_blocks = dealii::Utilities::MPI::sum(n_owned_row_blocks, partitioner->get_mpi_communicator());
    n_rows = dealii::Utilities::MPI::sum(row_blocks->total_size(), partitioner->get_mpi_communicator());

    AssertDimension(row_blocks->size(), dsp.n_rows());
    AssertDimension(col_blocks->size(), dsp.n_cols());

    setup_ghosts(dsp);
    matrix_is_ghosted = false;

    // do not reallocate import_data directly, but only upon request. It
    // is only used as temporary storage for compress() and
    // update_ghost_values, and we might have matrices where we never
    // call these methods and hence do not need to have the storage.
    import_data.reset();

    Assert(!symmetric || (row_blocks->total_size() == col_blocks->total_size()),
           ExcMessage("Symmetric matrices should be square"));

    Assert(!symmetric || (*row_blocks == *col_blocks),
           ExcMessage(
             "Symmetric matrices should have same rows and columns blocks"));

    // calculate the size we need to store all blocks
    size_type new_size = 0;
    data_start.resize(sp.n_nonzero_elements()+1, 0);
    auto data_start_it = ++data_start.begin();
    for (size_type r = 0; r < sp.n_rows(); ++r)
      {
        const auto & M = row_blocks->block_size(r);
        for (auto it = sp.begin(r); it != sp.end(r); ++it, ++data_start_it)
          {
            new_size += M * internal::padded_size<NumberType>(
                              col_blocks->block_size(it->column()));
            *data_start_it = new_size;
          }
      }

    Assert (data_start_it == data_start.end(), ExcInternalError());

    values.resize_fast(new_size);
    this->operator=(NumberType());
  }


  template <typename NumberType>
  BlockCSRMatrix<NumberType>::BlockCSRMatrix()
    : symmetric(false),
      n_row_blocks(0),
      n_owned_row_blocks(0),
      n_owned_rows(0),
      owned_row_start(0),
      n_rows(0),
      block_partitioner(false),
      matrix_is_ghosted(false),
      n_import_data(0)
  {
  }



  template <typename NumberType>
  BlockCSRMatrix<NumberType>::~BlockCSRMatrix()
  {
  }


  //--------------------------
  //----- element access -----
  //--------------------------
  template <typename NumberType>
  const NumberType &
  BlockCSRMatrix<NumberType>::operator()(const size_type i,
                                         const size_type j) const
  {
    Assert(i < this->m(), ExcMessage("AssertRange"));
    Assert(j < this->n(), ExcMessage("AssertRange"));

    // need to translate global row into local one
    Assert (!block_partitioner || partitioner->n_mpi_processes() == 1, ExcNotImplemented());
    const size_type i_ = block_partitioner ? i : partitioner->global_to_local(i);

    const std::pair<unsigned int, unsigned int> row_pair = row_blocks->global_to_local(i_);
    const std::pair<unsigned int, unsigned int> col_pair = col_blocks->global_to_local(j);

    const auto index = sp(row_pair.first, col_pair.first);

    if (index != SparsityPattern::invalid_entry)
      {
        return *(&values[0] + data_start[index] +
                 BlockCSRMatrix<NumberType>::local_index(
                   row_pair.second,
                   col_pair.second,
                   row_blocks->block_size(row_pair.first),
                   col_blocks->block_size(col_pair.first)));
      }
    else
      {
        Assert(false,
               ExcInvalidIndex(i,j));
        // return something to make compiler happy
        return values[0];
      }
  }



  template <typename NumberType>
  NumberType
  BlockCSRMatrix<NumberType>::el(const size_type i, const size_type j) const
  {
    AssertIndexRange(i, this->m());
    AssertIndexRange(j, this->n());

    Assert (!block_partitioner || partitioner->n_mpi_processes() == 1, ExcNotImplemented());
    const size_type i_ = block_partitioner ? i : partitioner->global_to_local(i);
    return local_el(i_, j);
  }



  template <typename NumberType>
  NumberType
  BlockCSRMatrix<NumberType>::local_el(const unsigned int i_, const size_type j) const
  {
    AssertIndexRange(i_, n_owned_rows);
    Assert(j < this->n(), ExcMessage("AssertRange"));

    const std::pair<unsigned int, unsigned int> row_pair = row_blocks->global_to_local(i_);
    const std::pair<unsigned int, unsigned int> col_pair = col_blocks->global_to_local(j);

    const auto index = sp(row_pair.first, col_pair.first);

    if (index != SparsityPattern::invalid_entry)
      {
        return *(&values[0] + data_start[index] +
                 BlockCSRMatrix<NumberType>::local_index(
                   row_pair.second,
                   col_pair.second,
                   row_blocks->block_size(row_pair.first),
                   col_blocks->block_size(col_pair.first)));
      }
    else
      {
        return NumberType();
      }
  }



  template <typename NumberType>
  NumberType &
  BlockCSRMatrix<NumberType>::operator()(const size_type i, const size_type j)
  {
    Assert(i < this->m(), ExcMessage("AssertRange"));
    Assert(j < this->n(), ExcMessage("AssertRange"));

    Assert (!block_partitioner || partitioner->n_mpi_processes() == 1, ExcNotImplemented());
    const size_type i_ = block_partitioner ? i : partitioner->global_to_local(i);

    const std::pair<unsigned int, unsigned int> row_pair = row_blocks->global_to_local(i_);
    const std::pair<unsigned int, unsigned int> col_pair = col_blocks->global_to_local(j);

    const auto index = sp(row_pair.first, col_pair.first);
    if (index != SparsityPattern::invalid_entry)
      {
        return *(&values[0] + data_start[index] +
                 BlockCSRMatrix<NumberType>::local_index(
                   row_pair.second,
                   col_pair.second,
                   row_blocks->block_size(row_pair.first),
                   col_blocks->block_size(col_pair.first)));
      }
    else
      {
        Assert(false,
               ExcInvalidIndex(i,j));
        // return something to make compiler happy
        return values[0];
      }
  }


  //--------------------------
  //----- mmult, Tmmult ------
  //--------------------------

  template <typename NumberType>
  void BlockCSRMatrix<NumberType>::mmult(BlockCSRMatrix<NumberType> &C,
                                         const BlockCSRMatrix<NumberType> &B,
                                         const bool rebuild_sparsity_C) const
  {
    Assert(!sp.empty(), ExcNotInitialized());
    Assert(!B.sp.empty(), ExcNotInitialized());
    Assert(!C.sp.empty(), ExcNotInitialized());

    const auto &sp_A = sp;
    const auto &sp_B = B.sp;

    Assert(sp_A.n_cols() == B.n_row_blocks,
           ExcDimensionMismatch(sp_A.n_cols(), B.n_row_blocks));

    // clear previous content of C
    if (rebuild_sparsity_C == true)
      {
        C.clear();

        // create a sparsity pattern for the matrix C.
        DynamicSparsityPattern dsp;
        Assert(false, ExcNotImplemented());
        (void)sp_A;
        (void)sp_B;
        // dsp.compute_mmult_pattern(sp_A, sp_B);

        // reinit matrix C from that information
        // C.reinit(dsp, this->get_row_blocks(), B.get_col_blocks());
      }

    Assert(C.sp.n_rows() == sp_A.n_rows(),
           ExcDimensionMismatch(C.sp.n_rows(), sp_A.n_rows()));
    Assert(C.sp.n_cols() == sp_B.n_cols(),
           ExcDimensionMismatch(C.sp.n_cols(), sp_B.n_cols()));

    // now compute the actual entries: a matrix-matrix product involves three
    // nested loops. One over the rows of A, for each row we then loop over all
    // the columns, and then we need to multiply each element with all the
    // elements in that row in B.
    // C_{ij} = A_{ik} B_{kj}

    C = value_type(0.);

    // FIXME: can overlap computation and communication if we
    // limit iteration to owned blocks
    B.update_ghost_values();

    // outer loop over i:
    for (unsigned int i = 0; i < this->n_owned_row_blocks; ++i)
      {
        auto A_ik = this->begin_local(i);
        const auto last_A_ik = this->end_local(i);

        // FIXME: keep a temp vector<NumberType> storage to sum
        // multiplications from A_{ik} B_{kj}

        const auto first_C_ij = C.begin_local(i);
        const auto last_C_ij = C.end_local(i);

        // middle loop over k
        while (A_ik != last_A_ik)
          {
            const auto k = A_ik->column();
            auto B_kj = B.begin(k);
            const auto last_B_kj = B.end(k);

            auto C_ij = first_C_ij;

            // inner loop over j
            while (B_kj != last_B_kj)
              {
                const auto j = B_kj->column();
                // advance C iterator to have correct column.
                // we may have sparsity of C which is not enough to store
                // the product exactly. For example, A and C have the same sparsity.
                // Therefore we search for C_ij only until its column is less than j:
                while (C_ij != last_C_ij && C_ij->column() < j)
                  ++C_ij;

                // we assume that if the sparsity of C is not enough to store
                // full product, at least we can store  C_ij = A_ij B_jj
                Assert((j != k) || (C_ij->column() == j), ExcInternalError());

                if (C_ij == last_C_ij || C_ij->column() != j)
                  {
                    ++B_kj;
                    continue;
                  }

                // actual multiplication: C_ij += A_ik * B_kj
                // when doing row-blocks, we essentially do (B^T A^T)^T = (A*B)
                const types::blas_int mm   = B.get_col_blocks()->block_size(j);
                const types::blas_int mm_p = internal::padded_size<NumberType>(mm); // stride of left
                const types::blas_int nn   = row_blocks->block_size(i);
                const types::blas_int kk   = col_blocks->block_size(k);
                const types::blas_int kk_p = internal::padded_size<NumberType>(kk); // stride of right
                const auto left            = B_kj->data();
                const auto right           = A_ik->data();
                const NumberType alpha     = 1.;
                const NumberType beta      = 1.; // adding

                gemm("N",
                    "N",
                    &mm,
                    &nn,
                    &kk,
                    &alpha,
                    left,
                    &mm_p,
                    right,
                    &kk_p,
                    &beta,
                    C_ij->data(),
                    &mm_p);

                ++B_kj;
              } // loop over j

            ++A_ik;
          } // loop over k
      }     // loop over i
  }



  template <typename NumberType>
  void BlockCSRMatrix<NumberType>::Tmmult(BlockCSRMatrix<NumberType> &C,
                                          const BlockCSRMatrix<NumberType> &B,
                                          const bool rebuild_sparsity_C) const
  {
    Assert(!sp.empty(), ExcNotInitialized());
    Assert(!B.sp.empty(), ExcNotInitialized());
    Assert(!C.sp.empty(), ExcNotInitialized());

    const auto &sp_A = sp;
    const auto &sp_B = B.sp;

    Assert(sp_A.n_rows() == sp_B.n_rows(),
           ExcDimensionMismatch(sp_A.n_rows(), sp_B.n_rows()));

    // clear previous content of C
    if (rebuild_sparsity_C == true)
      {
        C.clear();

        // create a sparsity pattern for the matrix.
        DynamicSparsityPattern dsp;
        Assert(false, ExcNotImplemented());
        (void)sp_A;
        (void)sp_B;
        // dsp.compute_Tmmult_pattern(sp_A, sp_B);

        // reinit matrix C from that information
        // C.reinit(dsp, this->get_col_blocks(), B.get_col_blocks(), C.is_symmetric());
      }

    const auto &sp_C = C.sp;

    Assert(C.n_row_blocks == sp_A.n_cols(),
           ExcDimensionMismatch(C.n_row_blocks, sp_A.n_cols()));
    Assert(sp_C.n_cols() == sp_B.n_cols(),
           ExcDimensionMismatch(sp_C.n_cols(), sp_B.n_cols()));

    // now compute the actual entries: a matrix-matrix product involves three
    // nested loops. One over the rows of A, for each row we then loop over all
    // the columns, and then we need to multiply each element with all the
    // elements in that row in B.
    // C_{kl} = A_{ik} B_{il}

    C = value_type(0.);
    const bool C_symmetric = C.is_symmetric();

    // outer most loop over i:
    for (unsigned int i = 0; i < this->n_owned_row_blocks; ++i)
      {
        auto A_ik = this->begin_local(i);
        const auto last_A_ik = this->end_local(i);

        const auto first_B_il = B.begin_local(i);
        const auto last_B_il = B.end_local(i);

        // middle loop over k:
        while (A_ik != last_A_ik)
          {
            const auto k = A_ik->column();
            auto C_kl = C.begin(k);
            const auto last_C_kl = C.end(k);

            // innermost loop over l:
            auto B_il = first_B_il;
            while (B_il != last_B_il)
              {
                const auto l = B_il->column();
                // find the right C_kl for given k and l.
                while (C_kl->column() != l)
                  {
                    ++C_kl;
                    (void)last_C_kl;
                    Assert(C_kl != last_C_kl,
                           ExcMessage("Could not find column " +
                                      std::to_string(l) + " in row " +
                                      std::to_string(k) + " of C matrix."));
                  }

                Assert(C_kl->column() == l, ExcInternalError());

                // now the product itself:
                if (!C_symmetric || (k >= l))
                  {
                    const NumberType    alpha = 1.;
                    const NumberType    beta  = 1.; // adding
                    // C_kl += A_ik^T * B_il
                    // when doing row-blocks, we essentially do (B^T A)^T = (A^T*B)
                    const types::blas_int mm   = B.get_col_blocks()->block_size(l);
                    const types::blas_int mm_p = internal::padded_size<NumberType>(mm);
                    const types::blas_int nn   = this->col_blocks->block_size(k);
                    const types::blas_int nn_p = internal::padded_size<NumberType>(nn);
                    const types::blas_int kk   = row_blocks->block_size(i);
                    gemm("N",
                         "T",
                         &mm,
                         &nn,
                         &kk,
                         &alpha,
                         B_il->data(),
                         &mm_p,
                         A_ik->data(),
                         &nn_p,
                         &beta,
                         C_kl->data(),
                         &mm_p);
                  }

                ++B_il;
              } // loop over l

            ++A_ik;
          } // loop over k
      }     // over i

    C.compress(VectorOperation::add);

    // fill in lower part
    if (C_symmetric)
      {
        for (types::global_dof_index i = 0; i < sp_C.n_rows(); ++i)
          {
            const types::blas_int n   = C.get_row_blocks()->block_size(i);
            const types::blas_int n_p = internal::padded_size<NumberType>(n);
            (void)n_p;
            const auto end = C.end(i);
            for (auto it = C.begin(i); it != end; ++it)
              if (i > it->column())
                {
                  const auto & j = it->column();
                  auto src = it->data();

                  // dst = srt^T
                  const types::blas_int m    = C.get_col_blocks()->block_size(j);
                  const types::blas_int m_p  = internal::padded_size<NumberType>(m);
                  (void)m_p;

                  // need to find C(j,i)
                  // FIXME: do binary search?
                  auto dst_it = C.begin(j);
                  while (dst_it->column() != i)
                    {
                      ++dst_it;
                      Assert(dst_it != C.end(j),
                             ExcMessage("Could not find column " +
                                        std::to_string(i) + " in row " +
                                        std::to_string(j) + " of C matrix."));

                    }

                  Assert (dst_it->column() == i, ExcInternalError());
                  auto dst = dst_it->data();

#ifdef DEAL_II_LAPACK_WITH_MKL
                  const NumberType one = 1.;
                  omatcopy('C', 'C', n, m, one, src, n_p, dst, m_p);
#else
                  for (types::blas_int ii = 0; ii < m; ++ii)
                    for (types::blas_int jj = 0; jj < n; ++jj)
                      *(dst + BlockCSRMatrix<NumberType>::local_index(ii,jj,m,n)) =
                        *(src + BlockCSRMatrix<NumberType>::local_index(jj,ii,n,m));
#endif
                }
          }
      }
  }



  template <typename NumberType>
  NumberType
  BlockCSRMatrix<NumberType>::Tr_Tmmult(const BlockCSRMatrix<NumberType> &B) const
  {
    Assert(!sp.empty(), ExcNotInitialized());
    Assert(!B.sp.empty(), ExcNotInitialized());


    Assert(sp.n_rows() == B.sp.n_rows(),
           ExcDimensionMismatch(sp.n_rows(), B.sp.n_rows()));

    // the result should be a square matrix:
    Assert(sp.n_cols() == B.sp.n_cols(),
           ExcDimensionMismatch(sp.n_cols(), B.sp.n_cols()));

    // we also require that column blocking is the same. Otherwise C_kk = A_ik
    // B_ik are not necessarily square matrices

    Assert(*col_blocks == (*B.get_col_blocks()),
          ExcMessage("Column blocking should be the same."));

    // now compute the actual entries: a matrix-matrix product involves three
    // nested loops. One over the rows of A, for each row we then loop over all
    // the columns, and then we need to multiply each element with all the
    // elements in that row in B.
    // C_{kl} = A_{ik} B_{il}
    //
    // In order to calculate the trace, we are only interested in k==l:
    // C_{kk} = A_{ik} B_{ik}

    // resize the vector to store diagonal blocks:

    // FIXME: we only need C to store local results, i.e. for
    // columns k that appear in locally owned parts of A and B.
    size_type C_size = 0;
    diagonal_data_start.resize(col_blocks->size() + 1);
    diagonal_data_start[0] = 0;
    {
      auto it = ++diagonal_data_start.begin();
      for (unsigned int i = 0; i < col_blocks->size(); ++i, ++it)
        {
          C_size +=
            col_blocks->block_size(i) *
            internal::padded_size<NumberType>(col_blocks->block_size(i));
          (*it) = C_size;
        }
      Assert (it == diagonal_data_start.end(), ExcInternalError());
    }

    Assert (C_size > 0, ExcInternalError());
    diagonal.resize_fast(C_size);
    diagonal.fill();

    // outer most loop over i:
    for (unsigned int i = 0; i < this->n_owned_row_blocks; ++i)
      {
        const auto last_A_ik = this->end_local(i);
        auto A_ik = this->begin_local(i);

        auto B_il = B.begin_local(i);
        const auto last_B_il = B.end_local(i);

        // we now sit on i-th row both in A and B, so we need to get all
        // matching non-zero columns. We can figure this out by iterating
        // both at the same time knowing that indices are ordered.

        // enter the while() loop only if both rows are non-empty
        // we shall break off the loop manually when we rich the end of rows
        const bool both_valid = (A_ik != last_A_ik && B_il != last_B_il);
        while (both_valid)
          {
            // if after incrementing below we hit the last one, break
            if (A_ik == last_A_ik)
              break;

            // if k > l, iterate l until we hit l>=k or end of row in B ==> done
            while (B_il != last_B_il && A_ik->column() > B_il->column())
              ++B_il;

            // if we hit the last one in Bi, break:
            if (B_il == last_B_il)
              break;

            Assert(B_il->column() >= A_ik->column(),
                   ExcInternalError());

            // if l > k, iterate k until we hit k>=l or end of row in A ==> done
            while (A_ik != last_A_ik && B_il->column() > A_ik->column())
              ++A_ik;

            // if we hit the last one in Ai, break
            if (A_ik == last_A_ik)
              break;

            Assert(A_ik->column() >= B_il->column(),
                   ExcInternalError());

            // at this point both should be valid iterators
            Assert(A_ik != last_A_ik && B_il != last_B_il,
                   ExcInternalError());

            // if columns match do element multiplication
            if (A_ik->column() == B_il->column())
              {
                // diagonal[A_ik->column()] += A_ik^T B_il
                // when doing row-blocks, we essentially do (B^T A)^T = (A^T*B)
                const types::blas_int mm   = B.get_col_blocks()->block_size(B_il->column());
                const types::blas_int mm_p = internal::padded_size<NumberType>(mm); // stride of B_il
                const types::blas_int nn   = this->col_blocks->block_size(A_ik->column());
                const types::blas_int nn_p = internal::padded_size<NumberType>(nn); // stride of A_ik
                const types::blas_int kk   = row_blocks->block_size(i);
                const NumberType      alpha = 1.;
                const NumberType      beta  = 1.; // adding
                gemm("N",
                     "T",
                     &mm,
                     &nn,
                     &kk,
                     &alpha,
                     B_il->data(),
                     &mm_p,
                     A_ik->data(),
                     &nn_p,
                     &beta,
                     diagonal.data() + diagonal_data_start[A_ik->column()],
                     &mm_p);

                // and move both iterators:
                ++B_il;
                ++A_ik;
                // if columns don't match, we will be moving iterators at the
                // beginning of the while loop
              }
          } // loop over k and l
      }     // over i

    // now go through diagonal and calculate trace
    value_type res = value_type(0.);
    {
      auto data = &diagonal[0];
      for (unsigned int i = 0; i < col_blocks->size(); ++i)
        {
          const auto & M = col_blocks->block_size(i);
          for (unsigned int ii = 0; ii < M; ++ii)
            res+= *(data + BlockCSRMatrix<NumberType>::local_index(
                                   ii, ii, M, M)
                                   );

          data+= M * internal::padded_size<NumberType>(M);
        }
    }

    return dealii::Utilities::MPI::sum(res, this->partitioner->get_mpi_communicator());
  }



  //--------------------------
  //----- ROWACCESSOR --------
  //--------------------------

  namespace BlockCSRMatrixIterators
  {

    template <typename NumberType>
    RowsAccessor<NumberType, true>::RowsAccessor(
      typename RowsAccessorBase<NumberType, true>::MatrixType *matrix,
      const std::vector<unsigned int> &local_rows)
      : RowsAccessorBase<NumberType, true>(matrix)
    {
      this->reinit(local_rows);
    }



    template <typename NumberType>
    RowsAccessor<NumberType, true>::RowsAccessor(
      typename RowsAccessorBase<NumberType, true>::MatrixType *matrix)
      : RowsAccessorBase<NumberType, true>(matrix)
    {
    }



    template <typename NumberType>
    RowsAccessor<NumberType, false>::RowsAccessor(
      typename RowsAccessorBase<NumberType, false>::MatrixType *matrix,
      const std::vector<unsigned int> &local_rows)
      : RowsAccessorBase<NumberType, false>(matrix)
    {
      this->reinit(local_rows);
    }



    template <typename NumberType>
    RowsAccessor<NumberType, false>::RowsAccessor(
      typename RowsAccessorBase<NumberType, false>::MatrixType *matrix)
      : RowsAccessorBase<NumberType, false>(matrix)
    {
    }



    template <typename NumberType, bool Constness>
    RowsAccessorBase<NumberType, Constness>::RowsAccessorBase(
      MatrixType *matrix)
      : matrix(matrix),
        block_col(numbers::invalid_unsigned_int),
        col_within_block(numbers::invalid_unsigned_int)
    {
    }



    template <typename NumberType, bool Constness>
    void RowsAccessorBase<NumberType, Constness>::reinit(
      const std::vector<unsigned int> &local_rows)
    {
      // first, go through all requested rows, get block numbers
      // and make them unique.
      std::vector<unsigned int> my_row_blocks(local_rows.size());
      auto it = my_row_blocks.begin();
      for (auto i : local_rows)
        {
          (*it) = matrix->row_blocks->global_to_local(i).first;
          ++it;
        }

      std::sort(my_row_blocks.begin(), my_row_blocks.end());
      my_row_blocks.erase(
        std::unique(my_row_blocks.begin(), my_row_blocks.end()),
        my_row_blocks.end());

      reinit_blocks(my_row_blocks);
    }



    template <typename NumberType, bool Constness>
    void RowsAccessorBase<NumberType, Constness>::reinit_blocks(
      const std::vector<unsigned int> &active_row_blocks_)
    {
      Assert(active_row_blocks_.size() > 0,
             ExcMessage("Provide at least one row block for access."));
      Assert(active_row_blocks_.size() ==
               std::set<unsigned int>(active_row_blocks_.begin(),
                                      active_row_blocks_.end())
                 .size(),
             ExcMessage("Provided row blocks are not unique"));

#ifdef DEBUG
      {
        std::vector<unsigned int> tmp(active_row_blocks_);
        std::sort(tmp.begin(), tmp.end());
        Assert(active_row_blocks_ == tmp,
               ExcMessage("Provided row blocks are not sorted"));
      }
#endif

      active_row_blocks.resize(active_row_blocks_.size());

      col_within_block = 0;

      // Find the minimum column block and set iterators
      block_col = std::numeric_limits<unsigned int>::max();
      for (unsigned int i = 0; i < active_row_blocks_.size(); ++i)
        {
          const auto &b = active_row_blocks_[i];
          DataType &d = active_row_blocks[i];
          d.block_size = matrix->row_blocks->block_size(b);
          d.block_start = matrix->row_blocks->block_start(b);
          d.it = matrix->begin_local(b);
          d.end = matrix->end_local(b);
          if (d.it != d.end)
            block_col = std::min(block_col, d.it->column());
        }

      // if all turned out to be empty
      Assert(block_col != std::numeric_limits<unsigned int>::max(),
             ExcNotImplemented());

      // set active flags and pointer to data:
      move_iterators();
    }



    //-----------------------------
    //---- RowsBlockAccessor ------
    //-----------------------------

    template <typename NumberType, bool Constness>
    RowsBlockAccessor<NumberType, Constness>::RowsBlockAccessor(
      MatrixType *matrix, const DoFInfo & dof_info)
      : matrix(matrix),
        col_block(numbers::invalid_dof_index),
        col_block_size(numbers::invalid_dof_index),
        col_block_size_p(numbers::invalid_dof_index),
        dof_info(dof_info)
    {
    }



    template <typename NumberType, bool Constness>
    void
    RowsBlockAccessor<NumberType, Constness>::clear()
    {
      col_block_size = numbers::invalid_dof_index;
      col_block_size_p = numbers::invalid_dof_index;
      col_block = numbers::invalid_dof_index;
      row_blocks.clear();
    }



    template <typename NumberType, bool Constness>
    types::global_dof_index RowsBlockAccessor<NumberType, Constness>::advance()
    {
      Assert(col_block != numbers::invalid_dof_index,
             ExcMessage("The accessor reached the last non-empty column blocks "
                        "and can not be advanced further"));

      // 1. increment all non-end iterators with column() < col_block + 1
      {
        const auto block_col_plus_one = col_block + 1;
        for (auto &el : row_blocks)
          if (el.it != el.end &&
              el.it->column() < block_col_plus_one)
            ++el.it;
      }

      // 2. get the minimum column among non-end iterators
      col_block = std::numeric_limits<types::global_dof_index>::max();
      for (auto &el : row_blocks)
        if (el.it != el.end)
          col_block = std::min(col_block, el.it->column());

      // 3. set pointer and active flags
      set_pointers_and_active_flags();
      // and return the column block
      return col_block;
    }



    template <typename NumberType, bool Constness>
    types::global_dof_index RowsBlockAccessor<NumberType, Constness>::reinit(
      const unsigned int cell,
      const unsigned int subcell)
    {
      // don't forget to clear internal data
      clear();

      Assert(dof_info.row_starts.size() > 0, ExcNotInitialized());
      AssertIndexRange(cell * VectorizedArray<NumberType>::n_array_elements +
                         subcell,
                       dof_info.row_starts.size() - 1);

      const auto &this_row_starts =
        dof_info
          .row_starts[cell * VectorizedArray<NumberType>::n_array_elements +
                      subcell];
      // we construct dof_info in such a way that row_starts+1 is valid
      // even when a subcell is the last filled
      const auto &next_row_starts =
        dof_info
          .row_starts[cell * VectorizedArray<NumberType>::n_array_elements +
                      subcell + 1];

      const auto total_blocks = next_row_starts.second - this_row_starts.second;

      row_blocks.reserve(total_blocks);

      // setup iterators, as well as figure out the column and whether or not
      // each RowBlock is active.
      col_block = std::numeric_limits<types::global_dof_index>::max();
      RowBlock block;
      unsigned int row_counter = 0;
      for (unsigned int b = 0; b < total_blocks; ++b)
        {
          AssertIndexRange(this_row_starts.second + b,
                           dof_info.block_indices.size());
          const auto &block_row =
            dof_info.block_indices[this_row_starts.second + b].first;
          const auto &n_rows =
            dof_info.block_indices[this_row_starts.second + b].second;

          Assert(n_rows > 0, ExcInternalError());

          block.block_row = block_row;
          block.it = matrix->begin_local(block_row);
          block.end = matrix->end_local(block_row);
          block.active = false;

          AssertIndexRange(this_row_starts.first + row_counter + n_rows - 1,
                           dof_info.dof_indices.size());
          block.dof_view.reinit(
            &dof_info.dof_indices[this_row_starts.first + row_counter], n_rows);

          if (block.it != block.end)
            col_block = std::min(col_block, block.it->column());

          row_blocks.push_back(block);
          row_counter+= n_rows;
        }

      set_pointers_and_active_flags();

      return col_block;
    }



    template <typename NumberType, bool Constness>
    void
    RowsBlockAccessor<NumberType, Constness>::set_pointers_and_active_flags()
    {
      // if there are no non-end iterators, the col_block
      // will be numeric_limits<types::global_dof_index>::max()
      // in this case reset data structures appropriately and return
      if (col_block >= matrix->col_blocks->size())
        {
          clear();
          return;
        }

      // if we have non-end iterators, set active flag
      // and pointer to row accordingly
      col_block_size = matrix->col_blocks->block_size(col_block);
      col_block_size_p = internal::padded_size<NumberType>(col_block_size);
      stride = col_block_size_p/VectorizedArray<NumberType>::n_array_elements;
      for (auto &block : row_blocks)
        {
          block.active =
            (block.it != block.end && block.it->column() == col_block);

          if (block.active)
            block.pointer =
              reinterpret_cast<vectorized_pointer>(block.it->data());
        }
    }

  } // namespace BlockCSRMatrixIterators



  //---------------------------------------
  //---- BCSR parallel data exchange ------
  //---------------------------------------


  template <typename NumberType>
  void
  BlockCSRMatrix<NumberType>::clear_mpi_requests ()
  {
#ifdef DEAL_II_WITH_MPI
    for (size_type j=0; j<compress_requests.size(); j++)
      {
        const int ierr = MPI_Request_free(&compress_requests[j]);
        AssertThrowMPI(ierr);
      }
    compress_requests.clear();
    for (size_type j=0; j<update_ghost_values_requests.size(); j++)
      {
        const int ierr = MPI_Request_free(&update_ghost_values_requests[j]);
        AssertThrowMPI(ierr);
      }
    update_ghost_values_requests.clear();
#endif
  }

  template <typename NumberType>
  void
  BlockCSRMatrix<NumberType>::compress(::dealii::VectorOperation::values operation)
  {
    compress_start (0, operation);
    compress_finish(operation);
  }



  template <typename NumberType>
  void
  BlockCSRMatrix<NumberType>::update_ghost_values() const
  {
    update_ghost_values_start ();
    update_ghost_values_finish ();
  }



  template <typename NumberType>
  void
  BlockCSRMatrix<NumberType>::compress_start(
        const unsigned int                communication_channel,
        ::dealii::VectorOperation::values operation)
  {
    (void)communication_channel;
    (void)operation;
    Assert(matrix_is_ghosted == false,
           ExcMessage("Cannot call compress() on a ghosted matrix"));

#ifdef DEAL_II_WITH_MPI
    // nothing to do for insert (only need to zero ghost entries in
    // compress_finish()). in debug mode we want to check consistency
    // of the inserted data, therefore the communication is still
    // initialized. Having different code in debug and optimized mode is
    // somewhat dangerous, but it really saves communication so it seems
    // still worthwhile.
#ifndef DEBUG
    if (operation == VectorOperation::insert)
      return;
#endif

    const dealii::Utilities::MPI::Partitioner &part = *partitioner;

    // nothing to do when we neither have import
    // nor ghost indices.
    if (part.n_ghost_indices() == 0 && part.n_import_indices() == 0)
      return;

    // make this function thread safe
    std::lock_guard<std::mutex> lock(mutex);

    const unsigned int n_import_targets = part.import_targets().size();
    const unsigned int n_ghost_targets = part.ghost_targets().size();

    // Need to send and receive the data. Use non-blocking communication,
    // where it is generally less overhead to first initiate the receive and
    // then actually send the data
    if (compress_requests.size() == 0)
      {
        // set channels in different range from update_ghost_values channels
        const unsigned int channel = communication_channel + 400;
        unsigned int current_index_start = 0;
        compress_requests.resize(n_import_targets + n_ghost_targets);

        // allocate import_data in case it is not set up yet
        if (import_data == nullptr && n_import_data > 0)
          import_data = std_cxx14::make_unique<NumberType[]>(n_import_data);

        // initiate the recieve
        AssertDimension(import_targets.size(), n_import_targets);
        for (unsigned int i = 0; i < n_import_targets; i++)
          {
            AssertThrow(
              static_cast<size_type>(import_targets[i]) *
                  sizeof(NumberType) <
                static_cast<size_type>(std::numeric_limits<int>::max()),
              ExcMessage(
                "Index overflow: Maximum message size in MPI is 2GB. "
                "The number of ghost entries times the size of 'Number' "
                "exceeds this value. This is not supported."));
            const int ierr =
              MPI_Recv_init(&import_data[current_index_start], // recieve buffer
                            import_targets[i] * sizeof(NumberType), // count
                            MPI_BYTE,                               // type
                            part.import_targets()[i].first,         // source
                            part.import_targets()[i].first +
                              part.n_mpi_processes() * channel, // tag
                            part.get_mpi_communicator(),        // communicator
                            &compress_requests[i]);             // request
            AssertThrowMPI(ierr);
            current_index_start += import_targets[i];
          }
        AssertDimension(current_index_start, n_import_data);

        // index of the first block in the first ghost row:
        // note that sp.rowstart is of size sp.n_rows() + 1
        AssertIndexRange(n_owned_row_blocks, sp.n_rows() + 1);
        AssertIndexRange(sp.rowstart[n_owned_row_blocks], data_start.size());
        current_index_start = data_start[sp.rowstart[n_owned_row_blocks]];
        AssertDimension(ghost_targets.size(), n_ghost_targets);
        for (unsigned int i = 0; i < n_ghost_targets; i++)
          {
            AssertThrow(
              static_cast<size_type>(ghost_targets[i]) *
                  sizeof(NumberType) <
                static_cast<size_type>(std::numeric_limits<int>::max()),
              ExcMessage(
                "Index overflow: Maximum message size in MPI is 2GB. "
                "The number of ghost entries times the size of 'Number' "
                "exceeds this value. This is not supported."));
            const int ierr = MPI_Send_init(
              &this->values[current_index_start],    // send buffer
              ghost_targets[i] * sizeof(NumberType), // count
              MPI_BYTE,                              // type
              part.ghost_targets()[i].first,         // destination
              part.this_mpi_process() + part.n_mpi_processes() * channel, // tag
              part.get_mpi_communicator(),               // communicator
              &compress_requests[n_import_targets + i]); // request
            AssertThrowMPI(ierr);
            current_index_start += ghost_targets[i];
          }
        AssertDimension(current_index_start, values.size());
      } // compress_requests.size() == 0

    AssertDimension(n_import_targets + n_ghost_targets,
                    compress_requests.size());
    if (compress_requests.size() > 0)
      {
        const int ierr =
          MPI_Startall(compress_requests.size(), &compress_requests[0]);
        AssertThrowMPI(ierr);
      }
#endif
  }



  template <typename NumberType>
  void
  BlockCSRMatrix<NumberType>::compress_finish(::dealii::VectorOperation::values operation)
  {
#ifdef DEAL_II_WITH_MPI
    // in optimized mode, no communication was started, so leave the
    // function directly (and only clear ghosts)
#ifndef DEBUG
    if (operation == VectorOperation::insert)
      {
        zero_out_ghosts();
        return;
      }
#endif

    const dealii::Utilities::MPI::Partitioner &part = *partitioner;

    // nothing to do when we neither have import nor ghost indices.
    if (part.n_ghost_indices() == 0 && part.n_import_indices() == 0)
      return;

    // make this function thread safe
    std::lock_guard<std::mutex> lock(mutex);

    const unsigned int n_import_targets = part.import_targets().size();
    const unsigned int n_ghost_targets = part.ghost_targets().size();

    if (operation != dealii::VectorOperation::insert)
      AssertDimension(n_ghost_targets + n_import_targets,
                      compress_requests.size());

    // first wait for the receive to complete
    if (compress_requests.size() > 0 && n_import_targets > 0)
      {
        const int ierr = MPI_Waitall(
          n_import_targets, &compress_requests[0], MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);


        // If the operation is no insertion, add the imported data to the
        // local values. For insert, nothing is done here (but in debug mode
        // we assert that the specified value is either zero or matches with
        // the ones already present
        if (operation != dealii::VectorOperation::insert)
          BlockCSRMatrix<NumberType>::read_write_import<false,true>(*this);
        else
          BlockCSRMatrix<NumberType>::read_write_import<false,false>(*this);
      }

    if (compress_requests.size() > 0 && n_ghost_targets > 0)
      {
        const int ierr = MPI_Waitall(n_ghost_targets,
                                     &compress_requests[n_import_targets],
                                     MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);
      }
    else
      AssertDimension(part.n_ghost_indices(), 0);

    zero_out_ghosts();
#else
    (void)operation;
#endif
  }


  namespace internal
  {
    // SFINAE tricks.
    // at compile time std::enable_if<True>::type won't result in substitution
    // failure and this function will be taken for this case. whereas
    // Addition==false will result in a substitution failure
    // (std::enable_if<false>::type does not exist) and therefore
    // the function below will be used.
    template <bool Addition, typename NumberType>
    inline typename std::enable_if<Addition>::type
    modifier(NumberType *import, NumberType *data, const unsigned int n = 1)
    {
      // addition for non-const data
      dealii::internal::VectorOperations::Vectorization_add_v<NumberType> vector_add(
        data, import);
      vector_add(0, n);
    }

    template <bool Addition, typename NumberType>
    inline typename std::enable_if<!Addition>::type
    modifier(NumberType *import, NumberType *data, const unsigned int n = 1)
    {
      // insertion to non-const data.
      // only check consistency (this code is not executed in release mode)
      for (unsigned int i = 0; i < n; ++i)
        if (*(import + i) == NumberType() ||
            std::abs(*(data + i) - *(import + i)) <=
              std::abs(*(data + i) * 1000. *
                       std::numeric_limits<NumberType>::epsilon()))
          {
          }
        else
          {
            std::stringstream ss;
            ss << "Called compress(VectorOperation::insert), but"
               << " the element received from a remote processor, value "
               << std::setprecision(16) << *(import + i)
               << ", does not match with the value " << std::setprecision(16)
               << *(data + i) << " on the owner processor";

            Assert(false, ExcMessage(ss.str()));
          }
    }

    // finally, this one will be used for Constness==True as
    // others lack `const` for the data (second) argument.
    template <bool Addition, typename NumberType>
    inline void modifier(NumberType *import,
                         const NumberType *data,
                         const unsigned int n = 1)
    {
      // setting import to data
      dealii::internal::VectorOperations::Vector_copy<NumberType, NumberType>
        vector_copy(data, import);
      vector_copy(0, n);
    }
  }

  template <typename NumberType>
  template <bool Constness, bool Addition>
  void BlockCSRMatrix<NumberType>::read_write_import(
    typename std::conditional<Constness,
                              const BlockCSRMatrix<NumberType> &,
                              BlockCSRMatrix<NumberType> &>::type mtr)
  {
    const dealii::Utilities::MPI::Partitioner &part = *mtr.partitioner;
    Assert(mtr.n_import_data != 0, ExcInternalError());
    auto read_write_position = mtr.import_data.get();
    auto i_target = part.import_targets().begin();
    Assert(part.import_indices().size() > 0, ExcInternalError());

    // for block partitioner things simplify quite a bit as we
    // can take the whole block row of data and add or set it
    if (mtr.block_partitioner)
      {
        for (auto &i_ind : part.import_indices())
          // i_ind = [begin,end)
          for (unsigned int r = i_ind.first; r < i_ind.second; ++r)
            {
              const auto M = mtr.row_blocks->block_size(r);
              auto it = mtr.begin_local(r);
              const auto end = mtr.end_local(r);
              for (; it != end; ++it)
                {
                  const auto N = mtr.col_blocks->block_size(it->column());
                  const auto N_p = internal::padded_size<NumberType>(N);
                  // now we can simply copy/add M*N from/to
                  // it->data() and read_write_position
                  internal::modifier<Addition>(
                    read_write_position, it->data(), M * N_p);

                  // go to the next slot of import data
                  read_write_position += M * N_p;
                }
            }

        // should be at the end of import data
        // Assert(read_write_position == mtr.import_data.end(), ExcInternalError());
        // and Partitioner pointer should point at the end

        return;
      }

    // for non-block partitioner, life is harder. Essentially we need
    // to translate row import indices into block. Within each block row
    // take a couple of rows according to import indices and put it where
    // appropriate.

    unsigned int total_count = 0;
    unsigned int current_block = numbers::invalid_unsigned_int;
    std::vector<unsigned int> block_indices;

    // function to write current collection of rows into import data
    auto write_into = [&]() {
      const auto m = block_indices.size();
      const auto block_start = mtr.row_blocks->block_start(current_block);

      Assert (m > 0, ExcInternalError());

      auto it = mtr.begin_local(current_block);
      const auto end = mtr.end_local(current_block);
      for (; it != end; ++it)
        {
          const auto N = mtr.col_blocks->block_size(it->column());
          const auto N_p = internal::padded_size<NumberType>(N);
          for (unsigned int i = 0; i < block_indices.size(); ++i)
            {
              Assert(block_indices[i] >= block_start, ExcInternalError());
              const unsigned int ii = block_indices[i] - block_start;
              Assert(ii < mtr.row_blocks->block_size(current_block),
                     ExcInternalError());
              auto rw = read_write_position + N_p * i;
              auto d = it->data() + N_p * ii;
              internal::modifier<Addition>(rw, d, N_p);
            }
          read_write_position += m * N_p;
        }
      block_indices.resize(0);
    };

    for (auto &i_ind : part.import_indices())
      // i_ind = [b,e)
      {
        for (unsigned int j = i_ind.first; j < i_ind.second; ++j)
          {
            const auto j_block = mtr.row_blocks->global_to_local(j).first;
            // if the row block change, ready to write the results
            // and reset the block_indices.
            // after we switch targets, current_block has wrong value, skip
            // comparison if total_count is zero
            if (current_block != numbers::invalid_unsigned_int && j_block != current_block)
              write_into();

            block_indices.push_back(j);
            current_block = j_block;
          }

        // increase the counter to know when we process the next target
        total_count += (i_ind.second - i_ind.first);

        // if we are done on this target, need to write data
        if (total_count == i_target->second)
          {
            total_count = 0;
            write_into();
            ++i_target;
            current_block = numbers::invalid_unsigned_int;
          }
      }
    // once we through, we should not have any remaining row block indices
    // to be processed
    Assert(block_indices.size() == 0, ExcInternalError());
    // should be at the end of import data
    // Assert(read_write_position == mtr.import_data.end(), ExcInternalError());
    // and Partitioner pointer should point at the end
    Assert(i_target == part.import_targets().end(), ExcInternalError());
  }



  template <typename NumberType>
  void
  BlockCSRMatrix<NumberType>::update_ghost_values_start(
      const unsigned int communication_channel) const
  {
#ifdef DEAL_II_WITH_MPI
    const dealii::Utilities::MPI::Partitioner &part = *partitioner;

    // nothing to do when we neither have import nor ghost indices.
    if (part.n_ghost_indices()==0 && part.n_import_indices()==0)
      return;

    // make this function thread safe
    std::lock_guard<std::mutex> lock(mutex);

    const unsigned int n_import_targets = part.import_targets().size();
    const unsigned int n_ghost_targets = part.ghost_targets().size();

    // Need to send and receive the data. Use non-blocking communication,
    // where it is generally less overhead to first initiate the receive and
    // then actually send the data
    if (update_ghost_values_requests.size() == 0)
      {
        update_ghost_values_requests.resize (n_import_targets+n_ghost_targets);

        AssertDimension (ghost_targets.size(), n_ghost_targets);

        // index of the first block in the first ghost row:
        // note that rowstart is of size  n_rows + 1.
        AssertIndexRange(n_owned_row_blocks, sp.n_rows() + 1);
        AssertIndexRange(sp.rowstart[n_owned_row_blocks], data_start.size());
        size_type current_index_start =
          data_start[sp.rowstart[n_owned_row_blocks]];
        // recall that part.ghost_targets()[i].first contains MPI rank and our
        // ghost_targets[i] gives the number of elements we expect.
        for (unsigned int i = 0; i < n_ghost_targets; i++)
          {
              // allow writing into ghost indices even though we are in a
              // const function
              const int ierr = MPI_Recv_init (const_cast<NumberType *>(&values[current_index_start]),  // where to write
                                              ghost_targets[i]*sizeof(NumberType),                     // count
                                              MPI_BYTE,                                                // type
                                              part.ghost_targets()[i].first,                           // source
                                              part.ghost_targets()[i].first +
                                              communication_channel*part.n_mpi_processes(),            // tag
                                              part.get_mpi_communicator(),                             // communicator
                                              &update_ghost_values_requests[i]);                       // request
              AssertThrowMPI (ierr);
              current_index_start += ghost_targets[i];
          }
        AssertDimension (current_index_start,
                         values.size());

        // allocate import_data in case it is not set up yet
        if (import_data == nullptr && n_import_data > 0)
          import_data = std_cxx14::make_unique<NumberType[]>(n_import_data);

        // initiate MPI_Send
        current_index_start = 0;
        AssertDimension(n_import_targets, import_targets.size());
        for (unsigned int i=0; i<n_import_targets; i++)
            {
              const int ierr = MPI_Send_init (&import_data[current_index_start],                  // buffer
                                              import_targets[i]*sizeof(NumberType),               // count
                                              MPI_BYTE,                                           // type
                                              part.import_targets()[i].first,                     // destination
                                              part.this_mpi_process() +
                                              communication_channel*part.n_mpi_processes(),       // tag
                                              part.get_mpi_communicator(),                        // communicator
                                              &update_ghost_values_requests[n_ghost_targets+i]);  // request
              AssertThrowMPI (ierr);
              current_index_start += import_targets[i];
            }
        AssertDimension(n_import_data, current_index_start);
      }

    // copy the data that is actually to be send to the import_data field
    // note that this has to be done in exactly the same was the target process
    // expects this to arrive into ghost part of the matrix. Recall that we block
    // ghost DoFs according to blocks they belong on their owning process.
    if (part.n_import_indices() > 0)
      {
        BlockCSRMatrix<NumberType>::read_write_import<true>(*this);
      }

    AssertDimension (n_import_targets+n_ghost_targets,
                     update_ghost_values_requests.size());
    if (update_ghost_values_requests.size() > 0)
      {
        const int ierr = MPI_Startall(update_ghost_values_requests.size(),
                                      &update_ghost_values_requests[0]);
        AssertThrowMPI(ierr);
      }
#else
    (void)communication_channel;
#endif
  }



  template <typename NumberType>
  void
  BlockCSRMatrix<NumberType>::update_ghost_values_finish() const
  {
#ifdef DEAL_II_WITH_MPI
      // wait for both sends and receives to complete, even though only
      // receives are really necessary. this gives (much) better performance
      AssertDimension (partitioner->ghost_targets().size() +
                       partitioner->import_targets().size(),
                       update_ghost_values_requests.size());
      if (update_ghost_values_requests.size() > 0)
        {
          // make this function thread safe
          std::lock_guard<std::mutex> lock(mutex);

          const int ierr = MPI_Waitall (update_ghost_values_requests.size(),
                                        &update_ghost_values_requests[0],
                                        MPI_STATUSES_IGNORE);
          AssertThrowMPI (ierr);
        }
#endif
      matrix_is_ghosted = true;
  }



  template <typename NumberType>
  void
  BlockCSRMatrix<NumberType>::zero_out_ghosts() const
  {
    const auto index_start = data_start[sp.rowstart[n_owned_row_blocks]];
    Assert(values.size() >= index_start, ExcInternalError());

    if (values.size() > index_start)
      std::fill_n(const_cast<NumberType *>(&values[index_start]),
                  values.size() - index_start,
                  NumberType());
    matrix_is_ghosted = false;
  }



  template <typename NumberType>
  NumberType BlockCSRMatrix<NumberType>::frobenius_norm() const
  {
    const auto N = data_start[sp.rowstart[n_owned_row_blocks]];
    NumberType res = 0.;
    for (unsigned int i = 0; i < N; ++i)
      res += values[i] * values[i];
    return std::sqrt(
      dealii::Utilities::MPI::sum(res, partitioner->get_mpi_communicator()));
  }



  template <typename NumberType>
  BlockCSRMatrix<NumberType>&
  BlockCSRMatrix<NumberType>::operator-=(const BlockCSRMatrix<NumberType> &B)
  {
    // FIXME:
    // Assert(sp == B.sp, ExcMessage("Sparsity patterns should match."));
    Assert(n_owned_row_blocks == B.n_owned_row_blocks,
           ExcMessage("Number of owned row blocks should match."));
    for (unsigned int r = 0; r < n_owned_row_blocks; ++r)
      AssertDimension(row_blocks->block_size(r), B.row_blocks->block_size(r));

    Assert(col_blocks == B.col_blocks,
           ExcMessage("Column blocks should match"));

    const dealii::types::blas_int N = data_start[sp.rowstart[n_owned_row_blocks]];
    constexpr dealii::types::blas_int one_int = 1;
    constexpr NumberType mult = -1.;

    axpy(&N, &mult, &B.values[0], &one_int, &values[0], &one_int);

    /*
    // FIXME: use values.data() in current deal.II
    dealii::internal::VectorOperations::functions<NumberType, NumberType, MemorySpace::Host>::
      subtract_vector(thread_loop_partitioner,
                      data_start[sp.rowstart[n_owned_row_blocks]],
                      &B.values[0],
                      &values[0]);
    */

    if (matrix_is_ghosted)
      update_ghost_values();

    return *this;
  }



  // explicit instantiations
  template class BlockCSRMatrix<double>;
  template class BlockCSRMatrixIterators::Iterator<double, true>;
  template class BlockCSRMatrixIterators::Iterator<double, false>;
  template class BlockCSRMatrixIterators::Accessor<double, true>;
  template class BlockCSRMatrixIterators::Accessor<double, false>;
  template class BlockCSRMatrixIterators::RowsAccessor<double, true>;
  template class BlockCSRMatrixIterators::RowsAccessor<double, false>;
  template class BlockCSRMatrixIterators::RowsAccessorBase<double, true>;
  template class BlockCSRMatrixIterators::RowsAccessorBase<double, false>;
  template class BlockCSRMatrixIterators::RowsBlockAccessor<double, true>;
  template class BlockCSRMatrixIterators::RowsBlockAccessor<double, false>;

  template void DoFInfo::initialize(
    const DoFHandler<2> &,
    const std::shared_ptr<const dealii::Utilities::MPI::Partitioner> &,
    const std::shared_ptr<const MatrixFree<2,double>> &,
    const std::shared_ptr<const BlockIndices> &,
    const std::vector<unsigned int> &);
  template void DoFInfo::initialize(
    const DoFHandler<3> &,
    const std::shared_ptr<const dealii::Utilities::MPI::Partitioner> &,
    const std::shared_ptr<const MatrixFree<3,double>> &,
    const std::shared_ptr<const BlockIndices> &,
    const std::vector<unsigned int> &);

} // namespace RealFAStDFT
