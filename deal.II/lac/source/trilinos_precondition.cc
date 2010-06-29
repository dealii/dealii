//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/trilinos_precondition.h>

#ifdef DEAL_II_USE_TRILINOS

#  include <lac/vector.h>
#  include <lac/sparse_matrix.h>
#  include <lac/trilinos_sparse_matrix.h>

#  include <Ifpack.h>
#  include <Ifpack_Chebyshev.h>
#  include <Teuchos_ParameterList.hpp>
#  include <Epetra_MultiVector.h>
#  include <ml_include.h>
#  include <ml_MultiLevelPreconditioner.h>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{

  PreconditionBase::PreconditionBase()
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
                    :
                    communicator (MPI_COMM_SELF)
#endif
  {}



  PreconditionBase::PreconditionBase(const PreconditionBase &base)
                    :
                    Subscriptor (),
                    preconditioner (base.preconditioner),
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
                    communicator (base.communicator),
#endif
		    vector_distributor (new Epetra_Map(*base.vector_distributor))
  {}



  PreconditionBase::~PreconditionBase()
  {}



  void PreconditionBase::clear ()
  {
    preconditioner.reset();
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    communicator = MPI_COMM_SELF;
#endif
    vector_distributor.reset();
  }



/* -------------------------- PreconditionJacobi -------------------------- */

  PreconditionJacobi::AdditionalData::
  AdditionalData (const double omega,
		  const double min_diagonal)
                  :
                  omega (omega),
		  min_diagonal (min_diagonal)
  {}



  void
  PreconditionJacobi::initialize (const SparseMatrix   &matrix,
				  const AdditionalData &additional_data)
  {
				// release memory before reallocation
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
			  ("point relaxation",
			   const_cast<Epetra_CrsMatrix*>(&matrix.trilinos_matrix()),
			   0));

    Ifpack_Preconditioner * ifpack = static_cast<Ifpack_Preconditioner*>
      (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
				     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: sweeps", 1);
    parameter_list.set ("relaxation: type", "Jacobi");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value",
			additional_data.min_diagonal);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



/* -------------------------- PreconditionSSOR -------------------------- */

  PreconditionSSOR::AdditionalData::
  AdditionalData (const double       omega,
		  const double       min_diagonal,
		  const unsigned int overlap)
                  :
                  omega        (omega),
		  min_diagonal (min_diagonal),
		  overlap      (overlap)
  {}



  void
  PreconditionSSOR::initialize (const SparseMatrix   &matrix,
				const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
			  ("point relaxation",
			   const_cast<Epetra_CrsMatrix*>(&matrix.trilinos_matrix()),
			   additional_data.overlap));

    Ifpack_Preconditioner * ifpack = static_cast<Ifpack_Preconditioner*>
      (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
				     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: sweeps", 1);
    parameter_list.set ("relaxation: type", "symmetric Gauss-Seidel");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value",
			additional_data.min_diagonal);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



/* -------------------------- PreconditionSOR -------------------------- */

  PreconditionSOR::AdditionalData::
  AdditionalData (const double       omega,
		  const double       min_diagonal,
		  const unsigned int overlap)
                  :
                  omega        (omega),
		  min_diagonal (min_diagonal),
		  overlap      (overlap)
  {}



  void
  PreconditionSOR::initialize (const SparseMatrix   &matrix,
			       const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
			  ("point relaxation",
			   const_cast<Epetra_CrsMatrix*>(&matrix.trilinos_matrix()),
			   additional_data.overlap));

    Ifpack_Preconditioner * ifpack = static_cast<Ifpack_Preconditioner*>
      (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
				     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: sweeps", 1);
    parameter_list.set ("relaxation: type", "Gauss-Seidel");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value",
			additional_data.min_diagonal);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



/* -------------------------- PreconditionIC -------------------------- */

  PreconditionIC::AdditionalData::
  AdditionalData (const unsigned int ic_fill,
		  const double       ic_atol,
		  const double       ic_rtol,
		  const unsigned int overlap)
                  :
                  ic_fill (ic_fill),
		  ic_atol (ic_atol),
		  ic_rtol (ic_rtol),
		  overlap (overlap)
  {}



  void
  PreconditionIC::initialize (const SparseMatrix   &matrix,
			      const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
			  ("IC",
			   const_cast<Epetra_CrsMatrix*>(&matrix.trilinos_matrix()),
			   additional_data.overlap));

    Ifpack_Preconditioner * ifpack = static_cast<Ifpack_Preconditioner*>
      (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
				     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("fact: level-of-fill",(int)additional_data.ic_fill);
    parameter_list.set ("fact: absolute threshold",additional_data.ic_atol);
    parameter_list.set ("fact: relative threshold",additional_data.ic_rtol);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



/* -------------------------- PreconditionILU -------------------------- */

  PreconditionILU::AdditionalData::
  AdditionalData (const unsigned int ilu_fill,
		  const double       ilu_atol,
		  const double       ilu_rtol,
		  const unsigned int overlap)
                  :
                  ilu_fill (ilu_fill),
		  ilu_atol (ilu_atol),
		  ilu_rtol (ilu_rtol),
		  overlap  (overlap)
  {}



  void
  PreconditionILU::initialize (const SparseMatrix   &matrix,
			       const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
			  ("ILU",
			   const_cast<Epetra_CrsMatrix*>(&matrix.trilinos_matrix()),
			   additional_data.overlap));

    Ifpack_Preconditioner * ifpack = static_cast<Ifpack_Preconditioner*>
      (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
				     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("fact: level-of-fill",(int)additional_data.ilu_fill);
    parameter_list.set ("fact: absolute threshold",additional_data.ilu_atol);
    parameter_list.set ("fact: relative threshold",additional_data.ilu_rtol);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



/* -------------------------- PreconditionILUT -------------------------- */

  PreconditionILUT::AdditionalData::
  AdditionalData (const double       ilut_drop,
		  const unsigned int ilut_fill,
		  const double       ilut_atol,
		  const double       ilut_rtol,
		  const unsigned int overlap)
                  :
                  ilut_drop (ilut_drop),
                  ilut_fill (ilut_fill),
		  ilut_atol (ilut_atol),
		  ilut_rtol (ilut_rtol),
		  overlap  (overlap)
  {}



  void
  PreconditionILUT::initialize (const SparseMatrix   &matrix,
				const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
			  ("ILUT",
			   const_cast<Epetra_CrsMatrix*>(&matrix.trilinos_matrix()),
			   additional_data.overlap));

    Ifpack_Preconditioner * ifpack = static_cast<Ifpack_Preconditioner*>
      (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
				     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("fact: drop value",additional_data.ilut_drop);
    parameter_list.set ("fact: level-of-fill",(int)additional_data.ilut_fill);
    parameter_list.set ("fact: absolute threshold",additional_data.ilut_atol);
    parameter_list.set ("fact: relative threshold",additional_data.ilut_rtol);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



/* ---------------------- PreconditionBlockDirect --------------------- */

  PreconditionBlockwiseDirect::AdditionalData::
  AdditionalData (const unsigned int overlap)
                  :
                  overlap  (overlap)
  {}



  void
  PreconditionBlockwiseDirect::initialize (const SparseMatrix   &matrix,
					   const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
			  ("Amesos",
			   const_cast<Epetra_CrsMatrix*>(&matrix.trilinos_matrix()),
			   additional_data.overlap));

    Ifpack_Preconditioner * ifpack = static_cast<Ifpack_Preconditioner*>
      (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
				     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



/* ---------------------- PreconditionBlockDirect --------------------- */

  PreconditionChebyshev::AdditionalData::
  AdditionalData (const unsigned int degree,
		  const double       max_eigenvalue,
		  const double       eigenvalue_ratio,
		  const double       min_eigenvalue,
		  const double       min_diagonal,
		  const bool         nonzero_starting)
                  :
                  degree  (degree),
		  max_eigenvalue (max_eigenvalue),
		  eigenvalue_ratio (eigenvalue_ratio),
		  min_eigenvalue (min_eigenvalue),
		  min_diagonal (min_diagonal),
		  nonzero_starting (nonzero_starting)
  {}



  void
  PreconditionChebyshev::initialize (const SparseMatrix   &matrix,
				     const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (new Ifpack_Chebyshev (&matrix.trilinos_matrix()));

    Ifpack_Chebyshev * ifpack = static_cast<Ifpack_Chebyshev*>
      (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
				     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("chebyshev: ratio eigenvalue",
			additional_data.eigenvalue_ratio);
    parameter_list.set ("chebyshev: min eigenvalue",
			additional_data.min_eigenvalue);
    parameter_list.set ("chebyshev: max eigenvalue",
			additional_data.max_eigenvalue);
    parameter_list.set ("chebyshev: degree",
			(int)additional_data.degree);
    parameter_list.set ("chebyshev: min diagonal value",
			additional_data.min_diagonal);
    parameter_list.set ("chebyshev: zero starting solution",
			!additional_data.nonzero_starting);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



/* -------------------------- PreconditionAMG -------------------------- */

  PreconditionAMG::AdditionalData::
  AdditionalData (const bool                             elliptic,
		  const bool                             higher_order_elements,
		  const unsigned int                     n_cycles,
		  const bool                             w_cycle,
		  const double                           aggregation_threshold,
		  const std::vector<std::vector<bool> > &constant_modes,
		  const unsigned int                     smoother_sweeps,
		  const unsigned int                     smoother_overlap,
		  const bool                             output_details)
                  :
                  elliptic (elliptic),
		  higher_order_elements (higher_order_elements),
		  n_cycles (n_cycles),
		  w_cycle (w_cycle),
		  aggregation_threshold (aggregation_threshold),
		  constant_modes (constant_modes),
		  smoother_sweeps (smoother_sweeps),
		  smoother_overlap (smoother_overlap),
		  output_details (output_details)
  {}



  void
  PreconditionAMG:: initialize (const SparseMatrix   &matrix,
				const AdditionalData &additional_data)
  {
    const unsigned int n_rows = matrix.m();

				        // Build the AMG preconditioner.
    Teuchos::ParameterList parameter_list;

    if (additional_data.elliptic == true)
      {
	ML_Epetra::SetDefaults("SA",parameter_list);
	parameter_list.set("smoother: type", "Chebyshev");

				   // uncoupled mode can give a lot of
				   // warnings or even fail when there
				   // are too many entries per row and
				   // aggreggation gets complicated, but
				   // MIS does not work if too few
				   // elements are located on one
				   // processor. work around these
				   // warnings by choosing the different
				   // strategies in different
				   // situations: for low order, always
				   // use the standard choice
				   // uncoupled. if higher order, right
				   // now we also just use Uncoupled,
				   // but we should be aware that maybe
				   // MIS might be needed
				   //
				   // TODO: Maybe there are any
				   // other/better options?
	if (additional_data.higher_order_elements)
	  {
	    //if (matrix.m()/matrix.matrix->Comm().NumProc() < 50000)
	    parameter_list.set("aggregation: type", "Uncoupled");
	    //else
	    //parameter_list.set("aggregation: type", "MIS");
	  }
      }
    else
      {
	ML_Epetra::SetDefaults("NSSA",parameter_list);
	parameter_list.set("aggregation: type", "Uncoupled");
	parameter_list.set("aggregation: block scaling", true);
      }

    parameter_list.set("smoother: sweeps",
		       static_cast<int>(additional_data.smoother_sweeps));
    parameter_list.set("cycle applications",
		       static_cast<int>(additional_data.n_cycles));
    if (additional_data.w_cycle == true)
      parameter_list.set("prec type", "MGW");
    else
      parameter_list.set("prec type", "MGV");

    parameter_list.set("smoother: Chebyshev alpha",10.);
    parameter_list.set("smoother: ifpack overlap",
		       static_cast<int>(additional_data.smoother_overlap));
    parameter_list.set("aggregation: threshold",
		       additional_data.aggregation_threshold);

    if (additional_data.output_details)
      parameter_list.set("ML output", 10);
    else
      parameter_list.set("ML output", 0);

    const Epetra_Map & domain_map = matrix.domain_partitioner();

    const unsigned int constant_modes_dimension =
      additional_data.constant_modes.size();
    Epetra_MultiVector distributed_constant_modes (domain_map,
						   constant_modes_dimension);

    if (constant_modes_dimension > 1)
      {
	const bool constant_modes_are_global =
	  additional_data.constant_modes[0].size() == n_rows;
	const unsigned int n_relevant_rows =
	  constant_modes_are_global ? n_rows : additional_data.constant_modes[0].size();
	const unsigned int my_size = domain_map.NumMyElements();
	if (constant_modes_are_global == false)
	  Assert (n_relevant_rows == my_size,
		  ExcDimensionMismatch(n_relevant_rows, my_size));
	Assert (n_rows ==
		static_cast<unsigned int>(distributed_constant_modes.GlobalLength()),
		ExcDimensionMismatch(n_rows,
				     distributed_constant_modes.GlobalLength()));

				        // Reshape null space as a
				        // contiguous vector of
				        // doubles so that Trilinos
				        // can read from it.
	for (unsigned int d=0; d<constant_modes_dimension; ++d)
	  for (unsigned int row=0; row<my_size; ++row)
	    {
	      int global_row_id = constant_modes_are_global ? domain_map.GID(row) : row;
	      distributed_constant_modes[d][row] = 
		additional_data.constant_modes[d][global_row_id];
	    }

	parameter_list.set("null space: type", "pre-computed");
	parameter_list.set("null space: dimension",
			   distributed_constant_modes.NumVectors());
	parameter_list.set("null space: vectors",
			   distributed_constant_modes.Values());
      }

    initialize (matrix, parameter_list);

    if (additional_data.output_details)
      {
	ML_Epetra::MultiLevelPreconditioner * multilevel_operator = 
	  dynamic_cast<ML_Epetra::MultiLevelPreconditioner*> (preconditioner.get());
	Assert (multilevel_operator != 0, 
		ExcMessage ("Preconditioner setup failed."));
	multilevel_operator->PrintUnused(0);
      }
  }



  void
  PreconditionAMG::initialize (const SparseMatrix           &matrix,
			       const Teuchos::ParameterList &ml_parameters)
  {
    preconditioner.reset ();
    preconditioner.reset (new ML_Epetra::MultiLevelPreconditioner
			  (matrix.trilinos_matrix(), ml_parameters));
  }



  template <typename number>
  void
  PreconditionAMG::
  initialize (const ::dealii::SparseMatrix<number> &deal_ii_sparse_matrix,
	      const AdditionalData                 &additional_data,
	      const double                          drop_tolerance,
	      const ::dealii::SparsityPattern      *use_this_sparsity)
  {
    preconditioner.reset();
    const unsigned int n_rows = deal_ii_sparse_matrix.m();

				        // Init Epetra Matrix using an
				        // equidistributed map; avoid
				        // storing the nonzero
				        // elements.
    vector_distributor.reset (new Epetra_Map(n_rows, 0, communicator));

    if (trilinos_matrix.get() == 0)
      trilinos_matrix.reset (new SparseMatrix());

    trilinos_matrix->reinit (*vector_distributor, *vector_distributor,
			     deal_ii_sparse_matrix, drop_tolerance, true,
			     use_this_sparsity);

    initialize (*trilinos_matrix, additional_data);
  }



  void PreconditionAMG::reinit ()
  {
    ML_Epetra::MultiLevelPreconditioner * multilevel_operator = 
      dynamic_cast<ML_Epetra::MultiLevelPreconditioner*> (preconditioner.get());
    multilevel_operator->ReComputePreconditioner();
  }



  void PreconditionAMG::clear ()
  {
    PreconditionBase::clear();
    trilinos_matrix.reset();
  }



  unsigned int PreconditionAMG::memory_consumption() const
  {
    unsigned int memory = sizeof(this);

				// todo: find a way to read out ML's data
				// sizes
    if (trilinos_matrix.get() != 0)
      memory += trilinos_matrix->memory_consumption();
    return memory;
  }




				// explicit instantiations
  template void PreconditionAMG::initialize (const ::dealii::SparseMatrix<double> &,
					     const AdditionalData &, const double,
					     const ::dealii::SparsityPattern*);
  template void PreconditionAMG::initialize (const ::dealii::SparseMatrix<float> &,
					     const AdditionalData &, const double,
					     const ::dealii::SparsityPattern*);

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
