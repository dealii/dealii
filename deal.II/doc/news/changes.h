/**
 * @page changes_after_7_0 Changes after Version 7.0

<p>
This is the list of changes made after the release of
deal.II version 7.0.0.
All entries are signed with the names of the author.
</p>


<a name="incompatible"></a>
<h3 style="color:red">Incompatibilities</h3>

<p style="color:red">
Following are a few modifications to the library that unfortunately
are incompatible with previous versions of the library, but which we
deem necessary for the future maintainability of the
library. Unfortunately, some of these changes will require
modifications to application programs. We apologize for the
inconvenience this causes.
</p>

<ol>
  <li>
</ol>


<a name="general"></a>
<h3>General</h3>

<ol>

<li> New: The new tutorial program step-46 shows how to couple different
models defined on subsets of
the domain, in this case Stokes flow around an elastic solid. The
trick here is that variables (here the flow velocity and pressure,
and the solid displacement) do not live on the entire domain, but
only on a part. The point of the program is how to represent this in
source code.
<br>
(Wolfgang Bangerth, 2011/04/30)

<li> Fixed: On Debian, the Trilinos packages use a different layout
of include files and library names. The <code>./configure</code>
script can now deal with this.
<br>
(Walter Landry, 2011/02/22)

<li> Improved: Linking the deal.II libraries on file systems that
are mounted remotely from a file server took painfully long. This
is now fixed by linking everything on the local file system
and only subsequently moving the file into its final location.
<br>
(Wolfgang Bangerth, 2011/01/28)

<li> Changed: Most classes in deal.II have a member function
<code>memory_consumption</code> that used to return an unsigned int.
However, on most 64-bit systems, unsigned int is still only 32-bit
wide, and consequently the return type does not provide enough
precision to return the size of very large objects. The return types
of all of these functions has been changed to std::size_t, which is
defined to be a type that can hold the sizes of all objects possible
on any system.
<br>
(Wolfgang Bangerth, 2011/01/22)

<li> Fixed: When using the <code>--enable-mpi</code> to
<code>./configure</code>, the script only tried <code>mpiCC</code>
as the MPI C++ compiler. However, on some systems, it is called
<code>mpicxx</code>. The script now tries that as well.
<br>
(Wolfgang Bangerth, 2011/01/22)

<li> Fixed: When using Trilinos and using the Intel C++ compiler,
we accidentally used invalid compiler flags that led to a warning
every time we compiled a file.
<br>
(Wolfgang Bangerth, 2011/01/22)

<li> Fixed: At the bottom of the page of tutorial programs we show a "plain"
version of the tutorial program. However, the script that generates this plain
version was broken and sometimes truncated the file. This
should be fixed now.
<br>
(Wolfgang Bangerth, 2011/01/18)

<li> Extended: Several missing instantiations of functions for triangulations and DoF handlers embedded in higher dimensional space have been added.
<br>
(Wolfgang Bangerth, 2011/01/15)
</ol>



<a name="specific"></a>
<h3>Specific improvements</h3>

<ol>
<li> Changed: DoFTools is now a namespace. It has long been a class that
had only public, static member functions, making the end result semantically
exactly equivalent to a namespace, which is also how it was used. This is
now also reflected in the actual code.
<br>
(Wolfgang Bangerth, 2011/04/27)

<li> New: The version of DoFTools::make_flux_sparsity_pattern that takes
the coupling masks is now also available for hp::DoFHandler objects.
<br>
(Wolfgang Bangerth, 2011/04/27)

<li> Fixed: If Triangulation::create_triangulation is called after an
hp::DoFHandler object is attached to the triangulation object, setting active
FE indices leads to a crash. The problem did not happen if the mesh was
refined before setting the FE indices. This is now fixed. In the process, the
Triangulation::RefinementListener::create_notification function was
introduced.
<br>
(Wolfgang Bangerth, 2011/04/22)

<li> Fixed: The function FEValuesViews::SymmetricTensor::divergence had a bug.
This is now fixed.
<br>
(Wolfgang Bangerth, Feifei Cheng, Venkat Vallala 2011/04/21)

<li> Fixed: Under some conditions, FEFaceValues applied to an FESystem element
that contained elements of type FE_Nothing would receive an erroneous
exception. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/04/17)

<li> New: There is now an operator* for the multiplication of a <code>SymmetricTensor@<2,dim@></code>
and a <code>Tensor@<1,dim@></code>.
<br>
(Wolfgang Bangerth, 2011/04/12)

<li> Fixed: Added some instantiations to make KellyErrorEstimator and SolutionTransfer
work in  codimension one. Fixed some dim in spacedim.
<br>
(Luca Heltai, 2011/04/11)

<li> Fixed: Added some instantiations to make anisotropic refinement work
in codimension one.
<br>
(Luca Heltai, 2011/03/31)

<li> Fixed: Corrections in the creation of the face and subface
interpolation matrices in the class FE_Nedelec.
<br>
(Markus Buerg, 2011/03/17)

<li> Fixed: In step-21, the inner iteration would sometimes not converge for
very coarse meshes because of numerical roundoff. This is now fixed by allowing
more than <code>rhs.size()</code> CG iterations if the number of degrees of freedom
is very small.
<br>
(Jichao Yin, WB, 2011/04/06)
</li>

<li> New: There is now a new function ConditionalOStream::get_stream().
<br>
(WB, 2011/03/09)
</li>

<li> Fixed: FESystem::get_unit_face_support_points would refuse to return
anything if one of the base elements did not have support points. This
condition has been relaxed: it now only doesn't return anything if this
base element has no support points and also has degrees of freedom on
the face.
<br>
(WB, 2011/03/07)
</li>

<li> Fixed: Objects of type FE_Nothing could be generated with multiple vector components
by passing an argument to the constructor. Yet, the FE_Nothing::get_name() function
always just returned the string <code>FE_Nothing@<dim@>()</code> independently of the
number of components. This is now fixed.
<br>
(WB, 2011/03/07)
</li>

<li> Fixed: PETScWrappers:MPI:SparseMatrix and apply_boundary_values() produced an error in debug mode about non-existant SparsityPattern entries. Reason: clear_rows() also eliminated the whole row in the PETSc-internal SparsityPattern, which resulted in an error in the next assembly process.
<br>
(Timo Heister, 2011/02/23)
</li>

<li> Fixed: It wasn't possible to use the FE_Nothing element inside an FESystem
object and hand the result over to an FEValues object. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/02/18)
</li>

<li> New: There is now a function DataOutBase::write_visit_record that does
the equivalent for VisIt that DataOutBase::write_pvtu_record does for ParaView:
generate a file that contains a list of all other VTK or VTU files of which the
current parallel simulation consists.
<br>
(Wolfgang Bangerth, 2011/02/16)
</li>

<li> New: There is now a function TrilinosWrappers::VectorBase::minimal_value.
<br>
(Wolfgang Bangerth, 2011/02/16)
</li>

<li> Fixed: TableBase::operator= could not be compiled if the type of the
elements of the table was <code>bool</code>. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/02/16)
</li>

<li> Improved: The GridGenerator::hyper_shell function generated meshes in 3d
that are valid but of poor quality upon refinement. There is now an additional
option to generate a coarse mesh of 96 cells that has a much better quality.
<br>
(Wolfgang Bangerth, 2011/02/12)
</li>

<li> Fixed: There are systems where the <code>libz</code> library is installed
but the <code>zlib.h</code> header file is not available. Since the latter
condition was not tested, this would result in compiler errors. This is now
fixed.
<br>
(Wolfgang Bangerth, 2011/02/09)
</li>

<li> Fixed: Prolongation and restriction matrices were not computed at all
for elements of type FE_DGQ if <code>dim@<spacedim</code>. Consequently,
consumers of this information, such as the SolutionTransfer class or
the DoFCellAccess::set_dof_values_by_interpolation function did not
work either and simply returned zero results. This is now fixed.
<br>
(M. Sebastian Pauletti, Wolfgang Bangerth, 2011/02/09)
</li>

<li> Fixed: When refining a mesh with codimension one, edges were refined using
the same manifold description as adjacent cells, but this ignored that a
boundary indicator might have been purposefully set for edges that are truly at
the boundary of the mesh. For such edges, the boundary indicator is now honored.
<br>
(M. Sebastian Pauletti, Wolfgang Bangerth, 2011/02/09)
</li>

<li> Fixed: The functions VectorTools::compute_mean_value and
VectorTools::integrate_difference now also work for distributed
triangulations of type parallel::distributed::Triangulation.
<br>
(Wolfgang Bangerth, 2011/02/07)
</li>

<li> Changed: If the <code>libz</code> library was detected during library
configuration, the function DataOutBase::write_vtu now writes data in compressed
format, saving a good fraction of disk space (80-90% for big output files).
<br>
(Wolfgang Bangerth, 2011/01/28)
</li>

<li> New: Trilinos and PETSc vectors now have a function has_ghost_elements().
<br>
(Timo Heister, 2011/01/26)
</li>

<li> Changed: The TrilinosWrappers::MPI::BlockVector::compress function now takes an
argument (with a default value) in exactly the same way as the
TrilinosWrappers::MPI::Vector::compress function already did.
<br>
(Wolfgang Bangerth, 2011/01/21)
</li>

<li> Fixed: When calling DataOut::build_patches with a mapping, requesting more
than one subdivision, and when <code>dim@<spacedim</code>, then some cells
were not properly mapped. This is now fixed.
<br>
(Wolfgang Bangerth, 2011/01/18)
</li>

<li> New: Restructured the internals of PETScWrappers::Precondition* to allow a
PETSc PC object to exist without a solver. New: use Precondition*::vmult() to
apply the preconditioner once. Preconditioners now have a default constructor
and an initialize() function and are no longer initialized in the solver call,
but in the constructor or initialize().
<br>
(Timo Heister, 2011/01/17)
</li>

<li> Fixed: Boundary conditions in the step-23 tutorial program are now
applied correctly. Matrix columns get eliminated with the used method
and introduce some contribution to the right hand side coming from
inhomogeneous boundary values. The old implementation did not reset the
matrix columns before applying new boundary values.<br>
(Martin Stoll, Martin Kronbichler, 2011/01/14)
</li>

<li> Extended: <code>base/tensor.h</code> has an additional collection of
contractions between three tensors (<i>ie</i>. <code>contract3</code>).
This can be useful for writing matrix/vector assembly in a more compact
form than before.<br>
(Toby D. Young, 2011/01/12)
</ol>


*/
