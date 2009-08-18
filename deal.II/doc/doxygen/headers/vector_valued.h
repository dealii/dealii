//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------


/**
 * @defgroup vector_valued Handling vector valued problems
 *
 *
 * Vector-valued problems are partial differential equations in which the
 * solution variable is not a scalar function, but a vector-valued function or
 * a set of functions. This includes, for example:
 * <ul>
 *   <li>The elasticity equation discussed in @ref step_8 "step-8",
 *       @ref step_17 "step-17", and @ref step_18 "step-18" in which the
 *       solution is the vector-valued displacement at each point.
 *   <li>The mixed Laplace equation and its extensions discussed in
 *       @ref step_20 "step-20", and @ref step_21 "step-21" in which the
 *       solution is the scalar pressure and the vector-valued velocity
 *       at each point.
 *   <li>The Stokes equation and its extensions discussed in
 *       @ref step_22 "step-22", and @ref step_31 "step-31" in which again
 *       the solution is the scalar pressure and the vector-valued velocity
 *       at each point.
 *   <li>Complex-valued solutions consisting of real and imaginary parts, as
 *       discussed for example in @ref step_29 "step-29".
 * </ul>
 * 
 * This page gives an overview of how to implement such vector-valued problems
 * efficiently in deal.II.
 *
 * <table class="tutorial" width="50%">
 * <tr><th><b>%Table of contents</b></th></tr>
 * <tr><td width="100%" valign="top">
 * <ol>
 *  <li> @ref VVPhilosophy "A philosophical view of vector-valued problems"
 *  <li> @ref VVFEs "Describing finite element spaces"
 *  <li> @ref VVAssembling "Assembling linear systems"
 *  <li> @ref VVAlternative "An alternative approach"
 *  <li> @ref VVBlockSolvers "Block solvers"
 *  <li> @ref VVExtracting "Extracting data from solutions"
 * </ol> </td> </tr> </table>
 *
 *
 *
 * @anchor VVPhilosophy
 * <h3>A philosophical view of vector-valued problems</h3>
 *
 * The way one deals systematically with vector-valued problems is not
 * fundamentally different from scalar problems in that one needs to come up
 * first with a weak (variational) formulation of the problem that takes into
 * account all the solution variables. To understand how to do that, let us
 * consider a simple example, the mixed Laplace equations discussed in
 * @ref step_20 "step-20":
@f{eqnarray*}
  \textbf{u} + \nabla p &=& 0,
  \\
  -\textrm{div}\; \textbf{u} &=& f,
@f}
 *
 * Here, we have <code>dim+1</code> solution components: the scalar pressure
 * $p$ and the vector-valued velocity $\textbf u$ with <code>dim</code> vector
 * components.
 *
 * A systematic way to get a weak or variational form for this and other
 * vector problems is to first consider it as a problem where the operators
 * and solution variables are written in vector and matrix form. For the
 * example, this would read as follows:
@f{eqnarray*}
  \left(
  \begin{array}{cc} \mathbf 1 & \nabla \\ -\nabla^T & 0 \end{array}
  \right)
  \left(
  \begin{array}{c} \mathbf u \\ p \end{array}
  \right)
  =
  \left(
  \begin{array}{c} \mathbf 0 \\ f \end{array}
  \right)
@f}
 *
 * This makes it clear that the solution
@f{eqnarray*}
  U =
  \left(
  \begin{array}{c} \mathbf u \\ p \end{array}
  \right)
@f}
 * indeed has <code>dim+1</code> components. We note that we could change the
 * ordering of the solution components $\textbf u$ and $p$ inside $U$ if we
 * also change columns of the matrix operator.
 *
 * Next, we need to think about test functions $V$. We want to multiply both
 * sides of the equation with them, then integrate over $\Omega$. The result
 * should be a scalar equality. We can achieve this by choosing $V$ also
 * vector valued as
@f{eqnarray*}
  V =
  \left(
  \begin{array}{c} \mathbf v \\ q \end{array}
  \right).
@f}
 *
 * It is convenient to multiply the matrix-vector equation by the test
 * function from the left, since this way we automatically get the correct
 * matrix later on (in the linear system, the matrix is also multiplied from
 * the right with the solution variable, not from the left), whereas if we
 * multiplied from the right then the matrix so assembled is the transpose of
 * the one we really want.
 *
 * With this in mind, let us multiply by $V$ and integrate to get the
 * following equation which has to hold for all test functions $V$:
@f{eqnarray*}
  \int_\Omega
  \left(
  \begin{array}{c} \mathbf v \\ q \end{array}
  \right)^T
  \left(
  \begin{array}{cc} \mathbf 1 & \nabla \\ -\nabla^T & 0 \end{array}
  \right)
  \left(
  \begin{array}{c} \mathbf u \\ p \end{array}
  \right)
  \ dx
  =
  \int_\Omega
  \left(
  \begin{array}{c} \mathbf v \\ q \end{array}
  \right)^T
  \left(
  \begin{array}{c} \mathbf 0 \\ f \end{array}
  \right)
  \ dx,
@f}
 * or equivalently:
@f{eqnarray*}
  (\mathbf v, \mathbf u)
  +
  (\mathbf v, \nabla p)
  -
  (q, \mathrm{div}\ \mathbf u)
  =
  (q,f),
@f}
 * where we have defined the scalar product
 * $(\mathbf v, \mathbf u) = \int_\Omega \mathbf v(x)
 * \cdot \mathbf u(x) \; dx$, and similarly if both parts of the scalar
 * product are scalar-valued functions (e.g. the pressure) rather
 * than vector-valued onces (like the velocity).
 *
 * We get the final form by integrating by part the second term:
@f{eqnarray*}
  (\mathbf v, \mathbf u)
  -
  (\mathrm{div}\ \mathbf v, p)
  -
  (q, \mathrm{div}\ \mathbf u)
  =
  (q,f) + (\mathbf n\cdot\mathbf v, p)_{\partial\Omega}.
@f}
 *
 * It is this form that we will later use in assembling the discrete weak form
 * into a matrix and a right hand side vector: the form in which we have
 * solution and test functions $U,V$ that each consist of a number of vector
 * components that we can extract.
 *
 *
 * @anchor VVFEs
 * <h3>Describing finite element spaces</h3>
 *
 * Once we have settled on this description, we need to find a way to describe
 * the vector-valued finite element space from which we draw solution and test
 * functions. This is where the FESystem class comes in: it composes
 * vector-valued finite element spaces from simpler one. For example, if we
 * were to attempt to use $Q_1$ elements for all <code>dim</code> components
 * of $\mathbf u$ and the one pressure component $p$, we could use the
 * following object:
 * @code
 *   FESystem<dim> finite_element (FE_Q<dim>(1), dim,
 *                                 FE_Q<dim>(1), 1);
 * @endcode
 *
 * This means that the final finite element will consist of <code>dim</code>
 * components made up of FE_Q elements of degree 1, and another one also of
 * degree 1. Of course, a simpler (and more efficient) way to achieve the same
 * is to use the following form instead:
 * @code
 *   FESystem<dim> finite_element (FE_Q<dim>(1), dim+1);
 * @endcode
 * Another way (also not efficient, but making it clear which components
 * belong to which element) would be
 * @code
 *   FESystem<dim> finite_element (FESystem<dim>(FE_Q<dim>(1), dim), 1,
 *                                 FE_Q<dim>(1),                     1);
 * @endcode
 * meaning that we first create a vector element with <code>dim</code>
 * components, each consisting of FE_Q elements of order 1. We then couple
 * one copy of this already vector-valued element to a single scalar
 * copy of the FE_Q element of order 1 which will describe the pressure
 * component.
 *
 * As it turns out these (equivalent) choices do not lead to a stable scheme
 * for the mixed Laplace equation. In @ref step_20 "step-20", we therefore use
 * a Raviart-Thomas element for the velocities. What exactly this means may be
 * of less concern for now except that the FE_RaviartThomas class describes
 * elements that already have <code>dim</code> components. For the pressure,
 * we use piecewise bi-/trilinear elements that may be discontinuous between
 * cells; this is done using the FE_DGQ class. The combined element will then
 * be described by
 * @code
 *   FESystem<dim> finite_element (FERaviartThomas<dim>(1), 1,
 *                                 FE_DGQ<dim>(1),          1);
 * @endcode
 * i.e. we combine a single copy of the Raviart-Thomas element with a single
 * copy of the element used for the pressure $p$.
 *
 * In this manner, we can combine the whole vector-valued element from its
 * individual components. However, it is worth pointing out that whichever
 * way you construct your finite element object, this object doesn't really
 * know what it represents. For example, for
 * @code
 *   FESystem<dim> finite_element (FE_Q<dim>(1), dim,
 *                                 FE_Q<dim>(1), 1);
 * @endcode
 * the constructed object really only knows that it has <code>dim+1</code>
 * vector components. It has no notion, however, which of these components
 * represent scalar fields (e.g. temperature, pressure, concentration) and/or
 * if any of its components are parts of vector fields (velocities,
 * displacements) or tensors (e.g. stresses). As a consequence, the FEValues
 * objects we use below to evaluate finite element shape functions at
 * quadrature points only knows that it has a finite element with a number
 * of vector components, but doesn't know how to group them. We will show
 * how to give these components a logical connection using the
 * FEValuesExtractors classes.
 *
 *
 * @anchor VVAssembling
 * <h3>Assembling linear systems</h3> 
 *
 * The next step is to assemble the linear system. How to do this for the
 * simple case of a scalar problem has been shown in many tutorial programs,
 * starting with @ref step_3 "step-3". Here we will show how to do it for
 * vector problems.
 *
 * How to do this is possibly best explained by showing an example
 * illustrating how the local contribution of a cell to the weak form of above
 * mixed Laplace equations could be assembled. This is essentially how @ref
 * step_20 "step-20" does it:
 * @code
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  ...
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);
      
      for (unsigned int q=0; q<n_q_points; ++q) 
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              local_matrix(i,j) += (fe_values[velocities].value (i, q) *
		                    fe_values[velocities].value (j, q)
                                    -
				    fe_values[velocities].divergence (i, q) *
				    fe_values[pressure].value (j, q)
                                    -
				    fe_values[pressure].value (i, q) *
				    fe_values[velocities].divergence (j, q)) *
                                    fe_values.JxW(q);

            local_rhs(i) += - fe_values[pressure].value (i, q)
                              rhs_values[q] *
                              fe_values.JxW(q);
          }
 * @endcode 
 *
 * So here's what is happening:
 * <ul>
 *   <li> The first thing we do is to declare "extractors" (see the
 *        FEValuesExtractors namespace). These are objects that don't
 *        do much except store which components of a vector-valued finite
 *        element constitute a single scalar component, or a tensor of
 *        rank 1 (i.e. what we call a "physical vector", always consisting
 *        of <code>dim</code> components). Here, we declare
 *        an object that represents the velocities consisting of
 *        <code>dim</code> components starting at component zero, and the
 *        extractor for the pressure, which is a scalar component at
 *        position <code>dim</code>.
 *
 *   <li> We then do our usual loop over all cells, shape functions, and
 *        quadrature points. In the innermost loop, we compute the local
 *        contribution of a pair of shape functions to the global matrix
 *        and right hand side vector. Recall that the cell contributions
 *        to the bilinear form (i.e. neglecting boundary terms) looked as
 *        follows, based on shape functions
 *        $V_i=\left(\begin{array}{c}\mathbf v_i q_i\end{array}\right),
 *         V_j=\left(\begin{array}{c}\mathbf v_j q_j\end{array}\right)$:
          @f{eqnarray*}
            (\mathbf v_i, \mathbf v_j)
	    -
	    (\mathrm{div}\ \mathbf v_i, q_j)
	    -
	    (q_i, \mathrm{div}\ \mathbf v_j)
          @f}
 *        whereas the implementation looked like this:
 *        @code
              local_matrix(i,j) += (fe_values[velocities].value (i, q) *
		                    fe_values[velocities].value (j, q)
                                    -
				    fe_values[velocities].divergence (i, q) *
				    fe_values[pressure].value (j, q)
                                    -
				    fe_values[pressure].value (i, q) *
				    fe_values[velocities].divergence (j, q)) *
                                    fe_values.JxW(q);
 *        @endcode
 *        The similarities are pretty obvious.
 *
 *   <li> Essentially, what happens in above code is this: when you do
 *        <code>fe_values[pressure]</code>, a so-called "view" is created, i.e.
 *        an object that unlike the full FEValues object represents not all
 *        components of a finite element, but only the one(s) represented by
 *        the extractor object <code>pressure</code> or
 *        <code>velocities</code>.
 *
 *   <li> These views can then be asked for information about these individual
 *        components. For example, when you write
 *        <code>fe_values[pressure].value(i,q)</code> you get the
 *        value of the pressure component of the $i$th shape function $V_i$ at
 *        the $q$th quadrature point. Because the extractor
 *        <code>pressure</code> represents a scalar component, the results of
 *        the operator <code>fe_values[pressure].value(i,q)</code> is a scalar
 *        number. On the other hand, the call
 *        <code>fe_values[velocities].value(i,q)</code> would produce the
 *        value of a whole set of <code>dim</code> components, which would
 *        be of type <code>Tensor@<1,dim@></code>.
 *
 *   <li> Other things that can be done with views is to ask for the gradient
 *        of a particular shape function's components described by an
 *        extractor. For example, <code>fe_values[pressure].gradient(i,q)</code>
 *        would represent the gradient of the scalar pressure component, which
 *        is of type <code>Tensor@<1,dim@></code>, whereas the gradient of the
 *        velocities components, <code>fe_values[velocities].value(i,q)</code>
 *        is a <code>Tensor@<2,dim@></code>, i.e. a matrix $G_{ij}$ that consits
 *        of entries $G_{ij}=\frac{\partial\phi_i}{\partial x_j}$. Finally,
 *        both scalar and vector views can be asked for the second derivatives
 *        ("Hessians") and vector views can be asked for the symmetric gradient,
 *        defined as $S_{ij}=\frac 12 \left[\frac{\partial\phi_i}{\partial x_j}
 *        + \frac{\partial\phi_j}{\partial x_i}\right]$ as well as the
 *        divergence $\sum_{d=0}^{dim-1} \frac{\partial\phi_d}{\partial x_d}$.
 * </ul>
 * Other examples of using extractors and views are shown in tutorial programs
 * @ref step_21 "step-21",
 * @ref step_22 "step-22",
 * @ref step_31 "step-31" and a few other programs.
 *
 *
 * @anchor VVAlternative
 * <h3>An alternative approach</h3>
 *
 * There are situations in which one can optimize the assembly of a matrix or
 * right hand side vector a bit, using knowledge of the finite element in
 * use. Consider, for example, the bilinear form of the elasticity equations
 * which we are concerned with first in @ref step_8 "step-8":
 *
@f[
  a({\mathbf u}, {\mathbf v}) =
  \left(
    \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
  \right)_\Omega
  +
  \sum_{i,j}
  \left(
    \mu \partial_i u_j, \partial_i v_j
  \right)_\Omega,
  +
  \sum_{i,j}
  \left(
    \mu \partial_i u_j, \partial_j v_i
  \right)_\Omega,
@f]
 *
 * Here, $\mathbf u$ is a vector function with <code>dim</code> components,
 * $\mathbf v$ the corresponding test function, and $\lambda,\mu$ are material
 * parameters. Given our discussions above, the obvious way to implement this
 * bilinear form would be as follows, using an extractor object that interprets
 * all <code>dim</code> components of the finite element as single vector,
 * rather than disjoint components:
 *
 * @code
      const FEValuesExtractors::Vector displacements (0);

      ...
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    const Tensor<2,dim> phi_i_grad
	      = fe_values[displacements].gradient (i,q_point);
	    const double phi_i_div
	      = fe_values[displacements].divergence (i,q_point);
	  
	    for (unsigned int j=0; j<dofs_per_cell; ++j) 
	      {
		const Tensor<2,dim> phi_j_grad
		  = fe_values[displacements].gradient (j,q_point);
		const double phi_j_div
		  = fe_values[displacements].divergence (j,q_point);

		cell_matrix(i,j) 
		  +=  (phi_i_div * phi_j_div *
		       lambda_values[q_point]
		       +
		       scalar_product(phi_i_grad, phi_j_grad) *
		       mu_values[q_point]
		       +
		       scalar_product(phi_i_grad, transpose(phi_j_grad)) *
		       mu_values[q_point])
		      *
		      fe_values.JxW(q_point);
	      }
	  }
 * @endcode 
 *
 * The scalar product between two tensors used in this bilinear form is
 * implemented as follows:
 *
 * @code
template <int dim>
double
scalar_product (const Tensor<2,dim> &u,
		const Tensor<2,dim> &v)
{
  double tmp = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      tmp += u[i][j] * v[i][j];
  return tmp;
}
 * @endcode 
 *
 * Now, this is not the code used in @ref step_8 "step-8". In fact,
 * if one used the above code over the one implemented in that program,
 * it would run about 8 per cent slower. It can be improved (bringing
 * down the penalty to about 4 per cent) by taking a close look at the
 * bilinear form. In fact, we can transform it as follows:
@f{eqnarray*} 
  a({\mathbf u}, {\mathbf v})
  &=&
  \left(
    \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
  \right)_\Omega
  +
  \sum_{i,j}
  \left(
    \mu \partial_i u_j, \partial_i v_j
  \right)_\Omega
  +
  \sum_{i,j}
  \left(
    \mu \partial_i u_j, \partial_j v_i
  \right)_\Omega
  \\
  &=&
  \left(
    \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
  \right)_\Omega
  +
  2
  \sum_{i,j}
  \left(
    \mu \partial_i u_j, \frac 12[\partial_i v_j + \partial_j v_i]
  \right)_\Omega  
  \\
  &=&
  \left(
    \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
  \right)_\Omega
  +
  2
  \sum_{i,j}
  \left(
    \mu \frac 12[\partial_i u_j + \partial_j u_i], \frac 12[\partial_i v_j + \partial_j v_i]
  \right)_\Omega
  \\
  &=&
  \left(
    \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
  \right)_\Omega
  +
  2
  \sum_{i,j}
  \left(
  \mu \varepsilon(\mathbf u), \varepsilon(\mathbf v)
  \right)_\Omega,
@f}
 * where $\varepsilon(\mathbf u) = \frac 12 \left([\nabla\mathbf u] +
 * [\nabla\mathbf u]^2\right)$ is the symmetrized gradient.
 * In the second to last step, we used that the scalar product between
 * an arbitrary tensor $\nabla\mathbf u$ and a symmetric tensor
 * $\frac 12[\partial_i v_j + \partial_j v_i]$ equals the scalar product
 * of the symmetric part of the former with the second tensor. Using the
 * techniques discussed above, the obvious way to implement this goes
 * like this:
 *
 * @code
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    const SymmetricTensor<2,dim> phi_i_symmgrad
	      = fe_values[displacements].symmetric_gradient (i,q_point);
	    const double phi_i_div
	      = fe_values[displacements].divergence (i,q_point);
	  
	    for (unsigned int j=0; j<dofs_per_cell; ++j) 
	      {
		const SymmetricTensor<2,dim> phi_j_symmgrad
		  = fe_values[displacements].symmetric_gradient (j,q_point);
		const double phi_j_div
		  = fe_values[displacements].divergence (j,q_point);

		cell_matrix(i,j) 
		  +=  (phi_i_div * phi_j_div *
		       lambda_values[q_point]
		       +
		       2 *
		       (phi_i_symmgrad * phi_j_symmgrad) *
		       mu_values[q_point])
		      *
		      fe_values.JxW(q_point));
	      }
	  }
 * @endcode 
 *
 * So if, again, this is not the code we use in @ref step_8 "step-8", what do
 * we do there? The answer rests on the finite element we use. There, we use the
 * following element:
 * @code
 *   FESystem<dim> finite_element (FE_Q<dim>(1), dim);
 * @endcode
 * In other words, the finite element we use consists of <code>dim</code> copies
 * of the same scalar element. This is what we call a @ref GlossPrimitive
 * "primitive" element: an element that may be vector-valued but where each
 * shape function has exactly one non-zero component. In other words: if the
 * $x$-component of a displacement shape function is nonzero, then the $y$-
 * and $z$-components must be zero and so on. What this means is that also
 * derived quantities based on shape functions inherit this sparsity property.
 * For example: the divergence
 * $\mathrm{div}\ \Phi(x,y,z)=\partial_x\varphi_x(x,y,z) +
 * \partial_y\varphi_y(x,y,z) + \partial_z\varphi_z(x,y,z)$
 * of a vector-valued shape function
 * $\Phi(x,y,z)=(\varphi_x(x,y,z), \varphi_y(x,y,z), \varphi_z(x,y,z))^T$ is,
 * in the present case, either
 * $\mathrm{div}\ \Phi(x,y,z)=\partial_x\varphi_x(x,y,z)$,
 * $\mathrm{div}\ \Phi(x,y,z)=\partial_y\varphi_y(x,y,z)$, or
 * $\mathrm{div}\ \Phi(x,y,z)=\partial_z\varphi_z(x,y,z)$, because exactly one
 * of the $\varphi_\ast$ is nonzero. Knowing this means that we can save a
 * number of computations that, if we were to do them, would only yield
 * zeros to add up.
 *
 * In a similar vein, if only one component of a shape function is nonzero,
 * then only one row of its gradient $\nabla\Phi$ is nonzero. What this means
 * for terms like $(\mu \nabla\Phi_i,\nabla\Phi_j)$, where the scalar product
 * between two tensors is defined as
 * $(\tau, \gamma)_\Omega=\int_\Omega \sum_{i,j=1}^d \tau_{ij} \gamma_{ij}$,
 * is that the term is only nonzero if both tensors have their nonzero
 * entries in the same row, which means that the two shape functions have
 * to have their single nonzero component in the same location.
 *
 * If we use this sort of knowledge, then we can in a first step avoid
 * computing gradient tensors if we can determine up front that their
 * scalar product will be nonzero, in a second step avoid
 * building the entire tensors and only get its nonzero components,
 * and in a final step simplify the scalar product by only considering
 * that index $i$ for the one nonzero row, rather than multiplying and
 * adding up zeros.
 *
 * The vehicle for all this is the ability to determine which vector
 * component is going to be nonzero. This information is provided by the
 * FiniteElement::system_to_component_index function. What can be done with
 * it, using the example above, is explained in detail in
 * @ref step_8 "step-8".
 *
 *
 * @anchor VVBlockSolvers
 * <h3>Block solvers</h3>
 *
 * Using techniques as shown above, it isn't particularly complicated to
 * assemble the linear system, i.e. matrix and right hand side, for a
 * vector-valued problem. However, then it also has to be solved. This is more
 * complicated. Naively, one could just consider the matrix as a whole. For
 * most problems, this matrix is not going to be definite (except for special
 * cases like the elasticity equations covered in @ref step_8 "step-8" and
 * @ref step_17 "step-17"). It will, often, also not be symmetric. This rather
 * general class of matrices presents problems for iterative solvers: the lack
 * of structural properties prevents the use of most efficient methods and
 * preconditioners. While it can be done, the solution process will therefore
 * most often be slower than necessary.
 *
 * The answer to this problem is to make use of the structure of the
 * problem. For example, for the mixed Laplace equations discussed above, the
 * operator has the form
@f{eqnarray*}
  \left(
  \begin{array}{cc} \mathbf 1 & \nabla \\ -\nabla^T & 0 \end{array}
  \right)
@f}
 *
 * It would be nice if this structure could be recovered in the linear system
 * as well. For example, after discretization, we would like to have a matrix
 * with the following block structure:
@f{eqnarray*}
  \left(
  \begin{array}{cc} M & B \\ B^T & 0 \end{array}
  \right),
@f}
 * where $M$ represents the mass matrix that results from discretizing the
 * identity operator $\mathbf 1$ and $B$ the equivalent of the gradient
 * operator.
 *
 * By default, this is not what happens, however. Rather, deal.II assigns
 * %numbers to degrees of freedom in a rather random manner. Consequently, if
 * you form a vector out of the values of degrees of freedom will not be
 * neatly ordered in a vector like
@f{eqnarray*}
  \left(
  \begin{array}{c} U \\ P \end{array}
  \right).
@f}
 * Rather, it will be a permutation of this, with %numbers of degrees of
 * freedom corresponding to velocities and pressures intermixed. Consequently,
 * the system matrix will also not have the nice structure mentioned above,
 * but with the same permutation or rows and columns.
 *
 * What is needed is to re-enumerate degrees of freedom so that velocities
 * come first and pressures last. This can be done using the
 * DoFRenumbering::component_wise function, as explained in @ref step_20
 * "step-20", @ref step_21 "step-21", @ref step_22 "step-22", and @ref step_31
 * "step-31". After this, at least the degrees of freedom are partitioned
 * properly.
 *
 * But then we still have to make use of it, i.e. we have to come up with a
 * solver that uses the structure. For example, in @ref step_20 "step-20", we
 * do a block elimination of the linear system
@f{eqnarray*}
  \left(
  \begin{array}{cc} M & B \\ B^T & 0 \end{array}
  \right)
  \left(
  \begin{array}{c} U \\ P \end{array}
  \right)
  =
  \left(
  \begin{array}{c} F \\ G \end{array}
  \right).  
@f}
 * What this system means, of course, is
@f{eqnarray*}
  MU + BP &=& F,\\
  B^TU  &=& G.
@f}
 *
 * So, if we multiply the first equation by $B^TM^{-1}$ and subtract the
 * second from the result, we get
@f{eqnarray*}
  B^TM^{-1}BP &=& B^TM^{-1}F-G.
@f}
 *
 * This is an equation that now only contains the pressure variables. If we
 * can solve it, we can in a second step solve for the velocities using
@f{eqnarray*}
  MU = F-BP.
@f}
 *
 * This has the advantage that the matrices $B^TM^{-1}B$ and $M$ that we have
 * to solve with are both symmetric and positive definite, as opposed to the
 * large whole matrix we had before.
 *
 * How a solver like this is implemented is explained in more detail in @ref
 * step_20 "step-20", @ref step_31 "step-31", and a few other tutorial
 * programs. What we would like to point out here is that we now need a way to
 * extract certain parts of a matrix or vector: if we are to multiply, say,
 * the $U$ part of the solution vector by the $M$ part of the global matrix,
 * then we need to have a way to access these parts of the whole.
 *
 * This is where the BlockVector, BlockSparseMatrix, and similar classes come
 * in. For all practical purposes, then can be used as regular vectors or
 * sparse matrices, i.e. they offer element access, provide the usual vector
 * operations and implement, for example, matrix-vector multiplications. In
 * other words, assembling matrices and right hand sides works in exactly the
 * same way as for the non-block versions. That said, internally they store
 * the elements of vectors and matrices in "blocks"; for example, instead of
 * using one large array, the BlockVector class stores it as a set of arrays
 * each of which we call a block. The advantage is that, while the whole thing
 * can be used as a vector, one can also access an individual block which
 * then, again, is a vector with all the vector operations.
 *
 * To show how to do this, let us consider the second equation $MU=F-BP$ to be
 * solved above. This can be achieved using the following sequence similar to
 * what we have in @ref step_20 "step-20":
 * @code
    Vector<double> tmp (solution.block(0).size());
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);
    

    SolverControl solver_control (solution.block(0).size(),
                                  1e-8*tmp.l2_norm());
    SolverCG<> cg (solver_control, vector_memory);
  
    cg.solve (system_matrix.block(0,0),
              solution.block(0),
	      tmp,
              PreconditionIdentity());        
 * @endcode
 *
 * What's happening here is that we allocate a temporary vector with as many
 * elements as the first block of the solution vector, i.e. the velocity
 * component $U$, has. We then set this temporary vector equal to the $(0,1)$
 * block of the matrix, i.e. $B$, times component 1 of the solution which is
 * the previously computed pressure $P$. The result is multiplied by $-1$, and
 * component 0 of the right hand side, $F$ is added to it. The temporary
 * vector now contains $F-BP$. The rest of the code snippet simply solves a
 * linear system with $F-BP$ as right hand side and the $(0,0)$ block of the
 * global matrix, i.e. $M$. Using block vectors and matrices in this way
 * therefore allows us to quite easily write rather complicated solvers making
 * use of the block structure of a linear system.
 *
 *
 *
 * @anchor VVExtracting
 * <h3>Extracting data from solutions</h3>
 *
 * Once one has computed a solution, it is often necessary to evaluate it at
 * quadrature points, for example to evaluate nonlinear residuals for the next
 * Newton iteration, to evaluate the finite element residual for error
 * estimators, or to compute the right hand side for the next time step in
 * a time dependent problem.
 *
 * The way this is done us to again use an FEValues object to evaluate
 * the shape functions at quadrature points, and with those also the
 * values of a finite element function. For the example of the mixed
 * Laplace problem above, consider the following code after solving:
 * @code
  std::vector<Vector<double> > local_solution_values (n_q_points,
                                                      Vector<double> (dim+1));
 
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      fe_values.get_function_values (solution,
                                     local_solution_values);
 * @endcode
 *
 * After this, the variable <code>local_solution_values</code> is a list
 * of vectors of a length equal to the number of quadrature points we
 * have initialized the FEValues object with; each of these vectors
 * has <code>dim+1</code> elements containing the values of the
 * <code>dim</code> velocities and the one pressure at a quadrature point.
 *
 * We can use these values to then construct other things like residuals.
 * However, the construct is a bit awkward. First, we have a
 * <code>std::vector</code> of <code>dealii::Vector</code>s, which always
 * looks strange. It is also inefficient because it implies dynamic memory
 * allocation for the outer vector as well as for all the inner vectors.
 * Secondly, maybe we are only interested in the velocities,
 * for example to solve an advection problem in a second stage (as, for
 * example, in @ref step_21 "step-21" or @ref step_31 "step-31"). In that
 * case, one would have to hand-extract these values like so:
 * @code
   for (unsigned int q=0; q<n_q_points; ++q)
     {
       Tensor<1,dim> velocity;
       for (unsigned int d=0; d<dim; ++d)
         velocity[d] = local_solution_values[q](d);
	 
       ... do something with this velocity ...				  
 * @endcode
 * Note how we convert from a dealii::Vector (which is simply a collection
 * of vector elements) into a <code>Tensor@<1,dim@></code> because the
 * velocity is a quantity characterized by <code>dim</code> elements that
 * have certain transformation properties under rotations of the coordinate
 * system.
 *
 * This code can be written more elegantly and efficiently using code like
 * the following:
 * @code
  std::vector<Tensor<1,dim> > local_velocity_values (n_q_points);

  const FEValuesExtractors::Vector velocities (0);
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      fe_values[velocities].get_function_values (solution,
                                                 local_velocity_values);
 * @endcode
 *
 * As a result, we here get the velocities right away, and in the
 * right data type (because we have described, using the extractor,
 * that the first <code>dim</code> components of the finite element
 * belong together, forming a tensor). The code is also more efficient:
 * it requires less dynamic memory allocation because the Tensor
 * class allocates its components as member variables rather than on
 * the heap, and we save cycles because we don't even bother computing
 * the values of the pressure variable at the quadrature points. On
 * the other hand, if we had been interested in only the pressure and
 * not the velocities, then the following code extracting scalar
 * values would have done:
 * @code
  std::vector<double> local_pressure_values (n_q_points);

  const FEValuesExtractors::Scalar pressure (dim);
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      fe_values[pressure].get_function_values (solution,
                                               local_pressure_values);
 * @endcode
 *
 * In similar cases, one sometimes needs the gradients or second
 * derivatives of the solution, or of individual scalar or vector
 * components. To get at those of all components of the solution,
 * the functions FEValuesBase::get_function_gradients and
 * FEValuesBase::get_function_hessians are the equivalent of the
 * function FEValuesBase::get_function_values used above.
 *
 * Likewise, to extract the gradients of scalar components,
 * FEValuesViews::Scalar::get_function_gradients and
 * FEValuesViews::Scalar::get_function_hessians do the job.
 * For vector-(tensor-)valued quantities, there are functions
 * FEValuesViews::Vector::get_function_gradients and
 * FEValuesViews::Vector::get_function_hessians, and in
 * addition
 * FEValuesViews::Vector::get_function_symmetric_gradients and
 * FEValuesViews::Vector::get_function_divergences.
 *
 * Moreover, there is a shortcut available in case when only the
 * Laplacians of the solution (which is the trace of the hessians) is
 * needed, usable for both scalar and vector-valued problems as
 * FEValuesViews::Scalar::get_function_laplacians and
 * FEValuesViews::Vector::get_function_laplacians.
 *
 * @ingroup feall feaccess
 */
 
