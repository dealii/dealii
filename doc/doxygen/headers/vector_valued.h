// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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



/**
 * @defgroup vector_valued Handling vector valued problems
 *
 *
 * Vector-valued problems are systems of partial differential
 * equations. These are problems where the
 * solution variable is not a scalar function, but a vector-valued function or
 * a set of functions. This includes, for example:
 * <ul>
 *   <li>The elasticity equation discussed in step-8,
 *       step-17, and step-18 in which the
 *       solution is the vector-valued displacement at each point.
 *   <li>The mixed Laplace equation and its extensions discussed in
 *       step-20, and step-21 in which the
 *       solution is the scalar pressure and the vector-valued velocity
 *       at each point.
 *   <li>The Stokes equation and its extensions discussed in
 *       step-22, and step-31 in which again
 *       the solution is the scalar pressure and the vector-valued velocity
 *       at each point.
 *   <li>Complex-valued solutions consisting of real and imaginary parts, as
 *       discussed for example in step-29.
 * </ul>
 *
 * This page gives an overview of how to implement such vector-valued problems
 * easily in deal.II. In particular, it explains the usage of the class
 * FESystem, which allows us to write code for systems of partial
 * differential very much like we write code for single equations.
 *
 * @dealiiVideoLecture{19,20,21}
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
 *  <li> @ref VVOutput "Generating graphical output"
 * </ol> </td> </tr> </table>
 *
 *
 *
 * @anchor VVPhilosophy
 * <h3>Examples of vector-valued problems</h3>
 *
 * The way one deals systematically with vector-valued problems is not
 * fundamentally different from scalar problems: first, we need a weak
 * (variational) formulation of the problem that takes into account
 * all the solution variables. After we did so, generating the system
 * matrix and solving the linear system follows the same outlines that
 * we are used to already.
 *
 * <h4>Linear elasticity</h4>
 *
 * Let us take for example the elasticity problem from step-8 and even
 * simplify it by choosing $\lambda = 0$ and $\mu = 1$ to highlight the important concepts. Therefore,
 * let consider the following weak formulation: find $\mathbf u \in
 * \mathbf V = H^1_0(\Omega; \mathbb R^3)$ such that for all $\mathbf
 * v\in V$ holds
 * @f[
 * a(u,v) \equiv 2\int_{\Omega} \mathbf D\mathbf u : \mathbf D\mathbf
 * v\,dx = \int_\Omega \mathbf f\cdot \mathbf v \,dx.
 * @f]
 * Here, <b>D</b> denotes the symmetric gradient defined by
 * $\mathbf Du = \tfrac12 (\nabla \mathbf u + (\nabla \mathbf u)^T)$
 * and the colon indicates double contraction of two tensors of rank 2
 * (the Frobenius inner product). This bilinear form looks indeed very
 * much like the bilinear form of the Poisson problem in step-3. The
 * only differences are
 * <ol>
 * <li>We replaced the gradient operator by the symmetric gradient;
 * this is actually not a significant difference, and everything said
 * here is true if your replace $\mathbf D$ by $\nabla$. Indeed, let
 * us do this to simplify the discussion:
 * @f[
 * a(u,v) \equiv \int_{\Omega} \nabla\mathbf u : \nabla\mathbf
 * v\,dx = \int_\Omega \mathbf f\cdot \mathbf v \,dx.
 * @f]
 * Note though, that this system is not very exciting, since we could
 * solve for the three components of <b>u</b> separately.
 *
 * <li> The trial and test functions are now from the space
 * $H^1_0(\Omega; \mathbb R^3)$, which can be viewed as three copies
 * of the scalar space $H^1_0(\Omega)$. And this is exactly, how we
 * are going to implement this space below, using FESystem.
 * </ol>
 *
 * But for now, let us look at the system a little more
 * closely. First, let us exploit that
 * <b>u</b>=(<i>u</i><sub>1</sub>,<i>u</i><sub>2</sub>,<i>u</i><sub>3</sub>)<sup>T</sup>
 * and <b>v</b> accordingly. Then, we can write the simplified equation in
 * coordinates as
 * @f[
 * a(u,v) = \int_\Omega \bigl(\nabla u_1\cdot \nabla v_1
 +\nabla u_2\cdot \nabla v_2+\nabla u_3\cdot \nabla v_3\bigr)\,dx
 = \int_\Omega \bigl(f_1v_1 + f_2 v_2 + f_3 v_3\bigr)\,dx.
 * @f]
 * We see, that this is just three copies of the bilinear form of the
 * Laplacian, one applied to each component (this is where the
 * formulation with the $\mathbf D$ is more exciting, and we want to
 * derive a framework that applies to that one as well). We can make
 * this weak form a system of differential equations again by choosing
 * special test functions: first, choose
 * <b>v</b>=(<i>v</i><sub>1</sub>,0,0)<sup>T</sup>, then
 * <b>v</b>=(0,<i>v</i><sub>2</sub>,0)<sup>T</sup>, and finally
 * <b>v</b>=(0,0,<i>v</i><sub>3</sub>)<sup>T</sup>. writing the outcomes below
 * each other, we obtain the system
 * @f[
 * \begin{matrix}
 * (\nabla u_1,\nabla v_1) &&& = (f_1, v_1)
 * \\
 * & (\nabla u_2,\nabla v_2) && = (f_2, v_2)
 * \\
 * && (\nabla u_3,\nabla v_3) & = (f_3, v_3)
 * \end{matrix}
 * @f]
 * where we used the standard inner product notation $(\mathbf
 * f,\mathbf g) =
 * \int_\Omega \mathbf f \cdot \mathbf g \,dx$. It is important for our understanding, that
 * we keep in mind that the latter form as a system of PDE is
 * completely equivalent to the original definition of the bilinear
 * form <i>a</i>(<i>u</i>,<i>v</i>), which does not immediately
 * exhibit this system structure. Let us close by writing the full
 * system of the elastic equation with symmetric gradient <b>D</b>:
 * @f[
 * \begin{matrix}
 * (\nabla u_1,\nabla v_1) + (\partial_1 u_1,\partial_1 v_1)
 * & (\partial_1 u_2,\partial_2 v_1)
 * & (\partial_1 u_3,\partial_3 v_1)
 * & = (f_1, v_1)
 * \\
 * (\partial_2 u_1,\partial_1 v_2)
 * & (\nabla u_2,\nabla v_2) + (\partial_2 u_2,\partial_2 v_2)
 * & (\partial_2 u_3,\partial_3 v_2)
 * & = (f_2, v_2)
 * \\
 * (\partial_3 u_1,\partial_1 v_3)
 * & (\partial_3 u_2,\partial_2 v_3)
 * & (\nabla u_3,\nabla v_3) + (\partial_3 u_3,\partial_3 v_3)
 * & = (f_3, v_3)
 * \end{matrix}.
 * @f]
 * Very formally, if we believe in operator valued matrices, we can
 * rewrite this in the form <b>v</b><sup>T</sup><b>Au</b> = <b>v</b><sup>T</sup><b>f</b> or
 * @f[
 * \begin{pmatrix}
 * v_1 \\ v_2 \\ v_3
 * \end{pmatrix}^T
 * \begin{pmatrix}
 * (\nabla \cdot,\nabla \cdot) + (\partial_1 \cdot,\partial_1 \cdot)
 * & (\partial_1 \cdot,\partial_2 \cdot)
 * & (\partial_1 \cdot,\partial_3 \cdot)
 * \\
 * (\partial_2 \cdot,\partial_1 \cdot)
 * & (\nabla \cdot,\nabla \cdot) + (\partial_2 \cdot,\partial_2 \cdot)
 * & (\partial_2 \cdot,\partial_3 \cdot)
 * \\
 * (\partial_3 \cdot,\partial_1 \cdot)
 * & (\partial_3 \cdot,\partial_2 \cdot)
 * & (\nabla \cdot,\nabla \cdot) + (\partial_3 \cdot,\partial_3 \cdot)
 * \end{pmatrix}
 * \begin{pmatrix}
 * u_1 \\ u_2 \\ u_3
 * \end{pmatrix}
 * =
 * \begin{pmatrix}
 * v_1 \\ v_2 \\ v_3
 * \end{pmatrix}^T
 * \begin{pmatrix} f_1 \\ f_2 \\ f_3\end{pmatrix}
 * @f]
 *
 * <h4>Mixed elliptic problems</h4>
 * Now, let us
 * consider a more complex example, the mixed Laplace equations discussed in
 * step-20 in three dimensions:
@f{eqnarray*}
  \textbf{u} + \nabla p &=& 0,
  \\
  -\textrm{div}\; \textbf{u} &=& f,
@f}
 *
 * Here, we have four solution components: the scalar pressure
 * $p \in L^2(\Omega)$ and the vector-valued velocity $\mathbf u \in
 * \mathbf V
 * = H^{\text{div}}_0(\Omega)$ with three vector
 * components. Note as important difference to the previous example,
 * that the vector space <b>V</b> is not just simply a copy of three
 * identical spaces/
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
 * indeed has four components. We note that we could change the
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
*
* We get the final form by integrating by part the second term:
@f{eqnarray*}
  (\mathbf v, \mathbf u)
  -
  (\mathrm{div}\ \mathbf v, p)
  -
  (q, \mathrm{div}\ \mathbf u)
  =
  (q,f) - (\mathbf n\cdot\mathbf v, p)_{\partial\Omega}.
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
 * Once we have settled on a bilinear form and a functional setting, we need to find a way to describe
 * the vector-valued finite element spaces from which we draw solution and test
 * functions. This is where the FESystem class comes in: it composes
 * vector-valued finite element spaces from simpler ones.
 * In the example of the elasticity problem, we need <code>dim</code>
 * copies of the same element, for instance
 * @code
 *   FESystem<dim> elasticity_element (FE_Q<dim>(1), dim);
 * @endcode
 * This will generate a vector valued space of dimension
 * <code>dim</code>, where each component is a
 * continuous bilinear element of type FE_Q. It will have <code>dim</code> times
 * as many basis functions as the corresponding FE_Q, and each of
 * these basis functions is a basis function of FE_Q, lifted into one
 * of the components of the vector.
 *
 * For the mixed Laplacian, the situation is more complex. First, we
 * have to settle on a pair of discrete spaces $\mathbf V_h \times Q_h
 * \subset H^{\text{div}}_0(\Omega) \times L^2_0(\Omega)$. One option
 * would be the stable Raviart-Thomas pair
 * @code
 *   FESystem<dim> rt_element (FE_RaviartThomas<dim>(1), 1,
 *                             FE_DGQ<dim>(1),          1);
 * @endcode
 * The first element in this system is already a vector valued
 * element of dimension <code>dim</code>, while the second is a
 * regular scalar element.
 *
 * Alternatively to using the stable Raviart-Thomas pair, we could
 * consider a stabilized formulation for the mixed Laplacian, for
 * instance the LDG method. There, we have the option of using the
 * same spaces for velocity components and pressure, namely
 * @code
 *   FESystem<dim> ldg_convoluted_element_1 (FE_DGQ<dim>(1), dim+1);
 * @endcode
 * This system just has <code>dim+1</code> equal copies of the same
 * discontinuous element, which not really reflects the structure of
 * the system. Therefore, we prefer
 * @code
 *   FESystem<dim> ldg_equal_element (FESystem<dim>(FE_DGQ<dim>(1), dim), 1,
 *                                    FE_DGQ<dim>(1),                     1);
 * @endcode
 * Here, we have a system of two elements, one vector-valued and one
 * scalar, very much like with the <code>rt_element</code>. Indeed, in
 * many codes, the two can be interchanged. This element also allows
 * us easily to switch to an LDG method with lower order approximation
 * in the velocity, namely
 * @code
 *   FESystem<dim> ldg_unequal_element (FESystem<dim>(FE_DGQ<dim>(1), dim), 1,
 *                                      FE_DGQ<dim>(2),                     1);
 * @endcode
 * It must be pointed out,
 * that this element is different from
 * @code
 *   FESystem<dim> ldg_convoluted_element_2 (FE_DGQ<dim>(1), dim,
 *                                           FE_DGQ<dim>(2), 1);
 * @endcode
 * While the constructor call is very similar to
 * <code>rt_element</code>, the result actually resembles more
 * <code>ldg_convoluted_element_1</code> in that this element produces
 * <code>dim+1</code> independent components.
 * A more detailed comparison of the resulting FESystem objects is below.
 *
 * <h4>Internal structure of FESystem</h4>
 *
 * FESystem has a few internal variables which reflect the internal
 * structure set up by the constructor. These can then also be used by
 * application programs to give structure to matrix assembling and
 * linear algebra. We give the names and values of these variables for
 * the examples above in the following table:
 * <table border="1">
 * <tr><th>System Element</th>
 * <th>FiniteElementData::n_blocks()</th>
 * <th>FiniteElementData::n_components()</th>
 * <th>FiniteElement::n_base_elements()</th>
 * </tr>
 * <tr><td><code>elasticity_element</code></td><td><code>dim</code></td><td><code>dim</code></td><td>1</td>
 * </tr>
 * <tr><td><code>rt_element</code></td><td>2</td><td><code>dim+1</code></td><td>2</td>
 * </tr>
 * <tr><td><code>ldg_equal_element</code></td><td>2</td><td><code>dim+1</code></td><td>2</td>
 * </tr>
 * <tr><td><code>ldg_convoluted_element_1</code></td><td><code>dim+1</code></td><td><code>dim+1</code></td><td>1</td>
 * </tr>
 * <tr><td><code>ldg_convoluted_element_2</code></td><td><code>dim+1</code></td><td><code>dim+1</code></td><td>2</td>
 * </tr>
 * </table>
 *
 * From this table, it is clear that the FESystem reflects a lot of
 * the structure of the system of differential equations in the cases
 * of the <code>rt_element</code> and the
 * <code>ldg_equal_element</code>, in that we have a vector valued and
 * a scalar variable. On the other hand, the convoluted elements do
 * not have this structure and we have to reconstruct it somehow when
 * assembling systems, as described below.
 *
 * At this point, it is important to note that nesting of two FESystem
 * object can give the whole FESystem a richer structure than just
 * concatenating them. This structure can be exploited by application
 * programs, but is not automatically so.
 *
 * @anchor VVAssembling
 * <h3>Assembling linear systems</h3>
 * The next step is to assemble the linear system. How to do this for the
 * simple case of a scalar problem has been shown in many tutorial programs,
 * starting with step-3. Here we will show how to do it for
 * vector problems. Corresponding to the different characterizations
 * of weak formulations above and the different system elements
 * created, we have a few options which are outlined below.
 *
 * The whole concept is probably best explained by showing an example
 * illustrating how the local contribution of a cell to the weak form of above
 * mixed Laplace equations could be assembled.
 *
 * <h4>A single FEValues and FEValuesExtractors</h4>
 * This is essentially how
 * step-20 does it:
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

            local_rhs(i) += fe_values[pressure].value (i, q)
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
 *        $V_i=\left(\begin{array}{c}\mathbf v_i \\ q_i\end{array}\right),
 *         V_j=\left(\begin{array}{c}\mathbf v_j \\ q_j\end{array}\right)$:
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
                                    fe_values[velocities].divergence (j, q)
                                   ) *
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
 *        velocities components, <code>fe_values[velocities].gradient(i,q)</code>
 *        is a <code>Tensor@<2,dim@></code>, i.e. a matrix $G_{ij}$ that consits
 *        of entries $G_{ij}=\frac{\partial\phi_i}{\partial x_j}$. Finally,
 *        both scalar and vector views can be asked for the second derivatives
 *        ("Hessians") and vector views can be asked for the symmetric gradient,
 *        defined as $S_{ij}=\frac 12 \left[\frac{\partial\phi_i}{\partial x_j}
 *        + \frac{\partial\phi_j}{\partial x_i}\right]$ as well as the
 *        divergence $\sum_{d=0}^{dim-1} \frac{\partial\phi_d}{\partial x_d}$.
 * </ul>
 * Other examples of using extractors and views are shown in tutorial programs
 * step-21,
 * step-22,
 * step-31 and several other programs.
 *
 * @note In the current context, when we talk about a vector (for example in
 * extracting the velocity components above), we mean the word in the sense
 * physics uses it: it has <code>spacedim</code> components that behave in
 * specific ways under coordinate system transformations. Examples include
 * velocity or displacement fields. This is opposed to how mathematics uses
 * the word "vector" (and how we use this word in other contexts in the
 * library, for example in the Vector class), where it really stands for a
 * collection of numbers. An example of this latter use of the word could be
 * the set of concentrations of chemical species in a flame; however, these
 * are really just a collection of scalar variables, since they do not change
 * if the coordinate system is rotated, unlike the components of a velocity
 * vector, and consequently, this FEValuesExtractors::Vector class should not
 * be used for this case.
 *
 *
 * @anchor VVAlternative
 * <h3>An alternative approach</h3>
 *
 * There are situations in which one can optimize the assembly of a matrix or
 * right hand side vector a bit, using knowledge of the finite element in
 * use. Consider, for example, the bilinear form of the elasticity equations
 * which we are concerned with first in step-8:
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
 * rather than disjoint scalar components:
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
                  +=  (lambda_values[q_point] *
                       phi_i_div * phi_j_div
                       +
                       mu_values[q_point] *
                       double_contract(phi_i_grad, phi_j_grad)
                       +
                       mu_values[q_point] *
                       double_contract(phi_i_grad, transpose(phi_j_grad))
                      ) *
                      fe_values.JxW(q_point);
              }
          }
 * @endcode
 *
 * Now, this is not the code used in step-8. In fact,
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
 * [\nabla\mathbf u]^T\right)$ is the symmetrized gradient.
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
                       mu_values[q_point]) *
                      fe_values.JxW(q_point));
              }
          }
 * @endcode
 *
 * So if, again, this is not the code we use in step-8, what do
 * we do there? The answer rests on the finite element we use. In step-8, we use the
 * following element:
 * @code
 *   FESystem<dim> finite_element (FE_Q<dim>(1), dim);
 * @endcode
 * In other words, the finite element we use consists of <code>dim</code> copies
 * of the same scalar element. This is what we call a @ref GlossPrimitive
 * "primitive" element: an element that may be vector-valued but where each
 * shape function has exactly one non-zero component. In other words: if the
 * $x$-component of a displacement shape function is nonzero, then the $y$-
 * and $z$-components must be zero and similarly for the other components.
 * What this means is that also
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
 * step-8.
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
 * cases like the elasticity equations covered in step-8 and
 * step-17). It will, often, also not be symmetric. This rather
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
 * "step-20", step-21, step-22, and @ref step_31
 * "step-31". After this, at least the degrees of freedom are partitioned
 * properly.
 *
 * But then we still have to make use of it, i.e. we have to come up with a
 * solver that uses the structure. For example, in step-20, we
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
 * step_20 "step-20", step-31, and a few other tutorial
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
 * what we have in step-20:
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
 * example, in step-21 or step-31). In that
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
 * For vector- (tensor-)valued quantities, there are functions
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
 *
 * @anchor VVOutput
 * <h3>Generating graphical output</h3>
 *
 * As mentioned above, an FESystem object may hold multiple vector components,
 * but it doesn't have a notion what they actually mean. As an example, take
 * the object
 * @code
 *   FESystem<dim> finite_element (FE_Q<dim>(1), dim+1);
 * @endcode
 * It has <code>dim+1</code> vector components, but what do they mean? Are they
 * the <code>dim</code> components of a velocity vector plus one pressure? Are
 * they the pressure plus the <code>dim</code> velocity components? Or are
 * they a collection of scalars?
 *
 * The point is that the FESystem class doesn't care. The <i>interpretation</i>
 * of what the components mean is up to the person who uses the element later,
 * for example in assembling a linear form, or in extracting data solution
 * components for a linearized system in the next Newton step. In almost
 * all cases, this interpretation happens at the place where it is needed.
 *
 * There is one case where one has to be explicit, however, and that is in
 * generating graphical output. The reason is that many file formats for
 * visualization want data that represents vectors (e.g. velocities,
 * displacements, etc) to be stored separately from scalars (pressures,
 * densities, etc), and there often is no way to group a bunch of scalars into
 * a vector field from within a visualization program.
 *
 * To achieve this, we need to let the DataOut class and friends know which
 * components of the FESystem form vectors (with <code>dim</code> components)
 * and which are scalars. This is shown, for example, in step-22 where we
 * generate output as follows:
 * @code
  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, solution_names,
                            DataOut<dim>::type_dof_data,
                            data_component_interpretation);
  data_out.build_patches ();
 * @endcode
 * In other words, we here create an array of <code>dim+1</code> elements in
 * which we store which elements of the finite element are vectors and which
 * are scalars; the array is filled with <code>dim</code> copies of
 * DataComponentInterpretation::component_is_part_of_vector and a single
 * trailing element of DataComponentInterpretation::component_is_scalar . The
 * array is then given as an extra argument to DataOut::add_data_vector to
 * explain how the data in the given solution vector is to be interpreted.
 * Visualization programs like Visit and Paraview will then offer to show
 * these <code>dim</code> components as vector fields, rather than as
 * individual scalar fields.
 *
 *
 * @ingroup feall feaccess
 */

