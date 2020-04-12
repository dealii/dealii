// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

/**
 * @defgroup auto_symb_diff Automatic and symbolic differentiation
 *
 * @brief A module dedicated to the implementation of functions and classes that relate
 * to automatic and symbolic differentiation.
 *
 * Below we provide a very brief introduction as to what automatic and symbolic differentiation are,
 * what variations of these computational/numerical schemes exist, and how they are integrated
 * within deal.II's framework. The purpose of all of these schemes is to automatically compute
 * the derivative of functions, or approximations of it, in cases where one does not want to
 * compute them by hand. Common examples are situations in the finite element context is where
 * one wants to solve a nonlinear problem that is given by requiring that some residual
 * $F(u,\nabla u)=0$ where $F$ is a complicated function that needs to be differentiated to
 * apply Newton's method; and situations where one is given a parameter dependent problem
 * ${\cal A}(q,u,\nabla u) = f$ and wants to form derivatives with regards to the parameters $q$, for example
 * to optimize an output functional with regards to $q$, or for a sensitivity analysis with
 * regards to $q$. One should think of $q$ as design parameters: say, the width
 * or shape of a wing, the stiffness coefficients of a material chosen to
 * build an object, the power sent to a device, the chemical composition of the
 * gases sent to a burner. In all of these cases, one should think of $F$ and $\cal A$ as <i>complicated</i>
 * and cumbersome to differentiate -- at least when doing it by hand. A relatively simple case of
 * a nonlinear problem that already highlights the tedium of computing derivatives by hand is shown in
 * step-15. However, in reality, one might, for example,
 * think about problems such as chemically reactive flows where the fluid equations have coefficients
 * such as the density and viscosity that depend strongly and nonlinearly on the chemical composition,
 * temperature, and pressure of the fluid at each point; and where the chemical species react with
 * each other based on reaction coefficients that also depend nonlinearly and in complicated
 * ways on the chemical composition, temperature, and pressure. In many cases, the exact formulas
 * for all of these coefficients can take several lines to write out, may include exponentials
 * and (harmonic or geometric) averages of several nonlinear terms, and/or may contain table
 * lookup of and interpolation between data points. Just getting these terms right is difficult
 * enough; computing derivatives of these terms is impractical in most applications and, in
 * reality, impossible to get right. Higher derivatives are even more impossible to do
 * without computer aid. Automatic or symbolic differentiation is a way out of this:
 * One only has to implement the function that computes these coefficients in terms
 * of their inputs only once, and gets the (correct!) derivatives without
 * further coding effort (though at a non-negligible computational cost either at run time, compile
 * time, or both).
 *
 *
 * @section auto_diff_1 Automatic differentiation
 *
 * <a href="https://en.wikipedia.org/wiki/Automatic_differentiation">Automatic differentiation </a>
 * (commonly also referred to as algorithmic differentiation),
 * is a numerical method that can be used to "automatically" compute the first, and perhaps higher-order,
 * derivatives of function(s) with respect to one or more input variables.
 * Although this comes at a certain computational cost, the benefits to using such a tool may be
 * significant. When used correctly the derivatives of often complicated functions can be computed
 * to a very high accuracy. Although the exact accuracy achievable by these frameworks largely
 * depends on their underlying mathematical formulation, some implementations compute with a precision
 * on the order of machine accuracy. Note that this is different to classical numerical
 * differentiation (using, for example, a finite difference approximation of a function by
 * evaluating it at different points),
 * which has an accuracy that depends on both the perturbation size as well as the chosen
 * finite-difference scheme; the error of these methods is measurably larger than
 * well-formulated automatic differentiation approaches.
 *
 * Three practical examples of auto-differentiation use within a finite-element context
 * would then be
 * - the quick prototyping of a new nonlinear formulation without the need to hand-compute
 *   the linearization itself,
 * - automatic linearization of finite-element residuals additively formed within complex
 *   multiphysics frameworks, and
 * - verification of user-implementations of linearizations for both cell-based calculations
 *   (e.g. a residual) and those based at a continuum point (e.g. tangents for nonlinear
 *   constitutive laws).
 *
 * There are quite a number of implementations for auto-differentiable numbers. They primarily
 * fall into two broad categories, namely <em>source code transformation</em> and
 * <em>operator overloading</em>.
 * The first method generates new, compilable code based on some input function that, when executed,
 * returns the derivatives of the input function. The second exploits the capability of <tt>C++</tt>
 * operator definitions to be overloaded for custom class types. Therefore  a class that represents
 * such an auto-differentiable number can, following each mathematical operation performed on or
 * with it, in principle evaluate and keep track of its value as well as that of its directional
 * derivative(s).
 * As the libraries exclusively implementing the <em>source code transformation</em> approach
 * collectively describe highly specialized tools that are to be used as function preprocessors, they
 * have no direct support within deal.II itself. The latter, however, represent specialized number
 * types that can be supported through the use of template  metaprogramming in the appropriate context.
 * Given the examples presented above, this means that the FEValues class (and friends), as well as
 * the Tensor and SymmetricTensor classes should support calculations performed with these specialized
 * numbers.
 * (In theory an entire program could be made differentiable. This could be useful in, for example,
 * the sensitivity analysis of solutions with respect to input parameters. However, to date this has
 * not been tested.)
 *
 * Implementations of specialized frameworks based on <em>operator overloading</em> typically fall into
 * one of three categories. In each, some customized data classes representing the floating point value
 * of an evaluated function and its derivative(s) by
 * -# exploiting <em>dual</em>/<em>complex-step</em>/<em>hyper-dual</em> formulations (occasionally
 *    called <em>tapeless</em> methods),
 * -# those utilizing <em>taping</em> strategies, and
 * -# those using compile-time optimization through <em>expression templates</em>.
 *
 * To provide some tentative insight into how these various implementations might look like in practice, we
 * offer the following generic summary of these approaches:
 * -# The first two <em>tapeless</em> approaches listed above (dual numbers and complex-step method) use some
 *    variation of a truncated Taylor series, along with a particular choice of definition for the perturbation
 *    parameter, to compute function derivatives using a finite-difference based approach. The "dual" number
 *    constitutes the accumulated directional derivatives computed simultaneously as the function values are
 *    evaluated; in the complex-step approach, the imaginary value effectively serves this purpose. The choice of
 *    the perturbation parameter determines the numerical qualities of the scheme, such as the influence of the
 *    truncation of the Taylor scheme; dual numbers do not contain any higher-order terms in their first derivative,
 *    while for the complex-step method these existent higher-order terms are neglected. It can be shown that
 *    both of these methods are not subject to subtractive cancellation errors and that, within their
 *    finite-difference scheme, they are not numerically sensitive to the internal \step-size chosen for the
 *    numerical perturbation. The dual number approach thus produces exact first derivatives, while the
 *    complex-step approximation does not. The standard implementation of the dual numbers, however, cannot yield
 *    exact values for second derivatives. Hyper-dual numbers take a different view of this idea, with numbers
 *    being represented in a form similar to quaternions (i.e. carrying additional non-real components) and the
 *    derivatives being computed from a high-order truncation of the Taylor series all four components. The outcome
 *    is that, with the appropriate implementation, both first and second derivatives can be computed exactly.
 * -# With <em>taped</em> approaches, a specified subregion of code is selected as one for which all
 *    operations executed with active (marked) input variables are tracked and recorded in a data structure
 *    referred to as a tape. At the end of the taped region, the recorded function(s) may be reevaluated
 *    by "replaying" the tape with a different set of input variables instead of recomputing the function
 *    directly. Assuming that the taped region represents a smooth function, arbitrarily high-order
 *    derivatives of the function then can be computed by referring to the code path tracked and stored on
 *    the tape.
 *    (This could perhaps be achieved, for example, through evaluation of the function around the point
 *    of interest.) There exist strategies to deal with situations where the taped function is not
 *    smooth at the evaluated point, or if it is not analytic. Furthermore, one might need to consider the
 *    case of branched functions, where the tape is no longer sequential, but rather forks off on a different
 *    evaluation path to that due to the original recorded inputs.
 * -# Methods based on <a href="https://en.wikipedia.org/wiki/Expression_templates">expression templates</a>
 *    leverage the computational graph
 *    (in this case, a <a href="https://en.wikipedia.org/wiki/Directed_acyclic_graph">directed acyclic graph (DAG)</a>),
 *    constructed from the abstract syntax tree (AST), that resolves the function output from its input values.
 *    The outermost leaves on the tree represent the independent variables or constants, and are transformed by unary
 *    operators and connected by binary operators (in the most simple case). Therefore, the operations performed on
 *    the function inputs is known at compile time, and with that the associated derivative operation can also be defined
 *    at the same time using the well-known rules of computing the derivative of an operation (such as
 *    the associativity of derivatives under addition and subtraction, the product rule, and the chain
 *    rule). The compiled output type returned by this operator need not be generic, but can rather be
 *    specialized based on the specific inputs (possibly carrying a differential history) given to that specific
 *    operator on the vertex of the DAG. In this way, a compile-time optimized set of instructions can be generated
 *    for the very specialized individual operations used to evaluate each intermediate result of the dependent
 *    function.
 *
 * Each of these methods, of course, has its advantages and disadvantages, and one may be more appropriate
 * than another for a given problem that is to be solved. As the aforementioned implementational details
 * (and others not discussed) may be hidden from the user, it may still be important to understand the
 * implications, run-time cost,  and potential limitations, of using any one of these "black-box"
 * auto-differentiable numbers.
 *
 * In addition to the supplied linked articles, resources used to furnish the details supplied here include:
 *
 * @code{.bib}
 * @InProceedings{Fike2011a,
 *   author    = {Fike, Jeffrey A and Alonso, Juan J},
 *   title     = {The Development of Hyper-Dual Numbers for Exact Second-Derivative Calculations},
 *   booktitle = {49th {AIAA} Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition},
 *   year      = {2011},
 *   volume    = {886},
 *   pages     = {124},
 *   month     = {jan},
 *   publisher = {American Institute of Aeronautics and Astronautics},
 *   doi       = {10.2514/6.2011-886},
 * }
 * @endcode
 *
 * @code{.bib}
 * @Manual{Walther2009a,
 *   title     = {Getting Started with ADOL-C},
 *   author    = {Walther, Andrea and Griewank, Andreas},
 *   year      = {2009},
 *   booktitle = {Combinatorial scientific computing},
 *   doi       = {10.1.1.210.4834},
 *   pages     = {181--202}
 * }
 * @endcode
 *
 * ### Exploitation of the chain-rule
 *
 * In the most practical sense, any of the above categories exploit the chain-rule to compute the total
 * derivative of a composite function. To perform this action, they typically use one of two mechanisms to
 * compute derivatives, specifically
 * - <em>forward-mode</em> (or <em>forward accumulation</em>) auto-differentiation, or
 * - <em>reverse-mode</em> (or <em>reverse accumulation</em>) auto-differentiation.
 *
 * As a point of interest, the <em>optimal Jacobian accumulation</em>, which performs a minimal set of
 * computations, lies somewhere between these two limiting cases. Its computation for a general composite
 * function remains an open problem in graph theory.
 *
 * With the aid of the diagram below (it and some of the listed details courtesy of this
 * <a href="https://en.wikipedia.org/wiki/Automatic_differentiation">Wikipedia article</a>),
 * let us think about the represention of the calculation of the function
 * $f (\mathbf{x}) = \sin (x_{1}) + x_{1} x_{2}$ and its derivatives:
 *
 * <div class="twocolumn" style="width: 80%">
 *   <div class="parent">
 *     <div class="img" align="center">
 *       <img src="https://upload.wikimedia.org/wikipedia/commons/a/a4/ForwardAccumulationAutomaticDifferentiation.png"
 *            alt="Forward mode automatic differentiation"
 *            width="400">
 *     </div>
 *     <div class="text" align="center">
 *       Forward mode automatic differentiation
 *     </div>
 *   </div>
 *   <div class="parent">
 *     <div class="img" align="center">
 *       <img src="https://upload.wikimedia.org/wikipedia/commons/a/a0/ReverseaccumulationAD.png"
 *            alt="Reverse mode automatic differentiation"
 *            width="400">
 *     </div>
 *     <div class="text" align="center">
 *       Reverse mode automatic differentiation
 *     </div>
 *   </div>
 * </div>
 *
 * Specifically, we will briefly describe what forward and reverse auto-differentiation are.
 * Note that in the diagram, along the edges of the graph in text are the directional
 * derivative of function $w$ with respect to the $i$-th variable, represented by
 * the notation $\dot{w} = \dfrac{d w}{d x_{i}}$.
 * The specific computations used to render the function value and its directional derivatives
 * for this example are tabulated in the
 * <a href="https://en.wikipedia.org/wiki/Automatic_differentiation">source article</a>.
 * For a second illustrative example, we refer the interested reader to
 * <a href="http://www.columbia.edu/~ahd2125/post/2015/12/5/">this article</a>.
 *
 * Consider first that any composite function $f(x)$, here represented as having two
 * independent variables, can be dissected into a composition of its elementary functions
 * @f[
 *   f (\mathbf{x})
 *   = f_{0} \circ f_{1} \circ f_{2} \circ \ldots \circ f_{n} (\mathbf{x})
 *   \quad .
 * @f]
 * As was previously mentioned, if each of the primitive operations $f_{n}$ is smooth and
 * differentiable, then the chain-rule can be universally employed to compute the total derivative of $f$,
 * namely $\dfrac{d f(x)}{d \mathbf{x}}$. What distinguishes the "forward" from the "reverse" mode
 * is how the chain-rule is evaluated, but ultimately both compute the total derivative
 * @f[
 *   \dfrac{d f (\mathbf{x})}{d \mathbf{x}}
 *   = \dfrac{d f_{0}}{d f_{1}} \dfrac{d f_{1}}{d f_{2}} \dfrac{d f_{2}}{d f_{3}} \ldots \dfrac{d f_{n} (\mathbf{x})}{d \mathbf{x}}
 *   \quad .
 * @f]
 *
 * In forward-mode, the chain-rule is computed naturally from the "inside out". The independent
 * variables are therefore fixed, and each sub-function $f'_{i} \vert_{f'_{i+1}}$ is computed
 * recursively and its result returned as inputs to the parent function. Encapsulating and fixing
 * the order of operations using parentheses, this means that we compute
 * @f[
 *   \dfrac{d f (\mathbf{x})}{d \mathbf{x}}
 *   = \dfrac{d f_{0}}{d f_{1}} \left( \dfrac{d f_{1}}{d f_{2}} \left(\dfrac{d f_{2}}{d f_{3}} \left(\ldots \left( \dfrac{d f_{n} (\mathbf{x})}{d \mathbf{x}} \right)\right)\right)\right)
 *   \quad .
 * @f]
 * The computational complexity of a forward-sweep is proportional to that of the input function.
 * However, for each directional derivative that is to be computed one sweep of the computational
 * graph is required.
 *
 * In reverse-mode, the chain-rule is computed somewhat unnaturally from the "outside in". The
 * values of the dependent variables first get computed and fixed, and then the preceding
 * differential operations are evaluated and multiplied in succession with the previous results
 * from left to right. Again, if we encapsulate and fix the order of operations using parentheses,
 * this implies that the reverse calculation is performed by
 * @f[
 * \dfrac{d f (\mathbf{x})}{d \mathbf{x}}
 *   = \left( \left( \left( \left( \left( \dfrac{d f_{0}}{d f_{1}} \right) \dfrac{d f_{1}}{d f_{2}} \right) \dfrac{d f_{2}}{d f_{3}} \right) \ldots \right) \dfrac{d f_{n} (\mathbf{x})}{d \mathbf{x}} \right)
 *   \quad .
 * @f]
 * The intermediate values $\dfrac{d f_{i-1}}{d f_{i}}$ are known as <em>adjoints</em>, which must be
 * computed and stored as the computational graph is traversed. However, for each dependent scalar function
 * one sweep of the computational graph renders all directional derivatives at once.
 *
 * Overall, the efficiency of each mode is determined by the number of independent (input) variables
 * and dependent (output) variables. If the outputs greatly exceed the inputs in number, then
 * forward-mode can be shown to be more efficient than reverse-mode. The converse is true when the
 * number of input variables greatly exceeds that of the output variables. This point may be used to
 * help inform which number type is most suitable for which set of operations are to be performed
 * using automatic differentiation. For example, in many applications for which second derivatives
 * are to be computed it is appropriate to combine both reverse- and forward-modes. The former would
 * then typically be used to calculate the first derivatives, and the latter the second derivatives.
 *
 * @subsection auto_diff_1_1 Supported automatic differentiation libraries
 *
 * We currently have validated implementations for the following number types
 * and combinations:
 *
 *  - Taped ADOL-C (n-differentiable, in theory, but internal drivers for up to second-order
 *    derivatives will be implemented)
 *  - Tapeless ADOL-C (once differentiable)
 *  - Forward-mode Sacado with dynamic memory allocation using expression templates (once differentiable)
 *  - Nested forward-mode Sacado using expression templates (twice differentiable)
 *  - Reverse-mode Sacado (once differentiable)
 *  - Nested reverse and dynamically-allocated forward-mode Sacado (twice differentiable)
 *
 * Note that in the above, "dynamic memory allocation" refers to the fact that the number of
 * independent variables need not be specified at compile time.
 *
 * The <a href="https://projects.coin-or.org/ADOL-C/browser/trunk/ADOL-C/doc/adolc-manual.pdf?format=raw">ADOL-C user manual</a>
 *
 * @code{.bib}
 * @Manual{Walther2009a,
 *   title     = {Getting Started with ADOL-C},
 *   author    = {Walther, Andrea and Griewank, Andreas},
 *   year      = {2009},
 *   booktitle = {Combinatorial scientific computing},
 *   doi       = {10.1.1.210.4834},
 *   pages     = {181--202},
 *   url       = {https://projects.coin-or.org/ADOL-C/browser/trunk/ADOL-C/doc/adolc-manual.pdf}
 * }
 * @endcode
 *
 * provides the principle insights into their taped and tapeless implementations, and how ADOL-C
 * can be incorporated into a user code.
 * Some further useful resources for understanding the implementation of ADOL-C, and possibilities
 * for how it may be used within a numerical code, include:
 *
 * @code{.bib}
 * @Article{Griewank1996a,
 *   author    = {Griewank, Andreas and Juedes, David and Utke, Jean},
 *   title     = {Algorithm 755: {ADOL-C}: a package for the automatic differentiation of algorithms written in {C/C++}},
 *   journal   = {ACM Transactions on Mathematical Software (TOMS)},
 *   year      = {1996},
 *   volume    = {22},
 *   number    = {2},
 *   pages     = {131--167},
 *   doi       = {10.1145/229473.229474},
 *   publisher = {ACM}
 * }
 * @endcode
 * @code{.bib}
 * @InCollection{Bischof2008a,
 *   author =    {Bischof, Christian and Guertler, Niels and Kowarz, Andreas and Walther, Andrea},
 *   title =     {Parallel reverse mode automatic differentiation for OpenMP programs with ADOL-C},
 *   booktitle = {Advances in Automatic Differentiation},
 *   publisher = {Springer},
 *   year =      {2008},
 *   pages =     {163--173}
 * }
 * @endcode
 * @code{.bib}
 * @InBook{Kulshreshtha2012a,
 *   chapter   = {Computing Derivatives in a Meshless Simulation Using Permutations in {ADOL}-C},
 *   pages     = {321--331},
 *   title     = {Recent Advances in Algorithmic Differentiation},
 *   publisher = {Springer Berlin Heidelberg},
 *   year      = {2012},
 *   author    = {Kshitij Kulshreshtha and Jan Marburger},
 *   editor    = {Forth S. and Hovland P. and Phipps E. and Utke J. and Walther A.},
 *   series    = {Lecture Notes in Computational Science and Engineering},
 *   doi       = {10.1007/978-3-642-30023-3_29},
 * }
 * @endcode
 * @code{.bib}
 * @InProceedings{Kulshreshtha2013a,
 *   author    = {Kulshreshtha, Kshitij and Koniaeva, Alina},
 *   title     = {Vectorizing the forward mode of ADOL-C on a GPU using CUDA},
 *   booktitle = {13th European AD Workshop},
 *   year      = {2013},
 *   month     = jun
 * }
 * @endcode
 *
 * Similarly, a selection of useful resources for understanding the implementation of Sacado
 * number types (in particular, how expression templating is employed and exploited) include:
 *
 * @code{.bib}
 * @InCollection{Bartlett2006a,
 *   author        = {Bartlett, R. A. and Gay, D. M. and Phipps, E. T.},
 *   title         = {Automatic Differentiation of C++ Codes for Large-Scale Scientific Computing},
 *   booktitle     = {International Conference on Computational Science {\textendash} {ICCS} 2006},
 *   publisher     = {Springer Berlin Heidelberg},
 *   year          = {2006},
 *   editor        = {Alexandrov, V.N. and van Albada, G.D. and Sloot, P.M.A. amd Dongarra, J.},
 *   pages         = {525--532},
 *   doi           = {10.1007/11758549_73},
 *   organization  = {Springer}
 * }
 * @endcode
 * @code{.bib}
 * @InBook{Gay2012a,
 *   chapter   = {Using expression graphs in optimization algorithms},
 *   pages     = {247--262},
 *   title     = {Mixed Integer Nonlinear Programming},
 *   publisher = {Springer New York},
 *   year      = {2012},
 *   author    = {Gay, D. M.},
 *   editor    = {Lee, J. and Leyffer, S.},
 *   isbn      = {978-1-4614-1927-3},
 *   doi       = {10.1007/978-1-4614-1927-3_8}
 * }
 * @endcode
 * @code{.bib}
 * @InBook{Phipps2012a,
 *   chapter     = {Efficient Expression Templates for Operator Overloading-based Automatic Differentiation},
 *   pages       = {309--319},
 *   title       = {Recent Advances in Algorithmic Differentiation},
 *   publisher   = {Springer},
 *   year        = {2012},
 *   author      = {Eric Phipps and Roger Pawlowski},
 *   editor      = {Forth S. and Hovland P. and Phipps E. and Utke J. and Walther A.},
 *   series      = {Lecture Notes in Computational Science and Engineering},
 *   volume      = {73},
 *   date        = {2012-05-15},
 *   doi         = {10.1007/978-3-642-30023-3_28},
 *   eprint      = {1205.3506v1},
 *   eprintclass = {cs.MS},
 *   eprinttype  = {arXiv}
 * }
 * @endcode
 *
 * The implementation of both forward- and reverse-mode Sacado numbers is quite intricate.
 * As of Trilinos 12.12, the implementation of math operations involves a lot of preprocessor
 * directives and macro programming. Accordingly, the code may be hard to follow and there
 * exists no meaningful companion documentation for these classes.
 * So, a useful resource for understanding the principle implementation of these numbers
 * can be found at
 * <a href="https://trilinos.org/docs/dev/packages/sacado/doc/html/classSacado_1_1Fad_1_1SimpleFad.html">this link for the Sacado::Fad::SimpleFad class</a>
 * that outlines a reference (although reportedly inefficient) implementation of a
 * forward-mode auto-differentiable number that does not use expression templates.
 * (Although not explicitly stated, it would appear that the Sacado::Fad::SimpleFad class
 * is implemented in the spirit of dual numbers.)
 *
 * @subsection auto_diff_1_2 How automatic differentiation is integrated into deal.II
 *
 * Since the interface to each automatic differentiation library is so vastly different,
 * a uniform internal interface to each number will be established in the near future.
 * The goal will be to allow some driver classes (that provide the core functionality,
 * and will later be introduced in the next section) to have a consistent mechanism to interact
 * with different auto-differentiation libraries. Specifically, they need to be able to correctly
 * initialize and finalize data that is to be interpreted as the dependent and independent
 * variables of a formula.
 *
 * A summary of the files that implement the interface to the supported auto-differentiable
 * numbers is as follows:
 *
 * - ad_drivers.h: Provides classes that act as drivers to the interface of internally supported
 *   automatic differentiation libraries. These are used internally as an intermediary to the
 *   helper classes that we provide.
 * - ad_helpers.h: Provides a set of classes to help perform automatic differentiation in a
 *   number of different contexts. These are detailed in @ref auto_diff_1_3.
 * - ad_number_types.h: Introduces an enumeration (called a type code) for the
 *   auto-differentiable number combinations that will be supported by the driver classes.
 *   The rationale behind the use of this somewhat restrictive mechanism is discussed below.
 * - ad_number_traits.h: Declare some internal classes that are to be specialized for
 *   each auto-differentiation library and/or number type. These are subsequently used to
 *   provide a uniform interface to the classes through the NumberTraits and ADNumberTraits
 *   classes which are extensively used throughout of drivers. We also provide some mechanisms
 *   to easily query select properties of these numbers, i.e. some type traits.
 * - adolc_math.h: Extension of the ADOL-C math operations that allow these numbers to be used
 *   consistently throughout the library.
 * - adolc_number_types.h: Implementation of the internal classes that define how we
 *   use ADOL-C numbers.
 * - adolc_product_types.h: Defines some product and scalar types that allow the use of
 *   ADOL-C numbers in conjunction with the Tensor and SymmetricTensor classes.
 * - sacado_math.h: Extension of the Sacado math operations that allow these numbers to be used
 *   consistently throughout the library.
 * - sacado_number_types.h: Implementation of the internal classes that define how we
 *   use the supported Sacado numbers.
 * - sacado_product_types.h: Defines some product and scalar types that allow the use of
 *   the supported Sacado numbers in conjunction with the Tensor and SymmetricTensor
 *   classes.
 *
 * By using type codes for each supported number type, we artificially limit the type
 * of auto-differentiable numbers that can be used within the library. This design choice
 * is due to the fact that its not trivial to ensure that each number type is correctly
 * initialized and that all combinations of nested (templated) types remain valid for all
 * operations performed by the library.
 * Furthermore, there are some lengthy functions within the library that are instantiated
 * for the supported number types and have internal checks that are only satisfied when a
 * auto-differentiable number, of which the library has knowledge, is used. This again
 * ensures that the integrity of all computations is maintained.
 * Finally, using a simple enumeration as a class template parameter ultimately makes it
 * really easy to switch between the type used in production code with little to no further
 * amendments required to user code.
 *
 * @subsubsection auto_diff_1_3 User interface to the automatic differentiation libraries
 *
 * The deal.II library offers a unified interface to the automatic differentiation libraries that
 * we support. To date, the helper classes have been developed for the following contexts:
 *
 * - Classes designed to operate at the quadrature point level (or any general continuum point):
 *   - ScalarFunction: Differentiation of a scalar-valued function. One typical use would be the
 *                     the development of constitutive laws directly from a strain energy function.
 *   - VectorFunction: Differentiation of a vector-valued function. This could be used to
 *                     linearize the kinematic variables of a constitutive law, or assist in solving
 *                     the evolution equations of local internal variables.
 * - Classes designed to operate at the cell level:
 *   - EnergyFunctional: Differentiation of a scalar-valued energy functional, such as might arise
 *                       from variational formulations.
 *   - ResidualLinearization: Differentiation of a vector-valued finite element residual, leading to
 *                            its consistent linearization.
 *
 * Naturally, it is also possible for users to manage the initialization and derivative
 * computations themselves.
 *
 * The most up-to-date examples of how this is done using ADOL-C can be found in
 * - their <a href="https://projects.coin-or.org/ADOL-C/browser/trunk/ADOL-C/doc/adolc-manual.pdf?format=raw">user manual</a>,
 * - their <a href="https://gitlab.com/adol-c/adol-c/tree/master/ADOL-C/examples">development repository</a>, and
 * - our <a href="https://github.com/dealii/dealii/tree/master/tests/adolc">test-suite</a>,
 *
 * while for Sacado, illustrative examples can be found in
 * - their <a href="https://github.com/trilinos/Trilinos/tree/master/packages/sacado/example">development repository</a>,
 * - a <a href="https://github.com/dealii/code-gallery/tree/master/Quasi_static_Finite_strain_Compressible_Elasticity">code-gallery example</a>, and
 * - our <a href="https://github.com/dealii/dealii/tree/master/tests/sacado">test-suite</a>.
 *
 *
 * @section symb_diff_1 Symbolic expressions and differentiation
 *
 * <a href="https://en.wikipedia.org/wiki/Symbolic_differentiation">Symbolic differentiation</a> is,
 * in terms of its design and usage, quite different to automatic differentiation.
 * Underlying any symbolic library is a computer algebra system (CAS) that implements a
 * language and collection of algorithms to manipulate symbolic (or "string-like") expressions.
 * This is most similar, from a philosophical point of view, to how algebraic operations would be
 * performed by hand.
 *
 * To help better distinguish between symbolic differentiation and numerical methods like automatic
 * differentiation, let's consider a very simple example.
 * Suppose that the function $f(x,y) = [2x+1]^{y}$, where $x$ and $y$ are variables that are independent
 * of one another.
 * By applying the chain-rule, the derivatives of this function are simply
 * $\dfrac{d f(x,y)}{d x} = 2y[2x+1]^{y-1}$ and
 * $\dfrac{d f(x,y)}{d y} = [2x+1]^{y} \ln(2x+1)$.
 * These are exactly the results that you get from a CAS after defining the symbolic variables
 * `x` and `y`, defining the symbolic expression `f = pow(2x+1, y)` and computing the
 * derivatives `diff(f, x)` and `diff(f, y)`.
 * At this point there is no assumption of what `x` and `y` represent; they may later be interpreted
 * as plain (scalar) numbers, complex numbers, or something else for which the power and natural
 * logarithm functions are well defined.
 * Obviously this means that there is also no assumption about which point to evaluate either
 * the expression or its derivatives.
 * One could readily take the expression for $\dfrac{d f(x, y)}{d x}$ and evaluate it
 * at $x=1, y=2.5$ and then later, with no recomputation of the derivative expression itself,
 * evaluate it at $x=3.25, y=-6$.
 * In fact, the interpretation of any symbolic variable or expression, and the inter-dependencies
 * between variables, may be defined or redefined at any point during their manipulation;
 * this leads to a degree of flexibility in computations that cannot be matched by
 * auto-differentiation.
 * For example, one could perform the permanent substitution
 * $g(x) = \dfrac{d f(x, y)}{d x} \vert_{y=1}$ and then recompute
 * $g(x)$ for several different values of $x$.
 * One could also post-factum express an interdependency between `x` and `y`, such as
 * $y \rightarrow y(x) := 2x$.
 * For such a case, this means that the initially computed derivatives
 * $\dfrac{d f(x, y)}{d x} \rightarrow \dfrac{\partial f(x, y(x))}{\partial x} = 2y(x) [2x+1]^{y(x)-1} = 4x[2x+1]^{2x-1}$ and
 * $\dfrac{d f(x, y)}{d y} \rightarrow \dfrac{\partial f(x, y(x))}{\partial y} = [2x+1]^{y(x)} \ln(2x+1) = [2x+1]^{2x} \ln(2x+1)$
 * truly represent partial derivatives rather than total derivatives.
 * Of course, if such an inter-dependency was explicitly defined before the derivatives
 * $\dfrac{d f(x, y(x))}{d x}$ and $\dfrac{d f(x, y(x))}{d y}$ are computed, then this
 * could correspond to the total derivative (which is the only result that auto-differentiation
 * is able to achieve for this example).
 *
 * Due to the sophisticated CAS that forms the foundation of symbolic operations, the types of
 * manipulations are not necessarily restricted to differentiation alone, but rather may span a
 * spectrum of manipulations relevant to discrete differential calculus, topics in pure
 * mathematics, and more.
 * The documentation for the <a href="https://www.sympy.org/en/index.html">SymPy</a> library gives
 * plenty of examples that highlight what a fully-fledged CAS is capable of.
 * Through the Differentiation::SD::Expression class, and the associated functions in the
 * Differentiation::SD namespace, we provide a wrapper to the high-performance
 * <a href="https://github.com/symengine/symengine">SymEngine</a> symbolic manipulation library
 * that has enriched operator overloading and a consistent interface that makes it easy and
 * "natural" to use.
 * In fact, this class can be used as a "drop-in" replacement for arithmetic types in many
 * situations, transforming the operations from being numeric to symbolic in nature; this is
 * made especially easy when classes are templated on the underlying number type.
 * Being focused on numerical simulation of PDEs, the functionality of the CAS that is exposed
 * within deal.II focuses on symbolic expression creation, manipulation, and differentiation.
 *
 * As a final note, it is important to recognize a major deficiency in deal.II's current implementation
 * of the interface to the supported symbolic library.
 * To date, convenience wrappers to SymEngine functionality is focused on manipulations that solely
 * involve dictionary-based (i.e., something reminiscent of "string-based") operations.
 * Although SymEngine performs these operations in an efficient manner, they are still known to be
 * computationally expensive, especially when the operations are performed on large expressions.
 * It should therefore be expected that the performance of the parts of code that perform
 * differentiation, symbolic substitution, etc., @b may be a limiting factor when using this in
 * production code.
 * In the future, deal.II will provide an interface to accelerate the evaluation of lengthy symbolic
 * expression through the @p BatchOptimizer class (which is already referenced in several places in
 * the documentation).
 * In particular, the @p BatchOptimizer will simultaneously optimize a collection of symbolic
 * expressions using methods such as common subexpression elimination (CSE), as well as by generating
 * high performance code-paths to evaluate these expressions through the use of a custom-generated
 * `std::function` or by compiling the expression using the LLVM JIT compiler.
 * Additionally, the level of functionality currently implemented effectively limits the use of
 * symbolic algebra to the traditional use case (i.e. scalar and tensor algebra, as might be useful to
 * define constitutive relations or complex functions for application as boundary conditions or
 * source terms).
 * In the future we will also implement classes to assist in performing assembly operations in
 * the same spirit as that which has been done in the Differentiation::AD namespace.
 *
 * A summary of the files that implement the interface to the supported symbolic differentiable
 * numbers is as follows:
 * - symengine_math.h: Implementation of math operations that allow the class that implements
 *   symbolic expressions to be used consistently throughout the library and in user code.
 *   It provides counterpart definitions for many of the math functions found in the standard
 *   namespace.
 * - symengine_number_traits.h: Provides some mechanisms to easily query select properties of
 *   symbolic numbers, i.e. some type traits.
 * - symengine_number_types.h: Implementation of the Differentiation::SD::Expression class that can
 *   be used to represent scalar symbolic variables, scalar symbolic expressions, and more.
 *   This Expression class has been given a full set of operators overloaded for all mathematical
 *   and logical operations that are supported by the SymEngine library and are considered useful
 *   within the context of numerical modeling.
 * - symengine_product_types.h: Defines some product and scalar types that allow the use of symbolic
 *   expressions in conjunction with the Tensor and SymmetricTensor classes.
 * - symengine_scalar_operations.h: Defines numerous operations that can be performed either on or
 *   with scalar symbolic expressions or variables.
 *   This includes (but is not limited to) the creation of scalar symbols, performing differentiation
 *   with respect to scalars, and symbolic substitution within scalar expressions.
 * - symengine_tensor_operations.h: Defines numerous operations that can be performed either on or
 *   with tensors of symbolic expressions or variables.
 *   This includes (but is not limited to) the creation of tensors of symbols, performing
 *   differentiation with respect to tensors of symbols, differentiation of tensors of symbols, and
 *   symbolic substitution within tensor expressions.
 * - symengine_types.h: Provides aliases for some types that are commonly used within the context of
 *   symbolic computations.
 * - symengine_utilities.h: Provides some utility functions that are useful within the context of
 *   symbolic computations.
 */
