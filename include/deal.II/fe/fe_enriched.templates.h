// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#ifndef dealii_fe_enriched_templates_h
#define dealii_fe_enriched_templates_h

#include <deal.II/fe/fe_enriched.h>

DEAL_II_NAMESPACE_OPEN

namespace ColorEnriched
{
  namespace internal
  {
    /**
     * Returns true if there is a connection between subdomains in the mesh
     * associated with @p hp::DoFHandler i.e., if the subdomains share at least
     * a vertex. The two subdomains are defined by predicates provided by
     * @p predicate_1 and @p predicate_2. A predicate is a function (or
     * object of a type with an operator()) which takes in a cell iterator and
     * gives a boolean. It is said to be active in a cell if it returns true.
     *
     * An example of a custom predicate is one that checks the distance from a
     * fixed point. Note that the operator() takes in a cell iterator. Using the
     * constructor, the fixed point and the distance can be chosen.
     * @code
     * <int dim>
     * struct predicate
     * {
     *     predicate(const Point<dim> p, const int radius)
     *     :p(p),radius(radius){}
     *
     *     template <class Iterator>
     *     bool operator () (const Iterator &i)
     *     {
     *         return ( (i->center() - p).norm() < radius);
     *     }
     *
     * private:
     *     Point<dim> p;
     *     int radius;
     *
     * };
     * @endcode
     * and then the function can be used as follows to find if the subdomains
     * are connected.
     * @code
     * find_connection_between_subdomains
     * (dof_handler,
     *  predicate<dim>(Point<dim>(0,0), 1)
     *  predicate<dim>(Point<dim>(2,2), 1));
     * @endcode
     *
     * @param[in] hp::DoFHandler object
     * @param[in] predicate_1 A function (or object of a type with an
     * operator()) defining the subdomain 1. The function takes in a cell and
     * returns a boolean.
     * @param[in] predicate_2 Same as @p predicate_1 but defines subdomain 2.
     * @return A boolean "true" if the subdomains share atleast a vertex.
     */
    template <int dim, int spacedim>
    bool
    find_connection_between_subdomains(
      const hp::DoFHandler<dim, spacedim> &    dof_handler,
      const predicate_function<dim, spacedim> &predicate_1,
      const predicate_function<dim, spacedim> &predicate_2);

    /**
     * Assign colors to subdomains using Graph coloring algorithm where each
     * subdomain is considered as a graph node. Subdomains which are
     * connected i.e share atleast a vertex have different color. Each subdomain
     * is defined using a predicate function of @p predicates.
     *
     * @param[in] dof_handler a hp::DoFHandler object
     * @param[in] predicates predicates defining the subdomains
     * @param[out] predicate_colors Colors (unsigned int) associated with each
     * subdomain.
     */
    template <int dim, int spacedim>
    unsigned int
    color_predicates(
      const hp::DoFHandler<dim, spacedim> &                 dof_handler,
      const std::vector<predicate_function<dim, spacedim>> &predicates,
      std::vector<unsigned int> &                           predicate_colors);

    /**
     * Used to construct data members @p cellwise_color_predicate_map and
     * @p fe_sets of Helper class. Inputs are hp::DoFHandler object,
     * vector of predicates and colors associated with them. Before calling
     * this function, colors can be assigned to predicates (i.e subdomains)
     * using the function color_predicates.
     *
     * Each active FE index has a set of colors associated with it.
     * A cell with an active FE index i has a set of colors given by
     * <code>fe_sets[i]</code>. An active FE index with color {a,b}
     * means that the cell has two active predicates (i.e they return true
     * for the cell) of color a and b.
     *
     * Eg: fe_sets = { {}, {1}, {2}, {1,2} } means
     * Cells with active FE index 0 have no predicates associated.
     * Cells with index 1 have a active predicate with color 1.
     * Cells with index 2 have a active predicate with color 2.
     * Cells with index 3 have active predicates with color 1 and color 2.
     *
     * A map of maps cellwise_color_predicate_map is used to associate
     * predicate colors in cells with predicate ids. For this purpose, each
     * cell is given a unique id which is stored in material id for now.
     * When the grid is refined, material id is inherited to the children, so
     * map which associates material id with color map will still be relevant.
     *
     * Now the color map can be explained with an example. If the cell with
     * material id 100 has active predicates 4 (color = 1) and 5 (color = 2),
     * the map will insert pairs (1, 4) and (2, 5) at key 100 (i.e unique id
     * of cell is mapped with a map which associates color with predicate id).
     *
     * @param[in] dof_handler hp::DoFHandler object
     * @param[in] predicates vector of predicates defining the subdomains.
     * <code>@p predicates[i]</code> returns true for a cell if it
     * belongs to subdomain with index i.
     * @param[in] predicate_colors vector of colors (unsigned int) associated
     * with each subdomain.
     * @param[out] cellwise_color_predicate_map A map of maps used to associate
     * predicate colors in cells with predicate ids.
     * @param[out] fe_sets a vector of color lists
     */
    template <int dim, int spacedim>
    void
    set_cellwise_color_set_and_fe_index(
      hp::DoFHandler<dim, spacedim> &                       dof_handler,
      const std::vector<predicate_function<dim, spacedim>> &predicates,
      const std::vector<unsigned int> &                     predicate_colors,
      std::map<unsigned int, std::map<unsigned int, unsigned int>>
        &                                  cellwise_color_predicate_map,
      std::vector<std::set<unsigned int>> &fe_sets);

    /**
     * A function that returns a vector of enrichment functions corresponding
     * to a color. The size of the vector is equal to total number of different
     * colors associated with predicates (i.e subdomains).
     *
     * Assume that a cell has a active predicates with ids 4 (color = 1) and
     * 5 (color = 2). cellwise_color_predicate_map has this information
     * provided we know the material id.
     *
     * Now a call to color_enrichment[1](cell) should in turn call
     * enrichments[4](cell).
     *
     * @param[in] num_colors number of colors for predicates
     * @param[in] enrichments vector of enrichment functions
     * @param[in] cellwise_color_predicate_map A map of maps used to associate
     * predicate colors in cells with predicate ids.
     * @param[out] color_enrichments A vector of functions that take in cell
     * and return a function pointer.
     */
    template <int dim, int spacedim>
    void
    make_colorwise_enrichment_functions(
      const unsigned int &                                    num_colors,
      const std::vector<std::shared_ptr<Function<spacedim>>> &enrichments,
      const std::map<unsigned int, std::map<unsigned int, unsigned int>>
        &cellwise_color_predicate_map,
      std::vector<std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>>
        &color_enrichments);


    /**
     * Creates a hp::FECollection object constructed using FE_Enriched
     * elements which itself is constructed using color enrichment functions
     * and is of size equal to number of colors.
     *
     * @param[in] fe_sets a vector of color lists
     * @param[in] color_enrichments A vector of functions that take in cell
     * and return a function pointer.
     * @param[in] fe_base base FiniteElement
     * @param[in] fe_enriched enriched FiniteElements
     * @param[in] fe_nothing a finite element with zero degrees of freedom
     * @param[out] fe_collection a collection of
     * finite elements
     */
    template <int dim, int spacedim>
    void
    make_fe_collection_from_colored_enrichments(
      const unsigned int &num_colors,
      const std::vector<std::set<unsigned int>>
        &fe_sets, // total list of color sets possible
      const std::vector<std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>>
        &color_enrichments,               // color wise enrichment functions
      const FE_Q<dim, spacedim> &fe_base, // basic FE element
      const FE_Q<dim, spacedim>
        &fe_enriched, // FE element multiplied by enrichment function
      const FE_Nothing<dim, spacedim> &fe_nothing,
      hp::FECollection<dim, spacedim> &fe_collection);
  } // namespace internal
} // namespace ColorEnriched

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_fe_enriched_templates_h
