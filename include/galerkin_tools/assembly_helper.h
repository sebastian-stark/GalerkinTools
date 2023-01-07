// --------------------------------------------------------------------------
// Copyright (C) 2020 by Sebastian Stark
//
// This file is part of the GalerkinTools library
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef GALERKINTOOLS_ASSEMBLYHELPER_H_
#define GALERKINTOOLS_ASSEMBLYHELPER_H_

#include <vector>
#include <utility>
#include <tuple>
#include <map>
#include <string>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/base/conditional_ostream.h>

#include <galerkin_tools/config.h>
#include <galerkin_tools/triangulation_system.h>
#include <galerkin_tools/dirichlet_constraint.h>
#include <galerkin_tools/dof_handler_system.h>
#include <galerkin_tools/fe_values_interface.h>
#include <galerkin_tools/total_potential.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * Special Function class allowing to pass the cell information when asking for the initial value of an unknown field at a finite element support point.
 * Objects of this type can be given to the constructor of an IndependentField; and when the initial solution vector is set up with
 * a call to AssemblyHelper::get_initial_fields_vector(), FunctionCell::cell will be appropriately set before each call to the Function::value(),
 * so that one can find out through FunctionCell::get_cell() on which cell the support point lies, for
 * which Function::value() currently has to provide with initial values.
 *
 * @tparam	dim			dimensionality of the object for which FunctionCell has to provide with initial values
 *
 * @tparam	spacedim	spatial dimension
 *
 */
template<unsigned int dim, unsigned int spacedim>
class FunctionCell : public dealii::Function<spacedim>
{
private:

	/**
	 * the current cell
	 */
	const TriaIterator<CellAccessor<dim, spacedim>>*
	cell = nullptr;

public:

	/**
	 * constructor
	 *
	 * @param[in]	n_components	Number of components of function
	 */
	FunctionCell(const unsigned int n_components = 1);

	/**
	 * Destructor, avoid that this class can be instantiated
	 */
	virtual
	~FunctionCell() = 0;

	/**
	 * @param[in]	cell	FunctionCell::cell
	 */
	void
	set_cell(const TriaIterator<CellAccessor<dim, spacedim>>& cell);

	/**
	 * @return	FunctionCell::cell
	 */
	const TriaIterator<CellAccessor<dim, spacedim>>&
	get_cell()
	const;
};

/**
 * The AssemblyHelper class puts it all together and provides with methods for assembly of the finite element system.
 *
 * Essentially, an AssemblyHelper object is a combination of a TriangulationSystem object (which includes the definition of the domain
 * portions and interface portions), a TotalPotential object, a corresponding DoFHandlerSystem object, and mapping
 * objects defining the mapping to be used on the domain and on the interface.
 *
 * Large part of the class rearranges the data involved in the problem definition in a way which is suitable for
 * efficient assembly.
 *
 * @todo	It might be worthwhile to allow for the usage of different mappings on different domain and interface portions.
 *
 * @todo	Routines for treatment of DG terms should be implemented.
 *
 * @tparam	spacedim	spatial dimension
 */
template<unsigned int spacedim>
class AssemblyHelper
{

private:

/** @name Problem definition
 *
 *  These members define the problem to be solved.
 */
///@{

	/**
	 * The total potential \f$\Pi\f$
	 */
	const TotalPotential<spacedim>
	total_potential;

	/**
	 * The triangulation system
	 */
	TriangulationSystem<spacedim>&
	tria_system;

	/**
	 * Mapping to be used for domain cells (should be consistent with AssemblyHelper::mapping_interface)
	 */
	const SmartPointer<const Mapping<spacedim, spacedim>>
	mapping_domain;

	/**
	 * Mapping to be used for interface cells (should be consistent with AssemblyHelper::mapping_domain)
	 */
	const SmartPointer<const Mapping<spacedim-1, spacedim>>
	mapping_interface;

///@}


/** @name Domain and interface portion indexing
 *
 *  These members define consecutive indices for domain and interface portions for internal use.
 */
///@{

	/**
	 * Map between the types::material_id%s associated with domain cells and a unique "internal material index",
	 * which is assigned upon construction of an AssemblyHelper object.
	 *
	 * Whereas the types::material_id%s need not be consecutive numbers, the internal material index is consecutive
	 * and starts from zero. This allows to reduce map lookups to a minimum. In particular, if a domain cell
	 * is visited during fe system assembly, the internal material index corresponding to the types::material_id is looked up and no further
	 * lookups are required afterwards.
	 */
	std::map<types::material_id, const unsigned int>
	material_id_to_internal_index_domain;

	/**
	 * Map between interface "subportion" (identified by types::material_id of interface,
	 * types::material_id of - side and types::material_id of + side) and a unique
	 *  "internal material index", which is assigned upon construction of an AssemblyHelper object.
	 *
	 * Whereas the types::material_id%s need not be consecutive numbers, the internal material index is consecutive
	 * and starts from zero. This allows to reduce map lookups to a minimum. In particular, if an interface cell
	 * is visited during fe system assembly, the internal material index corresponding to the types::material_id%s of
	 * the interface cell and the domain cells on both sides of
	 * the interface is looked up and no further lookups are required afterwards.
	 *
	 * If an interface subportion is actually a boundary, the third types::material_id in the tuple is
	 * numbers::invalid_material_id.
	 *
	 * It is advantageous to treat every combination of (interface types::material_id, domain types::material_id on - side,
	 * domain types::material_id on + side) separately  because the finite elements involved in the computation of
	 * interface related scalar functionals may be different on different subportions.
	 */
	std::map<
				std::tuple<	const types::material_id,
							const types::material_id,
							const types::material_id>,
							const unsigned int
			>
	material_ids_to_internal_index_interface;

///@}

// finite elements and dof handling

/** @name Finite elements and dof handling
 *
 *  These members comprise the finite elements to be used on the respective domain and interface portions as
 *  well as the DoFHandlerSystem (which contains all dof information).
 */
///@{

	/**
	 * hp::FECollection to be used on domain. This contains the finite element systems to be used for different domain portions
	 * identified by different types::material_id%s.
	 *
	 * The components of the finite element systems are sorted according to AssemblyHelper::global_component_indices_u_omega.
	 */
	hp::FECollection<spacedim, spacedim>
	fe_collection_domain;

	/**
	 * hp::FECollection for interface. This contains the finite element systems to be used for different interface portions
	 * identified by different types::material_id%s.
	 *
	 * The components of the finite element systems are sorted according to AssemblyHelper::global_component_indices_u_sigma.
	 */
	hp::FECollection<spacedim-1, spacedim>
	fe_collection_interface;

	/**
	 * Map between types::material_id%s on domain and corresponding fe index in AssemblyHelper::fe_collection_domain
	 */
	std::map<types::material_id, unsigned int>
	material_id_to_fe_system_id_domain;

	/**
	 * Map between types::material_id%s on interface and corresponding fe index in AssemblyHelper::fe_collection_interface
	 */
	std::map<types::material_id, unsigned int>
	material_id_to_fe_system_id_interface;

	/**
	 * The DoFHandlerSystem
	 */
	DoFHandlerSystem<spacedim>
	dof_handler_system;

///@}


/** @name Members re-organizing AssemblyHelper::total_potential
 *
 *  These members re-organize the information provided by AssemblyHelper::total_potential in a way
 *  suitable for efficient assembly of the finite element system.
 */
///@{

// data structures which re-organize the definition of the problem given by AssemblyHelper::total_potential
// in a way suitable for assembly of the finite element system

	/**
	 * vector of domain related independent fields \f$u^\Omega_\epsilon\f$ involved in the definition of the total potential
	 * AssemblyHelper::total_potential.
	 * These domain related independent fields are sorted by IndependentField::name
	 */
	std::vector<SmartPointer<const IndependentField<spacedim, spacedim>>>
	u_omega;

	/**
	 * vector of interface related independent fields \f$u^\Sigma_\eta\f$ involved in the definition of the total potential
	 * AssemblyHelper::total_potential.
	 * These interface related independent fields are sorted by IndependentField::name
	 */
	std::vector<SmartPointer<const IndependentField<spacedim-1, spacedim>>>
	u_sigma;

	/**
	 * vector of independent scalars \f$C_\iota\f$ involved in the definition of the problem.
	 * These independent scalars are sorted by IndependentField<0, spacedim>::name
	 */
	std::vector<SmartPointer<const IndependentField<0, spacedim>>>
	C;

	/**
	 * Map between domain related IndependentField objects and global component indices of the respective first components
	 * ("global component indices" refers to the indexing of components in AssemblyHelper::fe_collection_domain).
	 * The IndependentFields are distinguished by memory address. I.e., two different IndependentField objects with the
	 * same contents are not considered the same.
	 */
	std::map<SmartPointer<const IndependentField<spacedim, spacedim>>, const unsigned int>
	global_component_indices_u_omega;

	/**
	 * Map between interface related IndependentField objects and global component indices of the respective first components
	 * ("global component indices" refers to the indexing of components in AssemblyHelper::fe_collection_interface).
	 * The IndependentFields are distinguished by memory address. I.e., two different IndependentField objects with the
	 * same contents are not considered the same.
	 */
	std::map<SmartPointer<const IndependentField<spacedim-1, spacedim>>, const unsigned int>
	global_component_indices_u_sigma;

	/**
	 * Map between independent scalars and "global indices" for the independent scalars. This is essentially
	 * the inverse of the member AssemblyHelper::C (in that the global index corresponding to a certain independent scalar
	 * is the index into AssemblyHelper::C containing this independent scalar).
	 */
	std::map<SmartPointer<const IndependentField<0, spacedim>>, const unsigned int>
	global_indices_C;

	/**
	 * AssemblyHelper::scalar_functionals_domain[@p i] contains the domain related scalar functionals which have non-zero integrands
	 * \f$h^\Omega_\rho\f$ on the domain portion with internal index @p i. If the contribution of a certain domain cell associated
	 * with internal index @p i to the finite element system is to be assembled, it is looped over all scalar functionals contained in
	 * AssemblyHelper::scalar_functionals_domain[@p i].
	 */
	std::vector<std::vector<SmartPointer<const ScalarFunctional<spacedim, spacedim>>>>
	scalar_functionals_domain;

	/**
	 * AssemblyHelper::scalar_functionals_interface[@p i] contains the interface related scalar functionals which have non-zero integrands
	 * \f$h^\Sigma_\tau\f$ on the interface portion with internal index @p i. If the contribution of a certain interface cell associated
	 * with internal index @p i to the finite element system is to be assembled, it is looped over all scalar functionals contained in
	 * AssemblyHelper::scalar_functionals_interface[@p i].
	 */	std::vector<std::vector<SmartPointer<const ScalarFunctional<spacedim-1, spacedim>>>>
	scalar_functionals_interface;

	/**
	 * AssemblyHelper::scalar_functionals_domain[@p i][ AssemblyHelper::scalar_functionals_domain_nonprimitive[@p i][@p k] ] is the @p k-th domain related
	 * scalar functional entering non-primitively into the total potential (i.e., the scalar functional is part of at
	 * least one TotalPotentialContribution with the member TotalPotentialContribution::is_primitive set to @p false)
	 */
	std::vector<std::vector<unsigned int>>
	scalar_functionals_domain_nonprimitive;

	/**
	 * AssemblyHelper::scalar_functionals_interface[@p i][ AssemblyHelper::scalar_functionals_interface_nonprimitive[@p i][@p k] ] the @p k-th interface related
	 * scalar functional entering non-primitively into the total potential (i.e., the scalar functional is part of at
	 * least one TotalPotentialContribution with the member TotalPotentialContribution::is_primitive set to @p false)
	 */
	std::vector<std::vector<unsigned int>>
	scalar_functionals_interface_nonprimitive;

	/**
	 * Map assigning indices to all domain related scalar functionals entering into at least one non-primitive TotalPotentialContribution.
	 */
	std::map<SmartPointer<const ScalarFunctional<spacedim, spacedim>>, const unsigned int>
	scalar_functionals_domain_nonprimitive_indices;

	/**
	 * Map assigning indices to all interface related scalar functionals entering into at least one non-primitive TotalPotentialContribution.
	 * Note that the indices here do not start from 0, but rather from
	 * AssemblyHelper::scalar_functionals_domain_nonprimitive_indices .size().
	 */
	std::map<SmartPointer<const ScalarFunctional<spacedim-1, spacedim>>, const unsigned int>
	scalar_functionals_interface_nonprimitive_indices;

	/**
	 * Map assigning indices to all domain related scalar functionals entering into at least one primitive TotalPotentialContribution.
	 */
	std::map<SmartPointer<const ScalarFunctional<spacedim, spacedim>>, const unsigned int>
	scalar_functionals_domain_primitive_indices;

	/**
	 * Map assigning indices to all interface related scalar functionals entering into at least one primitive TotalPotentialContribution.
	 * Note that the indices here do not start from 0, but rather from
	 * AssemblyHelper::scalar_functionals_domain_primitive_indices .size().
	 */
	std::map<SmartPointer<const ScalarFunctional<spacedim-1, spacedim>>, const unsigned int>
	scalar_functionals_interface_primitive_indices;

	/**
	 * Map between domain related scalar functionals and a vector of pairs of
	 * ("contribution to total potential", "scalar functional index within contribution to total potential"),
	 * where "contribution to total potential" is an index into TotalPotential::total_potential_contributions
	 * and "scalar functional index within contribution to total potential" is an index into
	 * TotalPotentialContribution::H_omega. This member essentially describes how a domain related scalar
	 * functional factors into the TotalPotential
	 */
	std::map<SmartPointer<const ScalarFunctional<spacedim, spacedim>>, std::vector<std::pair<const unsigned int, const unsigned int>>>
	contributions_scalar_functionals_domain_total_potential;

	/**
	 * Map between interface related scalar functionals and a vector of pairs of
	 * ("contribution to total potential", "scalar functional index within contribution to total potential"),
	 * where "contribution to total potential" is an index into TotalPotential::total_potential_contributions
	 * and "scalar functional index within contribution to total potential" is an index into
	 * TotalPotentialContribution::H_sigma (but shifted by the size of TotalPotentialContribution::H_omega).
	 * This member essentially describes how an interface related scalar
	 * functional factors into the TotalPotential
	 */
	std::map<SmartPointer<const ScalarFunctional<spacedim-1, spacedim>>, std::vector<std::pair<const unsigned int, const unsigned int>>>
	contributions_scalar_functionals_interface_total_potential;

	/**
	 * Number of scalar functionals entering non-primitively into the total potential. This is the
	 * sum of the sizes of the maps AssemblyHelper::scalar_functionals_domain_nonprimitive_indices
	 * and AssemblyHelper::scalar_functionals_interface_nonprimitive_indices
	 */
	unsigned int
	n_scalar_functionals_nonprimitive;

	/**
	 * Number of scalar functionals entering primitively into the total potential. This is the
	 * the sum of the sizes of the maps AssemblyHelper::scalar_functionals_domain_primitive_indices
	 * and AssemblyHelper::scalar_functionals_interface_primitive_indices
	 */
	unsigned int
	n_scalar_functionals_primitive;

///@}


/** @name Members providing information about shape function indexing, etc.
 *
 *  These data structures are organized in a way that shape function related
 *  information needed again and again for the functionalities provided by
 *  AssemblyHelper (in particular assembly of the finite element system)
 *  can be retrieved efficiently.
 */
///@{

	/**
	 *  %Mapping between global components of domain related FESystem's (which are stored in
	 *  AssemblyHelper::fe_collection_domain) and local shape function indices.
	 *  AssemblyHelper::components_to_shapefuns_domain[@p i][@p j] contains the local shape
	 *  function indices of shape functions contributing to the global component @p j of the @p i-th
	 *  FESystem in AssemblyHelper::fe_collection_domain.
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	components_to_shapefuns_domain;

	/**
	 *  %Mapping between global components of interface related FESystem's (which are stored in
	 *  AssemblyHelper::fe_collection_interface) and local shape function indices.
	 *  AssemblyHelper::components_to_shapefuns_interface[@p i][@p j] contains the local shape
	 *  function indices of shape functions contributing to the global component @p j of the @p i-th
	 *  FESystem in AssemblyHelper::fe_collection_interface.
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	components_to_shapefuns_interface;

	/**
	 *  This member is similar to AssemblyHelper::components_to_shapefuns_domain. However, the information is
	 *  additionally broken down in cell faces. In particular,
	 *  AssemblyHelper::components_to_shapefuns_domain_facewise[@p i][@p j][@p k] contains the indices
	 *  of shape functions related to the global component @p j of the @p i-th FESystem in AssemblyHelper::fe_collection_domain
	 *  which have support on face @p k.
	 *
	 *  It is noted that not all shape functions included in AssemblyHelper::components_to_shapefuns_domain[@p i][@p j]
	 *  may be present in AssemblyHelper::components_to_shapefuns_domain_facewise[@p i][@p j] because there may
	 *  be shape functions which don't have support on a face or which are not associated with support points at all.
	 *
	 *  The main use of this data structure is to facilitate the assembly of Dirichlet type constraints imposed on domain related
	 *  independent fields on interfaces (such constraints can only be applied directly if the shape functions related to
	 *  the constrained independent field are associated with support points).
	 */
	std::vector<std::vector<std::vector<std::vector<unsigned int>>>>
	components_to_shapefuns_domain_facewise;

	/**
	 * This member contains information about the shape functions coupling for a certain domain related scalar functional.
	 * Essentially, this data structure transforms the cell related indexing of shape functions of deal.II into a scalar
	 * functional related shape function indexing. This is important because the assembly of the contributions of each
	 * cell to the finite element system happens scalar functional wise.
	 *
	 * AssemblyHelper::coupled_dof_indices_scalar_functionals_domain[@p u][@p v] contains all (cell related) dof indices
	 * which couple for the @p v-th scalar functional on the domain portion with internal index @p u (the ordering of
	 * the scalar functionals on the domain portion with internal index @p u is given by
	 * AssemblyHelper::scalar_functionals_domain[@p u]).
	 *
	 * Using this data structure it is e.g. possible to make sure that only those entries are included in the sparsity pattern
	 * which are really needed. As an example, consider the case that there are two domain related scalar functionals
	 * \f$H^\Omega_1[e^\Omega_1(u^\Omega_1, u^\Omega_2)]\f$ and \f$H^\Omega_2[e^\Omega_2(u^\Omega_2, u^\Omega_3)]\f$.
	 * For this case, the shape functions related to \f$u^\Omega_1\f$ will couple to those related to \f$u^\Omega_2\f$;
	 * and the shape functions related to \f$u^\Omega_2\f$ will couple to those related to \f$u^\Omega_3\f$; however there is no coupling
	 * between \f$u^\Omega_1\f$ and \f$u^\Omega_3\f$. The standard approach of assuming that all shape functions living on
	 * a certain cell couple would miss this subtlety and allocate entries in the finite element system matrix which
	 * are always zero. However, by using the scalar functional related shape function indexing provided by the member
	 * AssemblyHelper::coupled_dof_indices_scalar_functionals_domain this is avoided.
	 *
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	coupled_dof_indices_scalar_functionals_domain;

	/**
	 * Essentially the same as AssemblyHelper::coupled_dof_indices_scalar_functionals_domain. However, only dofs related to local independent field are included
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	coupled_dof_indices_scalar_functionals_domain_local;

	/**
	 * Essentially the same as AssemblyHelper::coupled_dof_indices_scalar_functionals_domain. However, only dofs related to nonlocal independent field are included
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	coupled_dof_indices_scalar_functionals_domain_nonlocal;

	/**
	 * This member contains information about the shape functions of an interface cell coupling for a certain
	 * interface related scalar functional. Essentially, this data structure transforms the cell related indexing of shape functions of deal.II into a scalar
	 * functional related shape function indexing. This is important because the assembly of the contributions of each
	 * cell to the finite element system happens scalar functional wise.
	 *
	 * AssemblyHelper::coupled_dof_indices_scalar_functionals_interface[@p u][@p v] contains all (cell related) dof indices
	 * of an interface cell which couple for the @p v-th scalar functional on the interface portion with internal index @p u (the ordering of
	 * the scalar functionals on the interface portion with internal index @p u is given by
	 * AssemblyHelper::scalar_functionals_interface[@p u]).
	 *
	 * See also AssemblyHelper::coupled_dof_indices_scalar_functionals_domain for further explanations regarding the relevance of
	 * this data structure.
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	coupled_dof_indices_scalar_functionals_interface;

	/**
	 * This member contains information about the shape functions of a domain cell on the minus side of the interface
	 * coupling for a certain interface related scalar functional. Essentially, this data structure transforms
	 * the cell related indexing of shape functions of deal.II into a scalar
	 * functional related shape function indexing. This is important because the assembly of the contributions of each
	 * cell to the finite element system happens scalar functional wise.
	 *
	 * AssemblyHelper::coupled_dof_indices_scalar_functionals_interface_minus[@p u][@p v] contains all (cell related) dof indices
	 * of a domain cell on the minus side of the interface which couple for the @p v-th scalar functional on the interface portion
	 * with internal index @p u (the ordering of the scalar functionals on the interface portion with internal index @p u is given by
	 * AssemblyHelper::scalar_functionals_interface[@p u]).
	 *
	 * See also AssemblyHelper::coupled_dof_indices_scalar_functionals_domain for further explanations regarding the relevance of
	 * this data structure.
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	coupled_dof_indices_scalar_functionals_interface_minus;

	/**
	 * This member contains information about the shape functions of a domain cell on the plus side of the interface
	 * coupling for a certain interface related scalar functional. Essentially, this data structure transforms
	 *  the cell related indexing of shape functions of deal.II into a scalar
	 * functional related shape function indexing. This is important because the assembly of the contributions of each
	 * cell to the finite element system happens scalar functional wise.
	 *
	 * AssemblyHelper::coupled_dof_indices_scalar_functionals_interface_plus[@p u][@p v] contains all (cell related) dof indices
	 * of a domain cell on the plus side of the interface which couple for the @p v-th scalar functional on the interface portion
	 * with internal index @p u (the ordering of the scalar functionals on the interface portion with internal index @p u is given by
	 * AssemblyHelper::scalar_functionals_interface[@p u]).
	 *
	 * See also AssemblyHelper::coupled_dof_indices_scalar_functionals_domain for further explanations regarding the relevance of
	 * this data structure.
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	coupled_dof_indices_scalar_functionals_interface_plus;

	/**
	 * This member contains information about the independent scalars coupling for a certain
	 * domain related scalar functional (see AssemblyHelper::C and AssemblyHelper::global_indices_C
	 * for the global indexing of the independent scalars). Essentially, this data structure transforms the global indexing of
	 * independent scalars into a scalar functional related independent scalar indexing.
	 * This is important because the assembly of the contributions of each cell to the finite element system happens scalar functional wise.
	 *
	 * AssemblyHelper::coupled_C_indices_scalar_functionals_domain[@p u][@p v] contains all global independent scalar indices
	 * coupling for the @p v-th scalar functional on the domain portion with internal index @p u (the ordering of
	 * the scalar functionals on the domain portion with internal index @p u is given by
	 * AssemblyHelper::scalar_functionals_domain[@p u]).
	 *
	 * See also AssemblyHelper::coupled_dof_indices_scalar_functionals_domain for further explanations regarding the relevance of
	 * this data structure.
	 *
	 * @todo	This data structure contains redundancy because the way how the independent scalars couple for
	 * 			a certain scalar functional does not depend on the domain portion. This redundancy has initially been introduced
	 * 			to allow for a consistent treatment of cell related dofs and independent scalar related dofs. However, its should be
	 * 			eliminated in later releases.
	 *
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	coupled_C_indices_scalar_functionals_domain;

	/**
	 * This member contains information about the independent scalars coupling for a certain
	 * interface related scalar functional (see AssemblyHelper::C and AssemblyHelper::global_indices_C
	 * for the indexing of the independent scalars). Essentially, this data structure transforms the global indexing of
	 * independent scalars into a scalar functional related independent scalar indexing.
	 * This is important because the assembly of the contributions of each cell to the finite element system happens scalar functional wise.
	 *
	 * AssemblyHelper::coupled_C_indices_scalar_functionals_interface[@p u][@p v] contains all global independent scalar indices
	 * coupling for the @p v-th scalar functional on the interface portion with internal index @p u (the ordering of
	 * the scalar functionals on the interface portion with internal index @p u is given by
	 * AssemblyHelper::scalar_functionals_interface[@p u]).
	 *
	 * See also AssemblyHelper::coupled_dof_indices_scalar_functionals_domain for further explanations regarding the relevance of
	 * this data structure.
	 *
	 * @todo	This data structure contains redundancy because the way how the independent scalars couple for
	 * 			a certain scalar functional does not depend on the interface portion. This redundancy has initially been introduced
	 * 			to allow for a consistent treatment of cell related dofs and independent scalar related dofs. However, its should be
	 * 			eliminated in later releases.
	 */
	std::vector<std::vector<std::vector<unsigned int>>>
	coupled_C_indices_scalar_functionals_interface;


	/**
	 * This member contains information about how domain related dependent fields \f$e^\Omega_\lambda\f$ are
	 * related to shape functions. In particular, it relates shape functions to the terms \f$a^\Omega_{\lambda\epsilon} u^\Omega_\epsilon\f$
	 * (see DependentField<spacedim, spacedim> for further information).
	 *
	 * AssemblyHelper::a_omega[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the domain portion with internal index @p u.
	 *
	 * Each AssemblyHelper::a_omega[@p u][@p v][@p k][@p l] corresponds to a single term \f$a^\Omega_{\lambda \epsilon} u^\Omega_\epsilon\f$,
	 * where	@p get<@p 0>(AssemblyHelper::a_omega[@p u][@p v][@p k][@p l]) is the coefficient \f$a^\Omega_{\lambda \epsilon}\f$,
	 * 			@p get<@p 1>(AssemblyHelper::a_omega[@p u][@p v][@p k][@p l]) is the global component corresponding to the independent field \f$u^\Omega_\epsilon\f$,
	 * and 		@p get<@p 2>(AssemblyHelper::a_omega[@p u][@p v][@p k][@p l]) contains all the (scalar functional related) shape function indices of a domain cell
	 * corresponding to shape functions which are non-zero in the relevant global component.
	 */
	std::vector<std::vector<std::vector<
							std::vector<  std::tuple<const double, const unsigned int, std::vector<unsigned int>> >
							>>>
	a_omega;

	/**
	 * This member contains information about how domain related dependent fields \f$e^\Omega_\lambda\f$ are
	 * related to shape functions. In particular, it relates shape functions to the terms \f$b^\Omega_{\lambda \epsilon i} \dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\f$
	 * (see DependentField<spacedim, spacedim> for further information).
	 *
	 * AssemblyHelper::b_omega[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the domain portion with internal index @p u.
	 *
	 * Each AssemblyHelper::b_omega[@p u][@p v][@p k][@p l] corresponds to a single term \f$b^\Omega_{\lambda \epsilon i} \dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\f$,
	 * where	@p get<@p 0>(AssemblyHelper::b_omega[@p u][@p v][@p k][@p l]) is the coefficient \f$b^\Omega_{\lambda \epsilon i}\f$,
	 * 			@p get<@p 1>(AssemblyHelper::b_omega[@p u][@p v][@p k][@p l]) is the global component corresponding to the independent field \f$u^\Omega_\epsilon\f$,
	 * 			@p get<@p 2>(AssemblyHelper::b_omega[@p u][@p v][@p k][@p l]) is the derivative \f$i\f$ of the independent field \f$u^\Omega_\epsilon\f$,
	 * and 		@p get<@p 3>(AssemblyHelper::b_omega[@p u][@p v][@p k][@p l]) contains all the (scalar functional related) shape function indices of a domain cell
	 * corresponding to shape functions which are non-zero in the relevant global component.
	 */
	std::vector<std::vector<std::vector<
							std::vector<  std::tuple<const double, const unsigned int, const unsigned int, std::vector<unsigned int>> >
							>>>
	b_omega;

	/**
	 * This member contains information about how domain related dependent fields \f$e^\Omega_\lambda\f$ are
	 * related to independent scalars. In particular, it relates scalar functional related independent scalar
	 * indices to terms \f$c^\Omega_{\lambda\iota} C_\iota\f$
	 * (see DependentField<spacedim, spacedim> for further information).
	 *
	 * AssemblyHelper::c_omega[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the domain portion with internal index @p u.
	 *
	 * Each AssemblyHelper::c_omega[@p u][@p v][@p k][@p l] corresponds to a single term \f$c^\Omega_{\lambda\iota} C_\iota\f$,
	 * where	@p get<@p 0>(AssemblyHelper::c_omega[@p u][@p v][@p k][@p l]) is the coefficient \f$c^\Omega_{\lambda\iota}\f$,
	 * and		@p get<@p 1>(AssemblyHelper::c_omega[@p u][@p v][@p k][@p l]) is the scalar functional related independent scalar index corresponding to the independent scalar \f$C_\iota\f$
	 */
	std::vector<std::vector<std::vector<
							std::vector<std::tuple<const double, unsigned int>>
							>>>
	c_omega;

	/**
	 * This member contains information about the constant terms \f$d^\Omega_\lambda\f$ of domain related dependent fields \f$e^\Omega_\lambda\f$.
	 * (see DependentField<spacedim, spacedim> for further information).
	 *
	 * AssemblyHelper::d_omega[@p u][@p v][@p k] contains \f$d^\Omega_\lambda\f$ of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the domain portion with internal index @p u.
	 * For computation of a dependent variable \f$d^\Omega_{\lambda}\f$ is added to the other contributions to the dependent variable.
	 */
	std::vector<std::vector<std::vector<double>>>
	d_omega;

	/**
	 * This member contains information about whether the dependent field is local.
	 *
	 * AssemblyHelper::e_omega_local[@p u][@p v][@p k] contains whether the @p k-th dependent field of the
	 * @p v-th scalar functional on the domain portion with internal index @p u is local.
	 */
	std::vector<std::vector<std::vector<bool>>>
	e_omega_local;

	/**
	 * This member contains information about how interface related dependent fields \f$e^\Sigma_\nu\f$ are
	 * related to shape functions. In particular, it relates shape functions to the terms \f$a^\Sigma_{\nu\eta} u^\Sigma_\eta\f$
	 * (see DependentField for further information).
	 *
	 * AssemblyHelper::a_sigma[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the interface (sub)portion with internal index @p u.
	 *
	 * Each AssemblyHelper::a_sigma[@p u][@p v][@p k][@p l] corresponds to a single term \f$a^\Sigma_{\nu\eta} u^\Sigma_\eta\f$,
	 * where	@p get<@p 0>(AssemblyHelper::a_sigma[@p u][@p v][@p k][@p l]) is the coefficient \f$a^\Sigma_{\nu \eta}\f$,
	 * 			@p get<@p 1>(AssemblyHelper::a_sigma[@p u][@p v][@p k][@p l]) is the global component corresponding to the independent field \f$u^\Sigma_\eta\f$,
	 * and 		@p get<@p 2>(AssemblyHelper::a_sigma[@p u][@p v][@p k][@p l]) contains all the (scalar functional related) shape function indices of an interface cell
	 * corresponding to shape functions which are non-zero in the relevant global component.
	 */
	std::vector<std::vector<std::vector<
							std::vector<  std::tuple<const double, const unsigned int, std::vector<unsigned int>> >
							>>>
	a_sigma;

	/**
	 * This member contains information about how interface related dependent fields \f$e^\Sigma_\nu\f$ are
	 * related to shape functions. In particular, it relates shape functions to the terms \f$b^\Sigma_{\nu \eta i} \dfrac{\partial u^\Sigma_\eta}{\partial x_i}\f$
	 * (see DependentField for further information).
	 *
	 * AssemblyHelper::b_sigma[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the interface (sub)portion with internal index @p u.
	 *
	 * Each AssemblyHelper::b_sigma[@p u][@p v][@p k][@p l] corresponds to a single term \f$b^\Sigma_{\nu \eta i} \dfrac{\partial u^\Sigma_\eta}{\partial x_i}\f$,
	 * where	@p get<@p 0>(AssemblyHelper::b_sigma[@p u][@p v][@p k][@p l]) is the coefficient \f$b^\Sigma_{\nu \eta i}\f$,
	 * 			@p get<@p 1>(AssemblyHelper::b_sigma[@p u][@p v][@p k][@p l]) is the global component corresponding to the independent field \f$u^\Sigma_\eta\f$,
	 * 			@p get<@p 2>(AssemblyHelper::b_sigma[@p u][@p v][@p k][@p l]) is the derivative \f$i\f$ of the independent field \f$u^\Sigma_\eta\f$,
	 * and		@p get<@p 3>(AssemblyHelper::b_sigma[@p u][@p v][@p k][@p l]) contains all the (scalar functional related) shape function indices of an interface cell
	 * corresponding to shape functions which are non-zero in the relevant global component.
	 */
	std::vector<std::vector<std::vector<
							std::vector<  std::tuple<const double, const unsigned int, const unsigned int, std::vector<unsigned int>> >
							>>>
	b_sigma;

	/**
	 * This member contains information about how interface related dependent fields \f$e^\Sigma_\nu\f$ are
	 * related to shape functions. In particular, it relates shape functions to the terms \f$a^-_{\nu\epsilon} (u^\Omega_\epsilon)^-\f$
	 * (see DependentField for further information).
	 *
	 * AssemblyHelper::a_minus[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the interface (sub)portion with internal index @p u.
	 *
	 * Each AssemblyHelper::a_minus[@p u][@p v][@p k][@p l] corresponds to a single term \f$a^-_{\nu\epsilon} (u^\Omega_\epsilon)^-\f$,
	 * where	@p get<@p 0>(AssemblyHelper::a_minus[@p u][@p v][@p k][@p l]) is the coefficient \f$a^-_{\nu\epsilon}\f$,
	 * 			@p get<@p 1>(AssemblyHelper::a_minus[@p u][@p v][@p k][@p l]) is the global component corresponding to the independent field \f$u^\Omega_\epsilon\f$,
	 * and 		@p get<@p 2>(AssemblyHelper::a_minus[@p u][@p v][@p k][@p l]) contains all the (scalar functional related) shape function indices of a domain cell on the minus side
	 * corresponding to shape functions which are non-zero in the relevant global component.
	 */
	std::vector<std::vector<std::vector<
							std::vector<  std::tuple<const double, const unsigned int, std::vector<unsigned int>> >
							>>>
	a_minus;

	/**
	 * This member contains information about how interface related dependent fields \f$e^\Sigma_\nu\f$ are
	 * related to shape functions. In particular, it relates shape functions to the terms \f$b^-_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^-\f$
	 * (see DependentField for further information).
	 *
	 * AssemblyHelper::b_minus[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the interface (sub)portion with internal index @p u.
	 *
	 * Each AssemblyHelper::b_minus[@p u][@p v][@p k][@p l] corresponds to a single term \f$b^-_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^-\f$,
	 * where	@p get<@p 0>(AssemblyHelper::b_minus[@p u][@p v][@p k][@p l]) is the coefficient \f$b^-_{\nu\epsilon i}\f$,
	 * 			@p get<@p 1>(AssemblyHelper::b_minus[@p u][@p v][@p k][@p l]) is the global component corresponding to the independent field \f$u^\Omega_\epsilon\f$,
	 * 			@p get<@p 2>(AssemblyHelper::b_minus[@p u][@p v][@p k][@p l]) is the derivative \f$i\f$ of the independent field \f$u^\Omega_\epsilon\f$,
	 * and		@p get<@p 3>(AssemblyHelper::b_minus[@p u][@p v][@p k][@p l]) contains all the (scalar functional related) shape function indices of a domain cell on the minus side
	 * corresponding to shape functions which are non-zero in the relevant global component.
	 */
	std::vector<std::vector<std::vector<
							std::vector<  std::tuple<const double, const unsigned int, const unsigned int, std::vector<unsigned int>> >
							>>>
	b_minus;

	/**
	 * This member contains information about how interface related dependent fields \f$e^\Sigma_\nu\f$ are
	 * related to shape functions. In particular, it relates shape functions to the terms \f$a^+_{\nu\epsilon} (u^\Omega_\epsilon)^+\f$
	 * (see DependentField for further information).
	 *
	 * AssemblyHelper::a_plus[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the interface (sub)portion with internal index @p u.
	 *
	 * Each AssemblyHelper::a_plus[@p u][@p v][@p k][@p l] corresponds to a single term \f$a^+_{\nu\epsilon} (u^\Omega_\epsilon)^+\f$,
	 * where	@p get<@p 0>(AssemblyHelper::a_plus[@p u][@p v][@p k][@p l]) is the coefficient \f$a^+_{\nu\epsilon}\f$,
	 * 			@p get<@p 1>(AssemblyHelper::a_plus[@p u][@p v][@p k][@p l]) is the global component corresponding to the independent field \f$u^\Omega_\epsilon\f$,
	 * and 		@p get<@p 2>(AssemblyHelper::a_plus[@p u][@p v][@p k][@p l]) contains all the (scalar functional related) shape function indices of a domain cell on the plus side
	 * corresponding to shape functions which are non-zero in the relevant global component.
	 */
	std::vector<std::vector<std::vector<
							std::vector<  std::tuple<const double, const unsigned int, std::vector<unsigned int>> >
							>>>
	a_plus;

	/**
	 * This member contains information about how interface related dependent fields \f$e^\Sigma_\nu\f$ are
	 * related to shape functions. In particular, it relates shape functions to the terms \f$b^+_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^+\f$
	 * (see DependentField for further information).
	 *
	 * AssemblyHelper::b_plus[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the interface (sub)portion with internal index @p u.
	 *
	 * Each AssemblyHelper::b_plus[@p u][@p v][@p k][@p l] corresponds to a single term \f$b^+_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^+\f$,
	 * where	@p get<@p 0>(AssemblyHelper::b_plus[@p u][@p v][@p k][@p l]) is the coefficient \f$b^+_{\nu\epsilon i}\f$,
	 * 			@p get<@p 1>(AssemblyHelper::b_plus[@p u][@p v][@p k][@p l]) is the global component corresponding to the independent field \f$u^\Omega_\epsilon\f$,
	 * 			@p get<@p 2>(AssemblyHelper::b_plus[@p u][@p v][@p k][@p l]) is the derivative \f$i\f$ of the independent field \f$u^\Omega_\epsilon\f$,
	 * and		@p get<@p 3>(AssemblyHelper::b_plus[@p u][@p v][@p k][@p l]) contains all the (scalar functional related) shape function indices of a domain cell on the plus side
	 * corresponding to shape functions which are non-zero in the relevant global component.
	 */
	std::vector<std::vector<std::vector<
							std::vector<  std::tuple<const double, const unsigned int, const unsigned int, std::vector<unsigned int>> >
							>>>
	b_plus;

	/**
	 * This member contains information about how interface related dependent fields \f$e^\Sigma_\nu\f$ are
	 * related to independent scalars. In particular, it relates scalar functional related independent scalar
	 * indices to terms \f$c^\Sigma_{\nu\iota} C_\iota\f$
	 * (see DependentField for further information).
	 *
	 * AssemblyHelper::c_sigma[@p u][@p v][@p k] is used for the evaluation of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the interface (sub)portion with internal index @p u.
	 *
	 * Each AssemblyHelper::c_sigma[@p u][@p v][@p k][@p l] corresponds to a single term \f$c^\Sigma_{\nu\iota} C_\iota\f$,
	 * where	@p get<@p 0>(AssemblyHelper::c_sigma[@p u][@p v][@p k][@p l]) is the coefficient \f$c^\Sigma_{\nu\iota}\f$,
	 * and		@p get<@p 1>(AssemblyHelper::c_sigma[@p u][@p v][@p k][@p l]) is the scalar functional related independent scalar index corresponding to the independent scalar \f$C_\iota\f$
	 */
	std::vector<std::vector<std::vector<
							std::vector<std::tuple<const double, unsigned int>>
							>>>
	c_sigma;

	/**
	 * This member contains information about the constant terms \f$d^\Sigma_\nu\f$ of domain related dependent fields \f$e^\Omega_\lambda\f$.
	 * (see DependentField for further information).
	 *
	 * AssemblyHelper::d_sigma[@p u][@p v][@p k] contains \f$d^\Sigma_\nu\f$ of the @p k-th dependent variable of the
	 * @p v-th scalar functional on the interface (sub)portion with internal index @p u.
	 * For computation of a dependent variable \f$d^\Sigma_{\nu}\f$ is added to the other contributions to the dependent variable.
	 */
	std::vector<std::vector<std::vector<double>>>
	d_sigma;

	/**
	 * This member contains information about whether the dependent field is local.
	 *
	 * AssemblyHelper::e_sigma_local[@p u][@p v][@p k] contains whether the @p k-th dependent field of the
	 * @p v-th scalar functional on the domain portion with internal index @p u is local.
	 */
	std::vector<std::vector<std::vector<bool>>>
	e_sigma_local;


///@}

/** @name Members containing the FEValues, FEFaceValues, and FESubfaceValues objects needed during assembly
 *
 *  These data structures make sure that the FEValues, FEFaceValues, and FESubfaceValues objects
 *  are re-used wherever possible and that only those objects are reinitialized for which
 *  this is really necessary when a new cell is visited.
 */
///@{

	/**
	 * FEValues objects needed for evaluation of scalar functionals on domain cells.
	 *
	 * AssemblyHelper::fe_values_domain[@p u][@p v] is used for the evaluation of the
	 * @p v-th scalar functional on the domain portion with internal index @p u.
	 *
	 * Using a shared pointer here allows to re-use FEValues objects if there are
	 * several scalar functionals defined on a certain domain portion and the
	 * integration of two or more of these scalar functionals is done based on the same
	 * quadrature rule.
	 */
	std::vector<std::vector<std::shared_ptr<FEValues<spacedim, spacedim>>>>
	fe_values_domain;

	/**
	 * FEValues objects needed for evaluation of scalar functionals on interface cells
	 *
	 * AssemblyHelper::fe_values_interface[@p u][@p v] is used for the evaluation of the
	 * @p v-th scalar functional on the interface (sub)portion with internal index @p u.
	 *
	 * @note	Using a shared pointers here allows to re-use FEValuesInterface objects,
	 * 			see also the note to AssemblyHelper::fe_values_domain
	 */
	std::vector<std::vector<std::shared_ptr<FEValuesInterface<spacedim>>>>
	fe_values_interface;

	/**
	 * AssemblyHelper::fe_values_domain_reinit[@p i] contains the FEValues objects for the domain portion with internal index @p i,
	 * which are actually needed when a particular cell of this domain portion is visited. By using this
	 * data structure, multiple initializations of FEValues objects with the same cell are avoided.
	 */
	std::vector<std::set<std::shared_ptr<FEValues<spacedim, spacedim>>>>
	fe_values_domain_reinit;

	/**
	 * AssemblyHelper::fe_values_interface_reinit[@p i] contains the FEValuesInterface objects for the interface (sub)portion with internal index @p i,
	 * which are actually needed when a particular cell of this interface (sub)portion is visited. By using this
	 * data structure, multiple initializations of FEValuesInterface objects with the same cell are avoided.
	 */
	std::vector<std::set<std::shared_ptr<FEValuesInterface<spacedim>>>>
	fe_values_interface_reinit;

	/**
	 *  The same as AssemblyHelper::fe_values_domain_reinit, but only contains those FEValues objects needed for
	 *  evaluation of scalar functionals entering non-primitively into the total potential
	 */
	std::vector<std::set<std::shared_ptr<FEValues<spacedim, spacedim>>>>
	fe_values_domain_reinit_nonprimitive;

	/**
	 * The same as AssemblyHelper::fe_values_interface_reinit, but only contains those FEValuesInterface objects needed for
	 * evaluation of scalar functionals entering non-primitively into the total potential
	 */
	std::vector<std::set<std::shared_ptr<FEValuesInterface<spacedim>>>>
	fe_values_interface_reinit_nonprimitive;

///@}


/** @name Miscellaneous members */
///@{

	/**
	 * AssemblyHelper::component_names_domain[@p i] is a name (name characterized by string and component) for the
	 * global component with global index @p i on the domain. This name is formed from the pair (IndependentField::name,
	 * @p component) of the underlying domain related independent field. The name is used to identify independent fields in
	 * output.
	 */
	std::vector<std::pair<std::string, unsigned int>>
	component_names_domain;

	/**
	 * AssemblyHelper::component_names_interface[@p i] is a name (name characterized by string and component) for the
	 * global component with global index @p i on the interface. This name is formed from the pair (IndependentField::name,
	 * @p component) of the underlying interface related independent field. The name is used to identify independent fields in
	 * output.
	 */
	std::vector<std::pair<std::string, unsigned int>>
	component_names_interface;

	/**
	 * A list of connections with which this object connects to the TriangulationSystem DoFHandlerSystem::tria_system
	 * to get notice when it changes.
	 */
	std::vector<boost::signals2::connection>
	tria_listeners;

	/**
	 * The id of this processor (which is always zero in sequential computations, and corresponds to the rank of the processor otherwise)
	 */
	const unsigned int
	this_proc;

	/**
	 * The number of participating processors
	 */
	const unsigned int
	n_procs;

	/**
	 * Stream used for standard output to screen (makes sure that output is printed only on one processor in parallel)
	 */
	ConditionalOStream
	pout;

	/**
	 * Absolute tolerance for checking proper alignment of corresponding quadrature points on interfaces (quadrature points on domain cell faces and corresponding interface cells)
	 */
	double
	quadrature_point_alignment_tol = 1e-8;

///@}


/** @name Member functions used during construction of an AssemblyHelper object
 *  These functions are essentially introduced to clean up tthe constructor of the class a bit.
 *  They are only used during construction of an AssemblyHelper object.
 */
///@{

	/**
	 * This function converts the dependent field definitions into a format suitable for
	 * assembly of the finite element system (essentially by defining how each dependent
	 * field is related to the shape functions and, possibly, the derivatives thereof).
	 * This member initializes the following member variables:
	 *  AssemblyHelper::coupled_dof_indices_scalar_functionals_domain,
	 * 	AssemblyHelper::coupled_dof_indices_scalar_functionals_interface,
	 * 	AssemblyHelper::coupled_dof_indices_scalar_functionals_interface_minus,
	 * 	AssemblyHelper::coupled_dof_indices_scalar_functionals_interface_plus,
	 * 	AssemblyHelper::coupled_C_indices_scalar_functionals_domain,
	 * 	AssemblyHelper::coupled_C_indices_scalar_functionals_interface,
	 * 	AssemblyHelper::a_omega,
	 * 	AssemblyHelper::b_omega,
	 * 	AssemblyHelper::c_omega,
	 * 	AssemblyHelper::d_omega,
	 * 	AssemblyHelper::a_sigma,
	 * 	AssemblyHelper::b_sigma,
	 * 	AssemblyHelper::a_minus,
	 * 	AssemblyHelper::b_minus,
	 * 	AssemblyHelper::a_plus,
	 * 	AssemblyHelper::b_plus,
	 * 	AssemblyHelper::c_sigma,
	 * 	AssemblyHelper::d_sigma.
	 */
	void
	convert_dependent_fields_to_shapefunctions();

	/**
	 * %Function initializing the hidden variables of the scalar functionals. This uses
	 * the functions ScalarFunctional::initial_vals_hidden and ScalarFunctional<spacedim, spacedim>::initial_vals_hidden,
	 * respectively.
	 *
	 * @warning		At present the transfer of hidden variables upon mesh refinement is not implemented.
	 * 				So, use hidden variables only if the mesh is not changed after the AssemblyHelper has been set up!
	 */
	void
	initialize_hidden_variables()
	const;

///@}

/** @name Member functions related to assembly of the finite element system
 *  These functions are all needed during assembly of the finite element system.
 */
///@{

//member functions related to assembly of the finite element system

	/**
	 * This distributes the dofs. This function is called automatically after construction of an AssemblyHelper
	 * object as well as after a change of the mesh. It makes sure that the active fe indices are updated before
	 * distributing the dofs.
	 */
	void
	distribute_dofs();

	/**
	 * Method to initialize the required FEValues objects for assembly of the contribution of a domain cell
	 * to the finite element system.
	 *
	 * @param[in]	cell			The domain cell
	 *
	 * @param[in]	internal_index	The internal index of the domain portion the domain cell belongs to.
	 * 								In principle, the function could determine the @p internal_index from
	 * 								@p cell. However, typically it is already known when the function is called
	 * 								and therefore it is passed here as an argument.
	 *
	 * @param[in]	nonprimitive	if @p true: initialize only those FEValues objects needed for evaluation of
	 * 								the scalar functionals entering non-primitively into the total potential
	 */
	void
	initialize_fe_values_domain(const typename DoFHandler<spacedim, spacedim>::active_cell_iterator& 	cell,
								const unsigned int 															internal_index,
								const bool 																	nonprimitive = false )
	const;

	/**
	 * Method to initialize the required FEValuesInterface objects for assembly of the contribution of an interface cell
	 * (and the neighboring domain cells) to the finite element system.
	 *
	 * @param[in]	interface_cell_domain_cells	The interface cell (and the neighboring domain cells)
	 *
	 * @param[in]	internal_index				The internal index of the interface (sub)portion the interface cell belongs to.
	 * 											In principle, the function could determine the @p internal_index from
	 * 											@p interface_cell_domain_cell. However, typically it is already known when the function is called
	 * 											and therefore it is passed here as an argument.
	 *
	 * @param[in]	nonprimitive				if @p true: initialize only those FEValuesInterface objects needed for evaluation of
	 * 											the scalar functionals entering non-primitively into the total potential
	 */
	void
	initialize_fe_values_interface(	const InterfaceCellDomainCellsDoF<spacedim>& 	interface_cell_domain_cells,
									const unsigned int								internal_index,
									const bool										nonprimitive = false )
	const;

	/**
	 * Method to compute the dependent fields on the domain (\f$e^\Omega_\lambda\f$) and the derivatives w.r.t. the relevant dofs
	 * at a quadrature point.
	 *
	 * It is assumed that all quantities except @p de_omega_dsol_T have the correct size when passed in and
	 * that the required FEValues objects are properly initialized.
	 *
	 * @param[in]	internal_index				internal index of domain portion
	 *
	 * @param[in]	scalar_functional_index		scalar functional index within domain portion
	 *
	 * @param[in]	q_point						quadrature point
	 *
	 * @param[in]	solution_u_omega			local solution vector for dofs related to domain related independent fields (in scalar functional related shape function indexing)
	 *
	 * @param[in]	solution_C					local solution vector for dofs related to independent scalars (in scalar functional related independent scalar indexing)
	 *
	 * @param[out]	e_omega						computed values of the dependent fields
	 *
	 * @param[out]	de_omega_dsol_T				derivatives of the dependent fields w.r.t. the local dofs
	 *											(each row in @p de_omega_dsol_T corresponds either to a domain related dof or to an independent scalar;
	 *											 the domain related dofs come first with the same indexing as for @p solution_u_omega, and the independent scalars follow with the
	 *											 same indexing as for @p solution_C; the size of @p de_omega_dsol_T is (@p solution_omega.size() + @p solution_C.size()) x @p e_omega.size())
	 *
	 * @param[in]	compute_derivative			indicates whether @p de_omega_dsol_T is to be computed or not
	 *
	 * @param[in]	ignore_constants			If @p true, the constant terms in the dependent fields are ignored (this is required for the computation of increments of dependent fields upon increments
	 * 											of the solution)
	 */
	void
	compute_e_omega(const unsigned int		internal_index,
					const unsigned int		scalar_functional_index,
					const unsigned int		q_point,
					const Vector<double>&	solution_u_omega,
					const Vector<double>&	solution_C,
					Vector<double>&			e_omega,
					FullMatrix<double>& 	de_omega_dsol_T,
					const bool				compute_derivative = true,
					const bool				ignore_constants = false)
	const;

	/**
	 * Method to compute the dependent fields on the domain (\f$e^\Sigma_\nu\f$) and the derivatives w.r.t. the relevant dofs
	 * at a quadrature point.
	 *
	 * It is assumed that all quantities except @p de_sigma_dsol_T have the correct size when passed in and
	 * that the required FEValuesInterface objects are properly initialized.
	 *
	 * The indexing of dofs used here is a bit awkward, but is necessary to treat duplicate dofs on the minus and the plus side.
	 * See Auxiliary::combine_dof_indices() for further information.
	 *
	 * @param[in]	internal_index								internal index of interface (sub)portion
	 *
	 * @param[in]	scalar_functional_index						scalar functional index within interface (sub)portion
	 *
	 * @param[in]	q_point										quadrature point
	 *
	 * @param[in]	solution_u_sigma							local solution vector for dofs related to interface related independent fields (in scalar functional related shape function indexing)
	 *
	 * @param[in]	solution_u_omega_minus						local solution vector for dofs related to domain related independent fields on minus side (in scalar functional related shape function indexing)
	 *
	 * @param[in]	solution_u_omega_plus						local solution vector for dofs related to domain related independent fields on plus side (in scalar functional related shape function indexing)
	 *
	 * @param[in]	solution_C									local solution vector for dofs related to independent scalars (in scalar functional related independent scalar indexing)
	 *
	 * @param[in]	dof_indices_interface_dof_indices_combined	@p dof_indices_interface_dof_indices_combined[@p i] is the row in @p de_sigma_dsol_T and @p dof_indices_global_combined
	 * 															corresponding to @p solution_u_sigma[@p i]
	 *
	 * @param[in]	dof_indices_minus_dof_indices_combined		@p dof_indices_minus_dof_indices_combined[@p i] is the row in @p de_sigma_dsol_T and @p dof_indices_global_combined
	 * 															corresponding to @p solution_u_omega_minus[@p i]
	 *
	 * @param[in]	dof_indices_plus_dof_indices_combined		@p dof_indices_plus_dof_indices_combined[@p i] is the row in @p de_sigma_dsol_T and @p dof_indices_global_combined
	 * 															corresponding to @p solution_u_omega_plus[@p i]
	 *
	 * @param[in]	dof_indices_C_dof_indices_combined			@p dof_indices_C_dof_indices_combined[@p i] is the row in @p de_sigma_dsol_T and @p dof_indices_global_combined
	 * 															corresponding to @p solution_C[@p i]
	 *
	 * @param[in]	dof_indices_global_combined					combined global dof indices, see also Auxiliary::combine_dof_indices(); this parameter is not really used
	 * 															inside the function except to determine the correct size of @p de_sigma_dsol_T
	 *
	 * @param[out]	e_sigma										computed values of the dependent fields
	 *
	 * @param[out]	de_sigma_dsol_T								derivatives of the dependent fields w.r.t. the local dofs
	 *															(each row in @p de_sigma_dsol_T corresponds either to an interface related dof, or to a domain related dof on the minus side,
	 *															or to a domain related dof on the plus side, or to an independent scalar;
	 *											 				the indexing is determined by @p dof_indices_interface_dof_indices_combined,
	 *											 				@p dof_indices_minus_dof_indices_combined, @p dof_indices_plus_dof_indices_combined,
	 *											 				@p dof_indices_C_dof_indices_combined and is consistent with the indexing in @p dof_indices_global_combined;
	 *											 				the size of @p de_sigma_dsol_T is @p dof_indices_global_combined.size() x @p e_sigma.size() )
	 *
	 * @param[in]	compute_derivative							indicates whether @p de_sigma_dsol_T is to be computed or not
	 *
	 * @param[in]	ignore_constants							If @p true, the constant terms in the dependent fields are ignored (this is required for the computation of increments of dependent fields upon increments
	 * 															of the solution)
	 */
	void
	compute_e_sigma(const unsigned int					internal_index,
					const unsigned int					scalar_functional_index,
					const unsigned int					q_point,
					const Vector<double>&				solution_u_sigma,
					const Vector<double>&				solution_u_omega_minus,
					const Vector<double>&				solution_u_omega_plus,
					const Vector<double>&				solution_C,
					const std::vector<unsigned int>&	dof_indices_interface_dof_indices_combined,
					const std::vector<unsigned int>&	dof_indices_minus_dof_indices_combined,
					const std::vector<unsigned int>&	dof_indices_plus_dof_indices_combined,
					const std::vector<unsigned int>&	dof_indices_C_dof_indices_combined,
					const std::vector<unsigned int>&	dof_indices_global_combined,
					Vector<double>&						e_sigma,
					FullMatrix<double>&					de_sigma_dsol_T,
					const bool							compute_derivative = true,
					const bool							ignore_constants = false)
	const;

	/**
	 * Method computing the values of the scalar functionals entering non-primitively into the total potential
	 *
	 * This method is required since the assembly of the finite element system cannot be done before the current values
	 * of the scalar functionals entering non-primitively into the total potential are known.
	 *
	 * @param[in]	solution								the global solution vector
	 *
	 * @param[in]	solution_ref_sets						a set of reference solution vectors (e.g. solutions at previous time steps),
	 * 														which may enter into the calculation of the scalar functionals
	 *
	 * @param[out]	nonprimitive_scalar_functional_values	Values of scalar functionals entering non-primitively into
	 * 														the total potential. The indexing is according to AssemblyHelper::scalar_functionals_domain_nonprimitive_indices
	 * 														and AssemblyHelper::scalar_functionals_interface_nonprimitive_indices
	 *
	 * @return												@p false if the evaluation of the scalar functionals was successful, and @p true
	 * 														if an error prevented the proper calculation of these quantities
	 *
	 * @tparam		VectorType								The type used for the solution vector
	 */
	template<class VectorType>
	bool
	get_nonprimitive_scalar_functional_values(	const VectorType&						solution,
												const std::vector<const VectorType*>	solution_ref_sets,
												Vector<double>& 						nonprimitive_scalar_functional_values )
	const;

	/**
	 * %Function returning nonprimitive index of a domain related scalar functional according to AssemblyHelper::scalar_functionals_domain_nonprimitive_indices
	 * as well as the primitive index of that scalar functional according to AssemblyHelper::scalar_functionals_domain_primitive_indices.
	 * Note that both indices may exist, because a certain scalar functional may enter different TotalPotentialContribution objects, with some being primitive and
	 * others non-primitive. If @p scalar_functional is not related to any non-primitive TotalPotentialContribution, the first of the
	 * returned indices will be @p -1; and if @p scalar_functional is not related to any primitive TotalPotentialContribution, the second
	 * of the returned indices is @p -1.
	 *
	 * @param[in]	scalar_functional	the domain related scalar functional for which the indices are requested
	 *
	 * @return	(non-primitive index, primitive index)
	 */
	std::pair<const int, const int>
	get_scalar_functional_indices(const ScalarFunctional<spacedim, spacedim>* scalar_functional)
	const;

	/**
	 * %Function returning nonprimitive index of an interface related scalar functional according to AssemblyHelper::scalar_functionals_interface_nonprimitive_indices
	 * as well as the primitive index of that scalar functional according to AssemblyHelper::scalar_functionals_interface_primitive_indices.
	 * Note that both indices may exist, because a certain scalar functional may enter different TotalPotentialContribution objects, with some being primitive and
	 * others non-primitive. If @p scalar_functional is not related to any non-primitive TotalPotentialContribution, the first of the
	 * returned indices will be @p -1; and if @p scalar_functional is not related to any primitive TotalPotentialContribution, the second
	 * of the returned indices is @p -1.
	 *
	 * @param[in]	scalar_functional	the interface related scalar functional for which the indices are requested
	 *
	 * @return	(non-primitive index, primitive index)
	 */
	std::pair<const int, const int>
	get_scalar_functional_indices(const ScalarFunctional<spacedim-1, spacedim>* scalar_functional)
	const;

	/**
	 * @param[out]	global_dof_indices_C	vector with the global dof indices corresponding to the independent scalars stored in AssemblyHelper::C
	 */
	void
	get_dof_indices_C(std::vector<unsigned int>& global_dof_indices_C)
	const;

	/**
	 * %Auxiliary function for creation of Dirichlet constraints
	 *
	 * @param[in]		domain_cell					The domain cell for which constraints are to be created
	 *
	 * @param[in]		face						The face of the domain cell on which constraints are to be created
	 *
	 * @param[in]		shapefuns					The shape functions to be constrained
	 *
	 * @param[in]		constraint					The constraint object describing the constraint
	 *
	 * @param[inout]	constraint_matrix			The constraint matrix into which the constraints go
	 *
	 * @param[in]		constraints_ignore			A constraint matrix with indices which should not be constrained again
	 */
	void
	make_dirichlet_constraints_recursion(const typename TriangulationSystem<spacedim>::DomainCell&	domain_cell,
										 const unsigned int											face,
										 const std::vector<unsigned int>&							shapefuns,
										 const DirichletConstraint<spacedim>&						constraint,
										 AffineConstraints<double>&									constraint_matrix,
										 const AffineConstraints<double>&							constraints_ignore)
	const;

///@}

public:

	/**
	 * The constructor of the class.
	 *
	 * @param[in]	total_potential			AssemblyHelper::total_potential (note: the total potential is copied over by the constructor; but as it contains pointers to the TotalPotentialContribution objects,
	 * 										changing the latter or the associated ScalarFunctional objects will affect the total potential)
	 *
	 * @param[in]	tria_system				AssemblyHelper::tria_system
	 *
	 * @param[in]	mapping_domain			AssemblyHelper::mapping_domain
	 *
	 * @param[in]	mapping_interface		AssemblyHelper::mapping_interface
	 *
	 * @param[in]	independent_scalars		Additional independent scalars to be included into the finite element system which are not appearing in the total potential
	 * 										(this may e.g. be constants appearing only in constraints)
	 */
	AssemblyHelper(	const TotalPotential<spacedim>&							total_potential,
					TriangulationSystem<spacedim>&							tria_system,
					const Mapping<spacedim, spacedim>&						mapping_domain,
					const Mapping<spacedim-1, spacedim>&					mapping_interface,
					const std::set<const IndependentField<0, spacedim>*>&	independent_scalars = std::set<const IndependentField<0, spacedim>*>());

	/**
	 * Destructor. The main task of the destructor is to release the memory allocated for hidden variables.
	 */
	~AssemblyHelper();

/** @name Methods for assembly of the finite element system */
///@{

	/**
	 * Method returning the global initial fields vector (i.e., a vector containing the initial dof values) as defined by IndependentField::initial_vals
	 * and IndependentField<0, spacedim>::initial_value.
	 * The vector has the format \f$\begin{pmatrix}
   	 *								\boldsymbol{u}^\mathrm{initial} \\
   	 *								\boldsymbol{0}
 	 *								\end{pmatrix}\f$.
 	 * Here, \f$\boldsymbol{u}^\mathrm{initial}\f$ contains the initial dof values, with the dofs of AssemblyHelper::dof_handler_system being first,
 	 * followed by the independent scalar dofs  (i.e., the size of \f$\boldsymbol{u}^\mathrm{initial}\f$ is
 	 *   AssemblyHelper::dof_handler_system.@link DoFHandlerSystem::n_dofs_domain() n_dofs_domain() @endlink +
 	 *   AssemblyHelper::dof_handler_system.@link DoFHandlerSystem::n_dofs_interface() n_dofs_interface() @endlink +
 	 *   AssemblyHelper::C@p .size()).
 	 * The \f$\boldsymbol{0}\f$ vector corresponds to the additional lines in the finite
 	 * element system associated with the scalar functionals and independent scalars entering the total potential non-primitively
 	 * (i.e., it has size AssemblyHelper::n_scalar_functionals_nonprimitive + AssemblyHelper::C@p .size()). Note, in this context,
 	 * that the independent scalars always enter the total potential non-primitively.
	 *
	 * @param[out]	initial_fields			initial fields (the vector must be passed in with the correct size,
	 * 										the appropriate size can be queried by AssemblyHelper::system_size())
	 *
	 * @param[in]	constraints				if a constraint matrix is passed,
	 * 										the initial fields will be made consistent with these constraints
	 *
	 * @tparam		VectorType				The type used for the initial fields vector
	 *
	 * @warning		Currently, this method only works if the underlying base elements are either associated with support points or only have a single vector component.
	 * 				In the former case, the values at the support points are assigned. In the latter case, a constant value based on the initial field value at the center
	 * 				is prescribed for each cell if the element has constant modes allowing for this procedure.
	 * 				Other cases are not currently implemented.
	 */
	template<class VectorType>
	void
	get_initial_fields_vector(	VectorType&							initial_fields,
								const AffineConstraints<double>*	constraints = nullptr)
	const;

	/**
	 * Method generating Dirichlet type constraints on domain related independent fields on interfaces.
	 *
	 * @param[out]	constraint_matrix		resulting constraint matrix with the Dirichlet type constraints
	 *
	 * @param[in]	dirichlet_constraints	constraints to be applied to the domain related fields
	 *
	 * @param[in]	constraints_ignore		All dofs which are already constrained in @p ignore_constrained will not again be constrained in @p constraint_matrix
	 *
	 * @warning		The constraint matrix will be cleared before the constraints are written.
	 * 				This means that if you have other constraints, it will be necessary to merge the constraint matrices.
	 * 				The @p local_lines property will be set to what DoFHandlerSystem::get_locally_relevant_dofs() of the object AssemblyHelper::dof_handler_system
	 * 				returns.

	 *
	 * @warning		Domain related independent fields can only be directly constrained using this method if they
	 * 				are discretized based on finite elements having support points! Attempting to constrain
	 * 				independent fields discretized in terms of elements not having support points with this function
	 * 				will have no effect.
	 */
	void
	make_dirichlet_constraints(	AffineConstraints<double>&									constraint_matrix,
								const std::vector<const DirichletConstraint<spacedim>*>&	dirichlet_constraints,
								const AffineConstraints<double>&							constraints_ignore = AffineConstraints<double>())
	const;

	/**
	 * An extended method generating Dirichlet type constraints on domain related independent fields on interfaces as well as point-wise Dirichlet type constraints.
	 *
	 * @param[out]	constraint_matrix		resulting constraint matrix with the Dirichlet type constraints
	 *
	 * @param[in]	dirichlet_constraints	constraints to be applied to the domain related fields
	 *
	 * @param[in]	point_constraints_omega	point constraints to be applied to the domain related fields
	 *
	 * @param[in]	point_constraints_sigma	point constraints to be applied to the interface related fields
	 *
	 * @param[in]	point_constraints_C		point constraints to be applied to the independent scalars
	 *
	 * @param[in]	constraints_ignore		All dofs which are already constrained in @p ignore_constrained will not again be constrained in @p constraint_matrix
	 *
	 * @warning		The constraint matrix will be cleared before the constraints are written.
	 * 				This means that if you have other constraints, it will be necessary to merge the constraint matrices.
	 * 				The @p local_lines property will be set to what DoFHandlerSystem::get_locally_relevant_dofs() of the object AssemblyHelper::dof_handler_system
	 * 				returns.

	 *
	 * @warning		Domain related independent fields can only be directly constrained using this method if they
	 * 				are discretized based on finite elements having support points! Attempting to constrain
	 * 				independent fields discretized in terms of elements not having support points with this function
	 * 				will have no effect.
	 */
	void
	make_dirichlet_constraints(	AffineConstraints<double>&											constraint_matrix,
								const std::vector<const DirichletConstraint<spacedim>*>&			dirichlet_constraints,
								const std::vector<const PointConstraint<spacedim, spacedim>*>&		point_constraints_omega,
								const std::vector<const PointConstraint<spacedim-1, spacedim>*>&	point_constraints_sigma,
								const std::vector<const PointConstraint<0, spacedim>*>&				point_constraints_C,
								const AffineConstraints<double>&									constraints_ignore = AffineConstraints<double>())
	const;

	/**
	 * Generate the current sparsity pattern by simulation of the assembly process,
	 * thereby taking into account the particular structure of the total potential function.
	 *
	 * The sparsity pattern must be initialized correctly outside this function. For the parallel case, don't forget to distribute the sparsity pattern
	 * after this function is called in order to make sure that every processor knows about its non-zero entries (entries in the sparsity pattern
	 * owned by a processor @p i may actually be written by a different processor @p j). See also the deal.II function SparsityTools::distribute_sparsity_pattern()
	 *
	 * @param[out]	dsp_K					the resulting sparsity pattern for the stretched system matrix \f$\boldsymbol{K}^\mathrm{s}\f$
	 *
	 * @param[in]	constraints				constraints to be taken into consideration when building the sparsity pattern
	 *
	 * @tparam		SparsityPatternType		the type of the sparsity pattern
	 */
	template<class SparsityPatternType>
	void
	generate_sparsity_pattern_by_simulation(SparsityPatternType&				dsp_K,
											const AffineConstraints<double>&	constraints)
	const;

	/**
	 * This method performs the actual assembly of the following (stretched) finite element system:
	 *	\f{equation*}
	 *	\boldsymbol{K}^\mathrm{s} \Delta \boldsymbol{\hat u}^\mathrm{s} = -\boldsymbol{f}^\mathrm{s}
	 *  \f}
	 *
	 * For further information about this finite element system, see the accompanying pdf file <a href='../notes/galerkin_tools.pdf'>galerkin_tools.pdf</a>. Note that during assembly constraints are
	 * directly incorporated into the system.
	 *
	 * @param[in]	solution				%Point of linearization \f$\boldsymbol{\hat u}^\mathrm{s} = \begin{pmatrix} \boldsymbol{\hat u} \\ \boldsymbol{\hat \lambda} \end{pmatrix}\f$,
	 * 										where the values of \f$\boldsymbol{\hat \lambda}\f$ are not used during assembly and are, therefore, arbitrary. The
	 * 										appropriate size of @p solution can be obtained from AssemblyHelper::system_size().
	 *
	 * @param[in]	solution_ref_sets		Sets of reference solution vectors, which can e.g. be the solution vectors of previous times steps. Make sure
	 * 										that you pass at least as many reference sets here as are required by the ScalarFunctional objects (see ScalarFunctional::get_h_sigma()
	 * 										and ScalarFunctional<spacedim, spacedim>::get_h_omega()) and the TotalPotentialContribution objects (see TotalPotentialContribution::get_potential_contribution())
	 *										to complete their computations. In general, the order of the sets of reference values in ScalarFunctional::get_h_sigma(),
	 * 										ScalarFunctional<spacedim, spacedim>::get_h_omega(), and TotalPotentialContribution::get_potential_contribution() corresponds to the
	 * 										order in @p solution_ref_sets. If e.g. a function ScalarFunctional::get_h_sigma() requires two reference sets of dependent variable values while @p solution_ref_sets
	 * 										contains more than two reference sets of dof values, the first two sets from @p solution_ref_sets will be used to obtain the reference sets of dependent variable values.
	 *
	 * @param[in]	constraints				Constraints to be applied to the finite element system. These must be the same as those used for
	 * 										generating the sparsity pattern for @p system_matrix with AssemblyHelper::generate_sparsity_pattern_by_simulation().
	 *
	 * @param[out]	potential_value			Value of the total potential (the correct calculation of the total potential value of course requires that all ScalarFunctional::get_h_sigma(),
	 * 										ScalarFunctional<spacedim, spacedim>::get_h_omega() and TotalPotentialContribution::get_potential_contribution() functions compute the respective
	 * 										values and not only compute the derivatives.
	 *
	 * @param[out]	f						The right hand side of the (stretched) finite element system \f$-\boldsymbol{f}^\mathrm{s}\f$ with @p constraints incorporated.
	 * 										This vector must be passed in with the correct size, which can be obtained from AssemblyHelper::system_size().
	 *
	 * @param[out]	K						The (stretched) system matrix \f$\boldsymbol{K}^\mathrm{s}\f$
	 * 										with @p constraints incorporated. The matrix must be passed in with the correct size, which can be obtained from AssemblyHelper::system_size(), and the
	 * 										appropriate sparsity pattern generated by AssemblyHelper::generate_sparsity_pattern_by_simulation() using @p constraints.
	 *
	 * @param[in]	requested_quantities	Tuple indicating which quantities are actually to be computed
	 * 										(e.g. (@p true, @p false, @p true) indicates that @p potential_value and @p K are to be computed)
	 *
	 * @param[out]	local_solution			Map between global dof indices and updated local solution in case that there are local independent fields
	 *
	 * @return								@p false if the assembly process was successful, and @p true
	 * 										if an error prevented proper assembly
	 *
	 * @note								Note that the finite element system returned will not ensure that the constrained dofs have the correct values after the solution of the finite element system
	 * 										I.e., after solution, the values of the constrained dofs must be computed from the unconstrained dofs using distribute_solution().
	 *
	 * @tparam		SolutionVectorType		The type used for @p solution and @p solution_ref_sets (in parallel this vector type must permit read access to ghosted entries while write access is not required)
	 *
	 * @tparam		RHSVectorType			The type used for @p f (in parallel this vector type must permit write access to ghosted entries while read access is not required)
	 *
	 * @tparam		MatrixType				The type used for @p K
	 */
	template<class SolutionVectorType, class RHSVectorType, class MatrixType>
	bool
	assemble_system(const SolutionVectorType&						solution,
					const std::vector<const SolutionVectorType*>	solution_ref_sets,
					const AffineConstraints<double>&				constraints,
					double&											potential_value,
					RHSVectorType&									f,
					MatrixType&										K,
					const std::tuple<bool,bool,bool>				requested_quantities = std::make_tuple(true, true, true),
					std::map<unsigned int, double>*					local_solution = nullptr)
	const;

	/**
	 * Method computing the values of the scalar functionals entering non-primitively into the total potential
	 *
	 * @param[in]	solution										the global solution vector
	 *
	 * @param[in]	solution_ref_sets								a set of reference solution vectors (e.g. solutions at previous time steps),
	 * 																which may enter into the calculation of the scalar functionals
	 *
	 * @param[out]	nonprimitive_scalar_functional_values_domain	Values of domain related scalar functionals entering non-primitively into
	 * 																the total potential.
	 *
	 * @param[out]	nonprimitive_scalar_functional_values_interface	Values of interface related scalar functionals entering non-primitively into
	 * 																the total potential.
	 *
	 * @return														@p false if the evaluation of the scalar functionals was successful, and @p true
	 * 																if an error prevented the proper calculation of these quantities
	 *
	 * @tparam		VectorType										The type used for the solution vector
	 */
	template<class VectorType>
	bool
	get_nonprimitive_scalar_functional_values(	const VectorType&													solution,
												const std::vector<const VectorType*>								solution_ref_sets,
												std::map<const ScalarFunctional<spacedim, spacedim>*, double>& 		nonprimitive_scalar_functional_values_domain,
												std::map<const ScalarFunctional<spacedim-1, spacedim>*, double>& 	nonprimitive_scalar_functional_values_interface)
	const;

	/**
	 * Method simply calling the ScalarFunctional<spacedim, spacedim>::get_h_omega() and ScalarFunctional::get_h_sigma() functions with the solution @p solution and
	 * the reference solution @p solution_ref_sets, without doing any assembly of the finite element system.
	 *
	 * @param[in]	solution								the global solution vector
	 *
	 * @param[in]	solution_ref_sets						a set of reference solution vectors (e.g. solutions at previous time steps),
	 * 														which may enter into the calculation of the scalar functionals
	 *
	 * @param[in]	scalar_functionals_domain_to_call		domain-related scalar functionals to be evaluated
	 *
	 * @param[in]	scalar_functionals_interface_to_call	interface-related scalar functionals to be evaluated
	 *
	 * @tparam		VectorType								The type used for the solution vector
	 */
	template<class VectorType>
	void
	call_scalar_functionals(const VectorType&												solution,
							const std::vector<const VectorType*>&							solution_ref_sets,
							const std::set<const ScalarFunctional<spacedim, spacedim>*>&	scalar_functionals_domain_to_call,
							const std::set<const ScalarFunctional<spacedim-1, spacedim>*>&	scalar_functionals_interface_to_call)
	const;

	/**
	 * Method determining the maximum permissible step length into a certain direction.
	 *
	 * @param[in]	solution										the current global solution vector
	 *
	 * @param[in]	solution_ref_sets								a set of reference solution vectors (e.g. solutions at previous time steps),
	 * 																which may enter into the calculation of the scalar functionals
	 *
	 * @param[in]	delta_solution									The direction into which the solution is updated
	 *
	 * @return														The maximum step length alpha such that solution + alpha*delta_solution is a permissible state
	 *
	 * @tparam		VectorType										The type used for the solution vector
	 */
	template<class VectorType>
	double
	get_maximum_step_length(const VectorType&						solution,
							const std::vector<const VectorType*>	solution_ref_sets,
							const VectorType&						delta_solution)
	const;


	/**
	 * This method compares the finite element system obtained with AssemblyHelper::assemble_system()
	 * with numerically computed equivalents.
	 *
	 * The numerically computed right hand side is based on a finite difference quotient of the total potential,
	 * and the numerically computed system matrix is based on a finite difference quotient of
	 * the right hand side of the finite element system. I.e., the numerically computed right hand side can only be
	 * "correct" (to within the accuracy of the finite difference approach) if the total potential is correctly implemented;
	 * and the numerically computed system matrix can only be "correct" (to within the accuracy of the finite difference approach)
	 * if the right hand side is correctly implemented. This fact can be used to check the implementation, which is the main purpose
	 * of this method.
	 *
	 * For comparison of the finite element systems, the dense system \f$\left( \boldsymbol{K} + \boldsymbol{L} \boldsymbol{\Pi} \boldsymbol{L}^\top \right) \Delta\boldsymbol{\hat u} = - \boldsymbol{f}\f$
	 * is used instead of the stretched system provided by AssemblyHelper::assemble_system(). This, however, means that this method is dealing internally with dense matrices. As a consequence,
	 * the method can only be used for very small test problems.
	 *
	 * @todo Presently, the method does not allow to take into account any constraints. This should be incorporated in future releases of the library.
	 *
	 * @param[in]	solution				%Point of linearization \f$\begin{pmatrix} \boldsymbol{\hat u} \\ \boldsymbol{\hat \lambda} \end{pmatrix}\f$,
	 * 										where the values of \f$\boldsymbol{\hat \lambda}\f$ are not used and are, therefore, arbitrary. The
	 * 										appropriate size of @p solution can be obtained from AssemblyHelper::system_size().
	 *
	 * @param[in]	solution_ref_sets		Sets of reference solution vectors, which can e.g. be the solution vectors of previous times steps. Make sure
	 * 										that you pass at least as many reference sets here as are required by the ScalarFunctional objects (see ScalarFunctional::get_h_sigma()
	 * 										and ScalarFunctional<spacedim, spacedim>::get_h_omega()) and the TotalPotentialContribution objects (see TotalPotentialContribution::get_potential_contribution())
	 *										to complete their computations. In general, the order of the sets of reference values in ScalarFunctional::get_h_sigma(),
	 * 										ScalarFunctional<spacedim, spacedim>::get_h_omega(), and TotalPotentialContribution::get_potential_contribution() corresponds to the
	 * 										order in @p solution_ref_sets. If e.g. a function ScalarFunctional::get_h_sigma() requires two reference sets of dependent variable values while @p solution_ref_sets
	 * 										contains more than two reference sets of dof values, the first two sets from @p solution_ref_sets will be used to obtain the reference sets of dependent variable values.
	 *
	 * @param[in]	detailed_printout_file	A file to which detailed printout is written (if no file name is provided, the results of the comparison will just be written to screen)
	 *
	 * @param[in]	epsilon					Step width for finite difference computation
	 */
	template<class VectorType>
	void
	compare_derivatives_with_numerical_derivatives(	const VectorType& 						solution,
													const std::vector<const VectorType*>	solution_ref_sets,
													const std::string						detailed_printout_file = "",
													const double							epsilon = 1e-8)
	const;

///@}


/** @name Methods for writing output */
///@{

	/**
	 * %Function to write dof output to *.vtu files (one for the domain and one for the interface).
	 * This only writes the dofs on the domain and the interface to the file; the independent scalars
	 * are currently not included in the output.
	 *
	 * @param[in]	solution				The solution \f$\begin{pmatrix} \boldsymbol{\hat u} \\ \boldsymbol{\hat \lambda} \end{pmatrix}\f$,
	 * 										where the values of \f$\boldsymbol{\hat \lambda}\f$ are not used and are, therefore, arbitrary. The
	 * 										appropriate size of @p solution can be obtained from AssemblyHelper::system_size().

	 * @param[in]	file_name_domain		A file name for the output of domain related dofs (note that the extension .vtu will be appended automatically)
	 *
	 * @param[in]	file_name_interface		A file name for the output of interface related dofs (note that the extension .vtu will be appended automatically)
	 *
	 * @param[in]	file_index				An index which is appended to the file name (but in front of the extension)
	 *
	 * @param[in]	dp_domain				DataPostprocessor objects for extra output on domain cells
	 *
	 * @param[in]	dp_interface			DataPostprocessor objects for extra output on interface cells
	 *
	 * @param[in]	n_subdivisions			The number of subdivisions of the cell (to get a better representation in case of curved inner cells, higher order elements, etc.)
	 *
	 * @return	The file names actually used for writing domain and interface output (including the file name extension)
	 */
	template<class VectorType>
	std::pair<const std::string, const std::string>
	write_output_independent_fields(	const VectorType&																	solution,
										const std::string																	file_name_domain,
										const std::string																	file_name_interface,
										const unsigned int																	file_index = 0,
										const std::vector<dealii::SmartPointer<const dealii::DataPostprocessor<spacedim>>>&	dp_domain = std::vector<dealii::SmartPointer<const dealii::DataPostprocessor<spacedim>>>(),
										const std::vector<dealii::SmartPointer<const dealii::DataPostprocessor<spacedim>>>&	dp_interface = std::vector<dealii::SmartPointer<const dealii::DataPostprocessor<spacedim>>>(),
										const unsigned int																	n_subdivisions = 1)
	const;

	/**
	 * Outputs the problem definition for diagnostic purposes to screen (e.g., how the
	 * total potential looks like, which finite elements are used, which quadrature rules are
	 * used, etc.).
	 *
	 * @todo		Currently this function does not print information about independent scalars and the constant terms in
	 * 				the dependent fields. This must be changed.
	 *
	 * @param[in] 	detailed_printout_shapefuns		if @p true print out a detailed summary about how dependent fields are
	 * 												related to shape functions
	 */
	void
	print_assembly_helper_definition(const bool detailed_printout_shapefuns = true)
	const;

///@}

/** @name Methods for comparing solutions
 *  These methods are essentially intended for studies of the convergence behavior. Two types of methods are offered.
 *  The first allows comparison with a different numerical solution obtained on a differently refined mesh; and the
 *  second allows comparison with an analytical solution.
 */
///@{

	/**
	 * %Function computing the "distance" of the solution vector @p solution of this AssemblyHelper
	 * to the solution @p other_solution of another AssemblyHelper object (the AssemblyHelpers must be the same
	 * apart from the mesh refinement, in particular they must be based on the same coarse mesh).
	 *
	 * Note that the values of the independent scalars are currently not taken into account in this method.
	 *
	 * The solution of the other AssemblyHelper is interpolated to the mesh
	 * of this AssemblyHelper, then both solutions are subtracted and finally the norm of the resulting difference is computed numerically on the mesh of
	 * this AssemblyHelper. This is done for the domain related and the interface related part separately.
	 *
	 * @todo Hanging node constraints are currently not taken care of after interpolation of the solution. Also not all VectorTools::NormType norms are implemented yet.
	 *
	 * @param[in]	solution					The solution \f$\begin{pmatrix} \boldsymbol{\hat u} \\ \boldsymbol{\hat \lambda} \end{pmatrix}\f$,
	 * 											where the values of \f$\boldsymbol{\hat \lambda}\f$ are not used and are, therefore, arbitrary. The
	 * 											appropriate size of @p solution can be obtained from AssemblyHelper::system_size().
	 *
	 * @param[in]	other_solution				The solution \f$\begin{pmatrix} \boldsymbol{\hat u} \\ \boldsymbol{\hat \lambda} \end{pmatrix}\f$
	 * 											of the other AssemblyHelper,
	 * 											where the values of \f$\boldsymbol{\hat \lambda}\f$ are not used and are, therefore, arbitrary. The
	 * 											appropriate size of @p solution can be obtained from AssemblyHelper::system_size().
	 *
	 * @param[in]	other_assembly_helper		The other AssemblyHelper object
	 *
	 * @param[in]	quadrature_domain			Quadrature scheme to be used on the domain for the computation of the norm
	 *
	 * @param[in]	quadrature_interface		Quadrature scheme to be used on the interface for the computation of the norm
	 *
	 * @param[in]	norm_type					Type of the norm (note: currently only VectorTools::NormType::@p L2_norm, VectorTools::NormType::@p Linfty_norm, VectorTools::NormType::@p H1_seminorm and VectorTools::NormType::@p W1infty_seminorm are implemented)
	 *
	 * @param[in]	component_mask_domain		Domain related solution components to be included in the calculation. If the ComponentMask
	 * 											is empty, all components will be included
	 *
	 * @param[in]	component_mask_interface	Interface related solution components to be included in the calculation. If the ComponentMask
	 * 											is empty, all components will be included
	 *
	 * @param[in]	exponent					Exponent of the norm if required. Currently this is unused because no norms with variable exponent are implemented.
	 *
	 * @param[in]	scaling_domain				Scaling factors to be used for errors of individual solution components on domain, if empty, scaling factors are set to 1.0
	 *
	 * @param[in]	scaling_interface			Scaling factors to be used for errors of individual solution components on interface, if empty, scaling factors are set to 1.0
	 *
	 * @return									The value of the norm computed on the domain and the interface, respectively
	 */
	template<class VectorType>
	std::pair<const double, const double>
	compute_distance_to_other_solution( const VectorType&				solution,
										const VectorType&				other_solution,
										const AssemblyHelper<spacedim>&	other_assembly_helper,
										const Quadrature<spacedim>		quadrature_domain,
										const Quadrature<spacedim-1>	quadrature_interface,
										const VectorTools::NormType		norm_type = VectorTools::NormType::L2_norm,
										const ComponentMask				component_mask_domain = ComponentMask(),
										const ComponentMask				component_mask_interface = ComponentMask(),
										const double					exponent = 2.0,
										const Vector<double>			scaling_domain = dealii::Vector<double>(),
										const Vector<double>			scaling_interface = dealii::Vector<double>())
	const;

	/**
	 * %Function computing the "distance" of the solution vector @p solution
	 * to an exact solution.
	 *
	 * The exact and the numerical solution are subtracted and finally the norm of the resulting difference is computed numerically on the mesh of
	 * this AssemblyHelper. This is done for the domain related and the interface related part separately.
	 *
	 * Note that the values of the independent scalars are currently not taken into account in this method.
	 *
	 * @param[in]	solution					The solution \f$\begin{pmatrix} \boldsymbol{\hat u} \\ \boldsymbol{\hat \lambda} \end{pmatrix}\f$,
	 * 											where the values of \f$\boldsymbol{\hat \lambda}\f$ are not used and are, therefore, arbitrary. The
	 * 											appropriate size of @p solution can be obtained from AssemblyHelper::system_size().
	 *
	 * @param[in]	exact_solution_domain		Exact solution on domain (use AssemblyHelper::get_u_omega_global_component_indices() to obtain information
	 * 											about the component indexing)
	 *
	 * @param[in]	exact_solution_interface	Exact solution on interface (use AssemblyHelper::get_u_sigma_global_component_indices() to obtain information
	 * 											about the component indexing)
	 *
	 * @param[in]	quadrature_domain			Quadrature scheme to be used on the domain for the computation of the norm
	 *
	 * @param[in]	quadrature_interface		Quadrature scheme to be used on the interface for the computation of the norm
	 *
	 * @param[in]	norm_type					Type of the norm
	 *
	 * @param[in]	component_mask_domain		Domain related solution components to be included in the calculation. If the ComponentMask
	 * 											is empty, all components will be included.
	 *
	 * @param[in]	component_mask_interface	Interface related solution components to be included in the calculation. If the ComponentMask
	 * 											is empty, all components will be included.
	 *
	 * @param[in]	exponent					Exponent of the norm if required
	 *
	 * @return									The value of the norm computed on the domain and the interface, respectively
	 */
	template<class VectorType>
	std::pair<const double, const double>
	compute_distance_to_exact_solution(	const VectorType&				solution,
										const Function<spacedim>&		exact_solution_domain,
										const Function<spacedim>&		exact_solution_interface,
										const Quadrature<spacedim>		quadrature_domain,
										const Quadrature<spacedim-1>	quadrature_interface,
										const VectorTools::NormType		norm_type = VectorTools::NormType::L2_norm,
										const ComponentMask				component_mask_domain = ComponentMask(),
										const ComponentMask				component_mask_interface = ComponentMask(),
										const double					exponent = 2.0)
	const;

///@}


/** @name Methods for querying information about the AssemblyHelper object. */
///@{

	/**
	 * %Function returning the TriangulationSystem
	 *
	 * @return	the TriangulationSystem AssemblyHelper::tria_system
	 */
	const TriangulationSystem<spacedim>&
	get_triangulation_system()
	const;

	/**
	 * %Function returning the TriangulationSystem
	 *
	 * @return	the TriangulationSystem AssemblyHelper::tria_system
	 */
	TriangulationSystem<spacedim>&
	get_triangulation_system();


	/**
	 * %Function returning the DoFHandlerSystem
	 *
	 * @return	the DoFHandlerSystem AssemblyHelper::dof_handler_system
	 */
	const DoFHandlerSystem<spacedim>&
	get_dof_handler_system()
	const;

	/**
	 * %Function returning the DoFHandlerSystem. This function does return a non-const reference and can,
	 * therefore, for example be used for reordering of dof indices if desired.
	 *
	 * @return	the DoFHandlerSystem AssemblyHelper::dof_handler_system
	 */
	DoFHandlerSystem<spacedim>&
	get_dof_handler_system();

	/**
	 * @return 		map between domain related IndependentField objects and the global component indices (in the FESystem) of the first components of the respective
	 * 				IndependentField objects
	 */
	std::map<const IndependentField<spacedim, spacedim>*, const unsigned int>
	get_u_omega_global_component_indices()
	const;

	/**
	 * @return 		map between interface related IndependentField objects and the global component indices (in the FESystem) of the first components of the respective
	 * 				IndependentField objects
	 */
	std::map<const IndependentField<spacedim-1, spacedim>*, const unsigned int>
	get_u_sigma_global_component_indices()
	const;

	/**
	 * Compute the global component index (in the domain related FESystem) of the first component of @p u_omega
	 *
	 * @param[in]	u_omega		The domain related independent field
	 *
	 * @return					global component index
	 */
	unsigned int
	get_u_omega_global_component_index(const IndependentField<spacedim, spacedim>& u_omega)
	const;

	/**
	 * Compute the global component index (in the interface related FESystem) of the first component of @p u_sigma
	 *
	 * @param[in]	u_sigma		The interface related independent field
	 *
	 * @return					global component index
	 */
	unsigned int
	get_u_sigma_global_component_index(const IndependentField<spacedim-1, spacedim>& u_sigma)
	const;

	/**
	 * @return		The size of the stretched finite element system
	 */
	unsigned int
	system_size()
	const;

	/**
	 * @return 	AssemblyHelper::n_scalar_functionals_nonprimitive + AssemblyHelper::C.size()
	 */
	unsigned int
	get_n_stretched_rows()
	const;

	/**
	 * @return	The number of independent scalars attached to the AssemblyHelper
	 */
	unsigned int
	get_n_C()
	const;

	/**
	 * Method returning the global dof index of an independent scalar.
	 *
	 * @param[in]	independent_scalar	independent scalar
	 *
	 * @return							the global dof index corresponding to @p independent_scalar
	 */
	unsigned int
	get_global_dof_index_C(const IndependentField<0,spacedim>* independent_scalar)
	const;

	/**
	 * @return		The set of locally relevant indices (including the stretched rows)
	 */
	const IndexSet
	get_locally_relevant_indices()
	const;

	/**
	 * @return		The set of locally owned dofs (including the stretched rows; ownership of the stretched rows is assigned to the last processor)
	 */
	const IndexSet
	get_locally_owned_indices()
	const;

	/**
	 * @return		Same as AssemblyHelper::get_locally_relevant_indices(), but block-wise (the first index contains the fe related indices, and the second block the rest)
	 */
	const std::vector<IndexSet>
	get_locally_relevant_indices_blocks()
	const;

	/**
	 * @return		Same as AssemblyHelper::get_locally_owned_indices(), but block-wise (the first index set contains the fe related indices, and the second index set the rest)
	 */
	const std::vector<IndexSet>
	get_locally_owned_indices_blocks()
	const;

	/**
	 *
	 * @param[in]	u_omega		independent field
	 *
	 * @param[in]	component	component
	 *
	 * @param[in]	p			point
	 *
	 * @param[in]	ignore_dofs	dofs to ignore during search - this allows for finding duplicate dofs
	 *
	 * @return		global dof index of @p component of @p u_omega at point @p point . This requires that the FE used
	 * 				for the independent field has support points defined.
	 */
	unsigned int
	get_dof_index_at_point_omega(	const IndependentField<spacedim, spacedim>* u_omega,
									const unsigned int							component,
									const Point<spacedim>						p,
									const std::set<unsigned int>				ignore_dofs = std::set<unsigned int>())
	const;

	/**
	 *
	 * @param[in]	u_sigma		independent field
	 *
	 * @param[in]	component	component
	 *
	 * @param[in]	p			point
	 *
	 * @param[in]	ignore_dofs	dofs to ignore during search - this allows for finding duplicate dofs
	 *
	 * @return		global dof index of @p component of @p u_sigma at point @p point . This requires that the FE used
	 * 				for the independent field has support points defined.
	 */
	unsigned int
	get_dof_index_at_point_sigma(	const IndependentField<spacedim-1, spacedim>*	u_sigma,
									const unsigned int								component,
									const Point<spacedim>							p,
									const std::set<unsigned int>					ignore_dofs = std::set<unsigned int>())
	const;

	/**
	 * %Function for printing out information about a dof
	 *
	 * @param[in] 	dof_index	The dof index
	 */
	void
	print_dof_information(const unsigned int dof_index)
	const;


///@}

	/** @name Miscellaneous memberfunctions */
	///@{

		/**
		 * @param[in]	quadrature_point_alignment_tol	Sets AssemblyHelper::quadrature_point_alignment_tol
		 */
		void
		set_quadrature_point_alignment_tol(const double quadrature_point_alignment_tol);

///@}

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_ASSEMBLYHELPER_H_ */
