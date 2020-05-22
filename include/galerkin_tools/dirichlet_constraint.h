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

#ifndef GALERKINTOOLS_DIRICHLETCONSTRAINT_H_
#define GALERKINTOOLS_DIRICHLETCONSTRAINT_H_

#include <vector>
#include <utility>

#include <deal.II/base/function.h>

#include <galerkin_tools/config.h>
#include <galerkin_tools/independent_field.h>
#include <galerkin_tools/triangulation_system.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * Class defining a Dirichlet type interface condition for a domain related
 * independent variable \f$u^\Omega_\epsilon\f$.
 *
 * The Dirichlet type condition has the form
 * \f$ \left(u^\Omega_\epsilon\right)^+ = b^\Omega_\epsilon + c^\Omega_\epsilon C_\iota\f$ or \f$ \left(u^\Omega_\epsilon\right)^- = b^\Omega_\epsilon + c^\Omega_\epsilon C_\iota\f$,
 * where \f$b^\Omega_\epsilon\f$ and \f$c^\Omega_\epsilon\f$ are prescribed functions, which may depend on
 * the spatial position on the interface; and + or - indicates the side of the interface on which the Dirichlet type
 * condition applies; and \f$C_\iota\f$ is a independent scalar (including this term allows to implement b.c.s of the form \f$u^\Omega_\epsilon = \mathrm{const.}\f$)
 *
 * @todo	It would be desirable to also allow for Dirichlet type conditions on interface related fields
 * 			(these conditions would then be imposed on codim 2 objects). This is not implemented yet.
 *
 * The DirichletConstraint class inherits from Subscriptor in order to be
 * able to check that DirichletConstraint objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam	spacedim	Spatial dimension of the problem
 */
template<unsigned int spacedim>
class DirichletConstraint : public Subscriptor
{

public:

	/**
	 * types::material_id%s determining the portions of the interface on which the constraint is applied
	 */
	const std::set<types::material_id>
	domain_of_constraint;

	/**
	 * IndependentField object to be constrained
	 * (defines together with DirichletConstraint::component \f$u^\Omega_\epsilon\f$)
	 */
	const SmartPointer<const IndependentField<spacedim, spacedim>>
	independent_field;

	/**
	 * component of IndependentField object to be constrained
	 * (defines together with DirichletConstraint::independent_field \f$u^\Omega_\epsilon\f$)
	 */
	const unsigned int
	component;

	/**
	 * Side of the interface (::InterfaceSide::@p plus or ::InterfaceSide::@p minus) on which constraint applies;
	 * this information is relevant if the independent field is discontinuous across the interface.
	 */
	const InterfaceSide
	side;

	/**
	 * A Function (or, rather, an instance of a derived class) specifying the constraint inhomogeneity
	 * \f$b^\Omega_\epsilon\f$ in dependence on the spatial position.
	 * The Function must have a single component and implement the method Function::value().
	 * If a @p nullptr is stored here, the constraint will be assumed homogeneous.
	 */
	const SmartPointer<const Function<spacedim>>
	constraint_inhomogeneity;

	/**
	 * \f$C_\iota\f$
	 */
	const SmartPointer<const IndependentField<0, spacedim>>
	independent_scalar;

	/**
	 * A Function (or, rather, an instance of a derived class) specifying
	 * \f$c^\Omega_\epsilon\f$ in dependence on the spatial position.
	 * The Function must have a single component and implement the method Function::value().
	 * If a @p nullptr is stored here, the function will be taken as uniformly one.
	 */
	const SmartPointer<const Function<spacedim>>
	coefficient_c;

	/**
	 * The constructor of the class
	 *
	 * @param[in]	independent_field			DirichletConstraint::independent_field
	 *
	 * @param[in]	component					DirichletConstraint::component
	 *
	 * @param[in]	side						DirichletConstraint::side
	 *
	 * @param[in]	domain_of_constraint		DirichletConstraint::domain_of_constraint
	 *
	 * @param[in]	constraint_inhomogeneity	DirichletConstraint::constraint_inhomogeneity
	 *
	 * @param[in]	independent_scalar			DirichletConstraint::independent_scalar
	 *
	 * @param[in]	coefficient_c				DirichletConstraint::coefficient_c
	 */
	DirichletConstraint(const IndependentField<spacedim, spacedim>& independent_field,
						const unsigned int							component,
						const InterfaceSide							side,
						const std::set<types::material_id> 			domain_of_constraint,
						const Function<spacedim>* const				constraint_inhomogeneity = nullptr,
						const IndependentField<0, spacedim>*		independent_scalar = nullptr,
						const Function<spacedim>* const				coefficient_c = nullptr);

	/**
	 * The destructor of DirichletConstraint essentially checks before destruction that the
	 * DirichletConstraint object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	virtual
	~DirichletConstraint();

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_DIRICHLETCONSTRAINT_H_ */
