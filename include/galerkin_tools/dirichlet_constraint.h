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
private:

	/**
	 * Bool indicating whether constraint is currently active
	 */
	bool
	constraint_is_active = true;

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

	/**
	 * Sets DirichletConstraint::constraint_is_active
	 *
	 * @param[in]	constraint_is_active		Value to be assigned to DirichletConstraint::constraint_is_active
	 */
	void set_constraint_is_active(const bool constraint_is_active);

	/**
	 * Return DirichletConstraint::constraint_is_active
	 */
	bool get_constraint_is_active()
	const;

	/**
	 * Set the time at which the constraint is evaluated
	 *
	 * @param[in]	time	The time at which the constraint is evaluated
	 */
	void
	set_time(const double time)
	const;

};

/**
 * Class defining a Dirichlet type condition for an interface related field at a single point.
 *
 * The Dirichlet type condition has the form
 * \f$u^\Sigma_\eta(\boldsymbol{X}) = b^\Sigma_\eta\f$,
 * where \f$b^\Sigma_\eta\f$ is a prescribed value.
 * The PointConstraint class inherits from Subscriptor in order to be
 * able to check that PointConstraint objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam dim		The dimension of the object on which the independent field is defined.
 * 					Note that also dim=0 is admissible in order to constrain an independent scalar.
 * 					In this case, however, PointConstraint::X is ignored
 *
 * @tparam spacedim	The spatial dimension of the problem
 */
template<unsigned int dim, unsigned int spacedim>
class PointConstraint : public Subscriptor
{
private:

	/**
	 * Bool indicating whether constraint is currently active
	 */
	bool
	constraint_is_active = true;

	/**
	 * point \f$\boldsymbol{X}\f$
	 */
	Point<spacedim>
	X;


public:

	/**
	 * IndependentField object to be constrained
	 * (defines together with PointConstraint::component \f$u^\Sigma_\eta\f$)
	 */
	const SmartPointer<const IndependentField<dim, spacedim>>
	independent_field;

	/**
	 * component of IndependentField object to be constrained
	 * (defines together with PointConstraint::independent_field \f$u^\Sigma_\eta\f$)
	 */
	const unsigned int
	component;

	/**
	 * A Function (or, rather, an instance of a derived class) specifying the constraint inhomogeneity
	 * \f$b^\Sigma_\eta\f$. The position variable of the function object is ignored. However, the Function::set_time() or PointConstraint::set_time() functionality can be used to change the constraint
	 * inhomogeneity with time.
	 * The Function must have a single component and implement the method Function::value().
	 * If a @p nullptr is stored here, the constraint will be assumed homogeneous.
	 */
	const SmartPointer<const Function<spacedim>>
	constraint_inhomogeneity;

	/**
	 * The constructor of the class
	 *
	 * @param[in]	independent_field			PointConstraint::independent_field
	 *
	 * @param[in]	component					PointConstraint::component
	 *
	 * @param[in]	X							PointConstraint::X
	 *
	 * @param[in]	constraint_inhomogeneity	PointConstraint::constraint_inhomogeneity
	 */
	PointConstraint(const IndependentField<dim, spacedim>& 	independent_field,
					const unsigned int						component,
					const Point<spacedim>					X,
					const Function<spacedim>* const			constraint_inhomogeneity = nullptr);

	/**
	 * The destructor of PointConstraint essentially checks before destruction that the
	 * PointConstraint object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	virtual
	~PointConstraint();

	/**
	 * Sets PointConstraint::constraint_is_active
	 *
	 * @param[in]	constraint_is_active		Value to be assigned to PointConstraint::constraint_is_active
	 */
	void
	set_constraint_is_active(const bool constraint_is_active);

	/**
	 * Sets PointConstraint::X
	 *
	 * @param[in]	X		Value to be assigned to PointConstraint::X
	 */
	void
	set_X(const Point<spacedim> X);

	/**
	 * Returns PointConstraint::X
	 *
	 */
	Point<spacedim>
	get_X()
	const;

	/**
	 * Return PointConstraint::constraint_is_active
	 */
	bool
	get_constraint_is_active()
	const;

	/**
	 * Set the time at which the constraint is evaluated
	 *
	 * @param[in]	time	The time at which the constraint is evaluated
	 */
	void
	set_time(const double time)
	const;

};

/**
 * Class defining a Dirichlet type condition for a domain related field at a single point.
 *
 * The Dirichlet type condition has the form
 * \f$u^\Omega_\epsilon(\boldsymbol{X}) = b^\Omega_\epsilon\f$,
 * where \f$b^\Omega_\epsilon\f$ is a prescribed value.
 * The PointConstraint<spacedim,spacedim> class inherits from Subscriptor in order to be
 * able to check that PointConstraint<spacedim,spacedim> objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam	spacedim	Spatial dimension of the problem
 */
template<unsigned int spacedim>
class PointConstraint<spacedim, spacedim> : public Subscriptor
{
private:

	/**
	 * Bool indicating whether constraint is currently active
	 */
	bool
	constraint_is_active = true;

	/**
	 * point \f$\boldsymbol{X}\f$
	 */
	Point<spacedim>
	X;


public:

	/**
	 * IndependentField object to be constrained
	 * (defines together with PointConstraint<spacedim,spacedim>::component \f$u^\Omega_\epsilon\f$)
	 */
	const SmartPointer<const IndependentField<spacedim, spacedim>>
	independent_field;

	/**
	 * component of IndependentField object to be constrained
	 * (defines together with PointConstraint<spacedim,spacedim>::independent_field \f$u^\Omega_\epsilon\f$)
	 */
	const unsigned int
	component;

	/**
	 * A Function (or, rather, an instance of a derived class) specifying the constraint inhomogeneity
	 * \f$b^\Omega_\epsilon\f$. The position variable of the function object is ignored. However, the Function::set_time() or PointConstraint<spacedim, spacedim>::set_time() functionality can be used to change the constraint
	 * inhomogeneity with time.
	 * The Function must have a single component and implement the method Function::value().
	 * If a @p nullptr is stored here, the constraint will be assumed homogeneous.
	 */
	const SmartPointer<const Function<spacedim>>
	constraint_inhomogeneity;

	/**
	 * The constructor of the class
	 *
	 * @param[in]	independent_field			PointConstraint<spacedim,spacedim>::independent_field
	 *
	 * @param[in]	component					PointConstraint<spacedim,spacedim>::component
	 *
	 * @param[in]	X							PointConstraint<spacedim,spacedim>::X
	 *
	 * @param[in]	constraint_inhomogeneity	PointConstraint<spacedim,spacedim>::constraint_inhomogeneity
	 */
	PointConstraint(const IndependentField<spacedim, spacedim>& independent_field,
					const unsigned int							component,
					const Point<spacedim>						X,
					const Function<spacedim>* const				constraint_inhomogeneity = nullptr);

	/**
	 * The destructor of PointConstraint<spacedim,spacedim> essentially checks before destruction that the
	 * PointConstraint<spacedim,spacedim> object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	virtual
	~PointConstraint();

	/**
	 * Sets PointConstraint<spacedim,spacedim>::constraint_is_active
	 *
	 * @param[in]	constraint_is_active		Value to be assigned to PointConstraint<spacedim,spacedim>::constraint_is_active
	 */
	void
	set_constraint_is_active(const bool constraint_is_active);

	/**
	 * Sets PointConstraint<spacedim,spacedim>::X
	 *
	 * @param[in]	X		Value to be assigned to PointConstraint<spacedim,spacedim>::X
	 */
	void
	set_X(const Point<spacedim> X);

	/**
	 * Returns PointConstraint::X
	 *
	 */
	Point<spacedim>
	get_X()
	const;

	/**
	 * Return PointConstraint<spacedim,spacedim>::constraint_is_active
	 */
	bool
	get_constraint_is_active()
	const;

	/**
	 * Set the time at which the constraint is evaluated
	 *
	 * @param[in]	time	The time at which the constraint is evaluated
	 */
	void
	set_time(const double time)
	const;

};


GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_DIRICHLETCONSTRAINT_H_ */
