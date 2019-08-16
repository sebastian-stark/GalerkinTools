// --------------------------------------------------------------------------
// Copyright (C) 2019 by Sebastian Stark
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

#ifndef GALERKINTOOLS_INDEPENDENTFIELD_H_
#define GALERKINTOOLS_INDEPENDENTFIELD_H_

#include <iostream>
#include <vector>
#include <string>

#include <deal.II/fe/fe_system.h>
#include <deal.II/base/function.h>

#include <galerkin_tools/config.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * This class is used to define "independent fields", which are the unknowns of the
 * problem to be solved. Currently domain related and interface related fields are
 * considered, which live on the domain \f$\Omega\f$ and the interface \f$\Sigma\f$,
 * respectively.
 *
 * Domain related fields are denoted by \f$u^\Omega_\epsilon(\boldsymbol{X})\f$,
 * where \f$\epsilon \in E=\left\{1 \hdots N^\mathrm{u,\Omega}\right\}\f$ enumerates the
 * independent fields, and \f$\boldsymbol{X} \in \mathcal{R}^n\f$ is the location in
 * Euclidean \f$n\f$-space.
 *
 * Interface related fields are denoted by \f$u^\Sigma_\eta(\boldsymbol{X})\f$,
 * where \f$\eta \in H=\left\{1 \hdots N^\mathrm{u,\Sigma}\right\}\f$.
 *
 * Each \f$u^\Omega_\epsilon\f$ and \f$u^\Sigma_\eta\f$, respectively, is considered
 * a scalar valued field. Vector and tensor valued independent fields may be
 * represented by combining several scalar valued fields.
 *
 * In order to allow for the use of vector valued finite elements, each IndependentField
 * object may represent a \f$K\f$-vector of scalar valued independent fields. The
 * restriction is, however, that either all \f$K\f$ vector components are discretized
 * with the same (scalar valued) finite element, or that the complete \f$K\f$-vector
 * is discretized by a single (vector valued) finite element.
 *
 * The IndependentField class inherits from Subscriptor in order to be
 * able to check that IndependentField objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam spacedim This template parameter represents the spatial dimension \f$n\f$
 * of the space wherein the problem is formulated.
 *
 * @tparam dim This template parameter represents the dimensionality of the object on
 * which the \f$K\f$-vector of independent fields represented by the IndependentField
 * object lives. E.g., an IndependentField <3,3> represents a field living on the domain
 * in Euclidean \f$3\f$-space, while an IndependentField <2,3> represents a field living
 * on an interface or boundary embedded in Euclidean \f$3\f$-space.
 *
 */
template<unsigned int dim, unsigned int spacedim>
class IndependentField : public Subscriptor
{

public:

	/**
	 * A string with a unique name for the IndependentField object, which is used to identify
	 * it in output and to ensure a well defined ordering of IndependentField objects.
	 */
	const std::string
	name;

	/**
	 * The finite element used for the discretization of the \f$K\f$-vector of independent fields
	 * represented by the IndependentField object. The finite element stored here may be either
	 * scalar valued (in which case IndependentField::n_components copies of IndependentField::fe
	 * are used to discretize all the vector components) or vector valued (in which case its number
	 * of components must equal IndependentField::n_components).
	 */
	const std::unique_ptr<const FiniteElement<dim, spacedim>>
	fe;

	/**
	 * The number of vector components associated with the IndependentField object. In case that
	 * the finite element IndependentField::fe is scalar valued, IndependentField::n_components
	 * corresponds to the multiplicity of the finite element. In case that the finite element
	 * IndependentField::fe is vector valued, the number of components of IndependentField::fe
	 * must correspond to IndependentField::n_components.
	 */
	const unsigned int
	n_components;

	/**
	 * A set of types::material_id indicating on which domain/interface portions the
	 * independent fields represented by the IndependentField object are non-zero.
	 * On all cells of the triangulation having a types::material_id listed in
	 * IndependentField::non_zero_regions, the finite element IndependentField::fe will be used
	 * for discretization. In contrast, on a cell with a material id not listed in
	 * IndependentField::non_zero_regions, the independent fields are assumed to be zero and are
	 * "discretized" by IndependentField::n_components copies of FE_Nothing.
	 */
	const std::set<types::material_id>
	non_zero_regions;

	/**
	 * A Function (or, rather, an instance of a derived class) specifying the initial values of the
	 * independent fields represented by IndependentField as a function of position in space.
	 * The Function must implement the method Function::value(), which
	 * provides with the initial values for all IndependentField::n_components components.
	 */
	const SmartPointer<const Function<spacedim>>
	initial_vals;

	/**
	 * The constructor of the class for the case that IndependentField::n_components copies of
	 * a scalar valued finite element are used for discretization of the
	 * \f$K\f$-vector of independent fields represented by the IndependentField
	 * object.
	 *
	 * @param[in]	name				IndependentField::name
	 *
	 * @param[in]	fe					IndependentField::fe, must be a scalar valued finite element
	 *
	 * @param[in]	n_components		IndependentField::n_components
	 *
	 * @param[in]	non_zero_regions	IndependentField::non_zero_regions
	 *
	 * @param[in]	initial_vals		IndependentField::initial_vals, if a null pointer is
	 * 									provided, the initial values of the independent fields
	 * 									are assumed to be zero.
	 */
	IndependentField(	const std::string 					name,
						const FiniteElement<dim, spacedim>& fe,
						const unsigned int 					n_components,
						const std::set<types::material_id>	non_zero_regions,
						const Function<spacedim> *const		initial_vals = nullptr);

	/**
	 * The constructor of the class for the case that a vector valued finite
	 * element is used for discretization of the \f$K\f$-vector of independent
	 * fields represented by the IndependentField object. With this constructor,
	 * IndependentField::n_components will be deduced from the number of components of the
	 * IndependentField::fe object.
	 *
	 * @param[in]	name				IndependentField::name
	 *
	 * @param[in]	fe					IndependentField::fe
	 *
	 * @param[in]	non_zero_regions	IndependentField::non_zero_regions
	 *
	 * @param[in]	initial_vals		IndependentField::initial_vals, if a null pointer is
	 * 									provided, the initial values of the independent fields
	 * 									are assumed to be zero.
	 */
	IndependentField(	const std::string					name,
						const FiniteElement<dim, spacedim>&	fe,
						const std::set<types::material_id>	non_zero_regions,
						const Function<spacedim> *const		initial_vals = nullptr);

	/**
	 * The destructor of IndependentField essentially checks before destruction that the
	 * IndependentField object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	~IndependentField();

};

/**
 * This  class is used to define an independent scalar \f$C_\iota\f$, where
 * \f$\iota \in I=\left\{1 \hdots N^\mathrm{C}\right\}\f$ enumerates the
 * independent scalars (which may be interpreted as
 * independent fields living on a single point in space).
 */
template<unsigned int spacedim>
class IndependentField<0, spacedim> : public Subscriptor
{

public:

	/**
	 * A string with a unique name for the IndependentField<0, spacedim> object, which is used to identify
	 * it in output and to ensure a well defined ordering of IndependentField<0, spacedim> objects.
	 */
	const std::string
	name;

	/**
	 * The initial value of the independent field.
	 */
	const double
	initial_value;

	/**
	 * The number of vector components associated with the IndependentField<0, spacedim> object. Here,
	 * only IndependentField<0, spacedim>::n_components = 1 is allowed, because no vector valued independent
	 * scalars are considered at present (although they may be represented by several IndependentField<0, spacedim>
	 * objects).
	 */
	const unsigned int
	n_components = 1;

	/**
	 * The constructor of the class.
	 *
	 * @param[in]	name			IndependentField<0, spacedim>::name
	 *
	 * @param[in]	initial_value	IndependentField<0, spacedim>::initial_value
	 */
	IndependentField(	const std::string	name,
						const double		initial_value = 0.0);

	/**
	 * The destructor of IndependentField<0, spacedim> essentially checks before destruction that the
	 * IndependentField<0, spacedim> object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	~IndependentField();

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_INDEPENDENTFIELD_H_ */
