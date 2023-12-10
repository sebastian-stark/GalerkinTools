// --------------------------------------------------------------------------
// Copyright (C) 2023 by Sebastian Stark
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

#ifndef INCLUDE_GALERKIN_TOOLS_FE_DGQ_ANISOTROPIC_H_
#define INCLUDE_GALERKIN_TOOLS_FE_DGQ_ANISOTROPIC_H_

#include <deal.II/base/config.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Implementation of scalar, discontinuous tensor product elements based on support points, which may be distributed anisotropically.
 */
template <int dim, int spacedim = dim>
class FE_DGQAnisotropic : public FE_Poly<dim, spacedim>
{
private:

	/**
	 * Degree in x-direction
	 */
	const unsigned int degree_x;

	/**
	 * Degree in y-direction
	 */
	const unsigned int degree_y;

	/**
	 * Degree in z-direction
	 */
	const unsigned int degree_z;

public:

	/**
	* Constructor for 3d anisotropic tensor product polynomials with support points generated from a set of points represented by deal.II quadrature rules.
	*
	* @param[in]	points_x	support points in x-direction
	*
	* @param[in]	points_y	support points in y-direction
	*
	* @param[in]	points_z	support points in z-direction
	*/
	FE_DGQAnisotropic(	const Quadrature<1>& points_x,
						const Quadrature<1>& points_y,
						const Quadrature<1>& points_z);

	/**
	* Constructor for 2d anisotropic tensor product polynomials with support points generated from a set of points represented by deal.II quadrature rules.
	*
	* @param[in]	points_x	support points in x-direction
	*
	* @param[in]	points_y	support points in y-direction
	*/
	FE_DGQAnisotropic(	const Quadrature<1>& points_x,
						const Quadrature<1>& points_y);

	/**
	* @copydoc FiniteElement::convert_generalized_support_point_values_to_dof_values()
	*/
	virtual void
	convert_generalized_support_point_values_to_dof_values(	const std::vector<Vector<double>>&	support_point_values,
															std::vector<double>&				nodal_values)
	const override;

	/**
	* @copydoc FiniteElement::clone()
	*/
	std::unique_ptr<FiniteElement<dim, spacedim>>
	clone()
	const
	override;

	/**
	* Return a string that uniquely identifies a finite element. This class
	* returns <tt>FE_DGQAnisotropic<dim>(degree_x, degree_y, degree_z)</tt>, with <tt>dim</tt>, <tt>degree_x</tt>, <tt>degree_y</tt> and
	* <tt>degree_z</tt> replaced by appropriate values.
	*/
	std::string
	get_name()
	const override;

	/**
	* @copydoc FiniteElement::compare_for_domination()
	*/
	FiniteElementDomination::Domination
	compare_for_domination(	const FiniteElement<dim, spacedim>&	fe_other,
							const unsigned int					codim)
	const
	override;

	/**
	* @copydoc FiniteElement::hp_vertex_dof_identities()
	*/
	std::vector<std::pair<unsigned int, unsigned int>>
	hp_vertex_dof_identities(const FiniteElement<dim, spacedim>& fe_other)
	const
	override;

	/**
	* @copydoc FiniteElement::hp_line_dof_identities()
	*/
	std::vector<std::pair<unsigned int, unsigned int>>
	hp_line_dof_identities(const FiniteElement<dim, spacedim>& fe_other)
	const
	override;

	/**
	* @copydoc FiniteElement::hp_quad_dof_identities()
	*/
	std::vector<std::pair<unsigned int, unsigned int>>
	hp_quad_dof_identities(	const FiniteElement<dim, spacedim>&	fe_other,
							const unsigned int					face_no = 0)
	const
	override;

};

DEAL_II_NAMESPACE_CLOSE

#endif /* INCLUDE_GALERKIN_TOOLS_FE_DGQ_ANISOTROPIC_H_ */
