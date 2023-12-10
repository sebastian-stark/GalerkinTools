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

#include <galerkin_tools/fe_dgq_anisotropic.h>

#include <vector>
#include <algorithm>

using namespace dealii;

DEAL_II_NAMESPACE_OPEN

namespace
{
	/**
	 * Helper function to set up the dpo vector of FE_DGQAnisotropic for the 3d case.
	 */
	std::vector<unsigned int>
	get_dpo_vector_fe_dgq_anisotropic(	const Quadrature<1>& points_x,
										const Quadrature<1>& points_y,
										const Quadrature<1>& points_z)
	{
		return {0, 0, 0, points_x.size() * points_y.size() * points_z.size()};
	}

	/**
	 * Helper function to set up the dpo vector of FE_DGQAnisotropic for the 2d case.
	 */
	std::vector<unsigned int>
	get_dpo_vector_fe_dgq_anisotropic(	const Quadrature<1>& points_x,
										const Quadrature<1>& points_y)
	{
	  return {0, 0, 0, points_x.size() * points_y.size()};
	}
} // namespace


template <int dim, int spacedim>
FE_DGQAnisotropic<dim, spacedim>::FE_DGQAnisotropic(const Quadrature<1>& points_x,
													const Quadrature<1>& points_y,
													const Quadrature<1>& points_z)
:
dealii::FE_Poly<dim, spacedim>(	dealii::AnisotropicPolynomials<dim>({Polynomials::generate_complete_Lagrange_basis(points_x.get_points()), Polynomials::generate_complete_Lagrange_basis(points_y.get_points()) , Polynomials::generate_complete_Lagrange_basis(points_z.get_points())}),
								FiniteElementData<dim>(get_dpo_vector_fe_dgq_anisotropic(points_x, points_y, points_z),
								1,
								std::max({points_x.size()-1, points_y.size()-1, points_z.size()-1}),
								FiniteElementData<dim>::Conformity::L2),
								std::vector<bool>(1, false),
								std::vector<ComponentMask>(1, std::vector<bool>(1, true))),
degree_x(points_x.size()),
degree_y(points_y.size()),
degree_z(points_z.size())
{
	Assert(dim == 3, ExcMessage("This constructor can only be used for dim == 3"));

	const QAnisotropic<dim> quadrature(points_x, points_y, points_z);
	this->unit_support_points = quadrature.get_points();
}


template <int dim, int spacedim>
FE_DGQAnisotropic<dim, spacedim>::FE_DGQAnisotropic(const Quadrature<1>& points_x,
													const Quadrature<1>& points_y)
:
dealii::FE_Poly<dim, spacedim>(	dealii::AnisotropicPolynomials<dim>({Polynomials::generate_complete_Lagrange_basis(points_x.get_points()), Polynomials::generate_complete_Lagrange_basis(points_y.get_points())}),
								FiniteElementData<dim>(get_dpo_vector_fe_dgq_anisotropic(points_x, points_y),
								1,
								std::max({points_x.size()-1, points_y.size()-1}),
								FiniteElementData<dim>::Conformity::L2),
								std::vector<bool>(1, false),
								std::vector<ComponentMask>(1, std::vector<bool>(1, true))),
degree_x(points_x.size()),
degree_y(points_y.size()),
degree_z(0)
{
	Assert(dim == 3, ExcMessage("This constructor can only be used for dim == 2"));

	const QAnisotropic<dim> quadrature(points_x, points_y);
	this->unit_support_points = quadrature.get_points();
}


template <int dim, int spacedim>
void
FE_DGQAnisotropic<dim, spacedim>::convert_generalized_support_point_values_to_dof_values(	const std::vector<Vector<double>>& support_point_values,
																							std::vector<double> &              nodal_values)
const
{
  AssertDimension(support_point_values.size(),
                  this->get_unit_support_points().size());
  AssertDimension(support_point_values.size(), nodal_values.size());
  AssertDimension(this->dofs_per_cell, nodal_values.size());

  for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
    {
      AssertDimension(support_point_values[i].size(), 1);

      nodal_values[i] = support_point_values[i](0);
    }
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_DGQAnisotropic<dim, spacedim>::clone() const
{
  return std::make_unique<FE_DGQAnisotropic<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::string
FE_DGQAnisotropic<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_DGQAnisopropic<" << dim << ">(" << this->degree_x << ", " << this->degree_y << ", " << this->degree_z << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_DGQAnisotropic<dim, spacedim>::compare_for_domination(	const FiniteElement<dim, spacedim>& /*fe_other*/,
															const unsigned int                  /*codim*/)
const
{

	Assert(false, ExcNotImplemented());
	// not implemented yet, maybe possible to copy over from FE_DGQ?
	return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGQAnisotropic<dim, spacedim>::hp_vertex_dof_identities(const FiniteElement<dim, spacedim>& /*fe_other*/)
const
{
	Assert(false, ExcNotImplemented());
	return {{0, 0}};
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGQAnisotropic<dim, spacedim>::hp_line_dof_identities(const FiniteElement<dim, spacedim>& /*fe_other*/)
const
{
	Assert(false, ExcNotImplemented());
	return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGQAnisotropic<dim, spacedim>::hp_quad_dof_identities(	const FiniteElement<dim, spacedim>& /*fe_other*/,
															const unsigned int                  /*face_no*/)
const
{
	Assert(false, ExcNotImplemented());
	return std::vector<std::pair<unsigned int, unsigned int>>();
}


// explicit instantiations
template class FE_DGQAnisotropic<2,2>;
template class FE_DGQAnisotropic<2,3>;
template class FE_DGQAnisotropic<3,3>;

DEAL_II_NAMESPACE_CLOSE
