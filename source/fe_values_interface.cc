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

#include <galerkin_tools/fe_values_interface.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
FEValuesInterface<spacedim>::FEValuesInterface(	const Mapping<spacedim-1, spacedim>&		mapping_interface,
												const Mapping<spacedim, spacedim>&			mapping_domain,
												const FiniteElement<spacedim-1, spacedim >&	fe_interface,
												const FiniteElement<spacedim, spacedim >&	fe_domain_minus,
												const FiniteElement<spacedim, spacedim >&	fe_domain_plus,
												const Quadrature<spacedim-1>&				quadrature,
												const UpdateFlags							update_flags_interface,
												const UpdateFlags							update_flags_domain_minus,
												const UpdateFlags							update_flags_domain_plus)
:
initialized_interface_refinement_case(InterfaceRefinementCase::equally_fine),
fe_values_interface(mapping_interface, fe_interface, quadrature, update_flags_interface),
fe_face_values_domain_minus(mapping_domain, fe_domain_minus, quadrature, update_flags_domain_minus),
fe_face_values_domain_plus(mapping_domain, fe_domain_plus, quadrature, update_flags_domain_plus),
fe_subface_values_domain_minus(mapping_domain, fe_domain_minus, quadrature, update_flags_domain_minus),
fe_subface_values_domain_plus(mapping_domain, fe_domain_plus, quadrature, update_flags_domain_plus),
n_quadrature_points(quadrature.size())
{
}

template<unsigned int spacedim>
void
FEValuesInterface<spacedim>::reinit(const InterfaceCellDomainCellsDoF<spacedim>& interface_cell_domain_cells)
{

	fe_values_interface.reinit(interface_cell_domain_cells.interface_cell);

	initialized_interface_refinement_case = interface_cell_domain_cells.refinement_case;

	switch(initialized_interface_refinement_case)
	{
		case InterfaceRefinementCase::at_boundary:
			fe_face_values_domain_minus.reinit(interface_cell_domain_cells.domain_cell_minus, interface_cell_domain_cells.face_minus);
			initialized_fe_values_domain_minus = &fe_face_values_domain_minus;
			break;
		case InterfaceRefinementCase::equally_fine:
			fe_face_values_domain_minus.reinit(interface_cell_domain_cells.domain_cell_minus, interface_cell_domain_cells.face_minus);
			fe_face_values_domain_plus.reinit(interface_cell_domain_cells.domain_cell_plus, interface_cell_domain_cells.face_plus);
			initialized_fe_values_domain_minus = &fe_face_values_domain_minus;
			initialized_fe_values_domain_plus = &fe_face_values_domain_plus;
			break;
		case InterfaceRefinementCase::minus_is_finer:
			fe_face_values_domain_minus.reinit(interface_cell_domain_cells.domain_cell_minus, interface_cell_domain_cells.face_minus);
			fe_subface_values_domain_plus.reinit(interface_cell_domain_cells.domain_cell_plus, interface_cell_domain_cells.face_plus, interface_cell_domain_cells.subface);
			initialized_fe_values_domain_minus = &fe_face_values_domain_minus;
			initialized_fe_values_domain_plus = &fe_subface_values_domain_plus;
			break;
		case InterfaceRefinementCase::plus_is_finer:
			fe_subface_values_domain_minus.reinit(interface_cell_domain_cells.domain_cell_minus, interface_cell_domain_cells.face_minus, interface_cell_domain_cells.subface);
			fe_face_values_domain_plus.reinit(interface_cell_domain_cells.domain_cell_plus, interface_cell_domain_cells.face_plus);
			initialized_fe_values_domain_minus = &fe_subface_values_domain_minus;
			initialized_fe_values_domain_plus = &fe_face_values_domain_plus;

			break;
	}
}

template<unsigned int spacedim>
void
FEValuesInterface<spacedim>::reinit(const InterfaceCellDomainCells<spacedim>& interface_cell_domain_cells)
{

	fe_values_interface.reinit(interface_cell_domain_cells.interface_cell);

	switch(interface_cell_domain_cells.refinement_case)
	{
		case InterfaceRefinementCase::at_boundary:
			fe_face_values_domain_minus.reinit(interface_cell_domain_cells.domain_cell_minus, interface_cell_domain_cells.face_minus);
			break;
		case InterfaceRefinementCase::equally_fine:
			fe_face_values_domain_minus.reinit(interface_cell_domain_cells.domain_cell_minus, interface_cell_domain_cells.face_minus);
			fe_face_values_domain_plus.reinit(interface_cell_domain_cells.domain_cell_plus, interface_cell_domain_cells.face_plus);
			break;
		case InterfaceRefinementCase::minus_is_finer:
			fe_face_values_domain_minus.reinit(interface_cell_domain_cells.domain_cell_minus, interface_cell_domain_cells.face_minus);
			fe_subface_values_domain_plus.reinit(interface_cell_domain_cells.domain_cell_plus, interface_cell_domain_cells.face_plus, interface_cell_domain_cells.subface);
			break;
		case InterfaceRefinementCase::plus_is_finer:
			fe_subface_values_domain_minus.reinit(interface_cell_domain_cells.domain_cell_minus, interface_cell_domain_cells.face_minus, interface_cell_domain_cells.subface);
			fe_face_values_domain_plus.reinit(interface_cell_domain_cells.domain_cell_plus, interface_cell_domain_cells.face_plus);
			break;
	}
}

template<unsigned int spacedim>
const Quadrature<spacedim-1>&
FEValuesInterface<spacedim>::get_quadrature()
const
{
	return fe_values_interface.get_quadrature();
}

template<unsigned int spacedim>
const
FEValues<spacedim-1,spacedim>&
FEValuesInterface<spacedim>::get_fe_values_interface()
const
{
	return fe_values_interface;
}

template<unsigned int spacedim>
const
FEFaceValuesBase<spacedim,spacedim>&
FEValuesInterface<spacedim>::get_fe_values_domain(const InterfaceSide& interface_side)
const
{
	if(interface_side == InterfaceSide::minus)
	{
		Assert(initialized_fe_values_domain_minus != nullptr, ExcMessage("You are asking for an fe values object which has not been initialized to a cell!"));
		return *initialized_fe_values_domain_minus;
	}
	else
	{
		//it's ok to return a nullpointer for the plus side if we are at the boundary
		Assert( (initialized_fe_values_domain_plus != nullptr) || (initialized_interface_refinement_case == InterfaceRefinementCase::at_boundary),
				ExcMessage("You are asking for an fe values object which has not been initialized to a cell!"));
		return *initialized_fe_values_domain_plus;
	}
}

template class FEValuesInterface<2>;
template class FEValuesInterface<3>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
