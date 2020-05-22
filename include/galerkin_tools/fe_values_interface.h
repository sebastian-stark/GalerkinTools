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

#ifndef GALERKINTOOLS_FEVALUESINTERFACE_H_
#define GALERKINTOOLS_FEVALUESINTERFACE_H_

#include <deal.II/fe/fe_values.h>

#include <galerkin_tools/config.h>
#include <galerkin_tools/dof_handler_system.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * This class collects the different FEValues, FEFaceValues, and FESubfaceValues objects needed for
 * cells at the interface. The situation there is a bit tricky, because, the domain cells adjacent
 * to the interface cell may be refined differently, thus requiring either an FEFaceValues or FESubfaceValues
 * object. This class makes sure that always the correct objects are used.
 */
template<unsigned int spacedim>
class FEValuesInterface
{

private:

	/**
	 * The refinement case at the interface cell with which the FEValues objects are
	 * initialized currently.
	 */
	InterfaceRefinementCase
	initialized_interface_refinement_case;

	/**
	 * The FEFaceValues or FESubfaceValues object currently in use for the domain cell on the minus side
	 * of the interface
	 */
	const FEFaceValuesBase<spacedim, spacedim>*
	initialized_fe_values_domain_minus = nullptr;

	/**
	 * The FEFaceValues or FESubfaceValues object currently in use for the domain cell on the minus side
	 * of the interface
	 */
	const FEFaceValuesBase<spacedim, spacedim>*
	initialized_fe_values_domain_plus = nullptr;

	/**
	 * The FEValues object for the interface cell
	 */
	FEValues<spacedim-1, spacedim>
	fe_values_interface;

	/**
	 * The FEFaceValues object for the domain cell face on the minus side
	 * (this is used in case the minus side is refined to the same degree
	 *  as the interface cell)
	 */
	FEFaceValues<spacedim, spacedim>
	fe_face_values_domain_minus;

	/**
	 * The FEFaceValues object for the domain cell face on the plus side
	 * (this is used in case the plus side is refined to the same degree
	 *  as the interface cell)
	 */
	FEFaceValues<spacedim, spacedim>
	fe_face_values_domain_plus;

	/**
	 * The FESubfaceValues object for the domain cell face on the minus side
	 * (this is used in case the minus side is coarser than the interface cell)
	 */
	FESubfaceValues<spacedim, spacedim>
	fe_subface_values_domain_minus;

	/**
	 * The FESubfaceValues object for the domain cell face on the plus side
	 * (this is used in case the plus side is coarser than the interface cell)
	 */
	FESubfaceValues<spacedim, spacedim>
	fe_subface_values_domain_plus;

public:

	/**
	 * The number of quadrature points
	 */
	const unsigned int
	n_quadrature_points;

	/**
	 * The constructor of the class
	 *
	 * @param[in]	mapping_interface			The mapping to be used for the interface cells
	 *
	 * @param[in]	mapping_domain				The mapping to be used for the domain cells
	 *
	 * @param[in]	fe_interface				The finite element for the interface cells
	 *
	 * @param[in]	fe_domain_minus				The finite element for the domain cells on the minus side of the interface
	 *
	 * @param[in]	fe_domain_plus				The finite element for the domain cells on the plus side of the interface
	 *
	 * @param[in]	quadrature					The Quadrature rule to be used
	 *
	 * @param[in]	update_flags_interface		::UpdateFlags for the interface cells
	 *
	 * @param[in]	update_flags_domain_minus	::UpdateFlags for the domain cells on the minus side of interface
	 *
	 * @param[in]	update_flags_domain_plus	::UpdateFlags for the domain cells on the plus side of interface
	 */
	FEValuesInterface(	const Mapping<spacedim-1, spacedim>&		mapping_interface,
						const Mapping<spacedim, spacedim>&			mapping_domain,
						const FiniteElement<spacedim-1, spacedim >&	fe_interface,
						const FiniteElement<spacedim, spacedim >&	fe_domain_minus,
						const FiniteElement<spacedim, spacedim >&	fe_domain_plus,
						const Quadrature<spacedim-1>&				quadrature,
						const UpdateFlags							update_flags_interface,
						const UpdateFlags							update_flags_domain_minus,
						const UpdateFlags							update_flags_domain_plus);

	/**
	 * Initialize the object with an InterfaceCellDomainCellsDoF. This initializes the required
	 * FEValues, FEFaceValues, and FESubfaceValues objects.
	 *
	 * @param[in]	interface_cell_domain_cells	The InterfaceCellDomainCellsDoF
	 */
	void
	reinit(const InterfaceCellDomainCellsDoF<spacedim>& interface_cell_domain_cells);

	/**
	 * Initialize the object with an InterfaceCellDomainCells. This initializes the required
	 * FEValues, FEFaceValues, and FESubfaceValues objects.
	 * However, note that the @p interface_cell_domain_cells does not contain any dof information.
	 * Therefore, the FEValues, FEFaceValues, and FESubfaceValues don't have the information to
	 * compute any values which require knowledge of the dof values.
	 *
	 * @param[in]	interface_cell_domain_cells	The InterfaceCellDomainCells
	 */
	void
	reinit(const InterfaceCellDomainCells<spacedim>& interface_cell_domain_cells);

	/**
	 * @return	Quadrature scheme in use
	 */
	const Quadrature<spacedim-1>&
	get_quadrature()
	const;

	/**
	 * @return	The presently initialized FEValues object for interface cells
	 */
	const
	FEValues<spacedim-1,spacedim>&
	get_fe_values_interface()
	const;

	/**
	 * @param[in]	interface_side	Either ::InterfaceSide::@p minus or ::InterfaceSide::@p plus depending on
	 * 								the side of the interface for which the FEFaceValues or FESubfaceValues
	 * 								is required
	 *
	 * @return						The presently initialized FEFaceValues or FESubfaceValues object for domain
	 * 								cells on the minus or plus side of the interface
	 */
	const
	FEFaceValuesBase<spacedim,spacedim>&
	get_fe_values_domain(const InterfaceSide& interface_side)
	const;

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_FEVALUESINTERFACE_H_ */
