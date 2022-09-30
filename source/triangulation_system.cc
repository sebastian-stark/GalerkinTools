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

#include <fstream>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <galerkin_tools/triangulation_system.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
tuple<const typename InterfaceCellDomainCells<spacedim>::InterfaceCell, const typename InterfaceCellDomainCells<spacedim>::DomainCell, const unsigned int, const typename InterfaceCellDomainCells<spacedim>::DomainCell, const unsigned int, const InterfaceRefinementCase, const unsigned int>
InterfaceCellDomainCells<spacedim>::convert_constructor_inputs(	const DomainCell&		domain_cell,
																const unsigned int		face,
																const InterfaceCell&	interface_cell,
																const InterfaceSide		interface_side)
{

	InterfaceRefinementCase refinement_case = InterfaceRefinementCase::equally_fine;
	unsigned int face_plus;
	unsigned int face_minus;
	unsigned int subface;

	//check that interface cell and underlying domain face have same centers
	Assert(	(interface_cell->center()).distance(domain_cell->face(face)->center()) < 1e-14,
			ExcMessage("Internal error: interface and domain face centers not at same location. Either this is a bug or you forgot to supply appropriate manifolds with the interface definition.!"));

	//check whether underlying domain face is at boundary
	if(domain_cell->face(face)->at_boundary())
	{
		//if at boundary: underlying domain cell must be the minus side
		Assert( interface_side == InterfaceSide::minus,
				ExcMessage("If an interface cell is at the boundary, the underlying domain cell must always be on the minus side of the interface!"));
		refinement_case = InterfaceRefinementCase::at_boundary;
	}

	//underlying domain cell on minus side
	if(interface_side == InterfaceSide::minus)
	{
		face_minus = face;
		//not at boundary
		if(refinement_case != InterfaceRefinementCase::at_boundary)
		{
			Assert( ( domain_cell->level() == domain_cell->neighbor(face)->level() ) || ( domain_cell->level() == (domain_cell->neighbor(face)->level() + 1) ),
					ExcMessage("Domain cell on minus side must be either equally refined than plus side or once finer if interface cell is defined based on minus side!"));
			if(domain_cell->level() == domain_cell->neighbor(face)->level())
			{
				refinement_case = InterfaceRefinementCase::equally_fine;
				face_plus = domain_cell->neighbor_of_neighbor(face);
				subface = 0;
			}
			else
			{
				refinement_case = InterfaceRefinementCase::minus_is_finer;
				const auto face_no = domain_cell->neighbor_of_coarser_neighbor(face);
				face_plus = face_no.first;
				subface = face_no.second;
			}
			return make_tuple(interface_cell, domain_cell, face_minus, domain_cell->neighbor(face), face_plus, refinement_case, subface);
		}
		else
		{
			subface = 0;
			face_plus = 0;
			return make_tuple(interface_cell, domain_cell, face_minus, domain_cell, face_plus, refinement_case, subface);
		}
	}
	//underlying domain cell on plus side
	else
	{
		face_plus = face;
		Assert( ( domain_cell->neighbor(face)->level() == domain_cell->level() ) || ( ( domain_cell->neighbor(face)->level()+1) == domain_cell->level() ),
				ExcMessage("Domain cell on plus side must be either equally refined than minus side or once finer if interface cell is defined based on plus side!"));
		if(domain_cell->neighbor(face)->level() == domain_cell->level())
		{
			refinement_case = InterfaceRefinementCase::equally_fine;
			face_minus = domain_cell->neighbor_of_neighbor(face);
			subface = 0;
		}
		else
		{
			refinement_case = InterfaceRefinementCase::plus_is_finer;
			const auto face_no = domain_cell->neighbor_of_coarser_neighbor(face);
			face_minus=face_no.first;
			subface=face_no.second;
		}
		return make_tuple(interface_cell, domain_cell->neighbor(face), face_minus, domain_cell, face_plus, refinement_case, subface);
	}

}


template<unsigned int spacedim>
InterfaceCellDomainCells<spacedim>::InterfaceCellDomainCells(	const DomainCell&		domain_cell,
																const unsigned int		face,
																const InterfaceCell&	interface_cell,
																const InterfaceSide		interface_side)
:
InterfaceCellDomainCells(convert_constructor_inputs(domain_cell, face, interface_cell, interface_side))
{
}

template<unsigned int spacedim>
InterfaceCellDomainCells<spacedim>::~InterfaceCellDomainCells()
{
	Assert(	n_subscriptions()==0,
			ExcMessage("You are about to destroy an InterfaceCellDomainCells object, which is currently in use! Make sure that all InterfaceCellDomainCells objects live at least as long as the objects using them!"));
}


template<unsigned int spacedim>
tuple<const types::material_id, const types::material_id, const types::material_id>
InterfaceCellDomainCells<spacedim>::get_material_ids()
const
{
	if(refinement_case == InterfaceRefinementCase::at_boundary)
		return make_tuple(interface_cell->material_id(), domain_cell_minus->material_id(), numbers::invalid_material_id);
	else
		return make_tuple(interface_cell->material_id(), domain_cell_minus->material_id(), domain_cell_plus->material_id());
}

template<unsigned int spacedim>
InterfaceCellDomainCells<spacedim>::InterfaceCellDomainCells(const tuple<const InterfaceCell, const DomainCell, const unsigned int, const DomainCell, const unsigned int, const InterfaceRefinementCase, const unsigned int> input)
:
interface_cell(get<0>(input)),
domain_cell_minus(get<1>(input)),
face_minus(get<2>(input)),
domain_cell_plus(get<3>(input)),
face_plus(get<4>(input)),
refinement_case(get<5>(input)),
subface(get<6>(input))
{
}

template<unsigned int spacedim>
TriangulationSystem<spacedim>::TriangulationSystem(	Triangulation<spacedim, spacedim>& 	tria_domain,
													const bool							fix_vertex_positions)
:
fix_vertex_positions(fix_vertex_positions),
tria_domain(&tria_domain)
{
	if(dynamic_cast<const dealii::parallel::distributed::Triangulation<spacedim, spacedim>*>(&tria_domain) != nullptr)
	{
		tria_listeners.push_back(tria_domain.signals.pre_distributed_refinement.connect(boost::bind(&TriangulationSystem<spacedim>::pre_refinement_domain, this)));
		tria_listeners.push_back(tria_domain.signals.post_distributed_refinement.connect(boost::bind(&TriangulationSystem<spacedim>::post_refinement_domain, this)));
	}
	else
	{
		tria_listeners.push_back(tria_domain.signals.pre_refinement.connect(boost::bind(&TriangulationSystem<spacedim>::pre_refinement_domain, this)));
		tria_listeners.push_back(tria_domain.signals.post_refinement.connect(boost::bind(&TriangulationSystem<spacedim>::post_refinement_domain, this)));
	}

	tria_interface = make_unique<Triangulation<spacedim-1, spacedim>>();

}

template<unsigned int spacedim>
TriangulationSystem<spacedim>::~TriangulationSystem()
{
	for(auto &connection : tria_listeners)
		connection.disconnect();
	tria_listeners.clear();

	Assert(	n_subscriptions() == 0,
			ExcMessage("You are about to destroy a TriangulationSystem, which is currently in use! Make sure that all TriangulationSystem objects live at least as long as the objects using them!"));
}

template<unsigned int spacedim>
const Triangulation<spacedim, spacedim>&
TriangulationSystem<spacedim>::get_triangulation_domain()
const
{
	return *tria_domain;
}

template<unsigned int spacedim>
const Triangulation<spacedim-1, spacedim>&
TriangulationSystem<spacedim>::get_triangulation_interface()
const
{
	Assert(	closed,
			ExcMessage("It does not make sense to ask for the interface triangulation before TriangulationSystem::close() has been called, because this would return an empty triangulation!"));
	return *tria_interface;
}

template<unsigned int spacedim>
Triangulation<spacedim, spacedim>&
TriangulationSystem<spacedim>::get_triangulation_domain()
{
	return *tria_domain;
}

template<unsigned int spacedim>
Triangulation<spacedim-1, spacedim>&
TriangulationSystem<spacedim>::get_triangulation_interface()
{
	Assert(	closed,
			ExcMessage("It does not make sense to ask for the interface triangulation before TriangulationSystem::close() has been called, because this would return an empty triangulation!"));
	return *tria_interface;
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::close()
{
	Assert( !closed,
			ExcMessage("A TriangulationSystem cannot be closed twice!"));
	generate_tria_interface_from_tria_domain();
	closed = true;
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::add_interface_cell(	const DomainCell&			cell,
													const unsigned int			face,
													const types::material_id	material_id)
{
	Assert(!closed, ExcMessage("Cells cannot be added after TriangulationSystem::close() has been called!"))
	coarse_domain_faces_material_ids[make_pair(cell, face)] = material_id;
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::add_interface_cells(vector< tuple<const DomainCell, const unsigned int, const types::material_id> > cells)
{
	for(const auto& cell : cells)
		add_interface_cell(get<0>(cell), get<1>(cell), get<2>(cell));
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::set_interface_manifold(	const types::manifold_id				manifold_id,
														const Manifold<spacedim-1, spacedim>&	manifold )
{
	Assert(!closed, ExcMessage("Manifolds cannot be added after TriangulationSystem::close() has been called!"))
	interface_manifolds.insert(make_pair(manifold_id, &manifold));
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::set_interface_manifolds(	const map<dealii::types::manifold_id, reference_wrapper<const dealii::Manifold<spacedim-1,spacedim>>> manifolds)
{
	for(const auto manifold : manifolds)
		set_interface_manifold(manifold.first, manifold.second);
}

template<unsigned int spacedim>
typename std::vector< InterfaceCellDomainCells<spacedim> >::iterator
TriangulationSystem<spacedim>::interface_begin_active()
{
	return active_interface_cell_domain_cells.begin();
}

template<unsigned int spacedim>
typename std::vector< InterfaceCellDomainCells<spacedim> >::iterator
TriangulationSystem<spacedim>::interface_end_active()
{
	return active_interface_cell_domain_cells.end();
}

template<unsigned int spacedim>
const typename std::vector< InterfaceCellDomainCells<spacedim> >&
TriangulationSystem<spacedim>::interface_active_iterators()
const
{
	return active_interface_cell_domain_cells;
}


template<unsigned int spacedim>typename std::vector< InterfaceCellDomainCells<spacedim> >::iterator
TriangulationSystem<spacedim>::interface_begin_coarse()
{
	return coarse_interface_cell_domain_cells.begin();
}

template<unsigned int spacedim>
typename std::vector< InterfaceCellDomainCells<spacedim> >::iterator
TriangulationSystem<spacedim>::interface_end_coarse()
{
	return coarse_interface_cell_domain_cells.end();
}

template<unsigned int spacedim>
const typename std::vector< InterfaceCellDomainCells<spacedim> >&
TriangulationSystem<spacedim>::interface_coarse_iterators()
const
{
	return coarse_interface_cell_domain_cells;
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::write_triangulations_vtk(string file_name_domain,
														string file_name_interface)
const
{
	if(file_name_domain != "")
	{
		ofstream ofstream_domain(file_name_domain);
		GridOut grid_out_domain;
		grid_out_domain.write_vtk(*tria_domain, ofstream_domain);
		ofstream_domain.close();
	}

	if(file_name_interface != "")
	{
		ofstream ofstream_interface(file_name_interface);
		GridOut grid_out_interface;
		grid_out_interface.write_vtk(*tria_interface, ofstream_interface);
		ofstream_interface.close();
	}
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::refine_global(const unsigned int times)
{
	Assert(closed, ExcMessage("Refinement is only possible after close() has been called!"))
	tria_domain->refine_global(times);
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::execute_coarsening_and_refinement()
{
	Assert(closed, ExcMessage("Refinement is only possible after close() has been called!"))
	tria_domain->execute_coarsening_and_refinement();
}

template<unsigned int spacedim>
bool
TriangulationSystem<spacedim>::check_active_interface_cell_domain_cells_consistency(const double tol)
const
{

	for(const auto& cell : active_interface_cell_domain_cells)
	{
		const auto& center_interface = cell.interface_cell->center();
		Point<spacedim> center_domain;
		if( (cell.refinement_case == InterfaceRefinementCase::minus_is_finer)
			||
			(cell.refinement_case == InterfaceRefinementCase::at_boundary)
			||
			(cell.refinement_case == InterfaceRefinementCase::equally_fine))
			center_domain = cell.domain_cell_minus->face(cell.face_minus)->center();
		else
			center_domain = cell.domain_cell_plus->face(cell.face_plus)->center();
		if(center_domain.distance(center_interface)>tol)
			return false;
	}
	return true;
}

template<unsigned int spacedim>
pair<const unsigned int, const unsigned int>
TriangulationSystem<spacedim>::get_this_proc_n_procs()
const
{
	return make_pair(0, 1);
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::generate_active_interface_cells_domain_cells(const bool no_assert)
{
	active_interface_cell_domain_cells.clear();
	for(const auto& coarse_interface_cell_domain_cells_n : coarse_interface_cell_domain_cells)
		generate_active_interface_cells_domain_cells_recursion(	coarse_interface_cell_domain_cells_n.domain_cell_minus,
																coarse_interface_cell_domain_cells_n.face_minus,
																coarse_interface_cell_domain_cells_n.interface_cell,
																no_assert);
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::generate_active_interface_cells_domain_cells_recursion(	const DomainCell&		domain_cell,
																						const unsigned int&		face,
																						const InterfaceCell&	interface_cell,
																						const bool				no_assert)
{
	if((domain_cell->refinement_case() != RefinementCase<spacedim>::isotropic_refinement) && (domain_cell->refinement_case() != RefinementCase<spacedim>::no_refinement))
	{
		for(unsigned int domain_cell_child_n = 0; domain_cell_child_n < domain_cell->n_children(); ++domain_cell_child_n)
		{
			for(unsigned int face_n = 0; face_n < GeometryInfo<spacedim>::faces_per_cell; ++face_n)
			{
				if(domain_cell->child(domain_cell_child_n)->face(face_n) == domain_cell->face(face))
				{
					generate_active_interface_cells_domain_cells_recursion(	domain_cell->child(domain_cell_child_n),
																			face_n,
																			interface_cell,
																			no_assert);
					return;
				}
			}
		}
	}

	if(interface_cell->has_children())
	{
		// interface cell has children
		// (1) at boundary:     domain cell must have children
		// (2) not at boundary: domain cell or its neighbor must have children
		if(domain_cell->has_children())
		{
			//domain cell has children->go one level deeper
			for(unsigned int interface_cell_child_n = 0; interface_cell_child_n < interface_cell->n_children(); ++interface_cell_child_n)
			{
				const unsigned int child_cell = GeometryInfo<spacedim>::child_cell_on_face(domain_cell->refinement_case(), face, interface_cell_child_n, domain_cell->face_orientation(face), domain_cell->face_flip(face), domain_cell->face_rotation(face));
				generate_active_interface_cells_domain_cells_recursion(	domain_cell->child(child_cell),
																		face,
																		interface_cell->child(interface_cell_child_n),
																		no_assert);
			}
		}
		else
		{
			// domain cell has no children
			//->domain cell cannot be at boundary and neighbor must have children;
			//  either the neighbor of the domain cell as well as the interface cell
			//  are refined exactly once more or the mesh is invalid
			Assert(	!(domain_cell->face(face)->at_boundary()),
					ExcMessage("Boundary and domain mesh inconsistent") );
			Assert(	domain_cell->neighbor(face)->has_children(),
					ExcMessage("Interface and domain mesh inconsistent!"));
			for(unsigned int interface_cell_child_n = 0; interface_cell_child_n < interface_cell->n_children(); ++interface_cell_child_n)
			{
				const DomainCell& domain_cell_neighbor_child_n = domain_cell->neighbor_child_on_subface(face, interface_cell_child_n);
				const unsigned int face_neighbor = domain_cell->neighbor_of_neighbor(face);
				Assert(	( domain_cell_neighbor_child_n->is_active() && interface_cell->child(interface_cell_child_n)->is_active() ) || no_assert,
						ExcMessage("Interface and domain mesh inconsistent!"));
				if(fix_vertex_positions)
				{
					for(unsigned int v = 0; v < GeometryInfo<spacedim>::vertices_per_face; ++ v)
					{
						if(domain_cell_neighbor_child_n->face(face_neighbor)->vertex(v).distance(interface_cell->child(interface_cell_child_n)->vertex(v)) > 1e-14)
						{
							interface_cell->child(interface_cell_child_n)->vertex(v) = domain_cell_neighbor_child_n->face(face_neighbor)->vertex(v);
						}
					}
				}

				Assert(	(domain_cell_neighbor_child_n->face(face_neighbor)->center()).distance(interface_cell->child(interface_cell_child_n)->center())<1e-14,
						ExcMessage("Internal error: interface and domain face centers not at same location. Either this is a bug or you forgot to supply appropriate manifolds with the interface definition!"));
				active_interface_cell_domain_cells.push_back(InterfaceCellDomainCells<spacedim>(domain_cell_neighbor_child_n, face_neighbor, interface_cell->child(interface_cell_child_n), InterfaceSide::plus));
			}
		}
	}
	// if the interface cell has no further children, the underlying domain cell may still be further refined in the case of anisotropic refinement
	else if((domain_cell->refinement_case() != RefinementCase<spacedim>::isotropic_refinement) && (domain_cell->refinement_case() != RefinementCase<spacedim>::no_refinement))
	{
		Assert(!(domain_cell->face(face)->has_children()), ExcMessage("Internal error! This indicates a bug!"));
		const unsigned int child_cell = GeometryInfo<spacedim>::child_cell_on_face(domain_cell->refinement_case(), face, 0, domain_cell->face_orientation(face), domain_cell->face_flip(face), domain_cell->face_rotation(face));
		generate_active_interface_cells_domain_cells_recursion(	domain_cell->child(child_cell),
																face,
																interface_cell,
																no_assert);
	}
	else
	{
		//the interface cell has no children
		//->neither the domain cell on the + nor the domain cell on the - cell can have children (otherwise the mesh is invalid)
		Assert(	( domain_cell->is_active() && interface_cell->is_active() ) || no_assert ,
				ExcMessage("Interface and domain mesh inconsistent!"));
		if(!(domain_cell->face(face)->at_boundary()))
			Assert(	domain_cell->neighbor(face)->is_active()  || no_assert,
					ExcMessage("Interface and domain mesh inconsistent!"));
		if(fix_vertex_positions)
		{
			for(unsigned int v = 0; v < GeometryInfo<spacedim>::vertices_per_face; ++ v)
			{
				if(domain_cell->face(face)->vertex(v).distance(interface_cell->vertex(v)) > 1e-14)
				{
					interface_cell->vertex(v) = domain_cell->face(face)->vertex(v);
				}
			}
		}
		Assert(	(domain_cell->face(face)->center()).distance(interface_cell->center())<1e-14,
				ExcMessage("Internal error: interface and domain face centers not at same location.  Either this is a bug or you forgot to supply appropriate manifolds with the interface definition."));
		active_interface_cell_domain_cells.push_back(InterfaceCellDomainCells<spacedim>(domain_cell, face, interface_cell, InterfaceSide::minus));
	}
}

template<unsigned int spacedim>
bool
TriangulationSystem<spacedim>::check_material_ids_recursion(const DomainCell& domain_cell)
const
{
	const types::material_id material_id = domain_cell->material_id();
	if(domain_cell->has_children())
		for(unsigned int child_n = 0; child_n < domain_cell->n_children(); ++child_n)
			if(domain_cell->child(child_n)->material_id() != material_id)
				return false;
			else
				return check_material_ids_recursion(domain_cell->child(child_n));
	return true;
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::pre_refinement_domain()
{
	Assert(closed, ExcMessage("Refinement is only possible after close() has been called!"))

	if(interface_active_iterators().size()==0)
		return;

	for(const auto& cell : interface_active_iterators())
	{
		switch (cell.refinement_case)
		{
			case InterfaceRefinementCase::minus_is_finer:
			{
				if( (cell.domain_cell_minus->coarsen_flag_set()) && (!cell.domain_cell_plus->refine_flag_set()))
					cell.interface_cell->set_coarsen_flag();
				else if ( (cell.domain_cell_minus->refine_flag_set()) && (cell.domain_cell_plus->refine_flag_set()) )
					cell.interface_cell->set_refine_flag();
				break;
			}
			case InterfaceRefinementCase::plus_is_finer:
			{
				if( (cell.domain_cell_plus->coarsen_flag_set()) && (!cell.domain_cell_minus->refine_flag_set()))
					cell.interface_cell->set_coarsen_flag();
				if ( (cell.domain_cell_plus->refine_flag_set()) && (cell.domain_cell_minus->refine_flag_set()) )
					cell.interface_cell->set_refine_flag();
				break;
			}
			case InterfaceRefinementCase::equally_fine:
			{
				// experimental: handle anisotropic refinement
				if( (cell.domain_cell_plus->coarsen_flag_set()) && (cell.domain_cell_minus->coarsen_flag_set()))
				{
					const auto& parent = cell.domain_cell_minus->parent();
					const auto face_refinement_case = GeometryInfo<spacedim>::face_refinement_case(	parent->refinement_case(),
																									cell.face_minus,
																									parent->face_orientation(cell.face_minus),
																									parent->face_flip(cell.face_minus),
																									parent->face_rotation(cell.face_minus));
					Assert( (face_refinement_case == RefinementCase<spacedim-1>::no_refinement) || (face_refinement_case == RefinementCase<spacedim-1>::isotropic_refinement), ExcMessage("No anisotropic refinement of interfaces or boundaries allowed!") );
					if(face_refinement_case == RefinementCase<spacedim-1>::isotropic_refinement)
						cell.interface_cell->set_coarsen_flag();
				}
				else if ( (cell.domain_cell_plus->refine_flag_set()) || (cell.domain_cell_minus->refine_flag_set()) )
				{
					const auto face_refinement_case_minus = GeometryInfo<spacedim>::face_refinement_case(	cell.domain_cell_minus->refine_flag_set(),
																											cell.face_minus,
																											cell.domain_cell_minus->face_orientation(cell.face_minus),
																											cell.domain_cell_minus->face_flip(cell.face_minus),
																											cell.domain_cell_minus->face_rotation(cell.face_minus));
					const auto face_refinement_case_plus = GeometryInfo<spacedim>::face_refinement_case(	cell.domain_cell_plus->refine_flag_set(),
																											cell.face_plus,
																											cell.domain_cell_plus->face_orientation(cell.face_plus),
																											cell.domain_cell_plus->face_flip(cell.face_plus),
																											cell.domain_cell_plus->face_rotation(cell.face_plus));
					Assert( (face_refinement_case_minus == RefinementCase<spacedim-1>::no_refinement) || (face_refinement_case_minus == RefinementCase<spacedim-1>::isotropic_refinement), ExcMessage("No anisotropic refinement of interfaces or boundaries allowed!") );
					Assert( (face_refinement_case_plus == RefinementCase<spacedim-1>::no_refinement) || (face_refinement_case_plus == RefinementCase<spacedim-1>::isotropic_refinement), ExcMessage("No anisotropic refinement of interfaces or boundaries allowed!") );
					if( (face_refinement_case_minus == RefinementCase<spacedim-1>::isotropic_refinement) || (face_refinement_case_plus == RefinementCase<spacedim-1>::isotropic_refinement) )
						cell.interface_cell->set_refine_flag();
				}
				break;
			}
			case InterfaceRefinementCase::at_boundary:
			{
				// experimental: handle anisotropic refinement
				if(cell.domain_cell_minus->coarsen_flag_set() )
				{
					const auto& parent = cell.domain_cell_minus->parent();
					const auto face_refinement_case = GeometryInfo<spacedim>::face_refinement_case(	parent->refinement_case(),
																									cell.face_minus,
																									parent->face_orientation(cell.face_minus),
																									parent->face_flip(cell.face_minus),
																									parent->face_rotation(cell.face_minus));
					Assert( (face_refinement_case == RefinementCase<spacedim-1>::no_refinement) || (face_refinement_case == RefinementCase<spacedim-1>::isotropic_refinement), ExcMessage("No anisotropic refinement of interfaces or boundaries allowed!") );
					if(face_refinement_case == RefinementCase<spacedim-1>::isotropic_refinement)
						cell.interface_cell->set_coarsen_flag();
				}
				else if( cell.domain_cell_minus->refine_flag_set() )
				{
					const auto face_refinement_case = GeometryInfo<spacedim>::face_refinement_case(	cell.domain_cell_minus->refine_flag_set(),
																									cell.face_minus,
																									cell.domain_cell_minus->face_orientation(cell.face_minus),
																									cell.domain_cell_minus->face_flip(cell.face_minus),
																									cell.domain_cell_minus->face_rotation(cell.face_minus));
					Assert( (face_refinement_case == RefinementCase<spacedim-1>::no_refinement) || (face_refinement_case == RefinementCase<spacedim-1>::isotropic_refinement), ExcMessage("No anisotropic refinement of interfaces or boundaries allowed!") );
					if(face_refinement_case == RefinementCase<spacedim-1>::isotropic_refinement)
						cell.interface_cell->set_refine_flag();
				}
				break;
			}
		}
	}
	pre_refinement();
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::post_refinement_domain()
{
	if(interface_active_iterators().size()>0)
	{
		tria_interface->execute_coarsening_and_refinement();
		generate_active_interface_cells_domain_cells();
	}
	post_refinement();
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::push_tria_listeners_to_end()
{
	for(auto &connection : tria_listeners)
		connection.disconnect();
	tria_listeners.clear();
	tria_listeners.push_back(tria_domain->signals.pre_refinement.connect(boost::bind(&TriangulationSystem<spacedim>::pre_refinement_domain, this)));
	tria_listeners.push_back(tria_domain->signals.post_refinement.connect(boost::bind(&TriangulationSystem<spacedim>::post_refinement_domain, this)));
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::generate_tria_interface_from_tria_domain()
{
	//reset interface triangulation in case it is not empty (this can happen in the parallel case
	//when the entire interface triangulation is rebuilt after refinement)
	tria_interface->clear();

	//attach manifolds
	for(const auto& interface_manifold : interface_manifolds)
	{
		//only attach manifolds here if they are not transfinite interpolation manifolds, otherwise later
		if(!dynamic_cast<const TransfiniteInterpolationManifold<spacedim-1, spacedim>*>(interface_manifold.second))
			tria_interface->set_manifold(interface_manifold.first, *(interface_manifold.second));
	}

	//this essentially creates the interface triangulation (approach similar as in GridGenerator::GridGenerator::extract_boundary_mesh() of deal.II)

	//map between vertex indices in volume and vertex indices on interface (vertex indices on
	//interface are additionally broken down into different interface material id's, because
	//there are duplicate nodes where interface portions with different material id's touch)
	map<
		unsigned int,
		map	<
			 types::material_id,
			 const unsigned int
			>
		>
	vertex_volume_vertex_interface;

	//vector of vertex locations
	vector<Point<spacedim>> vertices;

	//vector of cell properties (nodal connectivity, material id, boundary id, manifold id)
	vector<CellData<spacedim-1>> cell_data;

	//used in 3D only: edge properties (material id, boundary id, manifold id)
	SubCellData subcell_data;

	//vector for checking whether a vertex has already been visited (in order to avoid duplicate vertices)
	//if vertex i has already been visited during mesh assembly for an interface portion with
	//material id material_id, this material_id is included in touched[i]
	vector< set<types::material_id> > touched(tria_domain->n_vertices());

	//Check that material id's of refined mesh are consistent with material id's of coarse mesh
	for(const auto& cell : tria_domain->cell_iterators_on_level(0))
	{
		(void)cell;	//avoid unused variable compiler warning in release mode
		Assert(	check_material_ids_recursion(cell),
				ExcMessage("Material id's of refined domain mesh are not consistent with coarse mesh. This is not admissible!") );
	}

	//only if there is an interface mesh
	if(coarse_domain_faces_material_ids.size()>0)
	{
		//loop over interface cells
		for(const auto& cell_definition : coarse_domain_faces_material_ids)
		{

			const DomainCell cell=cell_definition.first.first;
			const unsigned int face_no=cell_definition.first.second;
			const types::material_id material_id=cell_definition.second;
			const auto face=cell->face(face_no);

			//assert that line flip face is not set (otherwise, create_triangulation will not work)
			for(unsigned int line_n=0; line_n<GeometryInfo<spacedim>::lines_per_face; line_n++)
				Assert(	face->line_orientation(line_n),
						ExcMessage("A line of a face making up the interface has non-standard orientation. This is currently not supported!") );

			//interface mesh must be defined by coarse domain mesh
			Assert(	cell->level() == 0,
					ExcMessage("Interfaces may only be defined by faces of cells of the coarse mesh!"));

			//assemble cell data
			CellData<spacedim-1> cell_data_cell;
			for(unsigned int vertex_n=0; vertex_n<GeometryInfo<spacedim-1>::vertices_per_cell; ++vertex_n)
			{
				const unsigned int vertex_index=face->vertex_index(vertex_n);
				if(touched[vertex_index].find(material_id)==touched[vertex_index].end())
				{
					vertices.push_back(face->vertex(vertex_n));
					touched[vertex_index].insert(material_id);
					vertex_volume_vertex_interface[vertex_index].insert(make_pair(material_id, vertices.size()-1));
				}
				cell_data_cell.vertices[vertex_n]=vertex_volume_vertex_interface[vertex_index][material_id];
			}
			cell_data_cell.material_id=material_id;
			cell_data_cell.manifold_id=face->manifold_id();

			//in 3D: add lines to subcell_data
			if(spacedim == 3)
			{
				for(unsigned int line_n=0; line_n<GeometryInfo<spacedim>::lines_per_face; ++line_n)
				{
					bool edge_found = false;
					const unsigned int vertex_0 = vertex_volume_vertex_interface[face->line(line_n)->vertex_index(0)][material_id];
					const unsigned int vertex_1 = vertex_volume_vertex_interface[face->line(line_n)->vertex_index(1)][material_id];
					for(const auto& boundary_line : subcell_data.boundary_lines)
					{
						const unsigned int vertex_0_in_subcell_data=boundary_line.vertices[0];
						const unsigned int vertex_1_in_subcell_data=boundary_line.vertices[1];
						if( ((vertex_0_in_subcell_data==vertex_0) && (vertex_1_in_subcell_data==vertex_1)) ||
							((vertex_0_in_subcell_data==vertex_1) && (vertex_1_in_subcell_data==vertex_0)) )
						{
							edge_found=true;
							break;
						}
					}
					if(edge_found==true)
						continue;

					CellData<1> line;
					line.vertices[0] = vertex_0;
					line.vertices[1] = vertex_1;
					line.boundary_id = numbers::internal_face_boundary_id;
					line.manifold_id = face->line(line_n)->manifold_id();
					subcell_data.boundary_lines.push_back(line);
				}
			}
			cell_data.push_back(cell_data_cell);
		}

		//create interface triangulation
		tria_interface->create_triangulation(vertices, cell_data, subcell_data);

		//initialize transfinite interpolation manifolds if these are present
		for(const auto& interface_manifold : interface_manifolds)
		{
			const auto transfinite_manifold = dynamic_cast<const TransfiniteInterpolationManifold<spacedim-1, spacedim>*>(interface_manifold.second);
			if(transfinite_manifold)
			{
				const_cast<TransfiniteInterpolationManifold<spacedim-1, spacedim>*>(transfinite_manifold)->initialize(*tria_interface);
				tria_interface->set_manifold(interface_manifold.first, *(interface_manifold.second));
			}
		}


		//Mapping between faces of the coarse domain mesh and corresponding cells of the coarse interface mesh
		//it is assumed here that create_triangulation does not change the order of the elements;
		//only do this if it has not been done during a previous call to close()
		if(coarse_interface_cell_domain_cells.size() == 0)
		{
			auto cell_definition = coarse_domain_faces_material_ids.begin();
			for(const auto& cell_interface : tria_interface->cell_iterators_on_level(0))
			{
				const DomainCell cell_domain=(cell_definition->first).first;
				const unsigned int face_domain=(cell_definition->first).second;
				coarse_interface_cell_domain_cells.push_back( InterfaceCellDomainCells<spacedim>(cell_domain, face_domain, cell_interface, InterfaceSide::minus) );
				++cell_definition;
			}
		}


		// Refinement and mapping between active faces of the domain mesh and corresponding active cells of the interface mesh
		generate_active_interface_cells_domain_cells(true);

		bool changed;
		do
		{
			changed = false;
			for(const auto & active_interface_cell_domain_cells_n : active_interface_cell_domain_cells)
			{
				const DomainCell& domain_cell_minus = active_interface_cell_domain_cells_n.domain_cell_minus;
				const unsigned int face_minus = active_interface_cell_domain_cells_n.face_minus;
				const unsigned int face_plus = active_interface_cell_domain_cells_n.face_plus;
				if(domain_cell_minus->face(face_minus)->at_boundary())
				{
					if( domain_cell_minus->face(face_minus)->has_children())
					{
						(active_interface_cell_domain_cells_n.interface_cell)->set_refine_flag();
						changed = true;
					}
				}
				else
				{
					const DomainCell& domain_cell_plus = active_interface_cell_domain_cells_n.domain_cell_plus;
					{
						// Checking the domain cells here seems superfluous. But if an artificial cell is involved in parallel, it seems to be possible that an artificial cell has no children, but its faces do have children.
						if( ( (domain_cell_minus->has_children()) && (domain_cell_minus->face(face_minus)->has_children()) ) ||	( domain_cell_plus->has_children() && (domain_cell_plus->face(face_plus)->has_children()) ) )
						{
							(active_interface_cell_domain_cells_n.interface_cell)->set_refine_flag();
							changed = true;
						}
					}
				}
			}
			if(changed)
			{
				tria_interface->execute_coarsening_and_refinement();
				generate_active_interface_cells_domain_cells(true);
			}
		}
		while(changed == true);
	}
}

//implementation of the parallel version of the TriangulationSystem
#ifdef DEAL_II_WITH_P4EST
#ifdef DEAL_II_WITH_MPI
namespace parallel
{

template<unsigned int dim, unsigned int spacedim>
Triangulation<dim, spacedim>::Triangulation(MPI_Comm mpi_communicator)
:
dealii::parallel::DistributedTriangulationBase<dim, spacedim>(mpi_communicator, dealii::Triangulation<spacedim-1, spacedim>::MeshSmoothing::none, false)
{
}

template<unsigned int dim, unsigned int spacedim>
void
Triangulation<dim, spacedim>::finalize_subdomain_assignment()
{
	this->update_number_cache();
}

template<unsigned int dim, unsigned int spacedim>
void
Triangulation<dim, spacedim>::update_cell_relations()
{
	Assert(false, ExcNotImplemented());
}

template<unsigned int dim, unsigned int spacedim>
void
Triangulation<dim, spacedim>::save(const string& /*filename*/)
const
{
	Assert(false, ExcNotImplemented());
}

template<unsigned int dim, unsigned int spacedim>
void
Triangulation<dim, spacedim>::load(const string& /*filename*/)
{
	Assert(false, ExcNotImplemented());
}

template<unsigned int dim, unsigned int spacedim>
void
Triangulation<dim, spacedim>::load(	const string&	/*filename*/,
					const bool	/*autopartition*/)
{
	Assert(false, ExcNotImplemented());
}

template<unsigned int dim, unsigned int spacedim>
bool
Triangulation<dim, spacedim>::is_multilevel_hierarchy_constructed()
const
{
	return false;
}

template<unsigned int spacedim>
TriangulationSystem<spacedim>::TriangulationSystem(	dealii::parallel::distributed::Triangulation<spacedim, spacedim>&	tria_domain,
													const bool															fix_vertex_positions)
:
dealii::GalerkinTools::TriangulationSystem<spacedim>(tria_domain, fix_vertex_positions),
mpi_communicator(tria_domain.get_communicator())
{
	//reset the interface triangulation to use a parallel::distributed::Triangulation
	this->tria_interface.reset(new dealii::GalerkinTools::parallel::Triangulation<spacedim-1, spacedim>(mpi_communicator));
}

template<unsigned int spacedim>
const dealii::parallel::distributed::Triangulation<spacedim, spacedim>&
TriangulationSystem<spacedim>::get_triangulation_domain()
const
{
	return dynamic_cast<const dealii::parallel::distributed::Triangulation<spacedim, spacedim>&>(*this->tria_domain);
}

template<unsigned int spacedim>
const dealii::GalerkinTools::parallel::Triangulation<spacedim-1, spacedim>&
TriangulationSystem<spacedim>::get_triangulation_interface()
const
{
	Assert(	this->closed,
			ExcMessage("It does not make sense to ask for the interface triangulation before TriangulationSystem::close() has been called, because this would return an empty triangulation!"));
	return dynamic_cast<const dealii::GalerkinTools::parallel::Triangulation<spacedim-1, spacedim>&>(*this->tria_interface);
}

template<unsigned int spacedim>
dealii::parallel::distributed::Triangulation<spacedim, spacedim>&
TriangulationSystem<spacedim>::get_triangulation_domain()
{
	return dynamic_cast<dealii::parallel::distributed::Triangulation<spacedim, spacedim>&>(*this->tria_domain);
}

template<unsigned int spacedim>
dealii::GalerkinTools::parallel::Triangulation<spacedim-1, spacedim>&
TriangulationSystem<spacedim>::get_triangulation_interface()
{
	Assert(	this->closed,
			ExcMessage("It does not make sense to ask for the interface triangulation before TriangulationSystem::close() has been called, because this would return an empty triangulation!"));
	return dynamic_cast<dealii::GalerkinTools::parallel::Triangulation<spacedim-1, spacedim>&>(*this->tria_interface);
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::close()
{
	//call the parent' close()
	dealii::GalerkinTools::TriangulationSystem<spacedim>::close();

	//update subdomain id's of interface
	update_interface_subdomain_ids();
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::write_meshes_per_processor_as_vtu(	string file_name_domain,
																	string file_name_interface)
const
{
	if(file_name_domain != "")
	{
		GridOut grid_out_domain;
		grid_out_domain.write_mesh_per_processor_as_vtu(*(this->tria_domain), file_name_domain);
	}

	if(file_name_interface != "")
	{
		GridOut grid_out_interface;
		grid_out_interface.write_mesh_per_processor_as_vtu(*(this->tria_interface), file_name_interface);
	}
}

template<unsigned int spacedim>
pair<const unsigned int, const unsigned int>
TriangulationSystem<spacedim>::get_this_proc_n_procs()
const
{
	MPI_Comm communicator = get_triangulation_domain().get_communicator();

	//the total number of processes
	int n_procs;
	MPI_Comm_size(communicator, &n_procs);

	// Get the rank of the current process
	int this_proc;
	MPI_Comm_rank(communicator, &this_proc);

	return make_pair((unsigned int)this_proc, (unsigned int)n_procs);
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::pre_refinement_domain()
{
	Assert(this->closed, ExcMessage("Refinement is only possible after close() has been called!"))
	this->pre_refinement();
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::post_refinement_domain()
{
	//completely reassemble the interface triangulation
	this->generate_tria_interface_from_tria_domain();

	//update subdomain id's of interface
	update_interface_subdomain_ids();

	this->post_refinement();
}

template<unsigned int spacedim>
void
TriangulationSystem<spacedim>::update_interface_subdomain_ids()
const
{
	const auto tria_interface_cast = dynamic_cast<dealii::GalerkinTools::parallel::Triangulation<spacedim-1, spacedim>*>(this->tria_interface.get());
	vector<bool> vertex_of_own_cell(tria_interface_cast->n_vertices(), false);
	const auto own_subdomain_id = tria_interface_cast->locally_owned_subdomain();

	for(const auto& interface_cell_domain_cells : this->interface_active_iterators())
	{
		const auto& interface_cell = (typename dealii::GalerkinTools::TriangulationSystem<spacedim>::ActiveInterfaceCell)interface_cell_domain_cells.interface_cell;
		const auto& domain_cell_minus = (typename dealii::GalerkinTools::TriangulationSystem<spacedim>::ActiveDomainCell)interface_cell_domain_cells.domain_cell_minus;
		interface_cell->set_subdomain_id(domain_cell_minus->subdomain_id());

		if(interface_cell->subdomain_id() == own_subdomain_id)
			for(unsigned int vertex = 0; vertex < GeometryInfo<spacedim-1>::vertices_per_cell; ++vertex)
				vertex_of_own_cell[interface_cell->vertex_index(vertex)] = true;
	}

	//finally make all ghost interface cells artificial, which actually don't share any vertex with an interface cell owned by the current processor
	//(such cells may exist at this point e.g. if the partitioning of domain cells between processors coincidentally goes along an interface and
	// the interface is not defined by domain cells owned by the current processor)
	for(const auto& interface_cell_domain_cells : this->interface_active_iterators())
	{
		const auto& interface_cell = (typename dealii::GalerkinTools::TriangulationSystem<spacedim>::ActiveInterfaceCell)interface_cell_domain_cells.interface_cell;
		//we only have to check cells which are ghost cells
		if( (interface_cell->subdomain_id() != own_subdomain_id) && (interface_cell->subdomain_id() != numbers::artificial_subdomain_id))
		{
			bool is_really_ghost = false;
			for(unsigned int vertex = 0; vertex < GeometryInfo<spacedim-1>::vertices_per_cell; ++vertex)
			{
				if(vertex_of_own_cell[interface_cell->vertex_index(vertex)])
				{
					is_really_ghost = true;
					break;
				}
			}
			if(!is_really_ghost)
				interface_cell->set_subdomain_id(numbers::artificial_subdomain_id);
		}
	}

	Assert( tria_interface_cast != nullptr, ExcMessage("Internal error!"));
	tria_interface_cast->finalize_subdomain_assignment();
}

}
#endif // DEAL_II_WITH_MPI
#endif // DEAL_II_WITH_P4EST

template class TriangulationSystem<2>;
template class TriangulationSystem<3>;

#ifdef DEAL_II_WITH_P4EST
#ifdef DEAL_II_WITH_MPI
	template class parallel::TriangulationSystem<2>;
	template class parallel::TriangulationSystem<3>;
	template class parallel::Triangulation<1,2>;
	template class parallel::Triangulation<2,3>;
#endif // DEAL_II_WITH_MPI
#endif // DEAL_II_WITH_P4EST

template struct InterfaceCellDomainCells<2>;
template struct InterfaceCellDomainCells<3>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
