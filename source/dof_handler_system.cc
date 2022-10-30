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

#include <galerkin_tools/dof_handler_system.h>
#include <galerkin_tools/tools.h>

#include <deal.II/dofs/dof_tools.h>

#include <unordered_map>
#include <unordered_set>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
InterfaceCellDoFIterator<spacedim>::InterfaceCellDoFIterator(	const TriaIterator<CellAccessor<spacedim-1, spacedim>>&	interface_cell,
																const DoFHandlerSystem<spacedim>& 						dof_handler_system)
:
DoFHandler<spacedim-1, spacedim>::active_cell_iterator(*interface_cell, &dof_handler_system.get_dof_handler_interface()),
dof_handler_system(dof_handler_system)
{
}

template<unsigned int spacedim>
InterfaceCellDoFIterator<spacedim>::InterfaceCellDoFIterator(const DoFHandlerSystem<spacedim>& dof_handler_system)
:
DoFHandler<spacedim-1, spacedim>::active_cell_iterator(dof_handler_system.get_dof_handler_interface().end()),
dof_handler_system(dof_handler_system)
{
}

template<unsigned int spacedim>
void
InterfaceCellDoFIterator<spacedim>::get_dof_indices(vector<types::global_dof_index >& dof_indices)
const
{
	this->accessor.get_dof_indices(dof_indices);
	const unsigned int n_dofs_domain = dof_handler_system.n_dofs_domain();
	for(auto& dof_indices_n : dof_indices)
		//shift global interface dof indices of the dof handler to the really global ones
		dof_indices_n += n_dofs_domain;

	//renumber
	if(dof_handler_system.dof_renumbering != nullptr)
		dof_handler_system.dof_renumbering->convert_dof_indices(dof_indices);
}

template<unsigned int spacedim>
DomainCellDoFIterator<spacedim>::DomainCellDoFIterator(	const TriaIterator<CellAccessor<spacedim, spacedim>>& 	domain_cell,
														const DoFHandlerSystem<spacedim>& 						dof_handler_system)
:
DoFHandler<spacedim, spacedim>::active_cell_iterator(*domain_cell, &dof_handler_system.get_dof_handler_domain()),
dof_handler_system(dof_handler_system)
{
}

template<unsigned int spacedim>
DomainCellDoFIterator<spacedim>::DomainCellDoFIterator(const DoFHandlerSystem<spacedim>& dof_handler_system)
:
DoFHandler<spacedim, spacedim>::active_cell_iterator(dof_handler_system.get_dof_handler_domain().end()),
dof_handler_system(dof_handler_system)
{
}

template<unsigned int spacedim>
void
DomainCellDoFIterator<spacedim>::get_dof_indices(vector<types::global_dof_index >& dof_indices)
const
{
	this->accessor.get_dof_indices(dof_indices);
	//renumber
	if(dof_handler_system.dof_renumbering != nullptr)
		dof_handler_system.dof_renumbering->convert_dof_indices(dof_indices);
}

template<unsigned int spacedim>
InterfaceCellDomainCellsDoF<spacedim>::InterfaceCellDomainCellsDoF(	const InterfaceCellDomainCells<spacedim>& 	interface_cell_domain_cell,
																	const DoFHandlerSystem<spacedim>& 			dof_handler_system)
:
interface_cell(InterfaceCellDoFIterator<spacedim>(interface_cell_domain_cell.interface_cell, dof_handler_system)),
domain_cell_minus(DomainCellDoFIterator<spacedim>(interface_cell_domain_cell.domain_cell_minus, dof_handler_system)),
face_minus(interface_cell_domain_cell.face_minus),
domain_cell_plus(DomainCellDoFIterator<spacedim>(interface_cell_domain_cell.domain_cell_plus, dof_handler_system)),
face_plus(interface_cell_domain_cell.face_plus),
refinement_case(interface_cell_domain_cell.refinement_case),
subface(interface_cell_domain_cell.subface)
{
}

template<unsigned int spacedim>
InterfaceCellDomainCellsDoF<spacedim>::~InterfaceCellDomainCellsDoF()
{
	Assert(n_subscriptions() == 0, ExcMessage("You are about to destroy an InterfaceCellDomainCellDoF, which is currently in use! Make sure that all InterfaceCellDomainCellDoF objects live at least as long as the objects using them!"));
}

template<unsigned int spacedim>
tuple<const types::material_id, const types::material_id, const types::material_id>
InterfaceCellDomainCellsDoF<spacedim>::get_material_ids()
const
{
	if(refinement_case==InterfaceRefinementCase::at_boundary)
			return make_tuple(interface_cell->material_id(), domain_cell_minus->material_id(), numbers::invalid_material_id);
		else
	return make_tuple(interface_cell->material_id(), domain_cell_minus->material_id(), domain_cell_plus->material_id());
}

template<unsigned int spacedim>
void
InterfaceCellDomainCellsDoF<spacedim>::get_dof_indices_local_global_interface(	vector<unsigned int>& dof_indices_local_global,
																				vector<unsigned int>& dof_indices_local_global_minus,
																				vector<unsigned int>& dof_indices_local_global_plus)
const
{
	//get the mapping between local and global dof indices on interface
	get_dof_indices_local_global_interface(dof_indices_local_global);

	//get the mappings between local and global dof indices on + and - side, respectively
	dof_indices_local_global_minus.resize(domain_cell_minus->get_fe().dofs_per_cell);
	domain_cell_minus.get_dof_indices(dof_indices_local_global_minus);
	if(refinement_case != InterfaceRefinementCase::at_boundary)
	{
		dof_indices_local_global_plus.resize(domain_cell_plus->get_fe().dofs_per_cell);
		domain_cell_plus.get_dof_indices(dof_indices_local_global_plus);
	}
	else
		dof_indices_local_global_plus.resize(0);
}

template<unsigned int spacedim>
void
InterfaceCellDomainCellsDoF<spacedim>::get_dof_indices_local_global_interface(vector<unsigned int>& dof_indices_local_global )
const
{
	//get the mapping between local and global dof indices on interface
	dof_indices_local_global.resize(interface_cell->get_fe().dofs_per_cell);
	interface_cell.get_dof_indices(dof_indices_local_global);
}

template<unsigned int spacedim>
DoFHandlerSystem<spacedim>::DoFHandlerSystem(const TriangulationSystem<spacedim>& tria_system):
tria_system(&tria_system)
{
	//create the dof handlers
	dof_handler_domain=shared_ptr<DoFHandler<spacedim, spacedim>>(new DoFHandler<spacedim, spacedim>(tria_system.get_triangulation_domain()));
	dof_handler_interface=shared_ptr<DoFHandler<spacedim-1, spacedim>>(new DoFHandler<spacedim-1, spacedim>(tria_system.get_triangulation_interface()));

	//generate active_interface_cell_domain_cells
	update_interface_domain_relation();

	//connect to triangulation system in order to make sure that after refinement of the triangulation also update_interface_domain_relation()
	//is called (in order to update active_interface_cell_domain_cells)
	tria_listeners.push_back(tria_system.post_refinement.connect(0, boost::bind(&DoFHandlerSystem<spacedim>::update_interface_domain_relation, this)));
}

template<unsigned int spacedim>
DoFHandlerSystem<spacedim>::~DoFHandlerSystem()
{
	for(auto &connection : tria_listeners)
		connection.disconnect();
	tria_listeners.clear();
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::distribute_dofs(const hp::FECollection<spacedim, spacedim>&		fe_collection_domain,
											const hp::FECollection<spacedim-1, spacedim>&	fe_collection_interface,
											const unsigned int								n_additional_dofs)
{
	const auto this_proc_n_procs = tria_system->get_this_proc_n_procs();
	dof_handler_domain->distribute_dofs(fe_collection_domain);

	dof_handler_interface->distribute_dofs(fe_collection_interface);
	dof_indices_C.resize(n_additional_dofs);
	locally_owned_dof_indices_C.reserve(n_additional_dofs);
	const unsigned int offset = n_dofs_domain() + n_dofs_interface();
	for(unsigned int m = 0; m < n_additional_dofs; ++m)
	{
		//warning:	the library generally assumes that dof_indices_C is contiguous and that the C's come after the domain and interface related dofs
		//			->don't change this
		dof_indices_C[m] = offset + m;
		//currently all C's are owned by the last processor
		if(this_proc_n_procs.first == (this_proc_n_procs.second - 1))
			locally_owned_dof_indices_C.push_back(offset + m);
	}
	set_locally_owned_dofs_standard_numbering();
	set_locally_relevant_dofs_standard_numbering();
	set_locally_owned_dofs();
	set_locally_relevant_dofs();

	n_dofs_per_processor.clear();
	n_dofs_per_processor.resize(this_proc_n_procs.second, 0);
	n_dofs_per_processor[this_proc_n_procs.first] = locally_owned_dofs.n_elements();

	if(this_proc_n_procs.second > 1)
	{
#ifdef DEAL_II_WITH_MPI
		const auto tria_domain_ptr = dynamic_cast<const dealii::parallel::distributed::Triangulation<spacedim, spacedim>*>(&(dof_handler_domain->get_triangulation()));
		Assert(tria_domain_ptr != nullptr, ExcMessage("Internal error!"));

		unsigned int send_value = n_dofs_per_processor[this_proc_n_procs.first];
		int ierr =  MPI_Allgather(&send_value, 1, MPI_UNSIGNED, n_dofs_per_processor.data(), 1, MPI_UNSIGNED, tria_domain_ptr->get_communicator());
		AssertThrowMPI(ierr);
#else
		Assert(this_proc_n_procs.second == 1, ExcMessage("Internal error: deal.II not configured with MPI, but computation run with more than one processor"));
#endif //DEAL_II_WITH_MPI
	}
}

template<unsigned int spacedim>
typename std::vector< InterfaceCellDomainCellsDoF<spacedim> >::iterator
DoFHandlerSystem<spacedim>::interface_begin_active()
{
	return active_interface_cell_domain_cells.begin();
}

template<unsigned int spacedim>
typename std::vector< InterfaceCellDomainCellsDoF<spacedim> >::iterator
DoFHandlerSystem<spacedim>::interface_end_active()
{
	return active_interface_cell_domain_cells.end();
}

template<unsigned int spacedim>
const typename std::vector< InterfaceCellDomainCellsDoF<spacedim> >&
DoFHandlerSystem<spacedim>::interface_active_iterators()
const
{
	return active_interface_cell_domain_cells;
}

template<unsigned int spacedim>
DomainCellDoFIterator<spacedim>
DoFHandlerSystem<spacedim>::domain_begin_active()
const
{
	return DomainCellDoFIterator<spacedim>(dof_handler_domain->begin_active(), *this);
}

template<unsigned int spacedim>
DomainCellDoFIterator<spacedim>
DoFHandlerSystem<spacedim>::domain_end_active()
const
{
	return DomainCellDoFIterator<spacedim>(*this);
}

template<unsigned int spacedim>
IteratorRange<DomainCellDoFIterator<spacedim>>
DoFHandlerSystem<spacedim>::domain_active_iterators()
const
{
	return IteratorRange<DomainCellDoFIterator<spacedim>>(domain_begin_active(), domain_end_active());
}

template<unsigned int spacedim>
unsigned int
DoFHandlerSystem<spacedim>::n_dofs_domain()
const
{
	return dof_handler_domain->n_dofs();
}

template<unsigned int spacedim>
unsigned int
DoFHandlerSystem<spacedim>::n_dofs_interface()
const
{
	return dof_handler_interface->n_dofs();
}

template<unsigned int spacedim>
unsigned int
DoFHandlerSystem<spacedim>::n_dofs_additional()
const
{
	return dof_indices_C.size();
}

template<unsigned int spacedim>
unsigned int
DoFHandlerSystem<spacedim>::n_dofs()
const
{
	return n_dofs_domain() + n_dofs_interface() + n_dofs_additional();
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::get_dof_indices(vector<unsigned int>& dof_indices)
const
{
	dof_indices = dof_indices_C;
	if(dof_renumbering != nullptr)
		dof_renumbering->convert_dof_indices(dof_indices);
}

template<unsigned int spacedim>
unsigned int
DoFHandlerSystem<spacedim>::get_dof_index(const unsigned int& dof_index)
const
{
	Assert(dof_index < dof_indices_C.size(), ExcMessage("You asked for the global index of a dof not related to a mesh which does not exist!"));
	vector<unsigned int> global_dof_index = {dof_indices_C[dof_index]};
	if(dof_renumbering != nullptr)
		dof_renumbering->convert_dof_indices(global_dof_index);

	return global_dof_index[0];
}


template<unsigned int spacedim>
const DoFHandler<spacedim, spacedim>&
DoFHandlerSystem<spacedim>::get_dof_handler_domain()
const
{
	return *(dof_handler_domain.get());
}

template<unsigned int spacedim>
const DoFHandler<spacedim-1, spacedim>&
DoFHandlerSystem<spacedim>::get_dof_handler_interface()
const
{
	return *(dof_handler_interface.get());
}

template<unsigned int spacedim>
DoFHandler<spacedim, spacedim>&
DoFHandlerSystem<spacedim>::get_dof_handler_domain()
{
	return *(dof_handler_domain.get());
}

template<unsigned int spacedim>
DoFHandler<spacedim-1, spacedim>&
DoFHandlerSystem<spacedim>::get_dof_handler_interface()
{
	return *(dof_handler_interface.get());
}

template<unsigned int spacedim>
const IndexSet&
DoFHandlerSystem<spacedim>::get_locally_owned_dofs()
const
{
	return this->locally_owned_dofs;
}

template<unsigned int spacedim>
const IndexSet&
DoFHandlerSystem<spacedim>::get_locally_relevant_dofs()
const
{
	return this->locally_relevant_dofs;
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::attach_dof_renumbering(const DoFRenumbering& dof_renumbering)
{
	this->dof_renumbering = &dof_renumbering;

#ifdef DEBUG
	//check that the renumbering scheme does not change the numbering range of the C's
	auto dof_indices_C_copy = dof_indices_C;
	dof_renumbering.convert_dof_indices(dof_indices_C_copy);
	sort(dof_indices_C_copy.begin(), dof_indices_C_copy.end());
	Assert(dof_indices_C_copy == dof_indices_C, ExcMessage("The dof renumbering scheme must not renumber the C's!"));
#endif // DEBUG

	set_locally_owned_dofs();
	set_locally_relevant_dofs();
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::make_hanging_node_constraints(AffineConstraints<double>& constraint_matrix)
const
{
	constraint_matrix.clear();
	const auto this_procs_n_procs = tria_system->get_this_proc_n_procs();
	if(this_procs_n_procs.second > 1)
		constraint_matrix.reinit(locally_relevant_dofs_standard_numbering);

	//hanging node constraints on interface
	//in 2d there are no hanging nodes on interfaces - so this has only to be done in 3d
	AffineConstraints<double> interface_hanging_node_constraints;
	if( (spacedim > 2) && (get_dof_handler_interface().n_dofs() > 0))
	{
		if(this_procs_n_procs.second > 1)
		{
			IndexSet locally_relevant_dofs_interface;
			DoFTools::extract_locally_relevant_dofs(get_dof_handler_interface(), locally_relevant_dofs_interface);
			//there is the problem that the extracted index set cannot accommodate the shifted dofs -> construct an enlarged index set
			IndexSet locally_relevant_dofs_interface_enlarged(n_dofs_domain() + n_dofs_interface());
			locally_relevant_dofs_interface_enlarged.add_indices(locally_relevant_dofs_interface);
			locally_relevant_dofs_interface_enlarged.compress();
			interface_hanging_node_constraints.reinit(locally_relevant_dofs_interface_enlarged);
		}

		DoFTools::make_hanging_node_constraints(get_dof_handler_interface(), interface_hanging_node_constraints);
		interface_hanging_node_constraints.shift(n_dofs_domain());
	}
	interface_hanging_node_constraints.close();

	//hanging node constraints on domain
	AffineConstraints<double> domain_hanging_node_constraints;
	if(this_procs_n_procs.second > 1)
	{
		IndexSet locally_relevant_dofs_domain;
		DoFTools::extract_locally_relevant_dofs(get_dof_handler_domain(), locally_relevant_dofs_domain);
		domain_hanging_node_constraints.reinit(locally_relevant_dofs_domain);
	}
	DoFTools::make_hanging_node_constraints(get_dof_handler_domain(), domain_hanging_node_constraints);
	domain_hanging_node_constraints.close();

	//merge interface constraints into domain constraints
	constraint_matrix.merge(interface_hanging_node_constraints, AffineConstraints<double>::MergeConflictBehavior::no_conflicts_allowed, true);

	//merge domain constraints into domain constraints
	constraint_matrix.merge(domain_hanging_node_constraints, AffineConstraints<double>::MergeConflictBehavior::no_conflicts_allowed, true);

	//renumber
	if(dof_renumbering != nullptr)
		Auxiliary::renumber_constraints(constraint_matrix, *dof_renumbering, false);
}

template<unsigned int spacedim>
template<class VectorType>
void
DoFHandlerSystem<spacedim>::split_vector(	const VectorType&	in_vect,
											VectorType&			out_vect_domain,
											VectorType&			out_vect_interface,
											VectorType&			out_vect_C)
const
{
	if(n_dofs_domain() > 0)
		split_vector_implementation(in_vect, out_vect_domain, 0, n_dofs_domain());
	if(n_dofs_interface() > 0)
	split_vector_implementation(in_vect, out_vect_interface, n_dofs_domain(), n_dofs_domain() + n_dofs_interface());
	if(n_dofs_additional() > 0)
		split_vector_implementation(in_vect, out_vect_C, n_dofs_domain() + n_dofs_interface(), n_dofs_domain() + n_dofs_interface() + n_dofs_additional());
}

template<unsigned int spacedim>
const vector<unsigned int>&
DoFHandlerSystem<spacedim>::get_n_dofs_per_processor()
const
{
	return n_dofs_per_processor;
}

template<unsigned int spacedim>
unsigned int
DoFHandlerSystem<spacedim>::get_single_dof_index_component_interface(const unsigned int component)
const
{
	unsigned int dof_index = numbers::invalid_dof_index;
	for(const auto& interface_cell_domain_cells : interface_active_iterators())
	{
		const auto& interface_cell = interface_cell_domain_cells.interface_cell;
		if(interface_cell->is_locally_owned())
		{
			const unsigned int n_dofs = interface_cell->get_fe().dofs_per_cell;
			for(unsigned int sf = 0; sf < n_dofs; ++sf)
			{
				if( (interface_cell->get_fe().is_primitive(sf)) && (interface_cell->get_fe().system_to_component_index(sf).first == component) )
				{
					vector<unsigned int> dof_indices_local_global;
					dof_indices_local_global.resize(n_dofs);
					interface_cell.get_dof_indices(dof_indices_local_global);
					if(locally_owned_dofs.is_element(dof_indices_local_global[sf]))
						dof_index = dof_indices_local_global[sf];
					// TODO: Why is here no break?
				}
			}
		}
	}

	// communicate results in parallel and choose the dof_index found on the processor with the smallest rank
	const auto this_proc_n_procs = tria_system->get_this_proc_n_procs();
	if(this_proc_n_procs.second > 1)
	{
		vector<unsigned int> dof_indices(this_proc_n_procs.second);
		dof_indices[this_proc_n_procs.first] = dof_index;
#ifdef DEAL_II_WITH_MPI
		const auto tria_domain_ptr = dynamic_cast<const dealii::parallel::distributed::Triangulation<spacedim, spacedim>*>(&(dof_handler_domain->get_triangulation()));
		Assert(tria_domain_ptr != nullptr, ExcMessage("Internal error!"));
		unsigned int send_value = dof_indices[this_proc_n_procs.first];
		int ierr =  MPI_Allgather(&send_value, 1, MPI_UNSIGNED, dof_indices.data(), 1, MPI_UNSIGNED, tria_domain_ptr->get_communicator());
		AssertThrowMPI(ierr);
#else
		Assert(this_proc_n_procs.second == 1, ExcMessage("Internal error: deal.II not configured with MPI, but computation run with more than one processor"));
#endif //DEAL_II_WITH_MPI
		for(unsigned int m = 0; m < dof_indices.size(); ++m)
		{
			if(dof_indices[m] != numbers::invalid_dof_index)
			{
				dof_index = dof_indices[m];
				break;
			}
		}
	}

	if(locally_relevant_dofs.is_element(dof_index))
		return dof_index;
	else
		return numbers::invalid_dof_index;
}

template<unsigned int spacedim>
unsigned int
DoFHandlerSystem<spacedim>::get_single_dof_index_component_domain(const unsigned int component)
const
{
	unsigned int dof_index = numbers::invalid_dof_index;
	for(const auto& domain_cell : domain_active_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			const unsigned int n_dofs = domain_cell->get_fe().dofs_per_cell;
			for(unsigned int sf = 0; sf < n_dofs; ++sf)
			{
				if( (domain_cell->get_fe().is_primitive(sf)) && (domain_cell->get_fe().system_to_component_index(sf).first == component) )
				{
					vector<unsigned int> dof_indices_local_global;
					dof_indices_local_global.resize(n_dofs);
					domain_cell.get_dof_indices(dof_indices_local_global);
					if(locally_owned_dofs.is_element(dof_indices_local_global[sf]))
						dof_index = dof_indices_local_global[sf];
					// TODO: Why is here no break?
				}
			}
		}
	}

	// communicate results in parallel and choose the dof_index found on the processor with the smallest rank
	const auto this_proc_n_procs = tria_system->get_this_proc_n_procs();
	if(this_proc_n_procs.second > 1)
	{
		vector<unsigned int> dof_indices(this_proc_n_procs.second);
		dof_indices[this_proc_n_procs.first] = dof_index;
#ifdef DEAL_II_WITH_MPI
		const auto tria_domain_ptr = dynamic_cast<const dealii::parallel::distributed::Triangulation<spacedim, spacedim>*>(&(dof_handler_domain->get_triangulation()));
		Assert(tria_domain_ptr != nullptr, ExcMessage("Internal error!"));
		unsigned int send_value = dof_indices[this_proc_n_procs.first];
		int ierr =  MPI_Allgather(&send_value, 1, MPI_UNSIGNED, dof_indices.data(), 1, MPI_UNSIGNED, tria_domain_ptr->get_communicator());
		AssertThrowMPI(ierr);
#else
		Assert(this_proc_n_procs.second == 1, ExcMessage("Internal error: deal.II not configured with MPI, but computation run with more than one processor"));
#endif //DEAL_II_WITH_MPI
		for(unsigned int m = 0; m < dof_indices.size(); ++m)
		{
			if(dof_indices[m] != numbers::invalid_dof_index)
			{
				dof_index = dof_indices[m];
				break;
			}
		}
	}

	if(locally_relevant_dofs.is_element(dof_index))
		return dof_index;
	else
		return numbers::invalid_dof_index;
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::get_dof_indices_component_interface(const unsigned int	component,
																set<unsigned int>&	indices)
const
{
	vector<unsigned int> dof_indices_local_global;

	for(const auto& interface_cell_domain_cells : interface_active_iterators())
	{
		const auto& interface_cell = interface_cell_domain_cells.interface_cell;
		if(interface_cell->is_locally_owned())
		{
			const unsigned int n_dofs = interface_cell->get_fe().dofs_per_cell;
			dof_indices_local_global.resize(n_dofs);
			interface_cell.get_dof_indices(dof_indices_local_global);
			for(unsigned int sf = 0; sf < n_dofs; ++sf)
			{
				if( (interface_cell->get_fe().is_primitive(sf)) && (interface_cell->get_fe().system_to_component_index(sf).first == component) )
					indices.insert(dof_indices_local_global[sf]);
			}
		}
	}
	return;
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::get_dof_indices_component_domain(	const unsigned int	component,
																set<unsigned int>&	indices)
const
{
	vector<unsigned int> dof_indices_local_global;

	for(const auto& domain_cell : domain_active_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			const unsigned int n_dofs = domain_cell->get_fe().dofs_per_cell;
			dof_indices_local_global.resize(n_dofs);
			domain_cell.get_dof_indices(dof_indices_local_global);
			for(unsigned int sf = 0; sf < n_dofs; ++sf)
			{
				if( (domain_cell->get_fe().is_primitive(sf)) && (domain_cell->get_fe().system_to_component_index(sf).first == component) )
				{
					indices.insert(dof_indices_local_global[sf]);
				}
			}
		}
	}
	return;
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::get_map_dof_index_component_interface(std::vector<unsigned int>& components)
const
{
	components.resize(n_dofs_interface(), 0);

	vector<unsigned int> dof_indices_local_global;
	for(const auto& interface_cell_domain_cells : interface_active_iterators())
	{
		const auto& interface_cell = interface_cell_domain_cells.interface_cell;
		if(interface_cell->is_locally_owned())
		{
			const unsigned int n_dofs = interface_cell->get_fe().dofs_per_cell;
			dof_indices_local_global.resize(n_dofs);
			interface_cell.get_dof_indices(dof_indices_local_global);
			for(unsigned int sf = 0; sf < n_dofs; ++sf)
				components[dof_indices_local_global[sf]] = interface_cell->get_fe().system_to_component_index(sf).first;
		}
	}
	return;
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::get_map_dof_index_component_domain(std::vector<unsigned int>& components)
const
{
	components.resize(n_dofs_domain(), 0);

	vector<unsigned int> dof_indices_local_global;
	for(const auto& domain_cell : domain_active_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			const unsigned int n_dofs = domain_cell->get_fe().dofs_per_cell;
			dof_indices_local_global.resize(n_dofs);
			domain_cell.get_dof_indices(dof_indices_local_global);
			for(unsigned int sf = 0; sf < n_dofs; ++sf)
				components[dof_indices_local_global[sf]] = domain_cell->get_fe().system_to_component_index(sf).first;
		}
	}
	return;
}

template<unsigned int spacedim>
struct PointComparator
{
	bool operator()(const Point<spacedim>& p1, const Point<spacedim>& p2)
	const
	{
	    return make_tuple(spacedim < 1 ? 0.0 : p1[0], spacedim < 2 ? 0.0 : p1[1], spacedim < 3 ? 0.0 : p1[2]) < make_tuple(spacedim < 1 ? 0.0 : p2[0], spacedim < 2 ? 0.0 : p2[1], spacedim < 3 ? 0.0 : p2[2]);
    }
};

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::get_map_dof_index_support_point_index_interface(vector<unsigned int>&					support_points,
																			map<unsigned int, Point<spacedim>>&		map_support_point_index_support_point_location,
																			const Mapping<spacedim-1, spacedim>&	mapping)
const
{
	support_points.resize(n_dofs_interface());
	map<Point<spacedim>, unsigned int, PointComparator<spacedim>> map_support_point_location_support_point_index;
	vector<unsigned int> dof_indices_local_global;
	Point<spacedim> support_point;
	Point<spacedim-1> unit_support_point_interface;
	unsigned int support_point_counter = 0;
	for(const auto& interface_cell_domain_cells : interface_active_iterators())
	{
		const auto& interface_cell = interface_cell_domain_cells.interface_cell;
		if(interface_cell->is_locally_owned())
		{
			const unsigned int n_dofs = interface_cell->get_fe().dofs_per_cell;
			dof_indices_local_global.resize(n_dofs);
			interface_cell.get_dof_indices(dof_indices_local_global);
			for(unsigned int sf = 0; sf < n_dofs; ++sf)
			{
				unit_support_point_interface = interface_cell->get_fe().unit_support_point(sf);
				support_point = mapping.transform_unit_to_real_cell(interface_cell, unit_support_point_interface);
				if(map_support_point_location_support_point_index.find(support_point) == map_support_point_location_support_point_index.end())
				{
					map_support_point_location_support_point_index[support_point] = support_point_counter;
					map_support_point_index_support_point_location[support_point_counter] = support_point;
					support_point_counter++;
				}
				support_points[dof_indices_local_global[sf]] = map_support_point_location_support_point_index[support_point];
			}
		}
	}
	return;
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::get_map_dof_index_support_point_index_domain(	vector<unsigned int>&				support_points,
																			map<unsigned int, Point<spacedim>>&	map_support_point_index_support_point_location,
																			const Mapping<spacedim>&			mapping)
const
{
	support_points.resize(n_dofs_domain());
	map<Point<spacedim>, unsigned int, PointComparator<spacedim>> map_support_point_location_support_point_index;
	vector<unsigned int> dof_indices_local_global;
	Point<spacedim> support_point, unit_support_point_domain;
	unsigned int support_point_counter = 0;
	for(const auto& domain_cell : domain_active_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			const unsigned int n_dofs = domain_cell->get_fe().dofs_per_cell;
			dof_indices_local_global.resize(n_dofs);
			domain_cell.get_dof_indices(dof_indices_local_global);
			for(unsigned int sf = 0; sf < n_dofs; ++sf)
			{
				unit_support_point_domain = domain_cell->get_fe().unit_support_point(sf);
				support_point = mapping.transform_unit_to_real_cell(domain_cell, unit_support_point_domain);
				if(map_support_point_location_support_point_index.find(support_point) == map_support_point_location_support_point_index.end())
				{
					map_support_point_location_support_point_index[support_point] = support_point_counter;
					map_support_point_index_support_point_location[support_point_counter] = support_point;
					support_point_counter++;
				}
				support_points[dof_indices_local_global[sf]] = map_support_point_location_support_point_index[support_point];
			}
		}
	}
	return;
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::update_interface_domain_relation()
{
	active_interface_cell_domain_cells.clear();
	for(const auto& interface_cell_domain_cell : tria_system->interface_active_iterators())
		active_interface_cell_domain_cells.push_back(InterfaceCellDomainCellsDoF<spacedim>(interface_cell_domain_cell, *this ));
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::set_locally_owned_dofs()
{
	locally_owned_dofs.clear();
	locally_owned_dofs.set_size(n_dofs());
	if(dof_renumbering != nullptr)
		Auxiliary::add_to_index_set(*dof_renumbering, locally_owned_dofs_standard_numbering, locally_owned_dofs);
	else
		locally_owned_dofs = locally_owned_dofs_standard_numbering;
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::set_locally_relevant_dofs()
{
	locally_relevant_dofs.clear();
	locally_relevant_dofs.set_size(n_dofs());
	if(dof_renumbering != nullptr)
		Auxiliary::add_to_index_set(*dof_renumbering, locally_relevant_dofs_standard_numbering, locally_relevant_dofs);
	else
		locally_relevant_dofs = locally_relevant_dofs_standard_numbering;
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::set_locally_owned_dofs_standard_numbering()
{
	const unsigned int n_dofs_total = n_dofs();

	locally_owned_dofs_standard_numbering.clear();
	locally_owned_dofs_standard_numbering.set_size(n_dofs_total);

	const IndexSet& locally_owned_dofs_domain = dof_handler_domain->locally_owned_dofs();
	locally_owned_dofs_standard_numbering.add_indices(locally_owned_dofs_domain);

	const IndexSet& locally_owned_dofs_interface = dof_handler_interface->locally_owned_dofs();
	locally_owned_dofs_standard_numbering.add_indices(locally_owned_dofs_interface, n_dofs_domain());

	IndexSet indices_C(n_dofs_total);
	for(const auto& index : locally_owned_dof_indices_C)
		indices_C.add_index(index);
	indices_C.compress();
	locally_owned_dofs_standard_numbering.add_indices(indices_C);

	locally_owned_dofs_standard_numbering.compress();
}

template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::set_locally_relevant_dofs_standard_numbering()
{
	const unsigned int n_dofs_total = n_dofs();

	locally_relevant_dofs_standard_numbering.clear();
	locally_relevant_dofs_standard_numbering.set_size(n_dofs_total);

	IndexSet locally_relevant_dofs_domain(n_dofs_domain());
	DoFTools::extract_locally_relevant_dofs(*dof_handler_domain, locally_relevant_dofs_domain);
	locally_relevant_dofs_standard_numbering.add_indices(locally_relevant_dofs_domain);

	IndexSet locally_relevant_dofs_interface(n_dofs_interface());
	DoFTools::extract_locally_relevant_dofs(*dof_handler_interface, locally_relevant_dofs_interface);
	locally_relevant_dofs_standard_numbering.add_indices(locally_relevant_dofs_interface, n_dofs_domain());

	IndexSet indices_C(n_dofs_total);
	for(const auto& index : dof_indices_C)
		indices_C.add_index(index);
	indices_C.compress();
	locally_relevant_dofs_standard_numbering.add_indices(indices_C);

	locally_relevant_dofs_standard_numbering.compress();
}


#ifdef DEAL_II_WITH_MPI
template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::split_vector_implementation(const LinearAlgebra::distributed::Vector<double>&	in_vect,
														LinearAlgebra::distributed::Vector<double>&			out_vect,
														const unsigned int									window_begin,
														const unsigned int									window_end)
const
{

	Assert( ( (int)window_end - (int)window_begin >= 0), ExcMessage("The slice must have positive size!") );
	Assert( locally_owned_dofs.is_contiguous(), ExcMessage("The locally owned indices of the standard numbering must form a contiguous range for this function to work!") )

	//slice output vector index sets
	IndexSet out_vect_locally_owned_indices = locally_owned_dofs_standard_numbering.get_view(window_begin, window_end);
	IndexSet out_vect_locally_relevant_indices = locally_relevant_dofs_standard_numbering.get_view(window_begin, window_end);

	//figure out where to take the locally owned values in out_vect from (note: it is not possible to assume here that a contiguous range in @p out_vect
	//maps to a contiguous range in @p in_vect!)
	vector<unsigned int> source_indices_locally_owned;
	source_indices_locally_owned.reserve(out_vect_locally_owned_indices.n_elements());
	for(const auto& index : out_vect_locally_owned_indices)
		source_indices_locally_owned.push_back(index + window_begin);
	if(dof_renumbering != nullptr)
		dof_renumbering->convert_dof_indices(source_indices_locally_owned);

	//reinit the output vector
	out_vect.reinit(out_vect_locally_owned_indices, out_vect_locally_relevant_indices, in_vect.get_mpi_communicator());

	//copy over locally owned indices
	auto source_index_locally_owned = source_indices_locally_owned.begin();
	for(const auto& index : out_vect_locally_owned_indices)
	{
		out_vect[index] = in_vect[*source_index_locally_owned];
		++source_index_locally_owned;
	}

	//communicate ghosts
	out_vect.update_ghost_values();
}
#endif // DEAL_II_WITH_MPI


template<unsigned int spacedim>
void
DoFHandlerSystem<spacedim>::split_vector_implementation(const Vector<double>&	in_vect,
														Vector<double>&			out_vect,
														const unsigned int		window_begin,
														const unsigned int		window_end)
const
{

	Assert( ( (int)window_end - (int)window_begin >= 0), ExcMessage("The slice must have positive size!") );

	//output vector index set
	IndexSet out_vect_locally_owned_indices(window_end - window_begin);
	out_vect_locally_owned_indices.add_range(0, window_end - window_begin);

	//figure out where to take the values in out_vect from
	vector<unsigned int> source_indices_locally_owned;
	source_indices_locally_owned.reserve(out_vect_locally_owned_indices.n_elements());
	for(const auto& index : out_vect_locally_owned_indices)
		source_indices_locally_owned.push_back(index + window_begin);
	if(dof_renumbering != nullptr)
		dof_renumbering->convert_dof_indices(source_indices_locally_owned);

	//reinit the output vector
	out_vect.reinit(window_end - window_begin);

	//copy over locally owned indices
	auto source_index_locally_owned = source_indices_locally_owned.begin();
	for(const auto& index : out_vect_locally_owned_indices)
	{
		out_vect[index] = in_vect[*source_index_locally_owned];
		++source_index_locally_owned;
	}
}

template class InterfaceCellDomainCellsDoF<2>;
template class InterfaceCellDomainCellsDoF<3>;
template class DoFHandlerSystem<2>;
template class DoFHandlerSystem<3>;
template class InterfaceCellDoFIterator<2>;
template class InterfaceCellDoFIterator<3>;
template class DomainCellDoFIterator<2>;
template class DomainCellDoFIterator<3>;


template
void
DoFHandlerSystem<2>::split_vector<Vector<double>>(	const Vector<double>&,
													Vector<double>&,
													Vector<double>&,
													Vector<double>&)
const;

template
void
DoFHandlerSystem<3>::split_vector<Vector<double>>(	const Vector<double>&,
													Vector<double>&,
													Vector<double>&,
													Vector<double>&)
const;


//if there's no p4est, this functionality doesn't make sense
#ifdef DEAL_II_WITH_MPI

template
void
DoFHandlerSystem<2>::split_vector<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																				LinearAlgebra::distributed::Vector<double>&,
																				LinearAlgebra::distributed::Vector<double>&,
																				LinearAlgebra::distributed::Vector<double>&)
const;

template
void
DoFHandlerSystem<3>::split_vector<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																				LinearAlgebra::distributed::Vector<double>&,
																				LinearAlgebra::distributed::Vector<double>&,
																				LinearAlgebra::distributed::Vector<double>&)
const;

#endif // DEAL_II_WITH_MPI

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
