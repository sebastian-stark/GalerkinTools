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

#include <galerkin_tools/tools.h>

#include <map>
#include <vector>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

void
Auxiliary::convert_local_indices_to_global_indices(	const vector<unsigned int>& dof_indices_local,
													vector<unsigned int>& 		dof_indices_global,
													const vector<unsigned int>& dof_indices_local_global )
{
	dof_indices_global.resize(dof_indices_local.size());
	for(unsigned int dof_indices_local_n = 0; dof_indices_local_n < dof_indices_local.size(); ++dof_indices_local_n)
		dof_indices_global[dof_indices_local_n] = dof_indices_local_global[dof_indices_local[dof_indices_local_n]];
}

void
Auxiliary::combine_dof_indices(	const vector<unsigned int>& dof_indices_global_interface,
								const vector<unsigned int>& dof_indices_global_minus,
								const vector<unsigned int>& dof_indices_global_plus,
								const vector<unsigned int>& dof_indices_global_C,
								vector<unsigned int>& 		dof_indices_interface_dof_indices_combined,
								vector<unsigned int>& 		dof_indices_minus_dof_indices_combined,
								vector<unsigned int>& 		dof_indices_plus_dof_indices_combined,
								vector<unsigned int>& 		dof_indices_C_dof_indices_combined,
								vector<unsigned int>& 		dof_indices_global_combined)
{
	dof_indices_interface_dof_indices_combined.resize(dof_indices_global_interface.size());
	dof_indices_minus_dof_indices_combined.resize(dof_indices_global_minus.size());
	dof_indices_plus_dof_indices_combined.resize(dof_indices_global_plus.size());
	dof_indices_C_dof_indices_combined.resize(dof_indices_global_C.size());

	const unsigned int dof_indices_global_size_interface = dof_indices_global_interface.size();
	const unsigned int dof_indices_global_size_interface_minus = dof_indices_global_interface.size()+dof_indices_global_minus.size();


	//assemble dof_indices_plus_dof_indices_combined and corresponding entries in dof_indices_global_combined
	//this needs a bit more work in order to avoid that there are duplicate global indices in dof_indices_global_combined
	map<unsigned int, unsigned int> map_minus;
	for(unsigned int dof_indices_global_minus_n=0; dof_indices_global_minus_n<dof_indices_global_minus.size(); ++dof_indices_global_minus_n)
		map_minus[dof_indices_global_minus[dof_indices_global_minus_n]]=dof_indices_global_minus_n;
	unsigned int counter=0;
	for(unsigned int dof_indices_global_plus_n=0; dof_indices_global_plus_n<dof_indices_global_plus.size(); ++dof_indices_global_plus_n)
	{
		const auto& map_minus_it = map_minus.find(dof_indices_global_plus[dof_indices_global_plus_n]);
		if(map_minus_it != map_minus.end())
			dof_indices_plus_dof_indices_combined[dof_indices_global_plus_n]=dof_indices_global_size_interface+map_minus_it->second;
		else
		{
			dof_indices_plus_dof_indices_combined[dof_indices_global_plus_n]=dof_indices_global_size_interface_minus+counter;
			++counter;
		}
	}
	dof_indices_global_combined.reserve(dof_indices_global_size_interface_minus+counter+dof_indices_global_C.size());
	dof_indices_global_combined.resize(dof_indices_global_size_interface_minus+counter);
	for(unsigned int dof_indices_global_plus_n=0; dof_indices_global_plus_n<dof_indices_global_plus.size(); ++dof_indices_global_plus_n)
		dof_indices_global_combined[dof_indices_plus_dof_indices_combined[dof_indices_global_plus_n]]=dof_indices_global_plus[dof_indices_global_plus_n];

	//assemble dof_indices_interface_dof_indices_combined and corresponding entries in dof_indices_global_combined
	for(unsigned int dof_indices_global_interface_n=0; dof_indices_global_interface_n<dof_indices_global_interface.size(); ++dof_indices_global_interface_n)
	{
		dof_indices_interface_dof_indices_combined[dof_indices_global_interface_n]=dof_indices_global_interface_n;
		dof_indices_global_combined[dof_indices_interface_dof_indices_combined[dof_indices_global_interface_n]]=dof_indices_global_interface[dof_indices_global_interface_n];

	}

	//assemble dof_indices_minus_dof_indices_combined and corresponding entries in dof_indices_global_combined
	for(unsigned int dof_indices_global_minus_n=0; dof_indices_global_minus_n<dof_indices_global_minus.size(); dof_indices_global_minus_n++)
	{
		dof_indices_minus_dof_indices_combined[dof_indices_global_minus_n]=dof_indices_global_minus_n+dof_indices_global_size_interface;
		dof_indices_global_combined[dof_indices_minus_dof_indices_combined[dof_indices_global_minus_n]]=dof_indices_global_minus[dof_indices_global_minus_n];
	}

	//finally assemble dof_indices_C_dof_indices_combined and corresponding entries in dof_indices_global_combined
	counter = 0;
	for(const auto& dof_indices_global_C_n : dof_indices_global_C)
	{
		dof_indices_global_combined.push_back(dof_indices_global_C_n);
		dof_indices_C_dof_indices_combined[counter]=dof_indices_global_combined.size()-1;
		++counter;
	}
}


void
Auxiliary::add_to_index_set(const DoFRenumbering&	dof_renumbering,
							const IndexSet& 		in,
							IndexSet& 				out)
{
	if(in.n_intervals() > 0)
	{
		for(auto interval = in.begin_intervals(); interval != in.end_intervals(); ++interval)
		{
			const auto new_intervals = dof_renumbering.convert_range(*(interval->begin()), interval->last());
			for(const auto& new_interval : new_intervals)
				out.add_range(new_interval.first, new_interval.second + 1);
		}
		out.compress();
	}
}

void
Auxiliary::renumber_constraints(AffineConstraints<double>&	constraint_matrix,
								const DoFRenumbering&		dof_renumbering,
								const bool					close)
{

	//do a first round to figure out the temporary memory needed
	const unsigned int n_constraints = constraint_matrix.n_constraints();
	unsigned int n_entries = 0;
	unsigned int max_entries = 0;
	for(const auto& line : constraint_matrix.get_lines())
	{
		if(line.entries.size() > max_entries)
			max_entries = line.entries.size();
		n_entries += line.entries.size();
	}

	//now copy the constraints to intermediate data structures
	vector<unsigned int> constrained_dofs;
	constrained_dofs.reserve(n_constraints);
	vector<unsigned int> n_entries_in_constraint;
	n_entries_in_constraint.reserve(n_constraints);
	vector<unsigned int> entries_indices;
	entries_indices.reserve(n_entries);
	vector<double> entries_coefficients;
	entries_coefficients.reserve(n_entries);
	vector<double> inhomogeneities;
	inhomogeneities.reserve(n_constraints);
	for(const auto& line : constraint_matrix.get_lines())
	{
		constrained_dofs.push_back(line.index);
		n_entries_in_constraint.push_back(line.entries.size());
		for(const auto& entry : line.entries)
		{
			entries_indices.push_back(entry.first);
			entries_coefficients.push_back(entry.second);
		}
		inhomogeneities.push_back(line.inhomogeneity);
	}

	//do the renumbering
	dof_renumbering.convert_dof_indices(constrained_dofs);
	dof_renumbering.convert_dof_indices(entries_indices);
	const IndexSet& local_lines_original = constraint_matrix.get_local_lines();
	IndexSet local_lines_renumbered(local_lines_original.size());
	Auxiliary::add_to_index_set(dof_renumbering, local_lines_original, local_lines_renumbered);
	//clear the constraint matrix
	constraint_matrix.clear();

	//set up the renumbered constraint_matrix object
	constraint_matrix.reinit(local_lines_renumbered);
	vector<pair<unsigned int, double>> temp_entries;
	temp_entries.reserve(max_entries);
	auto entry_indices_iterator = entries_indices.begin();
	auto entry_coefficients_iterator = entries_coefficients.begin();
	for(unsigned int line = 0; line < n_constraints; ++line)
	{
		temp_entries.resize(n_entries_in_constraint[line]);
		for(unsigned int entry = 0; entry < n_entries_in_constraint[line]; ++entry)
		{
			temp_entries[entry] = make_pair(*entry_indices_iterator, *entry_coefficients_iterator);
			++entry_indices_iterator;
			++entry_coefficients_iterator;
		}
		constraint_matrix.add_line(constrained_dofs[line]);
		constraint_matrix.add_entries(constrained_dofs[line], temp_entries);
		constraint_matrix.set_inhomogeneity(constrained_dofs[line], inhomogeneities[line]);
	}

	if(close)
		constraint_matrix.close();

}

#ifdef DEAL_II_WITH_MPI
template<unsigned int spacedim>
void
Auxiliary::compute_dof_renumbering_contiguous(	const DoFHandlerSystem<spacedim>&	dof_handler_system,
												DoFRenumberingOffset&				dof_renumbering_offset)
{
	dof_renumbering_offset.clear();

	const auto& dof_handler_domain = dof_handler_system.get_dof_handler_domain();
	const auto& dof_handler_interface = dof_handler_system.get_dof_handler_interface();

	//determine dof_offsets_domain and dof_offsets_interface
	const auto tria_domain_ptr = dynamic_cast<const dealii::parallel::Triangulation<spacedim, spacedim>*>(&(dof_handler_domain.get_triangulation()));
	if( (tria_domain_ptr != nullptr) && (dof_handler_system.n_dofs_domain() != 0) && (dof_handler_system.n_dofs_interface() != 0))
	{
		const auto tria_interface_ptr = dynamic_cast<const dealii::parallel::Triangulation<spacedim-1, spacedim>*>(&(dof_handler_interface.get_triangulation()));

		Assert(dof_handler_domain.locally_owned_dofs().is_contiguous(), ExcMessage("You can call compute_offsets only if the locally owned domain related dofs are contiguous!"));
		Assert(dof_handler_interface.locally_owned_dofs().is_contiguous(), ExcMessage("You can call compute_offsets only if the locally owned interface related dofs are contiguous!"));

		//number of locally owned dofs
		const auto n_locally_owned_dofs_per_processor_domain = dof_handler_domain.n_locally_owned_dofs_per_processor();
		const auto n_locally_owned_dofs_per_processor_interface = dof_handler_interface.n_locally_owned_dofs_per_processor();

		//number of participating processors and number of this processor
		const unsigned int n_procs = n_locally_owned_dofs_per_processor_domain.size();
		const unsigned int this_proc = tria_domain_ptr->locally_owned_subdomain();

		//global dof indices of first dofs on respective subdomains
		vector<unsigned int> dof_start_domain(n_procs);
		vector<unsigned int> dof_start_interface(n_procs);
		if( dof_handler_domain.locally_owned_dofs().n_elements() > 0)
			dof_start_domain[this_proc] = dof_handler_domain.locally_owned_dofs().nth_index_in_set(0);
		else
			dof_start_domain[this_proc] = 0;
		if( dof_handler_interface.locally_owned_dofs().n_elements() > 0)
			dof_start_interface[this_proc] = dof_handler_interface.locally_owned_dofs().nth_index_in_set(0) + dof_handler_domain.n_dofs();
		else
			dof_start_interface[this_proc] = 0;

		unsigned int send_value = dof_start_domain[this_proc];
		int ierr =  MPI_Allgather(&send_value, 1, MPI_UNSIGNED, dof_start_domain.data(), 1, MPI_UNSIGNED, tria_domain_ptr->get_communicator());
		AssertThrowMPI(ierr);

		send_value = dof_start_interface[this_proc];
		ierr =  MPI_Allgather(&send_value, 1, MPI_UNSIGNED, dof_start_interface.data(), 1, MPI_UNSIGNED, tria_domain_ptr->get_communicator());
		AssertThrowMPI(ierr);

		//set up dof_offsets_domain and dof_offsets_interface
		unsigned int current_dof_index = 0;
		const auto& ghost_owners_domain = tria_domain_ptr->ghost_owners();
		const auto& ghost_owners_interface = tria_interface_ptr->ghost_owners();
		//keep track by an index set for which dofs we locally have renumbering rules
		IndexSet renumbered_ranges(dof_handler_system.n_dofs());
		//do a first round to determine the offsets of the current processor
		for(unsigned int proc = 0; proc < n_procs; ++proc)
		{
			if( (proc == this_proc) && (n_locally_owned_dofs_per_processor_domain[proc] > 0) )
			{
				dof_renumbering_offset.add_range(dof_start_domain[proc], dof_start_domain[proc] + n_locally_owned_dofs_per_processor_domain[proc] - 1, (int)current_dof_index - (int)dof_start_domain[proc]);
				renumbered_ranges.add_range(dof_start_domain[proc], dof_start_domain[proc] + n_locally_owned_dofs_per_processor_domain[proc]);
			}
			current_dof_index += n_locally_owned_dofs_per_processor_domain[proc];

			if( ( proc == this_proc ) && (n_locally_owned_dofs_per_processor_interface[proc] > 0) )
			{
				dof_renumbering_offset.add_range(dof_start_interface[proc], dof_start_interface[proc] + n_locally_owned_dofs_per_processor_interface[proc] - 1, (int)current_dof_index - (int)dof_start_interface[proc]);
				renumbered_ranges.add_range(dof_start_interface[proc], dof_start_interface[proc] + n_locally_owned_dofs_per_processor_interface[proc]);
				break;
			}
			current_dof_index += n_locally_owned_dofs_per_processor_interface[proc];
		}

		//add C's
		vector<unsigned int> dof_indices_C;
		dof_handler_system.get_dof_indices(dof_indices_C);
		if(dof_indices_C.size() > 0)
			dof_renumbering_offset.add_range(dof_indices_C[0], dof_indices_C.back(), 0);

		//do a second round to determine the offsets of the ghost owners
		current_dof_index = 0;
		for(unsigned int proc = 0; proc < n_procs; ++proc)
		{
			if( (ghost_owners_domain.find(proc) != ghost_owners_domain.end()) && (n_locally_owned_dofs_per_processor_domain[proc] > 0) )
			{
				dof_renumbering_offset.add_range(dof_start_domain[proc], dof_start_domain[proc] + n_locally_owned_dofs_per_processor_domain[proc] - 1, (int)current_dof_index - (int)dof_start_domain[proc]);
				renumbered_ranges.add_range(dof_start_domain[proc], dof_start_domain[proc] + n_locally_owned_dofs_per_processor_domain[proc]);
			}
			current_dof_index += n_locally_owned_dofs_per_processor_domain[proc];

			if( (ghost_owners_interface.find(proc) != ghost_owners_interface.end()) && (n_locally_owned_dofs_per_processor_interface[proc] > 0))
			{
				dof_renumbering_offset.add_range(dof_start_interface[proc], dof_start_interface[proc] + n_locally_owned_dofs_per_processor_interface[proc] - 1, (int)current_dof_index - (int)dof_start_interface[proc]);
				renumbered_ranges.add_range(dof_start_interface[proc], dof_start_interface[proc] + n_locally_owned_dofs_per_processor_interface[proc]);
			}
			current_dof_index += n_locally_owned_dofs_per_processor_interface[proc];
		}

		//Still, there can be the case that a ghost cell is itself adjacent to a cell of yet another processor; so there may be the necessity to renumber dof indices which are neither owned
		//by the current processor nor by the ghost owners. We treat this by explicitly looking for these dofs and including the required ranges.
		IndexSet remaining_indices = dof_handler_system.get_locally_relevant_dofs();
		renumbered_ranges.compress();
		remaining_indices.subtract_set(renumbered_ranges);
		remaining_indices.compress();
		if(remaining_indices.n_intervals() > 0)
		{
			current_dof_index = 0;
			unsigned int range_0_start, range_0_end, range_1_start, range_1_end;
			for(unsigned int proc = 0; proc < n_procs; ++proc)
			{
				//add the additionally required ranges
				if(n_locally_owned_dofs_per_processor_domain[proc] > 0)
				{
					range_0_start = dof_start_domain[proc];
					range_0_end =  dof_start_domain[proc] + n_locally_owned_dofs_per_processor_domain[proc] - 1;
					for(auto range_1 = remaining_indices.begin_intervals(); range_1 != remaining_indices.end_intervals(); ++range_1)
					{
						range_1_start = *(range_1->begin());
						range_1_end = range_1->last();
						if( !( (range_0_end < range_1_start) || (range_1_end < range_0_start) ) )
							dof_renumbering_offset.add_range(range_0_start, range_0_end, (int)current_dof_index - (int)dof_start_domain[proc]);

					}
					current_dof_index += n_locally_owned_dofs_per_processor_domain[proc];
				}

				if(n_locally_owned_dofs_per_processor_interface[proc])
				{
					range_0_start = dof_start_interface[proc];
					range_0_end = dof_start_interface[proc] + n_locally_owned_dofs_per_processor_interface[proc] - 1;
					for(auto range_1 = remaining_indices.begin_intervals(); range_1 != remaining_indices.end_intervals(); ++range_1)
					{
						range_1_start = *(range_1->begin());
						range_1_end = range_1->last();
						if( !( (range_0_end < range_1_start) || (range_1_end < range_0_start) ) )
							dof_renumbering_offset.add_range(range_0_start, range_0_end, (int)current_dof_index - (int)dof_start_interface[proc]);

					}
					current_dof_index += n_locally_owned_dofs_per_processor_interface[proc];
				}
			}
		}
	}
	else
	{
		dof_renumbering_offset.add_range(0, dof_handler_system.n_dofs() - 1, 0);
	}

}
#endif // DEAL_II_WITH_MPI


template<unsigned int spacedim>
void
Auxiliary::compute_map_dofs(const DoFHandlerSystem<spacedim>&	dhs_1,
							const DoFHandlerSystem<spacedim>&	dhs_2,
							vector<unsigned int>& 				map_dofs)
{
	Assert(dhs_1.n_dofs() == dhs_2.n_dofs(), ExcMessage("The number of dofs of the dof handler systems must match!"));
	map_dofs.resize(dhs_1.n_dofs(), 0);

	const IndexSet& locally_owned_indices = dhs_1.get_locally_owned_dofs();
	IndexSet identified_indices(locally_owned_indices.size());

	//use float precision here to eliminate problems with round-off
	map<tuple<float, float, float>, vector<unsigned int>> dofs_1_domain, dofs_2_domain, dofs_1_interface, dofs_2_interface;
	Point<3> p;
	vector<unsigned int> dof_indices;

	for(const auto& domain_cell : dhs_1.domain_active_iterators())
	{
		if( !domain_cell->is_artificial() )
		{
			for(unsigned int m = 0; m < spacedim; ++m)
				p[m] = domain_cell->center()[m];
			dof_indices.resize(domain_cell->get_fe().dofs_per_cell);
			domain_cell.get_dof_indices(dof_indices);
			dofs_1_domain[make_tuple(p[0], p[1], p[2])] = dof_indices;
		}
	}

	for(const auto& domain_cell : dhs_2.domain_active_iterators())
	{
		if( !domain_cell->is_artificial() )
		{
			for(unsigned int m = 0; m < spacedim; ++m)
				p[m] = domain_cell->center()[m];
			dof_indices.resize(domain_cell->get_fe().dofs_per_cell);
			domain_cell.get_dof_indices(dof_indices);
			dofs_2_domain[make_tuple(p[0], p[1], p[2])] = dof_indices;
		}
	}

	for(const auto& interface_cell_domain_cells : dhs_1.interface_active_iterators())
	{
		if( !interface_cell_domain_cells.interface_cell->is_artificial() )
		{
			for(unsigned int m = 0; m < spacedim; ++m)
				p[m] = interface_cell_domain_cells.interface_cell->center()[m];
			dof_indices.resize(interface_cell_domain_cells.interface_cell->get_fe().dofs_per_cell);
			interface_cell_domain_cells.get_dof_indices_local_global_interface(dof_indices);
			dofs_1_interface[make_tuple(p[0], p[1], p[2])] = dof_indices;
		}
	}

	for(const auto& interface_cell_domain_cells : dhs_2.interface_active_iterators())
	{
		if( !interface_cell_domain_cells.interface_cell->is_artificial() )
		{
			for(unsigned int m = 0; m < spacedim; ++m)
				p[m] = interface_cell_domain_cells.interface_cell->center()[m];
			dof_indices.resize(interface_cell_domain_cells.interface_cell->get_fe().dofs_per_cell);
			interface_cell_domain_cells.get_dof_indices_local_global_interface(dof_indices);
			dofs_2_interface[make_tuple(p[0], p[1], p[2])] = dof_indices;
		}
	}

	vector<unsigned int> dof_indices_C_1, dof_indices_C_2;
	dhs_1.get_dof_indices(dof_indices_C_1);
	dhs_2.get_dof_indices(dof_indices_C_2);

	for(const auto& cell_1 : dofs_1_domain)
	{
		const auto& dofs_1 = cell_1.second;
		const auto cell_2 = dofs_2_domain.find(cell_1.first);
		Assert(cell_2 != dofs_2_domain.end(), ExcMessage("Didn't find a corresponding cell in the second triangulation!"));
		const auto& dofs_2 = cell_2->second;
		Assert(dofs_1.size() == dofs_2.size(), ExcMessage("There are corresponding cells in the dof handler systems which don't have the same number of dofs defined on them."));
		for(unsigned int dof = 0; dof < dofs_1.size(); ++dof)
			if(locally_owned_indices.is_element(dofs_1[dof]))
			{
				map_dofs[dofs_1[dof]] = dofs_2[dof];
				identified_indices.add_index(dofs_1[dof]);
			}
	}

	for(const auto& cell_1 : dofs_1_interface)
	{
		const auto& dofs_1 = cell_1.second;
		const auto cell_2 = dofs_2_interface.find(cell_1.first);
		Assert(cell_2 != dofs_2_interface.end(), ExcMessage("Didn't find a corresponding cell in the second triangulation!"));
		const auto& dofs_2 = cell_2->second;
		Assert(dofs_1.size() == dofs_2.size(), ExcMessage("There are corresponding cells in the dof handler systems which don't have the same number of dofs defined on them."));
		for(unsigned int dof = 0; dof < dofs_1.size(); ++dof)
			if(locally_owned_indices.is_element(dofs_1[dof]))
			{
				map_dofs[dofs_1[dof]] = dofs_2[dof];
				identified_indices.add_index(dofs_1[dof]);
			}
	}

	Assert(dof_indices_C_1.size() == dof_indices_C_2.size(), ExcMessage("There is not the same number of constants in the dof handlers."));
	for(unsigned int dof = 0; dof < dof_indices_C_1.size(); ++dof)
		if(locally_owned_indices.is_element(dof_indices_C_1[dof]))
		{
			map_dofs[dof_indices_C_1[dof]] = dof_indices_C_2[dof];
			identified_indices.add_index(dof_indices_C_1[dof]);
		}

	identified_indices.compress();
	Assert(identified_indices == locally_owned_indices, ExcMessage("Could not identify correct amount of corresponding dofs!"));

	//communicate
#ifdef DEAL_II_WITH_MPI
	//add up contributions of different processors
	const dealii::parallel::Triangulation<spacedim, spacedim>* tria_domain_ptr;
	tria_domain_ptr = dynamic_cast<const dealii::parallel::Triangulation<spacedim, spacedim>*>(&(dhs_1.get_dof_handler_domain().get_triangulation()));
	if(tria_domain_ptr == nullptr)
		tria_domain_ptr = dynamic_cast<const dealii::parallel::Triangulation<spacedim, spacedim>*>(&(dhs_2.get_dof_handler_domain().get_triangulation()));
	if(tria_domain_ptr != nullptr)
	{
		int ierr = MPI_Allreduce(MPI_IN_PLACE, map_dofs.data(), map_dofs.size(), MPI_UNSIGNED, MPI_SUM, tria_domain_ptr->get_communicator());
		AssertThrowMPI(ierr);
	}
#endif //DEAL_II_WITH_MPI

}

void
Auxiliary::split_vector(const Vector<double>&	in,
						Vector<double>&			out_0,
						Vector<double>&			out_1,
						const unsigned int		size_1)
{
	const unsigned int size_0 = in.size() - size_1;
	out_0.reinit(size_0, true);
	out_1.reinit(size_1, true);
	for(unsigned int m = 0; m < size_0; ++m)
		out_0[m] = in[m];
	for(unsigned int m = size_0; m < in.size(); ++m)
		out_1[m - size_0] = in[m];
}

void
Auxiliary::split_matrix(const FullMatrix<double>&	in,
						FullMatrix<double>&			out_00,
						FullMatrix<double>&			out_01,
						FullMatrix<double>&			out_10,
						FullMatrix<double>&			out_11,
						const unsigned int			size_1)
{
	Assert(in.m() == in.n(), ExcMessage("This function does only work for square matrices!"));
	const unsigned int size_in = in.m();
	const unsigned int size_0 = size_in - size_1;
	out_00.reinit(size_0, size_0, true);
	out_01.reinit(size_0, size_1, true);
	out_10.reinit(size_1, size_0, true);
	out_11.reinit(size_1, size_1, true);
	for(unsigned int m = 0; m < size_0; ++m)
		for(unsigned int n = 0; n < size_0; ++n)
			out_00(m, n) = in(m, n);
	for(unsigned int m = 0; m < size_0; ++m)
		for(unsigned int n = size_0; n < size_in; ++n)
			out_01(m, n - size_0) = in(m, n);
	for(unsigned int m = size_0; m < size_in; ++m)
		for(unsigned int n = 0; n < size_0; ++n)
			out_10(m - size_0, n) = in(m, n);
	for(unsigned int m = size_0; m < size_in; ++m)
		for(unsigned int n = size_0; n < size_in; ++n)
			out_11(m - size_0, n - size_0) = in(m, n);
}

#ifdef DEAL_II_WITH_MPI
bool
Auxiliary::communicate_bool(const bool		local_bool,
							const MPI_Comm& mpi_communicator)
{
	//collect sizes
	int n_procs, this_proc;
	MPI_Comm_size(mpi_communicator, &n_procs);
	MPI_Comm_rank(mpi_communicator, &this_proc);

	if(n_procs > 1)
	{
		vector<int> local_bools(n_procs, 0);

		local_bools[this_proc] = local_bool;

		int send_value = local_bools[this_proc];
		int ierr =  MPI_Allgather(&send_value, 1, MPI_INT, local_bools.data(), 1, MPI_INT, mpi_communicator);
		AssertThrowMPI(ierr);
		for(const auto& local_bool_n : local_bools)
			if(local_bool_n)
				return true;
		return false;
	}
	else
		return local_bool;

}
#endif // DEAL_II_WITH_MPI


#ifdef DEAL_II_WITH_MPI
template
void
Auxiliary::compute_dof_renumbering_contiguous<2>(	const DoFHandlerSystem<2>&,
													DoFRenumberingOffset&);

template
void
Auxiliary::compute_dof_renumbering_contiguous<3>(	const DoFHandlerSystem<3>&,
													DoFRenumberingOffset&);
#endif // DEAL_II_WITH_MPI

template
void
Auxiliary::compute_map_dofs<2>(	const DoFHandlerSystem<2>&,
								const DoFHandlerSystem<2>&,
								vector<unsigned int>&);

template
void
Auxiliary::compute_map_dofs<3>(	const DoFHandlerSystem<3>&,
								const DoFHandlerSystem<3>&,
								vector<unsigned int>&);

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
