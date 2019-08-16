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

#include <galerkin_tools/two_block_sparsity_pattern.h>

#include <vector>

#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/affine_constraints.templates.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

TwoBlockSparsityPattern::TwoBlockSparsityPattern(	const IndexSet&		locally_relevant_indices,
													const unsigned int	block_0_size)
{
	reinit(locally_relevant_indices, block_0_size);
}

void
TwoBlockSparsityPattern::reinit(const IndexSet&		locally_relevant_indices,
								const unsigned int	block_0_size)
{
	finalized = false;
	SparsityPattern sp_dummy;
	sp_A.copy_from(sp_dummy);
	sp_B.copy_from(sp_dummy);
	sp_C.copy_from(sp_dummy);
	sp_D.copy_from(sp_dummy);

	this->block_0_size = block_0_size;
	block_1_size = locally_relevant_indices.size() - block_0_size;
	total_dimension = block_0_size + block_1_size;

	const IndexSet indices_block_0 = locally_relevant_indices.get_view(0, block_0_size);
	const IndexSet indices_block_1 = locally_relevant_indices.get_view(block_0_size, total_dimension);

	dsp_A.reinit(block_0_size, block_0_size, indices_block_0);
	dsp_B.reinit(block_0_size, block_1_size, indices_block_0);
	dsp_C.reinit(block_0_size, block_1_size, indices_block_0);
	dsp_D.reinit(block_1_size, block_1_size, indices_block_1);
}

template<unsigned int spacedim>
void
TwoBlockSparsityPattern::reinit(const AssemblyHelper<spacedim>& assembly_helper)
{
	const auto& dof_handler_system = assembly_helper.get_dof_handler_system();
	reinit(assembly_helper.get_locally_relevant_indices(), dof_handler_system.n_dofs_domain() + dof_handler_system.n_dofs_interface());
}

unsigned int
TwoBlockSparsityPattern::n_rows()
const
{
	return total_dimension;
}

unsigned int
TwoBlockSparsityPattern::n_cols()
const
{
	return total_dimension;
}

void
TwoBlockSparsityPattern::add(	const unsigned int i,
								const unsigned int j)
{
	Assert(!finalized, ExcMessage("The TwoBlockSparsityPattern cannot be changed after finalize() has been called"));

	if(i < block_0_size)
	{
		if(j < block_0_size)
			dsp_A.add(i, j);
		else
			dsp_B.add(i, j - block_0_size);
	}
	else
	{
		if(j < block_0_size)
			dsp_C.add(j, i - block_0_size);
		else
			dsp_D.add(i - block_0_size, j - block_0_size);
	}
}

template <typename ForwardIterator>
void
TwoBlockSparsityPattern::add_entries(	const unsigned int	row,
										ForwardIterator		begin,
										ForwardIterator 	end,
										const bool      	indices_are_sorted)
{
	Assert(!finalized, ExcMessage("The TwoBlockSparsityPattern cannot be changed after finalize() has been called"));

	(void)indices_are_sorted;
	Assert(indices_are_sorted, ExcMessage("The method add_entries() is not yet implemented for the case that the indices are not sorted!"));

	//search for begin iterator of second block and set up the column indices for the second block
	ForwardIterator first_block_end = end;
	static vector<unsigned int> second_block;
	second_block.resize(0);
	for(auto it = begin; it != end; ++it)
	{
		if(*it >= block_0_size)
		{
			first_block_end = it;
			for(auto it = first_block_end; it != end; ++it)
				second_block.push_back(*it - block_0_size);
			break;
		}
	}

	if( row < block_0_size )
	{
		dsp_A.add_entries(row, begin, first_block_end, true);
		dsp_B.add_entries(row, second_block.begin(), second_block.end(), true);
	}
	else
	{
		for(auto it = begin; it != first_block_end; ++it)
			dsp_C.add(*it, row - block_0_size);
		dsp_D.add_entries(row - block_0_size, second_block.begin(), second_block.end());
	}
}

void
TwoBlockSparsityPattern::finalize()
{
	sp_A.copy_from(dsp_A);
	sp_B.copy_from(dsp_B);
	sp_C.copy_from(dsp_C);
	sp_D.copy_from(dsp_D);
	finalized = true;
}

const SparsityPattern&
TwoBlockSparsityPattern::get_sp_A()
const
{
	Assert(finalized, ExcMessage("You can only ask for the sp_A after finalize() has been called!"));

	return sp_A;
}

const SparsityPattern&
TwoBlockSparsityPattern::get_sp_B()
const
{
	Assert(finalized, ExcMessage("You can only ask for the sp_B after finalize() has been called!"));

	return sp_B;
}

const SparsityPattern&
TwoBlockSparsityPattern::get_sp_C()
const
{
	Assert(finalized, ExcMessage("You can only ask for the sp_C after finalize() has been called!"));

	return sp_C;
}

const SparsityPattern&
TwoBlockSparsityPattern::get_sp_D()
const
{
	Assert(finalized, ExcMessage("You can only ask for the sp_D after finalize() has been called!"));

	return sp_D;
}

bool
TwoBlockSparsityPattern::exists(const unsigned int i,
								const unsigned int j)
const
{
	if(!finalized)
	{
		if(i < block_0_size)
		{
			if(j < block_0_size)
				return dsp_A.exists(i, j);
			else
				return dsp_B.exists(i, j - block_0_size);
		}
		else
		{
			if(j < block_0_size)
				return dsp_C.exists(j, i - block_0_size);
			else
				return dsp_D.exists(i - block_0_size, j - block_0_size);
		}
	}
	else
	{
		if(i < block_0_size)
		{
			if(j < block_0_size)
				return sp_A.exists(i, j);
			else
				return sp_B.exists(i, j - block_0_size);
		}
		else
		{
			if(j < block_0_size)
				return sp_C.exists(j, i - block_0_size);
			else
				return sp_D.exists(i - block_0_size, j - block_0_size);
		}
	}
}

#ifdef DEAL_II_WITH_MPI

void
TwoBlockSparsityPattern::distribute(const IndexSet& locally_owned_indices, const MPI_Comm mpi_communicator)
{
	Assert(!finalized, ExcMessage("distribute() can only be called before finalization of the TwoBlockSparsityPattern"));

	//collect sizes
	int n_procs, this_proc;
	MPI_Comm_size(mpi_communicator, &n_procs);
	MPI_Comm_rank(mpi_communicator, &this_proc);

	if(n_procs > 1)
	{
		vector<unsigned int> n_dofs_per_processor_0(n_procs, 0), n_dofs_per_processor_1(n_procs, 0);

		n_dofs_per_processor_0[this_proc] = locally_owned_indices.get_view(0, block_0_size).n_elements();
		n_dofs_per_processor_1[this_proc] = locally_owned_indices.get_view(block_0_size, total_dimension).n_elements();

		unsigned int send_value = n_dofs_per_processor_0[this_proc];
		int ierr =  MPI_Allgather(&send_value, 1, MPI_UNSIGNED, n_dofs_per_processor_0.data(), 1, MPI_UNSIGNED, mpi_communicator);
		AssertThrowMPI(ierr);

		send_value = n_dofs_per_processor_1[this_proc];
		ierr =  MPI_Allgather(&send_value, 1, MPI_UNSIGNED, n_dofs_per_processor_1.data(), 1, MPI_UNSIGNED, mpi_communicator);
		AssertThrowMPI(ierr);

		SparsityTools::distribute_sparsity_pattern(dsp_A, n_dofs_per_processor_0, mpi_communicator, dsp_A.row_index_set());
		SparsityTools::distribute_sparsity_pattern(dsp_B, n_dofs_per_processor_0, mpi_communicator, dsp_B.row_index_set());
		SparsityTools::distribute_sparsity_pattern(dsp_C, n_dofs_per_processor_0, mpi_communicator, dsp_C.row_index_set());
		SparsityTools::distribute_sparsity_pattern(dsp_D, n_dofs_per_processor_1, mpi_communicator, dsp_D.row_index_set());
	}

}

#endif // DEAL_II_WITH_MPI

template
void
TwoBlockSparsityPattern::reinit<2>(const AssemblyHelper<2>&);

template
void
TwoBlockSparsityPattern::reinit<3>(const AssemblyHelper<3>&);

DEAL_II_NAMESPACE_CLOSE
GALERKIN_TOOLS_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_OPEN

template
void
AffineConstraints<double>::add_entries_local_to_global<GalerkinTools::TwoBlockSparsityPattern>(	const vector<AffineConstraints<double>::size_type>&,
																								GalerkinTools::TwoBlockSparsityPattern&,
																								const bool,
																								const Table<2, bool>&,
																								integral_constant<bool, false>)
const;

template
void
AffineConstraints<double>::add_entries_local_to_global<GalerkinTools::TwoBlockSparsityPattern>(	const vector<AffineConstraints<double>::size_type>&,
																								const vector<AffineConstraints<double>::size_type>&,
																								GalerkinTools::TwoBlockSparsityPattern&,
																								const bool,
																								const Table<2, bool>&)
const;

DEAL_II_NAMESPACE_CLOSE
