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

#ifndef GALERKIN_TOOLS_TWO_BLOCK_SPARSITY_PATTERN_H_
#define GALERKIN_TOOLS_TWO_BLOCK_SPARSITY_PATTERN_H_


#include <vector>
#include <iostream>

#include <deal.II/base/index_set.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#ifdef DEAL_II_WITH_MPI
#include <mpi/mpi.h>
#endif

#include <galerkin_tools/config.h>
#include <galerkin_tools/assembly_helper.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * This class implements a (dynamic) sparsity pattern corresponding to TwoBlockMatrix
 */
class TwoBlockSparsityPattern
{

private:

	/**
	 * The dynamic sparsity pattern for A
	 */
	DynamicSparsityPattern
	dsp_A;

	/**
	 * The dynamic sparsity pattern for B
	 */
	DynamicSparsityPattern
	dsp_B;

	/**
	 * The dynamic sparsity pattern for C
	 */
	DynamicSparsityPattern
	dsp_C;

	/**
	 * The dynamic sparsity pattern for D
	 */
	DynamicSparsityPattern
	dsp_D;

	/**
	 * The sparsity pattern for A
	 */
	SparsityPattern
	sp_A;

	/**
	 * The sparsity pattern for B
	 */
	SparsityPattern
	sp_B;

	/**
	 * The sparsity pattern for B
	 */
	SparsityPattern
	sp_C;

	/**
	 * The sparsity pattern for B
	 */
	SparsityPattern
	sp_D;

	/**
	 * The dimension of A
	 */
	unsigned int
	block_0_size = 0;

	/**
	 * The dimension of D
	 */
	unsigned int
	block_1_size = 0;

	/**
	 * The total dimension
	 */
	unsigned int
	total_dimension = 0;

	/**
	 * bool indicating whether the sparsity pattern is finalized (after finalization not changes are possible)
	 */
	bool
	finalized = false;

public:

	/**
	 * Default Constructor
	 */
	TwoBlockSparsityPattern() = default;

	/**
	 * Constructor
	 *
	 * @param[in]	locally_relevant_indices	The locally relevant indices
	 *
	 * @param[in]	block_0_size				The dimension of A (the dimension of B is implicitly given by @p locally_relevant_indices and @p block_0_size)
	 *
	 */
	TwoBlockSparsityPattern(const IndexSet&		locally_relevant_indices,
							const unsigned int	block_0_size);

	/**
	 * Virtual destructor to make the class polymorphic and allow dynamic casts.
	 */
	virtual
	~TwoBlockSparsityPattern() = default;

	/**
	 * Reinitialize the TwoBlockSparsityPattern
	 *
	 * @param[in]	locally_relevant_indices	The locally relevant indices
	 *
	 * @param[in]	block_0_size				The dimension of A (the dimension of B is implicitly given by @p locally_relevant_indices and @p block_0_size)
	 *
	 */
	void
	reinit(	const IndexSet&		locally_relevant_indices,
			const unsigned int	block_0_size);

	/**
	 * Reinitialize the TwoBlockSparsityPattern in a way that it is compatible with the functions of an AssemblyHelper.
	 * The first block will contain the dofs related to finite elements, and the second block is related to the independent
	 * scalars and the additional rows due to stretching of the system matrix. This choice ensures that the operations
	 * of AssemblyHelper scale in parallel.
	 *
	 * @param[in]	assembly_helper		The AssemblyHelper object
	 */
	template<unsigned int spacedim>
	void
	reinit(const AssemblyHelper<spacedim>& assembly_helper);

	/**
	 * The total number of rows
	 */
	unsigned int
	n_rows()
	const;

	/**
	 * The total number of columns
	 */
	unsigned int
	n_cols()
	const;

	/**
	 * %Function adding an entry to the sparsity pattern
	 *
	 * @param[in]	i	row
	 *
	 * @param[in]	j	column
	 *
	 */
	void
	add(const unsigned int i,
		const unsigned int j);

	/**
	 * Add a number of entries to a row of the sparsity pattern
	 *
	 * @param[in]	row					row to add to
	 *
	 * @param[in]	begin				iterator to first column index to add
	 *
	 * @param[in]	end					iterator past the last column index to add
	 *
	 * @param[in]	indices_are_sorted	This only exists for compatibility with deal.II and does not affect what's going on in this function.
	 * 									The function only works if the indices are actually sorted.
	 */
	template <typename ForwardIterator>
	void
	add_entries(const unsigned int	row,
				ForwardIterator 	begin,
				ForwardIterator 	end,
				const bool      	indices_are_sorted = true);

	/**
	 * This functions copies the dynamic sparsity patterns generated into the sparsity pattern. It does also clear
	 * the dynamic sparsity patterns. After finalize is called, no changes to the
	 * sparsity pattern are allowed, unless reinit() is called again.
	 */
	void
	finalize();

	/**
	 * @return		The sparsity pattern of A
	 */
	const SparsityPattern&
	get_sp_A()
	const;

	/**
	 * @return		The sparsity pattern of B
	 */
	const SparsityPattern&
	get_sp_B()
	const;

	/**
	 * @return		The sparsity pattern of C
	 */
	const SparsityPattern&
	get_sp_C()
	const;

	/**
	 * @return		The sparsity pattern of D
	 */
	const SparsityPattern&
	get_sp_D()
	const;

	/**
	 * Check wheter an entry exists in the sparsity pattern
	 *
	 * @param[in]	i		row
	 *
	 * @param[in]	j		column
	 */
	bool
	exists(	const unsigned int i,
			const unsigned int j)
	const;

#ifdef DEAL_II_WITH_MPI
	/**
	 * This function makes sure that all processors know about all their entries
	 *
	 * @param[in]	locally_owned_indices	The locally owned indices on this processor
	 *
	 * @param[in]	mpi_communicator		The MPI communicator
	 */
	void
	distribute(	const IndexSet&	locally_owned_indices,
				const MPI_Comm	mpi_communicator);

#endif //DEAL_II_WITH_MPI

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKIN_TOOLS_TWO_BLOCK_SPARSITY_PATTERN_H_ */
