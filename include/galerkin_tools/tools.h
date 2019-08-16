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

#ifndef GALERKINTOOLS_TOOLS_H_
#define GALERKINTOOLS_TOOLS_H_

#include <vector>

#include <deal.II/lac/la_parallel_vector.h>

#include <galerkin_tools/config.h>
#include <galerkin_tools/dof_renumbering.h>
#include <galerkin_tools/dof_handler_system.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

namespace Auxiliary
{

/**
 * %Auxiliary method converting local dof indices into global ones
 *
 * @param[in]	dof_indices_local			local dof indices
 * @param[out]	dof_indices_global			global dof indices
 * @param[in]	dof_indices_local_global	mapping between local and global dof indices
 */
void
convert_local_indices_to_global_indices(	const std::vector<unsigned int>&	dof_indices_local,
											std::vector<unsigned int>&			dof_indices_global,
											const std::vector<unsigned int>&	dof_indices_local_global);

/**
 * Auxiliary method to get a combined dof indexing at an interface
 * cell (which comprises dof indices of the interface cell, the adjacent domain cells, and the relevant independent scalars).
 *
 * This method eliminates duplicate dof's on the + and - side of the interface (which may arise
 * in the case that an independent field is continuous across the interface). This is done in order to eliminate
 * the case that two local dof indices correspond to the same global dof index (some deal.II functions
 * would not work for this case).
 *
 * @param[in]	dof_indices_global_interface				mapping between local dof indices and global ones for interface cell
 *
 * @param[in]	dof_indices_global_minus					mapping between local dof indices and global ones for domain cell on - side
 *
 * @param[in]	dof_indices_global_plus						mapping between local dof indices and global ones for domain cell on + side
 *
 * @param[in]	dof_indices_global_C						mapping between local dof indices and global ones for independent scalars
 *
 * @param[out]	dof_indices_interface_dof_indices_combined	defined such that
 * 															@p dof_indices_global_combined[@p dof_indices_interface_dof_indices_combined[@p i]]=@p dof_indices_global_interface[@p i];
 * 															this does relate the local dof indexing of deal.II to the combined local dof indexing
 *
 * @param[out]	dof_indices_minus_dof_indices_combined		defined such that
 * 															@p dof_indices_global_combined[@p dof_indices_minus_dof_indices_combined[@p i]]=@p dof_indices_global_minus[@p i];
 * 															this does relate the local dof indexing of deal.II to the combined local dof indexing
 *
 * @param[out]	dof_indices_plus_dof_indices_combined		defined such that
 * 															@p dof_indices_global_combined[@p dof_indices_plus_dof_indices_combined[@p i]]=@p dof_indices_global_plus[@p i];
 * 															this does relate the local dof indexing of deal.II to the combined local dof indexing
 *
 * @param[out]	dof_indices_C_dof_indices_combined			defined such that
 * 															@p dof_indices_global_combined[@p dof_indices_C_dof_indices_combined[@p i]]=@p dof_indices_global_C[@p i];
 * 															this does relate the local dof indexing of the C's to the combined local dof indexing
 *
 * @param[out]	dof_indices_global_combined					mapping between combined local dof indices and global dof indices (the important property is that
 * 															there are no duplicate global dof indices in this vector)
 */
void
combine_dof_indices(	const std::vector<unsigned int>&	dof_indices_global_interface,
						const std::vector<unsigned int>&	dof_indices_global_minus,
						const std::vector<unsigned int>&	dof_indices_global_plus,
						const std::vector<unsigned int>&	dof_indices_global_C,
						std::vector<unsigned int>&			dof_indices_interface_dof_indices_combined,
						std::vector<unsigned int>&			dof_indices_minus_dof_indices_combined,
						std::vector<unsigned int>&			dof_indices_plus_dof_indices_combined,
						std::vector<unsigned int>&			dof_indices_C_dof_indices_combined,
						std::vector<unsigned int>&			dof_indices_global_combined);

/**
 * This function takes the indices from the index set @p in, renumbers them according to the DoFRenumbering object provided,
 * and adds them to the IndexSet @p out.
 * Of course, this requires that the size of the IndexSet @p out is such that it can accommodate the renumbered indices.
 *
 * It is advantageous if the incoming index set is compressed.
 *
 * @param[in]		dof_renumbering		The DoFRenumbering object
 *
 * @param[in]		in					The indices to be renumbered
 *
 * @param[inout]	out					The index set to which the renumbered indices are added (this index set need not be empty)
 */
void
add_to_index_set(	const DoFRenumbering& dof_renumbering,
					const IndexSet& in,
					IndexSet& out);

/**
 * %Function renumbering constraints according to the renumbering scheme provided by @p dof_renumbering.
 *
 * @param[inout]	constraint_matrix	The constraint matrix to be renumbered
 *
 * @param[in]		dof_renumbering		The renumbering scheme to be applied (the default doesn't apply a renumbering)
 *
 * @param[in]		close				If @p true, the constraint matrix will be closed after renumbering
 */
void
renumber_constraints(	AffineConstraints<double>&	constraint_matrix,
						const DoFRenumbering&		dof_renumbering = DoFRenumbering(),
						const bool					close = true);

/**
 * Return a DoFRenumberingOffset object for a DoFHandlerSystem such that the renumbered dofs are contiguous on each processor. Of course, the renumbering
 * has to be recomputed after every distribution of dofs. The implementation currently requires that the domain related and interface related dof handlers
 * of the DoFHandlerSystem are both associated with a contiguous dof numbering on each processor.
 *
 * @param[in]	dof_handler_system		The DoFHandlerSystem for which the renumbering is to be done
 *
 * @param[out]	dof_renumbering_offset	The resulting DoFRenumberingOffset object
 */
template<unsigned int spacedim>
void
compute_dof_renumbering_contiguous(	const DoFHandlerSystem<spacedim>&	dof_handler_system,
									DoFRenumberingOffset&				dof_renumbering_offset);

/**
 * Compute a map between the dofs of @p dhs_1 and the dofs of @p dhs_2. This function requires that the triangulation underlying the dof handlers is based on the same mesh, that
 * the partitioning of the triangulations matches, and that the same finite elements are used.
 *
 * The function exists mainly for debug purposes (in order to be able to compare results obtained with a parallel computation with those of an equivalent sequential computation).
 *
 * The function identifies corresponding cells by equal positions of the cell centers. Therefore, there shouldn't be overlapping cells in the mesh!
 *
 * This function does not scale because the entire map is computed on every processor.
 *
 * @param[in]	dhs_1		The first dof handler system
 *
 * @param[in]	dhs_2		The second dof handler system
 *
 * @param[out]	map_dofs	The map between the dofs (@p map_dofs[i] contains the index of @p dhs_2 corresponding to index @p i of @p dhs_1)
 */
template<unsigned int spacedim>
void
compute_map_dofs(	const DoFHandlerSystem<spacedim>&	dhs_1,
					const DoFHandlerSystem<spacedim>&	dhs_2,
					std::vector<unsigned int>&			map_dofs);

/**
 * %Function to split a vector into two blocks
 *
 * @param[in]	in		The vector to be splitted
 *
 * @param[out]	out_0	The first part of the split
 *
 * @param[out]	out_1	The second part of the split
 *
 * @param[in]	size_1	The size of @p out_1
 */
void
split_vector(	const Vector<double>&	in,
				Vector<double>&			out_0,
				Vector<double>&			out_1,
				const unsigned int		size_1);

/**
 * %Function to split a square matrix into four blocks
 *
 * @param[in]	in		The matrix to be splitted
 *
 * @param[out]	out_00	The top left part of the split
 *
 * @param[out]	out_01	The top right part of the split
 *
 * @param[out]	out_10	The bottom left part of the split
 *
 * @param[out]	out_11	The bottom right part of the split
 *
 * @param[in]	size_1	The dimension of @p out_11
 */
void
split_matrix(	const FullMatrix<double>&	in,
				FullMatrix<double>&			out_00,
				FullMatrix<double>&			out_01,
				FullMatrix<double>&			out_10,
				FullMatrix<double>&			out_11,
				const unsigned int			size_1);

#ifdef DEAL_II_WITH_MPI

/**
 * %Auxiliary function communicating a bool such that the result is @p true on all processors if the bool is @p true on one of the processors.
 *
 * @param[in]	local_bool			The local value
 *
 * @param[in]	mpi_communicator	The MPI communicator
 *
 * @return						The global bool
 */
bool
communicate_bool(	const bool		local_bool,
					const MPI_Comm& mpi_communicator);

#endif //DEAL_II_WITH_MPI

}


GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_TOOLS_H_ */
