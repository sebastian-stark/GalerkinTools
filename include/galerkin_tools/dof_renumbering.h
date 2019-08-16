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

#ifndef INCLUDE_GALERKIN_TOOLS_DOF_RENUMBERING_H_
#define INCLUDE_GALERKIN_TOOLS_DOF_RENUMBERING_H_

#include <vector>
#include <tuple>

#include <deal.II/base/index_set.h>
#include <deal.II/base/subscriptor.h>

#include <galerkin_tools/config.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * This class is intended to describe a renumbering of dof indices
 *
 * The class inherits from Subscriptor in order to be
 * able to check that DoFRenumbering objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * The default implementation does no renumbering.
 */
class DoFRenumbering : public Subscriptor
{
public:

	/**
	 * This function renumbers the incoming dof indices
	 *
	 * @param[inout]	dof_indices		The vector with the dof indices to be renumbered
	 */
	virtual
	void
	convert_dof_indices(std::vector<unsigned int>& dof_indices)
	const;

	/**
	 * Converts a contiguous range [@p range_begin, @p range_end] of indices into one or several contiguous ranges in the new numbering
	 *
	 * @param[in]	range_begin		First index
	 *
	 * @param[in]	range_end		Past the last index
	 *
	 * @return	A vector with one or several ranges.
	 */
	virtual
	std::vector<std::pair<const unsigned int, const unsigned int>>
	convert_range(	const unsigned int range_begin,
					const unsigned int range_end)
	const;

	/**
	 * The destructor of DoFRenumbering essentially checks before destruction that the
	 * DoFRenumbering object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	virtual
	~DoFRenumbering();

};

/**
 * This renumbering scheme is based on applying offsets to index ranges.
 * Presently no measures are taken to ensure that the overall index set is consistent!
 */
class DoFRenumberingOffset : public DoFRenumbering
{
private:

	/**
	 * This data structure defines a fixed shift of dof numbers from a certain interval [dof_start_i, dof_end_i] to
	 * the interval [dof_start_i + offset_i, dof_end_i + offset_i]. The tuples (dof_start_i, dof_end_i, offset_i) define this shift.
	 * Multiple intervals can be added and the order of adding the intervals influences how the intervals are searched during conversion
	 * of dof indices. In particular, if a certain dof index is to be converted, it is initially tried to find it in the interval added first,
	 * if it is not found there, the second interval is tried, and so forth. @todo It would be worthwhile to implement a more effective
	 * conversion scheme.
	 */
	std::vector<std::tuple<const unsigned int, const unsigned int, const int>>
	dof_offsets;

public:

	/**
	 * Add an element to DoFRenumberingOffset::dof_offsets
	 *
	 * @param[in]	dof_start	The first index of the original interval
	 *
	 * @param[in]	dof_end		The last index of the original interval (or, rather, the index past the last index)
	 *
	 * @param[in]	offset		The offset to be applied
	 */
	void
	add_range(const unsigned int dof_start, const unsigned int dof_end, const int offset);

	/**
	 * This function renumbers the incoming dof indices
	 *
	 * @param[inout]	dof_indices		The vector with the dof indices to be renumbered
	 */
	virtual
	void
	convert_dof_indices(std::vector<unsigned int>& dof_indices)
	const;

	/**
	 * Converts a contiguous range [@p range_begin, @p range_end] of indices into one or several contiguous ranges in the new numbering
	 *
	 * @param[in]	range_begin		First index
	 *
	 * @param[in]	range_end		Past the last index
	 *
	 * @return	A vector with one or several ranges.
	 */
	virtual
	std::vector<std::pair<const unsigned int, const unsigned int>>
	convert_range(	const unsigned int range_begin,
					const unsigned int range_end)
	const;

	/**
	 * @return DoFRenumberingOffset::dof_offsets
	 */
	const std::vector<std::tuple<const unsigned int, const unsigned int, const int>>&
	get_dof_offsets()
	const;

	/**
	 * %Function resetting the object
	 */
	void
	clear();

	/**
	 * print the renumbering scheme to screen
	 *
	 * @param[in]	out		The output stream
	 */
	void
	print(std::ostream &out)
	const;

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* INCLUDE_GALERKIN_TOOLS_DOF_RENUMBERING_H_ */
