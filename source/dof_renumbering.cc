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

#include <galerkin_tools/dof_renumbering.h>

#include <iostream>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

using namespace std;

DoFRenumbering::~DoFRenumbering()
{
	Assert(n_subscriptions() == 0, ExcMessage("You are about to destroy an DoFRenumbering, which is currently in use! Make sure that all DoFRenumbering objects live at least as long as the objects using them!"));
}

void
DoFRenumbering::convert_dof_indices(std::vector<unsigned int>& /*dof_indices*/)
const
{
}

vector<pair<const unsigned int, const unsigned int>>
DoFRenumbering::convert_range(	const unsigned int range_begin,
								const unsigned int range_end)
const
{
	vector<pair<const unsigned int, const unsigned int>> returned_ranges;
	returned_ranges.push_back(make_pair(range_begin, range_end));
	return returned_ranges;
}

void
DoFRenumberingOffset::add_range(const unsigned int dof_start, const unsigned int dof_end, const int offset)
{
	dof_offsets.push_back(make_tuple(dof_start, dof_end, offset));
}

void
DoFRenumberingOffset::convert_dof_indices(std::vector<unsigned int>& dof_indices)
const
{
	for(auto& dof_index : dof_indices)
	{
		for(unsigned int m = 0; m < dof_offsets.size(); ++m)
		{
			if( (dof_index >= get<0>(dof_offsets[m])) && (dof_index <= get<1>(dof_offsets[m])) )
			{
				dof_index += get<2>(dof_offsets[m]);
				break;
			}
			Assert( (m < dof_offsets.size() - 1) , ExcMessage("No rule for renumbering this dof has been found. Most likely, the rules for the renumbering are inconsistent."));
		}
	}
}

vector<pair<const unsigned int, const unsigned int>>
DoFRenumberingOffset::convert_range(	const unsigned int range_begin,
										const unsigned int range_end)
const
{
	vector<pair<const unsigned int, const unsigned int>> returned_ranges;

	unsigned int current_range_begin = range_begin;
	for(;;)
	{
		for(unsigned int m = 0; m < dof_offsets.size(); ++m)
		{
			if( (current_range_begin >= get<0>(dof_offsets[m])) && (current_range_begin <= get<1>(dof_offsets[m])) )
			{
				const unsigned int converted_range_begin = current_range_begin + get<2>(dof_offsets[m]);
				const unsigned int converted_range_end = (range_end <= get<1>(dof_offsets[m])) ? range_end + get<2>(dof_offsets[m]) : get<1>(dof_offsets[m]) + get<2>(dof_offsets[m]);
				returned_ranges.push_back(make_pair(converted_range_begin, converted_range_end));
				current_range_begin = get<1>(dof_offsets[m]) + 1;
				break;
			}
			Assert( (m < dof_offsets.size() - 1) , ExcMessage("No rule for renumbering this dof has been found. Most likely, the rules for the renumbering are inconsistent."));
		}
		if(current_range_begin > range_end)
			break;
	}

	return returned_ranges;
}

const std::vector<std::tuple<const unsigned int, const unsigned int, const int>>&
DoFRenumberingOffset::get_dof_offsets()
const
{
	return dof_offsets;
}

void
DoFRenumberingOffset::clear()
{
	dof_offsets.clear();
}

void
DoFRenumberingOffset::print(ostream &out)
const
{
	for(const auto& dof_offset : dof_offsets)
	{
		out << "[" << get<0>(dof_offset) << ", " << get<1>(dof_offset) << "]";
		out << " --> ";
		out << "[" << get<0>(dof_offset) + get<2>(dof_offset) << ", " << get<1>(dof_offset) + get<2>(dof_offset) << "]";
		out << endl;
	}

}


GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

