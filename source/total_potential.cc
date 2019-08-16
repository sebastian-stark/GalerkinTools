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

#include <galerkin_tools/total_potential.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
void
TotalPotential<spacedim>::add_total_potential_contribution(const TotalPotentialContribution<spacedim>& total_potential_contribution)
{
	total_potential_contributions.push_back(&total_potential_contribution);
	for(const auto& H_omega : total_potential_contribution.H_omega)
		if(H_omega->e_omega.size() > max_dependent_vars)
			max_dependent_vars=H_omega->e_omega.size();
	for(const auto& H_sigma : total_potential_contribution.H_sigma)
		if(H_sigma->e_sigma.size() > max_dependent_vars)
			max_dependent_vars=H_sigma->e_sigma.size();
}

template class TotalPotential<2>;
template class TotalPotential<3>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
