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

#include <galerkin_tools/total_potential_contribution.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
TotalPotentialContribution<spacedim>::TotalPotentialContribution(	const vector<const ScalarFunctional<spacedim ,spacedim>*>&		H_omega,
																	const vector<const ScalarFunctional<spacedim-1, spacedim>*>&	H_sigma,
																	const vector<const IndependentField<0, spacedim>*>&				C)
:
is_primitive(false),
H_omega(H_omega),
H_sigma(H_sigma),
C(C)
{
	Assert(	H_omega.size() + H_sigma.size() + C.size() > 0,
			ExcMessage("You are trying to construct a TotalPotentialContribution which does not depend on any scalar functionals or independent scalars!"));
}

template<unsigned int spacedim>
TotalPotentialContribution<spacedim>::TotalPotentialContribution(const ScalarFunctional<spacedim, spacedim>& H_omega)
:
is_primitive(true),
H_omega({&H_omega})
{
}

template<unsigned int spacedim>
TotalPotentialContribution<spacedim>::TotalPotentialContribution(const ScalarFunctional<spacedim-1, spacedim>& H_sigma)
:
is_primitive(true),
H_sigma({&H_sigma})
{
}

template<unsigned int spacedim>
TotalPotentialContribution<spacedim>::~TotalPotentialContribution()
{
	Assert(	n_subscriptions()==0,
			ExcMessage("You are about to destroy a TotalPotentialContribution, which is currently in use! Make sure that all TotalPotentialContribution objects live at least as long as the objects using them!"));
}

template<unsigned int spacedim>
bool
TotalPotentialContribution<spacedim>::get_potential_contribution(	const Vector<double>& 			/*H_HS_C*/,
																	const vector<Vector<double>>& 	/*C_ref*/,
																	double& 						/*F*/,
																	Vector<double>& 				/*F_1*/,
																	FullMatrix<double>& 			/*F_2*/,
																	const tuple<bool,bool,bool>& 	/*requestedQuantities*/)
const
{
	Assert(false, ExcMessage("The function get_potential_contribution() of the base class TotalPotentialContribution should not be called under any circumstances. "
							 "Did you accidentally call it for a TotalPotentialContribution associated with is_primitive == true or did you forget to"
							 "overwrite the function in a derived class?"));
	return true;
}

template class TotalPotentialContribution<2>;
template class TotalPotentialContribution<3>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
