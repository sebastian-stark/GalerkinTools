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

#include <galerkin_tools/linear_material.h>

#include <iostream>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
bool
LinearMaterialDomain<spacedim>::get_h_omega(Vector<double>& 				e_omega,
											const vector<Vector<double>>&	/*e_omega_ref_sets*/,
											Vector<double>&					/*hidden_vars*/,
											const Point<spacedim>&			/*x*/,
											double&							h_omega,
											Vector<double>&					h_omega_1,
											FullMatrix<double>&				h_omega_2,
											const tuple<bool, bool, bool>	requested_quantities)
const
{
	Assert(e_omega.size()==this->e_omega.size(),ExcMessage("Called get_h_omega with invalid size of e_omega vector!"));

	Vector<double> C_e_omega(this->e_omega.size());
	C.vmult(C_e_omega,e_omega);

	if(get<0>(requested_quantities) )
	{
		h_omega=0.5*(e_omega*C_e_omega);
		h_omega+=e_omega*y;
	}

	if(get<1>(requested_quantities))
	{
		if(h_omega_1.size()!=this->e_omega.size())
			h_omega_1.reinit(this->e_omega.size());
		h_omega_1=y;
		h_omega_1+=C_e_omega;
	}
	if(get<2>(requested_quantities))
	{
		if( (h_omega_2.size()[0]!=this->e_omega.size()) || (h_omega_2.size()[1]!=this->e_omega.size()) )
			h_omega_2.reinit(this->e_omega.size(), this->e_omega.size());
		h_omega_2=C;
	}
	return false;
}

template<unsigned int spacedim>
LinearMaterialDomain<spacedim>::LinearMaterialDomain(	const vector<GalerkinTools::DependentField<spacedim,spacedim>>	e_omega,
														const set<types::material_id>									domain_of_integration,
														const Quadrature<spacedim>										quadrature,
														const FullMatrix<double>										C,
														const Vector<double>											y,
														const string													name)
:
ScalarFunctional<spacedim, spacedim>(e_omega, domain_of_integration, quadrature, name, 0),
C(C),
y(y)
{
	Assert(	(C.size()[0]==e_omega.size()) && (C.size()[1]==e_omega.size()),
			ExcMessage("Matrix of linear material must be square and have same the dimension as the number of dependent fields!") );
	Assert(	(y.size()==e_omega.size()),
			ExcMessage("Vector of linear material must have same dimension as number of dependent fields!") );
}

template<unsigned int spacedim>
bool
LinearMaterialInterface<spacedim>::get_h_sigma(	Vector<double>& 				e_sigma,
												const vector<Vector<double>>&	/*e_sigma_ref_sets*/,
												Vector<double>& 				/*hidden_vars*/,
												const Point<spacedim>& 			/*x*/,
												const Tensor<1,spacedim>& 		/*n*/,
												double& 						h_sigma,
												Vector<double>& 				h_sigma_1,
												FullMatrix<double>& 			h_sigma_2,
												const tuple<bool, bool, bool>	requested_quantities)
const
{
	Assert(e_sigma.size()==this->e_sigma.size(),ExcMessage("Called get_h_sigma with invalid size of e vector!"));

	Vector<double> C_e_sigma(this->e_sigma.size());
		C.vmult(C_e_sigma,e_sigma);

	if(get<0>(requested_quantities))
	{
		h_sigma=0.5*(e_sigma*C_e_sigma);
		h_sigma+=e_sigma*y;
	}

	if(get<1>(requested_quantities))
	{
		if(h_sigma_1.size()!=this->e_sigma.size())
			h_sigma_1.reinit(this->e_sigma.size());
		h_sigma_1=y;
		h_sigma_1+=C_e_sigma;
	}
	if(get<2>(requested_quantities))
	{
		if( (h_sigma_2.size()[0]!=this->e_sigma.size()) || (h_sigma_2.size()[1]!=this->e_sigma.size()) )
			h_sigma_2.reinit(this->e_sigma.size(), this->e_sigma.size());
		h_sigma_2=C;
	}
	return false;
}

template<unsigned int spacedim>
LinearMaterialInterface<spacedim>::LinearMaterialInterface(	const vector<GalerkinTools::DependentField<spacedim-1,spacedim>>	e_sigma,
															const set<types::material_id>										domain_of_integration,
															const Quadrature<spacedim-1>										quadrature,
															const FullMatrix<double>											C,
															const Vector<double>												y,
															const string														name)
:
ScalarFunctional<spacedim-1, spacedim>(e_sigma, domain_of_integration, quadrature, name, 0),
C(C),
y(y)
{
	Assert(	(C.size()[0]==e_sigma.size()) && (C.size()[1]==e_sigma.size()),
			ExcMessage("Matrix of linear material must be square and have same the dimension as the number of dependent fields!") );
	Assert(	(y.size()==e_sigma.size()),
			ExcMessage("Vector of linear material must have same dimension as number of dependent fields!") );
}

template class LinearMaterialDomain<2>;
template class LinearMaterialDomain<3>;
template class LinearMaterialInterface<2>;
template class LinearMaterialInterface<3>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
