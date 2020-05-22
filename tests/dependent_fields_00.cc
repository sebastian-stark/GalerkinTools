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

#include <iostream>

#include <deal.II/fe/fe_q.h>

#include <galerkin_tools/dependent_field.h>

using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;

template<unsigned int spacedim>
void test()
{
	IndependentField<spacedim,spacedim> u_omega("u_omega", FE_Q<spacedim>(2), spacedim, {0});
	IndependentField<spacedim,spacedim> p_omega("p_omega", FE_Q<spacedim>(1), 1, {0});
	IndependentField<spacedim-1,spacedim> u_sigma("u_sigma", FE_Q<spacedim-1,spacedim>(2), spacedim, {0});
	IndependentField<spacedim-1,spacedim> p_sigma("p_sigma", FE_Q<spacedim-1,spacedim>(1), 1, {0});
	IndependentField<0,spacedim> C_1("C_1");
	IndependentField<0,spacedim> C_2("C_2");

	DependentField<spacedim,spacedim> e_omega("e_omega");
	DependentField<spacedim-1,spacedim> e_sigma("e_sigma");

	e_omega.add_term(1.1);
	e_omega.add_term(1.2, u_omega, 0);
	e_omega.add_term(1.3, u_omega, 1);
	e_omega.add_term(1.4, u_omega, 1, 1);
	e_omega.add_term(1.5, u_omega, 1, 0);
	e_omega.add_term(1.6, u_omega, 0, {0,1});
	e_omega.add_term(1.7, u_omega, 1, {1,1});
	e_omega.add_term(1.8, C_1);
	e_omega.add_term(1.9, C_2);
	e_omega.add_term(2.0, p_omega, 0);
	e_omega.add_term(2.1, p_omega, 0, 1);
	e_omega.add_term(2.2, p_omega, 0, {1,0});

	e_omega.print();

	e_sigma.add_term(1.1);
	e_sigma.add_term(1.2, u_sigma, 0);
	e_sigma.add_term(1.3, u_sigma, 1);
	e_sigma.add_term(1.4, u_sigma, 1, 1);
	e_sigma.add_term(1.5, u_sigma, 1, 0);
	e_sigma.add_term(1.6, u_sigma, 0, {0,1});
	e_sigma.add_term(1.7, u_sigma, 1, {1,1});
	e_sigma.add_term(1.8, C_1);
	e_sigma.add_term(1.9, C_2);
	e_sigma.add_term(2.0, p_sigma, 0);
	e_sigma.add_term(2.1, p_sigma, 0, 1);
	e_sigma.add_term(2.2, p_sigma, 0, {1,0});
	e_sigma.add_term(2.3, u_omega, 0, InterfaceSide::minus);
	e_sigma.add_term(2.4, u_omega, 1, InterfaceSide::minus);
	e_sigma.add_term(2.5, u_omega, 1, 1, InterfaceSide::minus);
	e_sigma.add_term(2.6, u_omega, 0, 0, InterfaceSide::minus);
	e_sigma.add_term(2.7, u_omega, 0, {0,1}, InterfaceSide::minus);
	e_sigma.add_term(2.8, u_omega, 1, {1,1}, InterfaceSide::minus);
	e_sigma.add_term(2.9, p_omega, 0, InterfaceSide::minus);
	e_sigma.add_term(3.0, p_omega, 0, 1, InterfaceSide::minus);
	e_sigma.add_term(3.1, p_omega, 0, {1,0}, InterfaceSide::minus);
	e_sigma.add_term(3.2, u_omega, 0, InterfaceSide::plus);
	e_sigma.add_term(3.3, u_omega, 1, InterfaceSide::plus);
	e_sigma.add_term(3.4, u_omega, 1, 1, InterfaceSide::plus);
	e_sigma.add_term(3.5, u_omega, 0, 0, InterfaceSide::plus);
	e_sigma.add_term(3.6, u_omega, 0, {0,1}, InterfaceSide::plus);
	e_sigma.add_term(3.7, u_omega, 1, {1,1}, InterfaceSide::plus);
	e_sigma.add_term(3.8, p_omega, 0, InterfaceSide::plus);
	e_sigma.add_term(3.9, p_omega, 0, 1, InterfaceSide::plus);
	e_sigma.add_term(4.0, p_omega, 0, {1,0}, InterfaceSide::plus);

	e_sigma.print();
}

int main()
{
	test<2>();
	test<3>();
}
