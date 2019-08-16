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

#include <iostream>
#include <math.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/grid/manifold_lib.h>

#include <galerkin_tools/assembly_helper.h>
#include <galerkin_tools/linear_material.h>

using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;

//initial values u, v, r, t
template<unsigned int spacedim>
class LocationSquare : public Function<spacedim>
{
public:
	LocationSquare(const unsigned int n_components = 1)
	:
	Function<spacedim>(n_components)
	{}

	double
	value(	const Point<spacedim>& location,
			const unsigned int		component)
	const
	{
		if(component == 0)
			return location.square();
		return 0.0;
	}
};

//initial values w
template<unsigned int spacedim>
class LocationComponentSquare : public Function<spacedim>
{
public:
	LocationComponentSquare(const unsigned int n_components)
	:
	Function<spacedim>(n_components)
	{}

	double
	value(	const Point<spacedim>&	location,
			const unsigned int		component)
	const
	{
		return location[component]*location[component];
	}
};

//initial values s
template<unsigned int spacedim>
class LocationComponentSquareZ : public Function<spacedim>
{
public:
	LocationComponentSquareZ(const unsigned int n_components)
	:
	Function<spacedim>(n_components)
	{}

	double
	value(	const Point<spacedim>&	location,
			const unsigned int		component)
	const
	{
		if(component == 0)
			return location.square();
		if(component == 1)
			return location[2]*location[2];
		return 0.0;
	}
};

template<unsigned int spacedim>
class TPC_1 : public TotalPotentialContribution<spacedim>
{
public:
	TPC_1(	vector<const ScalarFunctional<spacedim, spacedim>*>		H_omega,
			vector<const ScalarFunctional<spacedim-1, spacedim>*>	H_sigma,
			vector<const IndependentField<0, spacedim>*>				C)
	:
	TotalPotentialContribution<spacedim>(H_omega, H_sigma, C)
	{
	}

	virtual
	bool get_potential_contribution(	const Vector<double>&			H_omega_H_sigma_C,
										const vector<Vector<double>>&	/*C_ref_sets*/,
										double&							Pi,
										Vector<double>&					Pi_1,
										FullMatrix<double>&				Pi_2,
										const tuple<bool,bool,bool>&	requested_quantities)
	const
	{
		if(get<0>(requested_quantities))
			Pi = 0.5*(H_omega_H_sigma_C[0]*H_omega_H_sigma_C[0]
					+ H_omega_H_sigma_C[1]*H_omega_H_sigma_C[1]*2.0
					+ H_omega_H_sigma_C[2]*H_omega_H_sigma_C[0]
					+ H_omega_H_sigma_C[3]*H_omega_H_sigma_C[3]
					+ H_omega_H_sigma_C[4]*H_omega_H_sigma_C[4]
					+ H_omega_H_sigma_C[5]*H_omega_H_sigma_C[5]
					+ H_omega_H_sigma_C[5]*H_omega_H_sigma_C[0]);
		if(get<1>(requested_quantities))
		{
			Pi_1[0] = H_omega_H_sigma_C[0] + 0.5*H_omega_H_sigma_C[2] + 0.5*H_omega_H_sigma_C[5];
			Pi_1[1] = H_omega_H_sigma_C[1]*2.0;
			Pi_1[2] = H_omega_H_sigma_C[0]*0.5;
			Pi_1[3] = H_omega_H_sigma_C[3];
			Pi_1[4] = H_omega_H_sigma_C[4];
			Pi_1[5] = H_omega_H_sigma_C[5] + 0.5*H_omega_H_sigma_C[0];
		}
		if(get<2>(requested_quantities))
		{
			Pi_2(0,0) = 1.0;
			Pi_2(1,1) = 2.0;
			Pi_2(3,3) = 1.0;
			Pi_2(0,2) = 0.5;
			Pi_2(2,0) = 0.5;
			Pi_2(4,4) = 1.0;
			Pi_2(5,5) = 1.0;
			Pi_2(0,5) = 0.5;
			Pi_2(5,0) = 0.5;
		}
		return false;
	}

	virtual ~TPC_1() = default;
};

template<unsigned int spacedim>
class TPC_2 : public TotalPotentialContribution<spacedim>
{
public:
	TPC_2(	vector<const ScalarFunctional<spacedim,spacedim>*>		H_omega,
			vector<const ScalarFunctional<spacedim-1,spacedim>*>	H_sigma,
			vector<const IndependentField<0,spacedim>*>				C)
	:
	TotalPotentialContribution<spacedim>(H_omega, H_sigma, C)
	{
	}

	virtual
	bool get_potential_contribution(	const Vector<double>&			H_omega_H_sigma_C,
										const vector<Vector<double>>&	/*C_ref_sets*/,
										double&							Pi,
										Vector<double>&					Pi_1,
										FullMatrix<double>&				Pi_2,
										const tuple<bool,bool,bool>&	requested_quantities)
	const
	{
		if(get<0>(requested_quantities))
			Pi = 0.5*(H_omega_H_sigma_C[0] - H_omega_H_sigma_C[1])*(H_omega_H_sigma_C[0] - H_omega_H_sigma_C[1]) + H_omega_H_sigma_C[2] + H_omega_H_sigma_C[2]*H_omega_H_sigma_C[0];
		if(get<1>(requested_quantities))
		{
			Pi_1[0] = H_omega_H_sigma_C[0] - H_omega_H_sigma_C[1] + H_omega_H_sigma_C[2];
			Pi_1[1] = H_omega_H_sigma_C[1] - H_omega_H_sigma_C[0];
			Pi_1[2] = 1.0 + H_omega_H_sigma_C[0];
		}
		if(get<2>(requested_quantities))
		{
			Pi_2(0,0) = 1.0;
			Pi_2(1,1) = 1.0;
			Pi_2(0,1) = -1.0;
			Pi_2(1,0) = -1.0;
			Pi_2(2,0) = 1.0;
			Pi_2(0,2) = 1.0;
		}
		return false;
	}

	virtual ~TPC_2() = default;
};


template<unsigned int spacedim>
void
check()
{
	Triangulation<spacedim, spacedim> tria_domain;
	GridGenerator::subdivided_hyper_cube(tria_domain, 2, -1.0, 1.0);

	TriangulationSystem<spacedim> tria_system(tria_domain);

	for(const auto& domain_cell : tria_domain.cell_iterators_on_level(0))
	{
		if(domain_cell->center()[0] > 0.0)
			domain_cell->set_material_id(2);
		else if(domain_cell->center()[1] > 0.0)
			domain_cell->set_material_id(3);
		else
		{
			domain_cell->set_material_id(1);
			domain_cell->set_refine_flag();
		}

		for(unsigned int face=0; face<dealii::GeometryInfo<spacedim>::faces_per_cell; ++face)
		{
			if( (domain_cell->material_id() == 1) && (domain_cell->face(face)->center()[1] > -1e-14))
				tria_system.add_interface_cell(domain_cell, face, 1);
			else if((domain_cell->material_id() == 2) && (domain_cell->face(face)->center()[0] < 1e-14))
				tria_system.add_interface_cell(domain_cell, face, 2);
			else if((domain_cell->material_id() == 2) && (domain_cell->face(face)->center()[0] > 1.0-1e-14))
				tria_system.add_interface_cell(domain_cell, face, 3);
		}
	}
	tria_system.close();
	tria_domain.execute_coarsening_and_refinement();

	LocationSquare<spacedim> location_square;
	LocationComponentSquare<spacedim> location_component_square(spacedim);
	LocationComponentSquareZ<spacedim> location_component_squareZ(spacedim-1);

	//u
	IndependentField<spacedim, spacedim> u("u", FE_Q<spacedim>(2), 1, {1,3}, &location_square);
	//v
	IndependentField<spacedim, spacedim> v("v", FE_Q<spacedim>(2), 1, {2}, &location_square);
	//w
	IndependentField<spacedim, spacedim> w("w", FE_Q<spacedim>(2), spacedim, {3}, &location_component_square);
	//r
	IndependentField<spacedim-1, spacedim> r("r", FE_Q<spacedim-1, spacedim>(2), 1, {1}, &location_square);
	//s
	IndependentField<spacedim-1, spacedim> s("s", FE_Q<spacedim-1, spacedim>(2), spacedim-1, {2}, &location_component_squareZ);
	//t
	IndependentField<spacedim-1, spacedim> t("t", FE_Q<spacedim-1, spacedim>(2), 1, {3}, &location_square);
	//C0
	IndependentField<0, spacedim> C0("C0", 1.5);
	//C1
	IndependentField<0, spacedim> C1("C1", 0.5);
	//C2
	IndependentField<0, spacedim> C2("C2", 3.4);

	//H_omega_1
	vector<DependentField<spacedim, spacedim>> dependent_fields_H_omega_1;

	DependentField<spacedim, spacedim> u_val("u");
	u_val.add_term(1.0, u);
	u_val.add_term(1.3, C0);
	u_val.add_term(2.0, C1);
	u_val.add_term(1.5, C2);
	dependent_fields_H_omega_1.push_back(u_val);

	DependentField<spacedim, spacedim> u_grad_x("u,x");
	u_grad_x.add_term(1.0, u, 0, 0);
	dependent_fields_H_omega_1.push_back(u_grad_x);

	DependentField<spacedim, spacedim> u_grad_y("u,y");
	u_grad_y.add_term(1.0, u, 0, 1);
	dependent_fields_H_omega_1.push_back(u_grad_y);

	DependentField<spacedim, spacedim> u_grad_z("u,z");
	if(spacedim == 3)
		u_grad_z.add_term(1.0, u, 0, 2);
	dependent_fields_H_omega_1.push_back(u_grad_z);

	FullMatrix<double> K(4);
	Vector<double> y(4);
	for(unsigned int m = 0;m < spacedim+1; ++m)
		K(m, m) = 1.0;
	LinearMaterialDomain<spacedim> H_omega_1(dependent_fields_H_omega_1, {1}, QGauss<spacedim>(3), K, y, "H_omega_1");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_omega_1(H_omega_1);

	//H_omega_2
	vector<DependentField<spacedim, spacedim>> dependent_fields_H_omega_2;

	DependentField<spacedim, spacedim> v_val("v");
	v_val.add_term(1.0, v);
	v_val.add_term(1.2, C2);
	v_val.add_term(1.3, C0);
	dependent_fields_H_omega_2.push_back(v_val);

	K.reinit(1, 1);
	y.reinit(1);
	for(unsigned int m = 0; m < 1; ++m)
		K(m, m) = 1.0;
	LinearMaterialDomain<spacedim> H_omega_2(dependent_fields_H_omega_2, {2}, QGauss<spacedim>(3), K, y, "H_omega_2");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_omega_2(H_omega_2);

	//H_omega_3
	vector<DependentField<spacedim, spacedim>> dependent_fields_H_omega_3;

	DependentField<spacedim, spacedim> div_w("div_w");
	div_w.add_term(1.0/3.0, w, 0, 0);
	div_w.add_term(1.0/3.0, w, 1, 1);
	if(spacedim == 3)
		div_w.add_term(1.0/3.0, w, 2, 2);
	dependent_fields_H_omega_3.push_back(div_w);
	dependent_fields_H_omega_3.push_back(u_val);

	K.reinit(2, 2);
	y.reinit(2);
	for(unsigned int m = 0; m < 2; ++m)
		K(m, m) = 1.0;
	LinearMaterialDomain<spacedim> H_omega_3(dependent_fields_H_omega_3, {3,1}, QGauss<spacedim>(3), K, y, "H_omega_3");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_omega_3(H_omega_3);

	//H_sigma_1
	vector<DependentField<spacedim-1, spacedim>> dependent_fields_H_sigma_1;

	DependentField<spacedim-1, spacedim> u_x_minus("u,x(-)");
	u_x_minus.add_term(1.0, u, 0, 0, InterfaceSide::minus);
	u_x_minus.add_term(1.3, C0);
	u_x_minus.add_term(1.4, C2);
	u_x_minus.add_term(1.5, C1);
	dependent_fields_H_sigma_1.push_back(u_x_minus);

	DependentField<spacedim-1, spacedim> w_x_x_plus_m_w_y_x_minus("w_x,x(+) - w_y,x(-)");
	w_x_x_plus_m_w_y_x_minus.add_term(0.5, w, 0, 0, InterfaceSide::plus);
	w_x_x_plus_m_w_y_x_minus.add_term(-1.0, w, 1, 0, InterfaceSide::minus);
	dependent_fields_H_sigma_1.push_back(w_x_x_plus_m_w_y_x_minus);

	DependentField<spacedim-1, spacedim> r_val("r");
	r_val.add_term(1.0, r);
	dependent_fields_H_sigma_1.push_back(r_val);

	DependentField<spacedim-1, spacedim> r_grad_x("r,x");
	r_grad_x.add_term(1.0, r, 0, 0);
	dependent_fields_H_sigma_1.push_back(r_grad_x);

	DependentField<spacedim-1, spacedim> r_grad_y("r,y");
	r_grad_y.add_term(1.0, r, 0, 1);
	dependent_fields_H_sigma_1.push_back(r_grad_y);

	K.reinit(5,5);
	y.reinit(5);
	for(unsigned int m = 0; m < 5; ++m)
		K(m, m) = 1.0;
	LinearMaterialInterface<spacedim> H_sigma_1(dependent_fields_H_sigma_1, {1}, QGauss<spacedim-1>(3), K, y, "H_sigma_1");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_sigma_1(H_sigma_1);

	//H_sigma_2
	vector<DependentField<spacedim-1, spacedim>> dependent_fields_H_sigma_2;

	DependentField<spacedim-1, spacedim> s_y("s_y");
	s_y.add_term(1.0, s, 0);
	dependent_fields_H_sigma_2.push_back(s_y);

	DependentField<spacedim-1, spacedim> s_z("s_z");
	if(spacedim == 3)
		s_z.add_term(1.0, s, 1);
	dependent_fields_H_sigma_2.push_back(s_z);

	DependentField<spacedim-1, spacedim> div_s("div_s");
	div_s.add_term(1.0, s, 0, 1);
	if(spacedim == 3)
		div_s.add_term(1.0, s, 1, 2);
	dependent_fields_H_sigma_2.push_back(div_s);

	DependentField<spacedim-1, spacedim> v_minus_m_w_y_plus("v(-) - w_y(+)");
	v_minus_m_w_y_plus.add_term(1.0, v, 0, InterfaceSide::minus);
	v_minus_m_w_y_plus.add_term(-1.0, w, 1, InterfaceSide::plus);
	dependent_fields_H_sigma_2.push_back(v_minus_m_w_y_plus);

	K.reinit(4,4);
	y.reinit(4);
	for(unsigned int m=0; m<4; ++m)
		K(m, m) = 1.0;
	LinearMaterialInterface<spacedim> H_sigma_2(dependent_fields_H_sigma_2, {2}, QGauss<spacedim-1>(3), K, y, "H_sigma_2");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_sigma_2(H_sigma_2);

	//H_sigma_3
	vector<DependentField<spacedim-1, spacedim>> dependent_fields_H_sigma_3;

	DependentField<spacedim-1, spacedim> v_val_sigma("v");
	v_val_sigma.add_term(1.0, v, 0, InterfaceSide::minus);
	dependent_fields_H_sigma_3.push_back(v_val_sigma);

	DependentField<spacedim-1, spacedim> v_grad_x("v,x");
	v_grad_x.add_term(1.0, v, 0, 0, InterfaceSide::minus);
	v_grad_x.add_term(1.0, C2);
	v_grad_x.add_term(1.2, C1);
	v_grad_x.add_term(1.3, C0);
	dependent_fields_H_sigma_3.push_back(v_grad_x);

	DependentField<spacedim-1, spacedim> v_grad_y("v,y");
	v_grad_y.add_term(1.0, v, 0, 1, InterfaceSide::minus);
	dependent_fields_H_sigma_3.push_back(v_grad_y);

	DependentField<spacedim-1, spacedim> v_grad_z("v,z");
	if(spacedim == 3)
		v_grad_z.add_term(1.0, v, 0, 2, InterfaceSide::minus);
	dependent_fields_H_sigma_3.push_back(v_grad_z);

	DependentField<spacedim-1, spacedim> t_val("t");
	t_val.add_term(1.0, t);
	dependent_fields_H_sigma_3.push_back(t_val);

	K.reinit(5,5);
	y.reinit(5);
	for(unsigned int m = 0; m < 5; ++m)
		K(m,m)=1.0;
	LinearMaterialInterface<spacedim> H_sigma_3(dependent_fields_H_sigma_3, {3}, QGauss<spacedim-1>(3), K, y, "H_sigma_3");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_sigma_3(H_sigma_3);

	//nonlinear part
	vector<const ScalarFunctional<spacedim,spacedim>*> H_omega;
	H_omega.push_back(&H_omega_1);
	H_omega.push_back(&H_omega_2);
	vector<const ScalarFunctional<spacedim-1,spacedim>*> H_sigma;
	H_sigma.push_back(&H_sigma_3);
	H_sigma.push_back(&H_sigma_2);
	vector<const IndependentField<0,spacedim>*> C;
	C.push_back(&C0);
	C.push_back(&C1);
	TPC_1<spacedim> nonlinear_total_potential_contribution_1(H_omega, H_sigma, C);

	H_omega.clear();
	H_sigma.clear();
	C.clear();
	H_omega.push_back(&H_omega_1);
	H_sigma.push_back(&H_sigma_1);
	C.push_back(&C0);
	TPC_2<spacedim> nonlinear_total_potential_contribution_2(H_omega, H_sigma, C);


	TotalPotential<spacedim> total_potential;
	total_potential.add_total_potential_contribution(total_potential_contribution_H_omega_1);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_omega_2);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_omega_3);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_sigma_1);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_sigma_2);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_sigma_3);
	total_potential.add_total_potential_contribution(nonlinear_total_potential_contribution_1);
	total_potential.add_total_potential_contribution(nonlinear_total_potential_contribution_2);

	MappingQGeneric<spacedim, spacedim> mapping_domain(1);
	MappingQGeneric<spacedim-1, spacedim> mapping_interface(1);

	AssemblyHelper<spacedim> assembly_helper(total_potential, tria_system, mapping_domain, mapping_interface);

	Vector<double> solution(assembly_helper.system_size());
	Vector<double> solution_ref(assembly_helper.system_size());
	vector<const Vector<double>*> solution_ref_sets(1);
	solution_ref_sets[0]=&solution_ref;

	assembly_helper.get_initial_fields_vector(solution);

	assembly_helper.compare_derivatives_with_numerical_derivatives(solution, solution_ref_sets);
}

int main()
{
	cout << "### 2D-Case, Cube ###\n\n";
	check<2>();

	cout << "\n### 3D-Case, Cube ###\n\n";
	check<3>();
}
