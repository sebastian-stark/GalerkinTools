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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <galerkin_tools/assembly_helper.h>
#include <galerkin_tools/linear_material.h>

using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;

unsigned int initial_val_type = 0;

template<unsigned int spacedim>
class InitialVal : public Function<spacedim>
{
public:
	InitialVal()
	:
	Function<spacedim>(spacedim)
	{}

	double
	value(	const Point<spacedim>&	location,
			const unsigned int		component)
	const
	{
		if(spacedim == 2)
		{
			if(initial_val_type == 0)
			{
				const double eps_xx = 1.0;
				const double eps_yy = 0.0;
				const double gamma_xy = 0.0;
				const double ux = 1.0;
				const double uy = 0.0;
				if(component == 0)
					return eps_xx*location[0] + 0.5*gamma_xy*location[1] + ux;
				if(component == 1)
					return eps_yy*location[1] + 0.5*gamma_xy*location[0] + uy;
			}
			else
			{
				const double eps_xx = 2.0;
				const double eps_yy = 0.0;
				const double gamma_xy = 0.0;
				const double ux = 2.0;
				const double uy = 0.0;
				if(component == 0)
					return eps_xx*location[0] + 0.5*gamma_xy*location[1] + ux;
				if(component == 1)
					return eps_yy*location[1] + 0.5*gamma_xy*location[0] + uy;
			}
			return 0.0;
		}
		else if(spacedim == 3)
		{
			if(initial_val_type == 0)
			{
				const double eps_xx = 1.0;
				const double eps_yy = 0.0;
				const double eps_zz = 0.0;
				const double gamma_xy = 0.0;
				const double gamma_yz = 0.0;
				const double gamma_zx = 0.0;
				const double ux = 1.0;
				const double uy = 0.0;
				const double uz = 0.0;
				if(component == 0)
					return eps_xx*location[0] + 0.5*gamma_xy*location[1] + 0.5*gamma_zx*location[2] + ux;
				if(component == 1)
					return eps_yy*location[1] + 0.5*gamma_yz*location[2] + 0.5*gamma_xy*location[0] + uy;
				if(component == 2)
					return eps_zz*location[2] + 0.5*gamma_zx*location[0] + 0.5*gamma_yz*location[1] + uz;
			}
			else
			{
				const double eps_xx = 2.0;
				const double eps_yy = 0.0;
				const double eps_zz = 0.0;
				const double gamma_xy = 0.0;
				const double gamma_yz = 0.0;
				const double gamma_zx = 0.0;
				const double ux = 2.0;
				const double uy = 0.0;
				const double uz = 0.0;
				if(component == 0)
					return eps_xx*location[0] + 0.5*gamma_xy*location[1] + 0.5*gamma_zx*location[2] + ux;
				if(component == 1)
					return eps_yy*location[1] + 0.5*gamma_yz*location[2] + 0.5*gamma_xy*location[0] + uy;
				if(component == 2)
					return eps_zz*location[2] + 0.5*gamma_zx*location[0] + 0.5*gamma_yz*location[1] + uz;
			}
			return 0.0;
		}
		return 0.0;
	}
};

template<unsigned int spacedim>
void
check()
{
	Triangulation<spacedim, spacedim> tria_domain_1;
	GridGenerator::hyper_cube(tria_domain_1);
	TriangulationSystem<spacedim> tria_system_1(tria_domain_1);
	for(const auto& domain_cell : tria_domain_1.cell_iterators_on_level(0))
		for(unsigned int face = 0; face<dealii::GeometryInfo < spacedim>::faces_per_cell; ++face)
			if(domain_cell->face(face)->at_boundary())
				tria_system_1.add_interface_cell(domain_cell, face, 0);
	tria_system_1.close();

	Triangulation<spacedim, spacedim> tria_domain_2;
	GridGenerator::hyper_cube(tria_domain_2);
	TriangulationSystem<spacedim> tria_system_2(tria_domain_2);
	for(const auto& domain_cell : tria_domain_2.cell_iterators_on_level(0))
		for(unsigned int face = 0; face < dealii::GeometryInfo<spacedim>::faces_per_cell; ++face)
			if(domain_cell->face(face)->at_boundary())
				tria_system_2.add_interface_cell(domain_cell, face, 0);
	tria_system_2.close();

	InitialVal<spacedim> initial_val;
	IndependentField<spacedim, spacedim> u("u", FE_Q<spacedim>(1), spacedim, {0}, &initial_val);
	IndependentField<spacedim-1, spacedim> v("v", FE_Q<spacedim-1, spacedim>(1), spacedim, {0}, &initial_val);

	FullMatrix<double> C(3);
	Vector<double> y(3);

	vector<DependentField<spacedim, spacedim>> dependent_fields_u;
	DependentField<spacedim, spacedim> u_x("u_x");
	u_x.add_term(1.0, u, 0);
	DependentField<spacedim, spacedim> u_y("u_y");
	u_y.add_term(1.0, u, 1);
	DependentField<spacedim, spacedim> u_z("u_z");
	if(spacedim == 3)
		u_z.add_term(1.0, u, 2);
	dependent_fields_u.push_back(u_x);
	dependent_fields_u.push_back(u_y);
	dependent_fields_u.push_back(u_z);

	vector<DependentField<spacedim-1, spacedim>> dependent_fields_v;
	DependentField<spacedim-1, spacedim> v_x("v_x");
	v_x.add_term(1.0, v, 0);
	DependentField<spacedim-1, spacedim> v_y("v_y");
	v_y.add_term(1.0, v, 1);
	DependentField<spacedim-1, spacedim> v_z("v_z");
	if(spacedim == 3)
		v_z.add_term(1.0, v, 2);
	dependent_fields_v.push_back(v_x);
	dependent_fields_v.push_back(v_y);
	dependent_fields_v.push_back(v_z);

	LinearMaterialDomain<spacedim> material_u(dependent_fields_u, {0}, QGauss<spacedim>(3), C, y, "Material_u");
	LinearMaterialInterface<spacedim> material_v(dependent_fields_v, {0}, QGauss<spacedim-1>(3), C, y, "Material_v");

	TotalPotentialContribution<spacedim> total_potential_contribution_u(material_u);
	TotalPotentialContribution<spacedim> total_potential_contribution_v(material_v);
	TotalPotential<spacedim> total_potential;
	total_potential.add_total_potential_contribution(total_potential_contribution_u);
	total_potential.add_total_potential_contribution(total_potential_contribution_v);

	MappingQGeneric<spacedim, spacedim> mapping_domain(1);
	MappingQGeneric<spacedim-1, spacedim> mapping_interface(1);

	initial_val_type = 0;
	tria_domain_1.refine_global(2);

	AssemblyHelper<spacedim> assembly_helper_1(total_potential, tria_system_1, mapping_domain, mapping_interface);
	AffineConstraints<double> constraints;
	constraints.close();
	Vector<double> solution_1(assembly_helper_1.system_size());
	assembly_helper_1.get_initial_fields_vector(solution_1, &constraints);

	initial_val_type=1;
	tria_domain_2.refine_global(3);

	AssemblyHelper<spacedim> assembly_helper_2(total_potential, tria_system_2, mapping_domain, mapping_interface);
	Vector<double> solution_2(assembly_helper_2.system_size());
	assembly_helper_2.get_initial_fields_vector(solution_2, &constraints);

	const auto norms = assembly_helper_2.compute_distance_to_other_solution(solution_2, solution_1, assembly_helper_1, QGauss<spacedim>(3), QGauss<spacedim-1>(3));

	cout << "D=" << sqrt(norms.first*norms.first + norms.second*norms.second) << endl;
}

int main()
{
	cout << "### 2D-Case, Cube ###\n\n";
	check<2>();

	cout << "\n### 3D-Case, Cube ###\n\n";
	check<3>();
}
