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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <galerkin_tools/assembly_helper.h>
#include <galerkin_tools/linear_material.h>

using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;

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
			const double eps_xx = 1.0;
			const double eps_yy = 0.0;
			const double gamma_xy = 0.0;
			if(component == 0)
				return eps_xx*location[0] + 0.5*gamma_xy*location[1];
			if(component == 1)
				return eps_yy*location[1] + 0.5*gamma_xy*location[0];
		}
		else if(spacedim == 3)
		{
			const double eps_xx = 1.0;
			const double eps_yy = 0.0;
			const double eps_zz = 0.0;
			const double gamma_xy = 0.0;
			const double gamma_yz = 0.0;
			const double gamma_zx = 0.0;
			if(component == 0)
				return eps_xx*location[0] + 0.5*gamma_xy*location[1] + 0.5*gamma_zx*location[2];
			if(component == 1)
				return eps_yy*location[1] + 0.5*gamma_yz*location[2] + 0.5*gamma_xy*location[0];
			if(component == 2)
				return eps_zz*location[2] + 0.5*gamma_zx*location[0] + 0.5*gamma_yz*location[1];
		}
		return 0.0;
	}
};

template<unsigned int spacedim>
void
check()
{
	Triangulation<spacedim, spacedim> tria_domain;
	GridGenerator::hyper_cube(tria_domain);

	TriangulationSystem<spacedim> tria_system(tria_domain);

	for(const auto& domain_cell : tria_domain.cell_iterators_on_level(0))
		for(unsigned int face = 0; face < GeometryInfo<spacedim>::faces_per_cell; ++face)
			if(domain_cell->face(face)->at_boundary())
				tria_system.add_interface_cell(domain_cell, face, 0);

	tria_system.close();

	InitialVal<spacedim> initial_val;
	IndependentField<spacedim,   spacedim> u("u", FE_Q<spacedim>(1), spacedim, {0}, &initial_val);

	//Linear material
	FullMatrix<double> C(6);
	Vector<double> y(6);
	for(unsigned int m = 0; m < 6; ++m)
		C(m, m) = 1.0;

	vector<DependentField<spacedim, spacedim>> dependent_fields;
	DependentField<spacedim,   spacedim> eps_xx("eps_xx"); 		eps_xx.add_term(1.0, u, 0, 0);
	DependentField<spacedim,   spacedim> eps_yy("eps_yy");		eps_yy.add_term(1.0, u, 1, 1);
	DependentField<spacedim,   spacedim> eps_zz("eps_zz");
	if(spacedim == 3)
		eps_zz.add_term(1.0, u, 2, 2);
	DependentField<spacedim,   spacedim> gamma_xy("gamma_xy"); 	gamma_xy.add_term(1.0, u, 0, 1);	gamma_xy.add_term(1.0, u, 1, 0);
	DependentField<spacedim,   spacedim> gamma_yz("gamma_yz");
	if(spacedim == 3)
	{
		gamma_yz.add_term(1.0, u, 1, 2);
		gamma_yz.add_term(1.0, u, 2, 1);
	}
	DependentField<spacedim,   spacedim> gamma_zx("gamma_zx");
	if(spacedim == 3)
	{	gamma_zx.add_term(1.0, u, 2, 0);
		gamma_zx.add_term(1.0, u, 0, 2);
	}
	dependent_fields.push_back(eps_xx);
	dependent_fields.push_back(eps_yy);
	dependent_fields.push_back(eps_zz);
	dependent_fields.push_back(gamma_xy);
	dependent_fields.push_back(gamma_yz);
	dependent_fields.push_back(gamma_zx);

	LinearMaterialDomain<spacedim> linear_material_domain(dependent_fields, {0}, QGauss<spacedim>(3), C, y, "Material");

	TotalPotentialContribution<spacedim> total_potential_contribution(linear_material_domain);
	TotalPotential<spacedim> total_potential;
	total_potential.add_total_potential_contribution(total_potential_contribution);

	MappingQGeneric<spacedim, spacedim> mapping_domain(1);
	MappingQGeneric<spacedim-1, spacedim> mapping_interface(1);

	AssemblyHelper<spacedim> assembly_helper(total_potential, tria_system, mapping_domain, mapping_interface);

	Vector<double> solution(assembly_helper.system_size());
	assembly_helper.get_initial_fields_vector(solution);
	Vector<double> solution_ref(assembly_helper.system_size());
	vector<const Vector<double>*> solution_ref_sets(1);
	solution_ref_sets[0] = &solution_ref;

	assembly_helper.compare_derivatives_with_numerical_derivatives(solution, solution_ref_sets);
}

int main()
{
	cout << "### 2D-Case, Cube ###\n\n";
	check<2>();

	cout << "\n### 3D-Case, Cube ###\n\n";
	check<3>();
}


