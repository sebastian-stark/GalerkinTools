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
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

#include <galerkin_tools/triangulation_system.h>

using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;

template<unsigned int spacedim>
void
check()
{
	Triangulation<spacedim, spacedim> tria_domain;
	GridGenerator::hyper_ball(tria_domain);
	tria_domain.refine_global(2);

	TriangulationSystem<spacedim> tria_system(tria_domain);
	SphericalManifold<spacedim-1, spacedim> sphericalManifold;
	tria_system.set_interface_manifold(0, sphericalManifold);

	unsigned int counter = 0;
	for(const auto& domain_cell : tria_domain.cell_iterators_on_level(0))
		for(unsigned int face = 0; face < dealii::GeometryInfo<spacedim>::faces_per_cell; ++face)
			if(domain_cell->face(face)->at_boundary())
			{
				tria_system.add_interface_cell(domain_cell, face, counter);
				++counter;
			}

	tria_system.close();

	GridOut grid_out_domain;
	grid_out_domain.write_vtk(tria_system.get_triangulation_domain(), cout);

	GridOut grid_out_interface;
	grid_out_interface.write_vtk(tria_system.get_triangulation_interface(), cout);
}

int main()
{
	cout << "### 2D-Case, Cube ###\n\n";
	check<2>();

	cout << "\n### 3D-Case, Cube ###\n\n";
	check<3>();
}
