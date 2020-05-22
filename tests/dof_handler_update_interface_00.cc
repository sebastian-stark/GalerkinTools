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

#include <galerkin_tools/triangulation_system.h>
#include <galerkin_tools/dof_handler_system.h>


using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;

int main()
{
	const unsigned int spacedim = 3;

	Triangulation<spacedim, spacedim> tria_domain;
	GridGenerator::hyper_cube(tria_domain);

	TriangulationSystem<spacedim> tria_system(tria_domain);
	for(const auto& cell : tria_domain.active_cell_iterators_on_level(0))
		for(unsigned int face=0; face<dealii::GeometryInfo<spacedim>::faces_per_cell; ++face)
			if(cell->face(face)->at_boundary())
				tria_system.add_interface_cell(cell, face, 0);
	tria_system.close();

	DoFHandlerSystem<spacedim> dof_handler_system(tria_system);

	cout << "N of interface cells before refinement: " <<  tria_system.interface_active_iterators().size() << " / " << dof_handler_system.interface_active_iterators().size() << endl;
	tria_domain.refine_global(1);
	cout << "N of interface cells after refinement: " <<  tria_system.interface_active_iterators().size() << " / " << dof_handler_system.interface_active_iterators().size() << endl;
}
