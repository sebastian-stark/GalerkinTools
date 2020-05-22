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
#include <ctime>

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
	//triangulation of domain
	Triangulation<spacedim, spacedim> tria_domain;
	GridGenerator::subdivided_hyper_cube(tria_domain,2,-1.0,1.0);

	TriangulationSystem<spacedim> tria_system(tria_domain);
	for(const auto& domain_cell : tria_domain.cell_iterators_on_level(0))
	{
		for(unsigned int face = 0; face<dealii::GeometryInfo<spacedim>::faces_per_cell; ++face)
		{
			if( (domain_cell->center()[0] < 0.) && (domain_cell->face(face)->center()[0] > -1e-14 ) )
				tria_system.add_interface_cell(domain_cell, face, 1);
			else if( (domain_cell->center()[0] < 0.) && (domain_cell->face(face)->center()[0] < -1.+1e-14 ) )
				tria_system.add_interface_cell(domain_cell, face, 0);
			else if( (domain_cell->center()[0] > 0.) && (domain_cell->face(face)->center()[0] > 1.-1e-14 ) )
				tria_system.add_interface_cell(domain_cell, face, 2);
			else if( (domain_cell->center()[1] < 0.) && (domain_cell->face(face)->center()[1] > -1e-14 ) )
				tria_system.add_interface_cell(domain_cell, face, 4);
			else if( (domain_cell->center()[1] < 0.) && (domain_cell->face(face)->center()[1] < -1.+1e-14 ) )
				tria_system.add_interface_cell(domain_cell, face, 3);
			else if( (domain_cell->center()[1] > 0.) && (domain_cell->face(face)->center()[1] > 1.-1e-14 ) )
				tria_system.add_interface_cell(domain_cell, face, 5);
		}
	}

	srand(time(0));
	tria_system.close();

	const unsigned int min_cell_count = 10000;
	for(;;)
	{
		for(const auto& cell : tria_domain.active_cell_iterators())
		{
			const unsigned int ref_case = rand() % (10*(spacedim - 1));
			if(ref_case == 1)
				cell->set_refine_flag();
			else if(ref_case > 1)
				cell->set_coarsen_flag();
		}
		tria_system.execute_coarsening_and_refinement();
		if(tria_domain.n_active_cells() > min_cell_count)
			break;
	}

	const bool meshes_are_consistent = tria_system.check_active_interface_cell_domain_cells_consistency();
	if(meshes_are_consistent)
		cout << "MESHES ARE CONSISTENT" << endl;
	else
		cout << "MESHES ARE INCONSISTENT" << endl;
}

int main()
{
	cout << "### 2D-Case, Cube ###\n\n";
	check<2>();

	cout << "\n### 3D-Case, Cube ###\n\n";
	check<3>();
}
