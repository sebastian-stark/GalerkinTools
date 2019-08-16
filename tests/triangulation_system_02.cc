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

//Headers for standard library
#include <iostream>
#include <string>

//deal.II headers
#include <deal.II/grid/grid_generator.h>
#include "tests.h"

//galerkin_tools headers
#include <galerkin_tools/triangulation_system.h>

//namespaces
using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;


template<unsigned int spacedim>
void
check()
{
	//setup triangulation
	dealii::parallel::distributed::Triangulation<spacedim> tria_domain(MPI_COMM_WORLD);
	GridGenerator::subdivided_hyper_cube(tria_domain,2,-1.0,1.0);
	dealii::GalerkinTools::parallel::TriangulationSystem<spacedim> tria_system(tria_domain);
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
	tria_system.close();

	//random refinement
	srand(time(0));
	const unsigned int min_cell_count = 10000;
	for(;;)
	{
		for(const auto& cell : tria_domain.active_cell_iterators())
		{
			if(cell->is_locally_owned())
			{
				const unsigned int ref_case = rand() % (10*(spacedim - 1));
				if(ref_case == 1)
					cell->set_refine_flag();
				else if(ref_case > 1)
					cell->set_coarsen_flag();
			}
		}
		tria_domain.execute_coarsening_and_refinement();
		if(tria_domain.n_global_active_cells() > min_cell_count)
			break;

	}

	//check consistency of meshes
	const bool meshes_are_consistent = tria_system.check_active_interface_cell_domain_cells_consistency();
	if(meshes_are_consistent)
		deallog << "MESHES ARE CONSISTENT" << endl;
	else
		deallog << "MESHES ARE INCONSISTENT" << endl;
}

int main(int argc, char **argv)
{
#ifdef DEAL_II_WITH_P4EST
	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
	MPILogInitAll log;
	check<2>();
	check<3>();
#else // DEAL_II_WITH_P4EST
	(void)argc;
	(void)argv;
	Assert(false, ExcMessage("p4est is not installed and, therefore, you cannot use the parallel version of the GalerkinTools library!"));
#endif // DEAL_II_WITH_P4EST

}
