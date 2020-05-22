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

//Headers for standard library
#include <iostream>
#include <string>

//deal.II headers
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/dofs/dof_tools.h>
#include "tests.h"

//galerkin_tools headers
#include <galerkin_tools/triangulation_system.h>
#include <galerkin_tools/dof_handler_system.h>
#include <galerkin_tools/dof_renumbering.h>
#include <galerkin_tools/tools.h>

//namespaces
using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;


template<unsigned int spacedim>
void
check()
{
	//fe collections
	dealii::FESystem<spacedim, spacedim> fe_system_domain_0(dealii::FE_Q<spacedim, spacedim>(2), 1, dealii::FE_DGQ<spacedim, spacedim>(1), 2);
	dealii::FESystem<spacedim-1, spacedim> fe_system_interface_0(dealii::FE_Q<spacedim-1, spacedim>(2), 1, dealii::FE_DGQ<spacedim-1, spacedim>(1), 2);
	dealii::FESystem<spacedim, spacedim> fe_system_domain_1(dealii::FE_Q<spacedim, spacedim>(2), 1, dealii::FE_Nothing<spacedim, spacedim>(), 2);
	dealii::FESystem<spacedim-1, spacedim> fe_system_interface_1(dealii::FE_Q<spacedim-1, spacedim>(2), 1, dealii::FE_Nothing<spacedim-1, spacedim>(), 2);
	dealii::hp::FECollection<spacedim, spacedim> fe_collection_domain;
	dealii::hp::FECollection<spacedim-1, spacedim> fe_collection_interface;
	fe_collection_domain.push_back(fe_system_domain_0);
	fe_collection_interface.push_back(fe_system_interface_0);
	fe_collection_domain.push_back(fe_system_domain_1);
	fe_collection_interface.push_back(fe_system_interface_1);

	//setup parallel triangulation
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

	//refinements
	const unsigned int cycles = 3;
	for(unsigned int cycle = 0; cycle < cycles; ++cycle)
	{
		for(const auto& domain_cell : tria_domain.active_cell_iterators())
			if( (domain_cell->center()[0] < 0.) && (domain_cell->center()[1] < 0.))
				if( spacedim == 2 )
					domain_cell->set_refine_flag();
				else if( domain_cell->center()[2] < 0. )
					domain_cell->set_refine_flag();
		tria_domain.execute_coarsening_and_refinement();
	}

	DoFRenumberingOffset dof_renumbering_offsets;

	//dof handlers
	dealii::GalerkinTools::DoFHandlerSystem<spacedim> dof_handler_system(tria_system);

	//set active fe indices
	for(const auto& domain_cell : dof_handler_system.get_dof_handler_domain().active_cell_iterators())
		if(domain_cell->is_locally_owned())
		{
			if(domain_cell->center()[0] < 0.)
				domain_cell->set_active_fe_index(0);
			else
				domain_cell->set_active_fe_index(1);
		}
	for(const auto& interface_cell : dof_handler_system.get_dof_handler_interface().active_cell_iterators())
		if(interface_cell->is_locally_owned())
		{
			if(interface_cell->center()[0] < -1e-12)
				interface_cell->set_active_fe_index(0);
			else
				interface_cell->set_active_fe_index(1);
		}


	//distribute dofs
	dof_handler_system.distribute_dofs(fe_collection_domain, fe_collection_interface);

	//original constraints
	AffineConstraints<double> constraints;
	dof_handler_system.make_hanging_node_constraints(constraints);
	constraints.close();
	constraints.print(deallog.get_file_stream());
	deallog << endl << "###################" << endl << endl;

	//compute renumbering
	Auxiliary::compute_dof_renumbering_contiguous<spacedim>(dof_handler_system, dof_renumbering_offsets);
	dof_handler_system.attach_dof_renumbering(dof_renumbering_offsets);
	const auto& offsets = dof_renumbering_offsets.get_dof_offsets();
	for(const auto& offset : offsets)
		deallog << "[" << get<0>(offset) << ", " << get<1>(offset) << "] -> " << "[" << get<0>(offset) + get<2>(offset) << ", " << get<1>(offset) + get<2>(offset) << "]" << endl;
	deallog << endl;

	//renumbered constraints
	dof_handler_system.attach_dof_renumbering(dof_renumbering_offsets);
	dof_handler_system.make_hanging_node_constraints(constraints);
	constraints.close();
	constraints.print(deallog.get_file_stream());


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
