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
#include <math.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include "tests.h"

#include <galerkin_tools/assembly_helper.h>
#include <galerkin_tools/linear_material.h>
#include <galerkin_tools/dirichlet_constraint.h>
#include <galerkin_tools/solver_wrapper.h>
#include <galerkin_tools/two_block_sparsity_pattern.h>
#include <galerkin_tools/two_block_matrix.h>

using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;

//class defining an inhomogeneous dirichlet constraint with fixed prescribed value (which does not depend on the location on the interface/boundary where the constraint is applied)
template<unsigned int spacedim>
class InhomogeneousDirichletBC : public dealii::Function<spacedim>
{
	private:
		const double prescribed_value;

	public:

		//override the relevant method of dealii::Function (note: the Function does have only a single component and location is unused)
		double
		value(const Point<spacedim>& /*location*/, const unsigned int /*component = 0*/)
		const
		{
			return prescribed_value;
		}

		InhomogeneousDirichletBC(const double prescribed_value = 0.0)
		:
		Function<spacedim>(),
		prescribed_value(prescribed_value)
		{
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
					+ H_omega_H_sigma_C[1]*H_omega_H_sigma_C[1]
					+ H_omega_H_sigma_C[2]*H_omega_H_sigma_C[2]
					+ H_omega_H_sigma_C[3]*H_omega_H_sigma_C[3]);
		if(get<1>(requested_quantities))
		{
			Pi_1[0] = H_omega_H_sigma_C[0];
			Pi_1[1] = H_omega_H_sigma_C[1];
			Pi_1[2] = H_omega_H_sigma_C[2];
			Pi_1[3] = H_omega_H_sigma_C[3];
		}
		if(get<2>(requested_quantities))
		{
			Pi_2(0,0) = 1.0;
			Pi_2(1,1) = 1.0;
			Pi_2(2,2) = 1.0;
			Pi_2(3,3) = 1.0;
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
			Pi = 0.5*(H_omega_H_sigma_C[0]*H_omega_H_sigma_C[0]
					+ H_omega_H_sigma_C[1]*H_omega_H_sigma_C[1]);
		if(get<1>(requested_quantities))
		{
			Pi_1[0] = H_omega_H_sigma_C[0];
			Pi_1[1] = H_omega_H_sigma_C[1];
		}
		if(get<2>(requested_quantities))
		{
			Pi_2(0,0) = 1.0;
			Pi_2(1,1) = 1.0;
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

	//u
	IndependentField<spacedim, spacedim> u("u", FE_Q<spacedim>(2), 1, {1, 2, 3});
	//r
	IndependentField<spacedim-1, spacedim> r("r", FE_Q<spacedim-1, spacedim>(2), 1, {1, 2, 3});
	//C0
	IndependentField<0, spacedim> C0("C0", 0.0);
	//C1
	IndependentField<0, spacedim> C1("C1", 0.0);


	//H_omega_1
	vector<DependentField<spacedim, spacedim>> dependent_fields_H_omega_1;

	DependentField<spacedim, spacedim> u_val("u");
	u_val.add_term(1.0, u);
	dependent_fields_H_omega_1.push_back(u_val);

	DependentField<spacedim, spacedim> C0_val("C0");
	C0_val.add_term(1.0, C0);
	dependent_fields_H_omega_1.push_back(C0_val);

	FullMatrix<double> M(2);
	Vector<double> y(2);
	for(unsigned int m = 0; m < y.size(); ++m)
	{
		M(m, m) = 1.0;
		y(m) = 1.0;
	}
	LinearMaterialDomain<spacedim> H_omega_1(dependent_fields_H_omega_1, {1, 2, 3}, QGauss<spacedim>(3), M, y, "H_omega_1");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_omega_1(H_omega_1);

	//H_sigma_1
	vector<DependentField<spacedim-1, spacedim>> dependent_fields_H_sigma_1;

	DependentField<spacedim-1, spacedim> r_val("r");
	r_val.add_term(1.0, r);
	dependent_fields_H_sigma_1.push_back(r_val);

	DependentField<spacedim-1, spacedim> C1_val("C1");
	C1_val.add_term(1.0, C1);
	dependent_fields_H_sigma_1.push_back(C1_val);

	LinearMaterialInterface<spacedim> H_sigma_1(dependent_fields_H_sigma_1, {1, 2, 3}, QGauss<spacedim-1>(3), M, y, "H_sigma_1");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_sigma_1(H_sigma_1);

	//nonlinear part
	vector<const ScalarFunctional<spacedim,spacedim>*> H_omega;
	H_omega.push_back(&H_omega_1);
	vector<const ScalarFunctional<spacedim-1,spacedim>*> H_sigma;
	H_sigma.push_back(&H_sigma_1);
	vector<const IndependentField<0,spacedim>*> C;
	C.push_back(&C0);
	C.push_back(&C1);
	TPC_1<spacedim> nonlinear_total_potential_contribution_1(H_omega, H_sigma, C);

	H_omega.clear();
	H_omega.push_back(&H_omega_1);
	H_sigma.clear();
	C.clear();
	C.push_back(&C0);
	TPC_2<spacedim> nonlinear_total_potential_contribution_2(H_omega, H_sigma, C);

	TotalPotential<spacedim> total_potential;
	total_potential.add_total_potential_contribution(total_potential_contribution_H_omega_1);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_sigma_1);
	total_potential.add_total_potential_contribution(nonlinear_total_potential_contribution_1);
	total_potential.add_total_potential_contribution(nonlinear_total_potential_contribution_2);

	MappingQGeneric<spacedim, spacedim> mapping_domain(1);
	MappingQGeneric<spacedim-1, spacedim> mapping_interface(1);

	//constraints
	InhomogeneousDirichletBC<spacedim> bc_ux(1.0);
	DirichletConstraint<spacedim> dc1(u, 0, InterfaceSide::minus, {1, 2}, &bc_ux);

	//assembly helper and solution fields
	AssemblyHelper<spacedim> assembly_helper(total_potential, tria_system, mapping_domain, mapping_interface);
	Vector<double> solution(assembly_helper.system_size());
	Vector<double> solution_ref(assembly_helper.system_size());
	vector<const Vector<double>*> solution_ref_sets(1);
	solution_ref_sets[0]=&solution_ref;

	AffineConstraints<double> dirichlet_constraints, hanging_node_constraints, constraints;
	assembly_helper.get_dof_handler_system().make_hanging_node_constraints(hanging_node_constraints);
	assembly_helper.make_dirichlet_constraints(dirichlet_constraints, {&dc1}, hanging_node_constraints);
	hanging_node_constraints.close();
	dirichlet_constraints.close();
	constraints.close();
	constraints.merge(hanging_node_constraints);
	constraints.merge(dirichlet_constraints);


	//sequential assembly and solve
	DynamicSparsityPattern dsp_K;
	dsp_K.reinit(assembly_helper.system_size(), assembly_helper.system_size());
	assembly_helper.generate_sparsity_pattern_by_simulation(dsp_K, constraints);
	SparsityPattern sp_K;
	sp_K.copy_from(dsp_K);

	double potential_value = 0.0;
	Vector<double> f(assembly_helper.system_size());
	SparseMatrix<double> K(sp_K);
	assembly_helper.assemble_system(solution, solution_ref_sets, constraints, potential_value, f, K, make_tuple(true, true, true));
	Vector<double> x(assembly_helper.system_size());
	SolverWrapperUMFPACK solver_wrapper_umfpack;
	solver_wrapper_umfpack.solve(K, x, f, false);

	//block sizes
	auto index_set_blocks_ = assembly_helper.get_locally_owned_indices_blocks();
	const vector<unsigned int> block_sizes = {index_set_blocks_[0].size(), index_set_blocks_[1].size()};

	//parallel assembly and solve

	//the total number of processes
	int n_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	// Get the rank of the current process
	int current_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &current_proc);

	//first we have to do the partitioning
	int approximate_dofs_per_proc = block_sizes[0] / n_procs;
	if(approximate_dofs_per_proc < 1)
		approximate_dofs_per_proc = 1;

	IndexSet local_rows(block_sizes[0] + block_sizes[1]);
	int temp_proc = 0;
	for(unsigned int m = 0; m < block_sizes[0]; ++m)
	{
		if( ((int)m >= approximate_dofs_per_proc * (temp_proc + 1)) && (temp_proc + 1 < n_procs) )
			++temp_proc;
		if(current_proc == temp_proc)
			local_rows.add_index(m);
	}
	if(current_proc == n_procs - 1)
		for(int m = 0; m < (int)block_sizes[1]; ++m)
			local_rows.add_index(m + block_sizes[0]);

	TwoBlockSparsityPattern dsp_K_mpi;
	dsp_K_mpi.reinit(local_rows, block_sizes[0]);
	assembly_helper.generate_sparsity_pattern_by_simulation(dsp_K_mpi, constraints);
	dsp_K_mpi.distribute(local_rows, MPI_COMM_WORLD);
	dsp_K_mpi.finalize();

	dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix> K_mpi(dsp_K_mpi, local_rows, MPI_COMM_WORLD);
	LinearAlgebra::distributed::Vector<double> x_mpi(local_rows, MPI_COMM_WORLD);
	const vector<IndexSet> index_set_blocks = {local_rows.get_view(0, block_sizes[0]), local_rows.get_view(block_sizes[0], block_sizes[0] + block_sizes[1])};
	PETScWrappers::MPI::BlockVector f_mpi(index_set_blocks, MPI_COMM_WORLD);
	for(const auto& row : local_rows)
	{
		IndexSet column_indices(local_rows.size());
		for(auto entry = K.begin(row); entry != K.end(row); ++entry )
			K_mpi.set(entry->row(), entry->column(), K(entry->row(), entry->column()));
		f_mpi(row) = f(row);
	}
	K_mpi.compress(VectorOperation::insert);
	f_mpi.compress(VectorOperation::insert);


	deallog << "START PARALLEL MUMPS SOLVE"<< endl;
	SolverWrapperPETSc solver_wrapper_petsc;
	solver_wrapper_petsc.solve(K_mpi, x_mpi, f_mpi, false);
	deallog << "END PARALLEL MUMPS SOLVE"<< endl;


	double local_error = 0.0;
	for(const auto& row : local_rows)
		local_error += fabs(x_mpi(row) - x(row));

	double global_error;
	MPI_Reduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if(current_proc == 0)
	{
		if(global_error/K.m() < 1e-12)
			deallog << "OK" << endl;
		else
			deallog << "ERROR (error = " << global_error/K.m() << ")" << endl;
	}


}

int main(int argc, char **argv)
{
#ifdef DEAL_II_WITH_PETSC

	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
	MPILogInitAll mpi_init_log;

	deallog << "### 2D-Case, Cube ###" << endl << endl;
	check<2>();

	deallog << "\n### 3D-Case, Cube ###" << endl << endl;
	check<3>();

#else // DEAL_II_WITH_PETSC
	(void)argc;
	(void)argv;
	Assert(false, ExcMessage("PETSc is not installed!"));
#endif // DEAL_II_WITH_PETSC
}
