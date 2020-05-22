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
#include <deal.II/lac/block_vector.h>
#include "tests.h"

#include <galerkin_tools/assembly_helper.h>
#include <galerkin_tools/linear_material.h>
#include <galerkin_tools/dirichlet_constraint.h>
#include <galerkin_tools/solver_wrapper.h>

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
	SolverWrapperUMFPACK solver_wrapper_umfpack;

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

	//stretched system
	DynamicSparsityPattern dsp_stretched;
	dsp_stretched.reinit(assembly_helper.system_size(), assembly_helper.system_size());
	assembly_helper.generate_sparsity_pattern_by_simulation(dsp_stretched, constraints);
	SparsityPattern sp_stretched;
	sp_stretched.copy_from(dsp_stretched);
	double potential_value_stretched = 0.0;
	Vector<double> f_stretched(assembly_helper.system_size());
	SparseMatrix<double> K_stretched(sp_stretched);

	assembly_helper.assemble_system(solution, solution_ref_sets, constraints, potential_value_stretched, f_stretched, K_stretched);

	Vector<double> x_stretched(assembly_helper.system_size());
	solver_wrapper_umfpack.solve(K_stretched, x_stretched, f_stretched);

	//non-stretched system
	BlockSolverWrapperUMFPACK block_solver_wrapper_umfpack;

	TwoBlockSparsityPattern sp_block;
	sp_block.reinit(assembly_helper);
	assembly_helper.generate_sparsity_pattern_by_simulation(sp_block, constraints);
	sp_block.finalize();
	double potential_value_block = 0.0;
	const vector<unsigned int> block_sizes = {sp_block.get_sp_A().n_rows(), sp_block.get_sp_D().n_rows()};
	BlockVector<double> f_block(block_sizes);
	TwoBlockMatrix<SparseMatrix<double>> K_block;
	K_block.reinit(sp_block);

	assembly_helper.assemble_system(solution, solution_ref_sets, constraints, potential_value_block, f_block, K_block);

	Vector<double> x_block(assembly_helper.system_size());
	block_solver_wrapper_umfpack.solve(K_block, x_block, f_block);

	Vector<double> x_diff(assembly_helper.system_size());
	for(unsigned int i = 0; i < assembly_helper.system_size(); ++i)
		x_diff[i] = x_stretched[i] - x_block[i];
	if(x_diff.linfty_norm() < 1e-12)
		deallog << "OK" << endl;
	else
		deallog << "ERROR" << endl;
}

int main()
{
	initlog();

	deallog << "### 2D-Case, Cube ###\n\n";
	check<2>();

	deallog << "\n### 3D-Case, Cube ###\n\n";
	check<3>();
}
