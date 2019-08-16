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

#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <galerkin_tools/assembly_helper.h>
#include <galerkin_tools/solver_wrapper.h>
#include <galerkin_tools/linear_material.h>

using namespace std;
using namespace dealii;
using namespace dealii::GalerkinTools;


/**
 * Function for velocity boundary condition
 */
template<unsigned int spacedim>
class FunctionV : public dealii::Function<spacedim>
{

public:

	double
	value(const Point<spacedim>& location, const unsigned int /*component*/)
	const
	{
		if(spacedim==2)
			return sin(numbers::PI*location[0]);
		else
			return sin(numbers::PI*location[0]) * sin(numbers::PI*location[1]);
		return 0.0;
	}

	FunctionV()
	:
	Function<spacedim>()
	{
	}
};

/**
 * Function definining the interaction term at the interface between solid and fluid
 */
template<unsigned int spacedim>
class InterfaceTerm:public ScalarFunctional<spacedim-1, spacedim>
{

private:

	bool get_h_sigma(	const Vector<double>&			/*e_sigma*/,
						const vector<Vector<double>>&	/*e_sigma_ref_sets*/,
						Vector<double>&					/*hidden_vars*/,
						const Point<spacedim>&			/*x*/,
						const Tensor<1,spacedim>&		n,
						double&							h_sigma,
						Vector<double>&					h_sigma_1,
						FullMatrix<double>&				h_sigma_2,
						const tuple<bool,bool,bool>		requested_quantities)
	const
	{

		const double n1 = n[0];
		const double n2 = n[1];
		const double n3 = ( spacedim == 2 ? 0.0 : n[2] );

		if(get<0>(requested_quantities))
			h_sigma=0.0;

		if(get<1>(requested_quantities))
			h_sigma_1.reinit(9);

		if(get<2>(requested_quantities))
		{
			h_sigma_2.reinit(9,9);
			h_sigma_2(6,0) = n1;
			h_sigma_2(7,1) = n2;
			h_sigma_2(8,2) = n3;
			h_sigma_2(7,3) = n1;
			h_sigma_2(6,3) = n2;
			h_sigma_2(8,4) = n2;
			h_sigma_2(7,4) = n3;
			h_sigma_2(8,5) = n1;
			h_sigma_2(6,5) = n3;
		}
		return 0;
	}

public:
	InterfaceTerm(	const vector<DependentField<spacedim-1,spacedim>> e_sigma,
					const set<types::material_id> domain_of_integration,
					const Quadrature<spacedim-1> quadrature)
	:
	ScalarFunctional<spacedim-1, spacedim>(e_sigma, domain_of_integration, quadrature, "Interface")
	{
	}

};


int main()
{

/**************
 * parameters *
 **************/

	const unsigned int spacedim = 2;			//spatial dimension
	const double viscosity = 2.0;				//Viscosity
	const double lambda = 1.0;					//Lame parameter
	const double mu = 1.0;						//Lame parameter
    const unsigned int stokes_degree = 1;		//degree of pressure element in stokes problem (velocity is one order higer)
    const unsigned int elasticity_degree = 1;	//degree of element in elasticity problem

	MappingQGeneric<spacedim, spacedim> mapping_domain(1);		//FE mapping on domain
	MappingQGeneric<spacedim-1, spacedim> mapping_interface(1);	//FE mapping on interfaces


/********
 * grid *
 ********/

		Triangulation<spacedim> tria_domain(Triangulation<spacedim>::maximum_smoothing);
		GridGenerator::subdivided_hyper_cube (tria_domain, 8, -1.0, 1.0);
		TriangulationSystem<spacedim> tria_system(tria_domain);
		for(auto& cell : tria_domain.active_cell_iterators())
		{
			if(	((fabs(cell->center()[0]) < 0.25) && (cell->center()[spacedim-1] > 0.5))
				 ||
				((fabs(cell->center()[0]) >= 0.25) && (cell->center()[spacedim-1] > -0.5)))
		        cell->set_material_id(1);
		      else
		        cell->set_material_id(0);
		}
		for(const auto& cell : tria_domain.active_cell_iterators())
		{
			for(unsigned int f=0; f<GeometryInfo<spacedim>::faces_per_cell; ++f)
				if(cell->face(f)->at_boundary())
				{
					if(cell->face(f)->center()[spacedim-1] == 1)
						tria_system.add_interface_cell(cell, f, 1);
					else if(cell->face(f)->center()[spacedim-1]<-0.5)
						tria_system.add_interface_cell(cell, f, 0);
				}
				else if((cell->material_id()==0) && (cell->neighbor(f)->material_id()==1))
					tria_system.add_interface_cell(cell, f, 2);
		}
		tria_system.close();
		tria_domain.refine_global(2);

/**********************
 * independent fields *
 **********************/

	IndependentField<spacedim, spacedim> u("u", FE_Q<spacedim>(elasticity_degree), spacedim, {0});	//displacement
	IndependentField<spacedim, spacedim> v("v", FE_Q<spacedim>(stokes_degree+1), spacedim, {1});	//velocity
	IndependentField<spacedim, spacedim> p("p", FE_Q<spacedim>(stokes_degree), 1, {1});				//pressure

/********************
 * dependent fields *
 ********************/

	//strain
	DependentField<spacedim, spacedim> eps_xx("eps_xx");
	DependentField<spacedim, spacedim> eps_yy("eps_yy");
	DependentField<spacedim, spacedim> eps_zz("eps_zz");
	DependentField<spacedim, spacedim> gamma_xy("gamma_xy");
	DependentField<spacedim, spacedim> gamma_yz("gamma_yz");
	DependentField<spacedim, spacedim> gamma_zx("gamma_zx");
	eps_xx.add_term(1.0, u, 0, 0);
	eps_yy.add_term(1.0, u, 1, 1);
	gamma_xy.add_term(1.0, u, 0, 1);
	gamma_xy.add_term(1.0, u, 1, 0);
	if(spacedim == 3)
	{
		eps_zz.add_term(1.0, u, 2, 2);
		gamma_yz.add_term(1.0, u, 1, 2);
		gamma_yz.add_term(1.0, u, 2, 1);
		gamma_zx.add_term(1.0, u, 2, 0);
		gamma_zx.add_term(1.0, u, 0, 2);
	}

	//stretching
	DependentField<spacedim, spacedim> d_xx("d_xx");
	DependentField<spacedim, spacedim> d_yy("d_yy");
	DependentField<spacedim, spacedim> d_zz("d_zz");
	DependentField<spacedim, spacedim> d2_xy("2*d_xy");
	DependentField<spacedim, spacedim> d2_yz("2*d_yz");
	DependentField<spacedim, spacedim> d2_zx("2*d_zx");
	d_xx.add_term(1.0, v, 0, 0);
	d_yy.add_term(1.0, v, 1, 1);
	d2_xy.add_term(1.0, v, 0, 1);
	d2_xy.add_term(1.0, v, 1, 0);
	if(spacedim == 3)
	{
		d_zz.add_term(1.0, v, 2, 2);
		d2_yz.add_term(1.0, v, 1, 2);
		d2_yz.add_term(1.0, v, 2, 1);
		d2_zx.add_term(1.0, v, 2, 0);
		d2_zx.add_term(1.0, v, 0, 2);
	}

	//divergence of velocity field
	DependentField<spacedim, spacedim> div_v("div_v");
	div_v.add_term(1.0, v, 0, 0);
	div_v.add_term(1.0, v, 1, 1);
	if(spacedim == 3)
		div_v.add_term(1.0, v, 2, 2);

	//pressure field
	DependentField<spacedim, spacedim> p_("p");
	p_.add_term(1.0, p);

	//stress tensor on interface
	DependentField<spacedim-1, spacedim> sig_if_xx("sig_xx");
	DependentField<spacedim-1, spacedim> sig_if_yy("sig_yy");
	DependentField<spacedim-1, spacedim> sig_if_zz("sig_zz");
	DependentField<spacedim-1, spacedim> sig_if_xy("sig_xy");
	DependentField<spacedim-1, spacedim> sig_if_yz("sig_yz");
	DependentField<spacedim-1, spacedim> sig_if_zx("sig_zx");
	sig_if_xx.add_term(-2.0*viscosity, v, 0, 0, InterfaceSide::plus);
	sig_if_xx.add_term(1.0, p, 0, InterfaceSide::plus);
	sig_if_yy.add_term(-2.0*viscosity, v, 1, 1, InterfaceSide::plus);
	sig_if_yy.add_term(1.0, p, 0, InterfaceSide::plus);
	sig_if_xy.add_term(-viscosity, v, 0, 1, InterfaceSide::plus);
	sig_if_xy.add_term(-viscosity, v, 1, 0, InterfaceSide::plus);
	if(spacedim == 3)
	{
		sig_if_zz.add_term(-2.0*viscosity, v, 2, 2, InterfaceSide::plus);
		sig_if_zz.add_term(1.0, p, 0, InterfaceSide::plus);
		sig_if_yz.add_term(-viscosity, v, 1, 2, InterfaceSide::plus);
		sig_if_yz.add_term(-viscosity, v, 2, 1, InterfaceSide::plus);
		sig_if_zx.add_term(-viscosity, v, 0, 2, InterfaceSide::plus);
		sig_if_zx.add_term(-viscosity, v, 2, 0, InterfaceSide::plus);
	}

	//displacement on interface
	DependentField<spacedim-1, spacedim> u_if_x("u_x");
	DependentField<spacedim-1, spacedim> u_if_y("u_y");
	DependentField<spacedim-1, spacedim> u_if_z("u_z");
	u_if_x.add_term(1.0, u, 0, InterfaceSide::minus);
	u_if_y.add_term(1.0, u, 1, InterfaceSide::minus);
	if(spacedim == 3)
		u_if_z.add_term(1.0, u, 2, InterfaceSide::minus);

/************************
 * potential definition *
 ************************/

	//mechanical free energy contribution
	FullMatrix<double> C(6);
	C(0,0) = C(1,1) = C(2,2) = 2.0*mu + lambda;
	C(0,1) = C(1,0) = C(1,2) = C(2,1) = C(0,2) = C(2,0) = lambda;
	C(3,3) = C(4,4) = C(5,5) = mu;
	LinearMaterialDomain<spacedim> psi(	{eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_zx},
										{0},
										QGauss<spacedim>(elasticity_degree+2),
										C,
										Vector<double>(6),
										"psi");

	//dissipation in fluid
	FullMatrix<double> D(6);
	D(0,0) = D(1,1) = D(2,2) = 2.0*viscosity;
	D(3,3) = D(4,4) = D(5,5) = viscosity;
	LinearMaterialDomain<spacedim> delta(	{d_xx, d_yy, d_zz, d2_xy, d2_yz, d2_zx},
											{1},
											QGauss<spacedim>(stokes_degree+2),
											D,
											Vector<double>(6),
											"delta");

	//pressure term
	FullMatrix<double> C_p(2);
	C_p(0,1)=C_p(1,0)=-1.;
	LinearMaterialDomain<spacedim> pressure_term(	{p_, div_v},
													{1},
													QGauss<spacedim>(stokes_degree+2),
													C_p,
													Vector<double>(2),
													"pressure term");

	//interface term
	InterfaceTerm<spacedim> interface_term(	{sig_if_xx, sig_if_yy, sig_if_zz, sig_if_xy, sig_if_yz, sig_if_zx, u_if_x, u_if_y, u_if_z},
											{2},
											QGauss<spacedim-1>(max(stokes_degree+2, elasticity_degree + 2)));

	TotalPotentialContribution<spacedim> psi_tpc(psi);
	TotalPotentialContribution<spacedim> delta_tpc(delta);
	TotalPotentialContribution<spacedim> pressure_term_tpc(pressure_term);
	TotalPotentialContribution<spacedim> interface_term_tpc(interface_term);

 	TotalPotential<spacedim> total_potential;
 	total_potential.add_total_potential_contribution(psi_tpc);
 	total_potential.add_total_potential_contribution(delta_tpc);
	total_potential.add_total_potential_contribution(pressure_term_tpc);
	total_potential.add_total_potential_contribution(interface_term_tpc);

/******************
 * AssemblyHelper *
 ******************/

	AssemblyHelper<spacedim> assembly_helper(total_potential, tria_system, mapping_domain, mapping_interface);

/***************
 * constraints *
 ***************/

	vector<const DirichletConstraint<spacedim>*> dirichlet_constraints_vector;

	//u
	DirichletConstraint<spacedim> dc_1(u, 0, InterfaceSide::minus, {0});
	DirichletConstraint<spacedim> dc_2(u, 1, InterfaceSide::minus, {0});
	DirichletConstraint<spacedim> dc_3(u, spacedim-1, InterfaceSide::minus, {0});

	//v
	FunctionV<spacedim> inhomogeneous_bc_v;
	DirichletConstraint<spacedim> dc_4(v, 0, InterfaceSide::minus, {1});
	DirichletConstraint<spacedim> dc_5(v, 0, InterfaceSide::plus, {2});
	DirichletConstraint<spacedim> dc_6(v, 1, InterfaceSide::plus, {2});
	DirichletConstraint<spacedim> dc_7(v, spacedim-1, InterfaceSide::plus, {2});
	DirichletConstraint<spacedim> dc_8(v, 1, InterfaceSide::minus, {1}, &inhomogeneous_bc_v);
	DirichletConstraint<spacedim> dc_9(v, 1, InterfaceSide::minus, {1}, &inhomogeneous_bc_v);
	DirichletConstraint<spacedim> dc_10(v, 1, InterfaceSide::minus, {1});
	DirichletConstraint<spacedim> dc_11(v, spacedim-1, InterfaceSide::minus, {1}, &inhomogeneous_bc_v);

	dirichlet_constraints_vector.push_back(&dc_1);
	dirichlet_constraints_vector.push_back(&dc_2);
	if(spacedim == 3)
		dirichlet_constraints_vector.push_back(&dc_3);
	dirichlet_constraints_vector.push_back(&dc_4);
	dirichlet_constraints_vector.push_back(&dc_5);
	dirichlet_constraints_vector.push_back(&dc_6);
	if(spacedim == 3)
		dirichlet_constraints_vector.push_back(&dc_7);
	if(spacedim == 2)
	{
		dirichlet_constraints_vector.push_back(&dc_8);
		dirichlet_constraints_vector.push_back(&dc_9);
	}
	else
	{
		dirichlet_constraints_vector.push_back(&dc_10);
		dirichlet_constraints_vector.push_back(&dc_11);
	}

	AffineConstraints<double> hanging_node_constraints;
	assembly_helper.get_dof_handler_system().make_hanging_node_constraints(hanging_node_constraints);
	hanging_node_constraints.close();

	AffineConstraints<double> dirichlet_constraints;
	assembly_helper.make_dirichlet_constraints(	dirichlet_constraints,
												dirichlet_constraints_vector,
												hanging_node_constraints);
	dirichlet_constraints.close();

	AffineConstraints<double> constraints;
	constraints.close();

	constraints.merge(hanging_node_constraints);
	constraints.merge(dirichlet_constraints);

/*********
 * solve *
 *********/

	Vector<double> solution(assembly_helper.system_size());
	Vector<double> solution_ref(assembly_helper.system_size());
	vector<const Vector<double>*> solution_ref_sets(1);
	solution_ref_sets[0] = &solution_ref;
	assembly_helper.get_initial_fields_vector(solution);

	DynamicSparsityPattern dsp(assembly_helper.system_size(), assembly_helper.system_size());
	assembly_helper.generate_sparsity_pattern_by_simulation(dsp, constraints);
	SparsityPattern sp;
	sp.copy_from(dsp);

	double potential_value_seq;
	SparseMatrix<double> K(sp);
	Vector<double> f(assembly_helper.system_size());

	assembly_helper.assemble_system(solution, solution_ref_sets, constraints, potential_value_seq, f, K);

	SolverWrapperUMFPACK solver_wrapper;
	solver_wrapper.solve(K, solution, f);

	assembly_helper.write_output_independent_fields(solution, "output_domain", "");

}
