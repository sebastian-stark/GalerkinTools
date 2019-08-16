#include <iostream>
#include <math.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/block_vector.h>
#include "tests.h"

#include <galerkin_tools/assembly_helper.h>
#include <galerkin_tools/linear_material.h>
#include <galerkin_tools/dof_renumbering.h>
#include <galerkin_tools/tools.h>
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

//initial values u, v, r, t
template<unsigned int spacedim>
class LocationSquare : public Function<spacedim>
{
public:
	LocationSquare(const unsigned int n_components = 1)
	:
	Function<spacedim>(n_components)
	{}

	double
	value(	const Point<spacedim>& location,
			const unsigned int		component)
	const
	{
		if(component == 0)
			return location.square();
		return 0.0;
	}
};

//initial values w
template<unsigned int spacedim>
class LocationComponentSquare : public Function<spacedim>
{
public:
	LocationComponentSquare(const unsigned int n_components)
	:
	Function<spacedim>(n_components)
	{}

	double
	value(	const Point<spacedim>&	location,
			const unsigned int		component)
	const
	{
		return location[component]*location[component];
	}
};

//initial values s
template<unsigned int spacedim>
class LocationComponentSquareZ : public Function<spacedim>
{
public:
	LocationComponentSquareZ(const unsigned int n_components)
	:
	Function<spacedim>(n_components)
	{}

	double
	value(	const Point<spacedim>&	location,
			const unsigned int		component)
	const
	{
		if(component == 0)
			return location.square();
		if(component == 1)
			return location[2]*location[2];
		return 0.0;
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
					+ H_omega_H_sigma_C[1]*H_omega_H_sigma_C[1]*2.0
					+ H_omega_H_sigma_C[2]*H_omega_H_sigma_C[0]
					+ H_omega_H_sigma_C[3]*H_omega_H_sigma_C[3]
					+ H_omega_H_sigma_C[4]*H_omega_H_sigma_C[4]
					+ H_omega_H_sigma_C[5]*H_omega_H_sigma_C[5]
					+ H_omega_H_sigma_C[5]*H_omega_H_sigma_C[0]
					+ H_omega_H_sigma_C[6]*H_omega_H_sigma_C[6]);
		if(get<1>(requested_quantities))
		{
			Pi_1[0] = H_omega_H_sigma_C[0] + 0.5*H_omega_H_sigma_C[2] + 0.5*H_omega_H_sigma_C[5];
			Pi_1[1] = H_omega_H_sigma_C[1]*2.0;
			Pi_1[2] = H_omega_H_sigma_C[0]*0.5;
			Pi_1[3] = H_omega_H_sigma_C[3];
			Pi_1[4] = H_omega_H_sigma_C[4];
			Pi_1[5] = H_omega_H_sigma_C[5] + 0.5*H_omega_H_sigma_C[0];
			Pi_1[6] = H_omega_H_sigma_C[6];
		}
		if(get<2>(requested_quantities))
		{
			Pi_2(0,0) = 1.0;
			Pi_2(1,1) = 2.0;
			Pi_2(3,3) = 1.0;
			Pi_2(0,2) = 0.5;
			Pi_2(2,0) = 0.5;
			Pi_2(4,4) = 1.0;
			Pi_2(5,5) = 1.0;
			Pi_2(0,5) = 0.5;
			Pi_2(5,0) = 0.5;
			Pi_2(6,6) = 1.0;
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
			Pi = 0.5*(H_omega_H_sigma_C[0] - H_omega_H_sigma_C[1])*(H_omega_H_sigma_C[0] - H_omega_H_sigma_C[1]) + H_omega_H_sigma_C[2] + H_omega_H_sigma_C[2]*H_omega_H_sigma_C[0];
		if(get<1>(requested_quantities))
		{
			Pi_1[0] = H_omega_H_sigma_C[0] - H_omega_H_sigma_C[1] + H_omega_H_sigma_C[2];
			Pi_1[1] = H_omega_H_sigma_C[1] - H_omega_H_sigma_C[0];
			Pi_1[2] = 1.0 + H_omega_H_sigma_C[0];
		}
		if(get<2>(requested_quantities))
		{
			Pi_2(0,0) = 1.0;
			Pi_2(1,1) = 1.0;
			Pi_2(0,1) = -1.0;
			Pi_2(1,0) = -1.0;
			Pi_2(2,0) = 1.0;
			Pi_2(0,2) = 1.0;
		}
		return false;
	}

	virtual ~TPC_2() = default;
};


template<unsigned int spacedim>
void
check()
{

	//parallel triangulation
	dealii::parallel::distributed::Triangulation<spacedim> tria_domain(MPI_COMM_WORLD);
	GridGenerator::subdivided_hyper_cube(tria_domain, 2, -1.0, 1.0);

	dealii::GalerkinTools::parallel::TriangulationSystem<spacedim> tria_system(tria_domain);

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
	tria_domain.refine_global(1);

	//sequential triangulation
	Triangulation<spacedim> tria_domain_seq;
	GridGenerator::subdivided_hyper_cube(tria_domain_seq, 2, -1.0, 1.0);

	TriangulationSystem<spacedim> tria_system_seq(tria_domain_seq);

	for(const auto& domain_cell : tria_domain_seq.cell_iterators_on_level(0))
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
				tria_system_seq.add_interface_cell(domain_cell, face, 1);
			else if((domain_cell->material_id() == 2) && (domain_cell->face(face)->center()[0] < 1e-14))
				tria_system_seq.add_interface_cell(domain_cell, face, 2);
			else if((domain_cell->material_id() == 2) && (domain_cell->face(face)->center()[0] > 1.0-1e-14))
				tria_system_seq.add_interface_cell(domain_cell, face, 3);
		}
	}
	tria_system_seq.close();
	tria_domain_seq.execute_coarsening_and_refinement();
	tria_domain_seq.refine_global(1);

	LocationSquare<spacedim> location_square;
	LocationComponentSquare<spacedim> location_component_square(spacedim);
	LocationComponentSquareZ<spacedim> location_component_squareZ(spacedim-1);

	//u
	IndependentField<spacedim, spacedim> u("u", FE_Q<spacedim>(2), 1, {1,3}, &location_square);
	//v
	IndependentField<spacedim, spacedim> v("v", FE_Q<spacedim>(2), 1, {2}, &location_square);
	//w
	IndependentField<spacedim, spacedim> w("w", FE_Q<spacedim>(2), spacedim, {3}, &location_component_square);
	//r
	IndependentField<spacedim-1, spacedim> r("r", FE_Q<spacedim-1, spacedim>(2), 1, {1}, &location_square);
	//s
	IndependentField<spacedim-1, spacedim> s("s", FE_Q<spacedim-1, spacedim>(2), spacedim-1, {2}, &location_component_squareZ);
	//t
	IndependentField<spacedim-1, spacedim> t("t", FE_Q<spacedim-1, spacedim>(2), 1, {3}, &location_square);
	//C0
	IndependentField<0, spacedim> C0("C0", 1.5);
	//C1
	IndependentField<0, spacedim> C1("C1", 0.5);
	//C2
	IndependentField<0, spacedim> C2("C2", 3.4);

	//H_omega_1
	vector<DependentField<spacedim, spacedim>> dependent_fields_H_omega_1;

	DependentField<spacedim, spacedim> u_val("u");
	u_val.add_term(1.0, u);
	u_val.add_term(1.3, C0);
	u_val.add_term(2.0, C1);
	u_val.add_term(1.5, C2);
	dependent_fields_H_omega_1.push_back(u_val);

	DependentField<spacedim, spacedim> u_grad_x("u,x");
	u_grad_x.add_term(1.0, u, 0, 0);
	dependent_fields_H_omega_1.push_back(u_grad_x);

	DependentField<spacedim, spacedim> u_grad_y("u,y");
	u_grad_y.add_term(1.0, u, 0, 1);
	dependent_fields_H_omega_1.push_back(u_grad_y);

	DependentField<spacedim, spacedim> u_grad_z("u,z");
	if(spacedim == 3)
		u_grad_z.add_term(1.0, u, 0, 2);
	dependent_fields_H_omega_1.push_back(u_grad_z);

	FullMatrix<double> K(4);
	Vector<double> y(4);
	for(unsigned int m = 0;m < spacedim+1; ++m)
		K(m, m) = 1.0;
	LinearMaterialDomain<spacedim> H_omega_1(dependent_fields_H_omega_1, {1}, QGauss<spacedim>(3), K, y, "H_omega_1");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_omega_1(H_omega_1);

	//H_omega_2
	vector<DependentField<spacedim, spacedim>> dependent_fields_H_omega_2;

	DependentField<spacedim, spacedim> v_val("v");
	v_val.add_term(1.0, v);
	v_val.add_term(1.2, C2);
	v_val.add_term(1.3, C0);
	dependent_fields_H_omega_2.push_back(v_val);

	K.reinit(1, 1);
	y.reinit(1);
	for(unsigned int m = 0; m < 1; ++m)
		K(m, m) = 1.0;
	LinearMaterialDomain<spacedim> H_omega_2(dependent_fields_H_omega_2, {2}, QGauss<spacedim>(3), K, y, "H_omega_2");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_omega_2(H_omega_2);

	//H_omega_3
	vector<DependentField<spacedim, spacedim>> dependent_fields_H_omega_3;

	DependentField<spacedim, spacedim> div_w("div_w");
	div_w.add_term(1.0/3.0, w, 0, 0);
	div_w.add_term(1.0/3.0, w, 1, 1);
	if(spacedim == 3)
		div_w.add_term(1.0/3.0, w, 2, 2);
	dependent_fields_H_omega_3.push_back(div_w);
	dependent_fields_H_omega_3.push_back(u_val);

	K.reinit(2, 2);
	y.reinit(2);
	for(unsigned int m = 0; m < 2; ++m)
		K(m, m) = 1.0;
	LinearMaterialDomain<spacedim> H_omega_3(dependent_fields_H_omega_3, {3,1}, QGauss<spacedim>(3), K, y, "H_omega_3");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_omega_3(H_omega_3);

	//H_omega_4
	vector<DependentField<spacedim, spacedim>> dependent_fields_H_omega_4;

	DependentField<spacedim, spacedim> w_x("w_x");
	w_x.add_term(1.0, w, 0);

	DependentField<spacedim, spacedim> w_y("w_y");
	w_y.add_term(1.0, w, 1);

	DependentField<spacedim, spacedim> w_z("w_z");
	if(spacedim == 3)
		w_z.add_term(1.0, w, 2);

	dependent_fields_H_omega_4.push_back(u_val);
	dependent_fields_H_omega_4.push_back(v_val);
	dependent_fields_H_omega_4.push_back(w_x);
	dependent_fields_H_omega_4.push_back(w_y);
	dependent_fields_H_omega_4.push_back(w_z);

	K.reinit(5, 5);
	y.reinit(5);
	for(unsigned int m = 0; m < 5; ++m)
		K(m, m) = 1.0;
	LinearMaterialDomain<spacedim> H_omega_4(dependent_fields_H_omega_4, {1,2,3}, QGauss<spacedim>(3), K, y, "H_omega_4");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_omega_4(H_omega_4);

	//H_sigma_1
	vector<DependentField<spacedim-1, spacedim>> dependent_fields_H_sigma_1;

	DependentField<spacedim-1, spacedim> u_x_minus("u,x(-)");
	u_x_minus.add_term(1.0, u, 0, 0, InterfaceSide::minus);
	u_x_minus.add_term(1.3, C0);
	u_x_minus.add_term(1.4, C2);
	u_x_minus.add_term(1.5, C1);
	dependent_fields_H_sigma_1.push_back(u_x_minus);

	DependentField<spacedim-1, spacedim> w_x_x_plus_m_w_y_x_minus("w_x,x(+) - w_y,x(-)");
	w_x_x_plus_m_w_y_x_minus.add_term(0.5, w, 0, 0, InterfaceSide::plus);
	w_x_x_plus_m_w_y_x_minus.add_term(-1.0, w, 1, 0, InterfaceSide::minus);
	dependent_fields_H_sigma_1.push_back(w_x_x_plus_m_w_y_x_minus);

	DependentField<spacedim-1, spacedim> r_val("r");
	r_val.add_term(1.0, r);
	dependent_fields_H_sigma_1.push_back(r_val);

	DependentField<spacedim-1, spacedim> r_grad_x("r,x");
	r_grad_x.add_term(1.0, r, 0, 0);
	dependent_fields_H_sigma_1.push_back(r_grad_x);

	DependentField<spacedim-1, spacedim> r_grad_y("r,y");
	r_grad_y.add_term(1.0, r, 0, 1);
	dependent_fields_H_sigma_1.push_back(r_grad_y);

	K.reinit(5,5);
	y.reinit(5);
	for(unsigned int m = 0; m < 5; ++m)
		K(m, m) = 1.0;
	LinearMaterialInterface<spacedim> H_sigma_1(dependent_fields_H_sigma_1, {1}, QGauss<spacedim-1>(3), K, y, "H_sigma_1");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_sigma_1(H_sigma_1);

	//H_sigma_2
	vector<DependentField<spacedim-1, spacedim>> dependent_fields_H_sigma_2;

	DependentField<spacedim-1, spacedim> s_x("s_x");
	s_x.add_term(1.0, s, 0);

	DependentField<spacedim-1, spacedim> s_y("s_y");
	s_y.add_term(1.0, s, 0);
	dependent_fields_H_sigma_2.push_back(s_y);

	DependentField<spacedim-1, spacedim> s_z("s_z");
	if(spacedim == 3)
		s_z.add_term(1.0, s, 1);
	dependent_fields_H_sigma_2.push_back(s_z);

	DependentField<spacedim-1, spacedim> div_s("div_s");
	div_s.add_term(1.0, s, 0, 1);
	if(spacedim == 3)
		div_s.add_term(1.0, s, 1, 2);
	dependent_fields_H_sigma_2.push_back(div_s);

	DependentField<spacedim-1, spacedim> v_minus_m_w_y_plus("v(-) - w_y(+)");
	v_minus_m_w_y_plus.add_term(1.0, v, 0, InterfaceSide::minus);
	v_minus_m_w_y_plus.add_term(-1.0, w, 1, InterfaceSide::plus);
	dependent_fields_H_sigma_2.push_back(v_minus_m_w_y_plus);

	K.reinit(4,4);
	y.reinit(4);
	for(unsigned int m=0; m<4; ++m)
		K(m, m) = 1.0;
	LinearMaterialInterface<spacedim> H_sigma_2(dependent_fields_H_sigma_2, {2}, QGauss<spacedim-1>(3), K, y, "H_sigma_2");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_sigma_2(H_sigma_2);

	//H_sigma_3
	vector<DependentField<spacedim-1, spacedim>> dependent_fields_H_sigma_3;

	DependentField<spacedim-1, spacedim> v_val_sigma("v");
	v_val_sigma.add_term(1.0, v, 0, InterfaceSide::minus);
	dependent_fields_H_sigma_3.push_back(v_val_sigma);

	DependentField<spacedim-1, spacedim> v_grad_x("v,x");
	v_grad_x.add_term(1.0, v, 0, 0, InterfaceSide::minus);
	v_grad_x.add_term(1.0, C2);
	v_grad_x.add_term(1.2, C1);
	v_grad_x.add_term(1.3, C0);
	dependent_fields_H_sigma_3.push_back(v_grad_x);

	DependentField<spacedim-1, spacedim> v_grad_y("v,y");
	v_grad_y.add_term(1.0, v, 0, 1, InterfaceSide::minus);
	dependent_fields_H_sigma_3.push_back(v_grad_y);

	DependentField<spacedim-1, spacedim> v_grad_z("v,z");
	if(spacedim == 3)
		v_grad_z.add_term(1.0, v, 0, 2, InterfaceSide::minus);
	dependent_fields_H_sigma_3.push_back(v_grad_z);

	DependentField<spacedim-1, spacedim> t_val("t");
	t_val.add_term(1.0, t);
	dependent_fields_H_sigma_3.push_back(t_val);

	K.reinit(5,5);
	y.reinit(5);
	for(unsigned int m = 0; m < 5; ++m)
		K(m,m)=1.0;
	LinearMaterialInterface<spacedim> H_sigma_3(dependent_fields_H_sigma_3, {3}, QGauss<spacedim-1>(3), K, y, "H_sigma_3");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_sigma_3(H_sigma_3);

	//H_sigma_4
	vector<DependentField<spacedim-1, spacedim>> dependent_fields_H_sigma_4;
	dependent_fields_H_sigma_4.push_back(r_val);
	dependent_fields_H_sigma_4.push_back(s_x);
	dependent_fields_H_sigma_4.push_back(s_y);
	dependent_fields_H_sigma_4.push_back(s_z);
	dependent_fields_H_sigma_4.push_back(t_val);

	K.reinit(5, 5);
	y.reinit(5);
	for(unsigned int m = 0; m < 5; ++m)
		K(m, m) = 1.0;
	LinearMaterialInterface<spacedim> H_sigma_4(dependent_fields_H_sigma_4, {1,2,3}, QGauss<spacedim-1>(3), K, y, "H_sigma_4");
	TotalPotentialContribution<spacedim> total_potential_contribution_H_sigma_4(H_sigma_4);


	//nonlinear part
	vector<const ScalarFunctional<spacedim,spacedim>*> H_omega;
	H_omega.push_back(&H_omega_1);
	H_omega.push_back(&H_omega_2);
	vector<const ScalarFunctional<spacedim-1,spacedim>*> H_sigma;
	H_sigma.push_back(&H_sigma_3);
	H_sigma.push_back(&H_sigma_2);
	vector<const IndependentField<0,spacedim>*> C;
	C.push_back(&C0);
	C.push_back(&C1);
	C.push_back(&C2);
	TPC_1<spacedim> nonlinear_total_potential_contribution_1(H_omega, H_sigma, C);

	H_omega.clear();
	H_sigma.clear();
	C.clear();
	H_omega.push_back(&H_omega_1);
	H_sigma.push_back(&H_sigma_1);
	C.push_back(&C0);
	TPC_2<spacedim> nonlinear_total_potential_contribution_2(H_omega, H_sigma, C);

	TotalPotential<spacedim> total_potential;
	total_potential.add_total_potential_contribution(total_potential_contribution_H_omega_1);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_omega_2);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_omega_3);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_omega_4);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_sigma_1);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_sigma_2);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_sigma_3);
	total_potential.add_total_potential_contribution(total_potential_contribution_H_sigma_4);
	total_potential.add_total_potential_contribution(nonlinear_total_potential_contribution_1);
	total_potential.add_total_potential_contribution(nonlinear_total_potential_contribution_2);

	//constraints
	InhomogeneousDirichletBC<spacedim> bc_ux(1.0);
	DirichletConstraint<spacedim> dc1(u, 0, InterfaceSide::minus, {1, 2}, &bc_ux);

	MappingQGeneric<spacedim, spacedim> mapping_domain(1);
	MappingQGeneric<spacedim-1, spacedim> mapping_interface(1);

	//parallel assembly helper and renumbering
	DoFRenumberingOffset dof_renumbering;
	AssemblyHelper<spacedim> assembly_helper(total_potential, tria_system, mapping_domain, mapping_interface);
	Auxiliary::compute_dof_renumbering_contiguous(assembly_helper.get_dof_handler_system(), dof_renumbering);
	dof_renumbering.add_range(assembly_helper.system_size() - assembly_helper.get_n_stretched_rows(), assembly_helper.system_size() - 1, 0);
	assembly_helper.get_dof_handler_system().attach_dof_renumbering(dof_renumbering);

	//index sets
	const auto locally_owned_indices = assembly_helper.get_locally_owned_indices();
	const auto locally_relevant_indices = assembly_helper.get_locally_relevant_indices();

	//initial vectors
	LinearAlgebra::distributed::Vector<double> solution(locally_owned_indices, locally_relevant_indices, MPI_COMM_WORLD);
	LinearAlgebra::distributed::Vector<double> solution_ref(locally_owned_indices, locally_relevant_indices, MPI_COMM_WORLD);
	vector<const LinearAlgebra::distributed::Vector<double>*> solution_ref_sets(1);
	solution_ref_sets[0] = &solution_ref;
	assembly_helper.get_initial_fields_vector(solution);
	solution_ref.update_ghost_values();

	//assemble system
	AffineConstraints<double> dirichlet_constraints, hanging_node_constraints, constraints;
	constraints.reinit(locally_relevant_indices);
	assembly_helper.get_dof_handler_system().make_hanging_node_constraints(hanging_node_constraints);
	assembly_helper.make_dirichlet_constraints(dirichlet_constraints, {&dc1}, hanging_node_constraints);
	hanging_node_constraints.close();
	dirichlet_constraints.close();
	constraints.close();
	constraints.merge(hanging_node_constraints, AffineConstraints<double>::MergeConflictBehavior::no_conflicts_allowed, true);
	constraints.merge(dirichlet_constraints, AffineConstraints<double>::MergeConflictBehavior::no_conflicts_allowed, true);

	TwoBlockSparsityPattern sp_K;
	sp_K.reinit(assembly_helper);
	assembly_helper.generate_sparsity_pattern_by_simulation(sp_K, constraints);
	sp_K.distribute(locally_owned_indices, MPI_COMM_WORLD);
	sp_K.finalize();

	double potential_value;
	dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix> K_mpi(sp_K, locally_owned_indices, MPI_COMM_WORLD);
	PETScWrappers::MPI::BlockVector f_mpi(assembly_helper.get_locally_owned_indices_blocks(), MPI_COMM_WORLD);
	assembly_helper.assemble_system(solution, solution_ref_sets, constraints, potential_value, f_mpi, K_mpi);

	//solve parallel
	LinearAlgebra::distributed::Vector<double> x_mpi;
	x_mpi.reinit(locally_owned_indices, locally_relevant_indices, MPI_COMM_WORLD);
	SolverWrapperPETSc solver_wrapper_petsc;
	solver_wrapper_petsc.solve(K_mpi, x_mpi, f_mpi, false);
	constraints.distribute(x_mpi);

	//sequential for comparison
	AssemblyHelper<spacedim> assembly_helper_seq(total_potential, tria_system_seq, mapping_domain, mapping_interface);
	Vector<double> solution_seq(assembly_helper_seq.system_size());
	Vector<double> solution_ref_seq(assembly_helper_seq.system_size());
	vector<const Vector<double>*> solution_ref_sets_seq(1);
	solution_ref_sets_seq[0] = &solution_ref_seq;
	assembly_helper_seq.get_initial_fields_vector(solution_seq);

	AffineConstraints<double> dirichlet_constraints_seq, hanging_node_constraints_seq, constraints_seq;
	assembly_helper_seq.get_dof_handler_system().make_hanging_node_constraints(hanging_node_constraints_seq);
	assembly_helper_seq.make_dirichlet_constraints(dirichlet_constraints_seq, {&dc1}, hanging_node_constraints_seq);
	hanging_node_constraints_seq.close();
	dirichlet_constraints_seq.close();
	constraints_seq.close();
	constraints_seq.merge(hanging_node_constraints_seq);
	constraints_seq.merge(dirichlet_constraints_seq);

	TwoBlockSparsityPattern sp_K_seq;
	sp_K_seq.reinit(assembly_helper_seq);
	assembly_helper_seq.generate_sparsity_pattern_by_simulation(sp_K_seq, constraints_seq);
	sp_K_seq.finalize();

	double potential_value_seq;
	dealii::GalerkinTools::TwoBlockMatrix<SparseMatrix<double>> K_seq(sp_K_seq);
	const vector<unsigned int> block_sizes = {sp_K_seq.get_sp_A().n_rows(), sp_K_seq.get_sp_D().n_rows()};
	BlockVector<double> f_seq(block_sizes);
	assembly_helper_seq.assemble_system(solution_seq, solution_ref_sets_seq, constraints_seq, potential_value_seq, f_seq, K_seq);

	//solve sequential
	Vector<double> x_seq(assembly_helper_seq.system_size());
	BlockSolverWrapperUMFPACK solver_wrapper_umfpack;
	solver_wrapper_umfpack.solve(K_seq, x_seq, f_seq, false);
	constraints_seq.distribute(x_seq);

	//compare
	vector<unsigned int> map_dofs;
	Auxiliary::compute_map_dofs(assembly_helper.get_dof_handler_system(), assembly_helper_seq.get_dof_handler_system(), map_dofs);
	for(unsigned int m = 0; m < assembly_helper.get_n_stretched_rows(); ++m)
		map_dofs.push_back(map_dofs.size());

	double error_f = 0.0;
	for(const auto& index : locally_owned_indices)
		error_f += fabs(f_mpi[index] - f_seq[map_dofs[index]]);
	if(locally_owned_indices.n_elements() > 0)
		error_f = error_f / locally_owned_indices.n_elements();
	if(error_f < 1e-10)
		deallog << "OK" << endl;
	else
		deallog << "NOT OK" << endl;

	double error_K = 0.0;
	unsigned int element_count = 0;
	const auto& A_mpi = K_mpi.get_A();
	const auto& B_mpi = K_mpi.get_B();
	const auto& C_mpi = K_mpi.get_C();
	const auto& D_mpi = K_mpi.get_D();
	const auto& A_seq = K_seq.get_A();
	const auto& B_seq = K_seq.get_B();
	const auto& C_seq = K_seq.get_C();
	const auto& D_seq = K_seq.get_D();
	auto block_0_size =  A_mpi.m();
	for(const auto& row : locally_owned_indices.get_view(0, block_0_size))
	{
		for(auto element = A_mpi.begin(row); element != A_mpi.end(row); ++element)
		{
			error_K += fabs(element->value() - A_seq(map_dofs[element->row()], map_dofs[element->column()]));
			++element_count;
		}
		for(auto element = B_mpi.begin(row); element != B_mpi.end(row); ++element)
		{
			error_K += fabs(element->value() - B_seq(map_dofs[element->row()], element->column()));
			++element_count;
		}
		for(auto element = C_mpi.begin(row); element != C_mpi.end(row); ++element)
		{
			error_K += fabs(element->value() - C_seq(map_dofs[element->row()], element->column()));
			++element_count;
		}
	}
	for(const auto& row : locally_owned_indices.get_view(block_0_size, locally_owned_indices.size()))
	{
		for(auto element = D_mpi.begin(row); element != D_mpi.end(row); ++element)
		{
			error_K += fabs(element->value() - D_seq(element->row(), element->column()));
			++element_count;
		}
	}

	if(element_count > 0)
		error_K = error_K / element_count;
	if(error_K < 1e-10)
		deallog << "OK" << endl;
	else
		deallog << "NOT OK" << endl;

	double error_x = 0.0;
	element_count = 0;
	for(const auto& row : locally_owned_indices)
	{
		error_x += fabs(x_mpi[row] - x_seq[map_dofs[row]]);
		++element_count;
	}
	if(element_count > 0)
		error_x = error_x / element_count;
	if(error_x < 1e-8)
		deallog << "OK" << endl;
	else
		deallog << "NOT OK" << endl;
}

int main(int argc, char **argv)
{

#ifdef DEAL_II_WITH_P4EST
	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
	MPILogInitAll log;
	check<2>();
	check<3>();
#else // DEAL_II_WITH_P4EST
	(void)argc;
	(void)argv;
	Assert(false, ExcMessage("p4est is not installed and, therefore, you cannot use the parallel version of the GalerkinTools library!"));
#endif // DEAL_II_WITH_P4EST
}
