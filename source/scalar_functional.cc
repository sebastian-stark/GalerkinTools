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

#include <galerkin_tools/scalar_functional.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
ScalarFunctional<spacedim, spacedim>::ScalarFunctional(	vector<DependentField<spacedim,spacedim>> 	e_omega,
														const set<types::material_id> 				domain_of_integration,
														const Quadrature<spacedim> 					quadrature,
														const string 								name,
														const unsigned int 							n_ref_sets,
														const unsigned int 							n_hidden,
														const Function<spacedim> *const 			initial_vals_hidden):
domain_of_integration(domain_of_integration),
initial_vals_hidden(initial_vals_hidden),
e_omega(e_omega),
quadrature(quadrature),
n_hidden(n_hidden),
name(name),
n_ref_sets(n_ref_sets)

{
	Assert(e_omega.size() > 0, ExcMessage("No dependent variables defined for scalar functional!"));
	if(initial_vals_hidden != nullptr)
		Assert(initial_vals_hidden->n_components == n_hidden, ExcMessage("The function for the initial values of the hidden variables must have exactly the same number of components as there are hidden variables"));
}

template<unsigned int dim, unsigned int spacedim>
ScalarFunctional<dim, spacedim>::ScalarFunctional(	vector<DependentField<dim,spacedim>>	e_sigma,
													const set<types::material_id>			domain_of_integration,
													const Quadrature<dim>					quadrature,
													const string							name,
													const unsigned int						n_ref_sets,
													const unsigned int						n_hidden,
													const Function<spacedim> *const			initial_vals_hidden):
domain_of_integration(domain_of_integration),
initial_vals_hidden(initial_vals_hidden),
e_sigma(e_sigma),
quadrature(quadrature),
n_hidden(n_hidden),
name(name),
n_ref_sets(n_ref_sets)
{
	Assert(e_sigma.size() > 0, ExcMessage("No dependent variables defined for scalar functional!"));
	if(initial_vals_hidden != nullptr)
		Assert(initial_vals_hidden->n_components == n_hidden, ExcMessage("The function for the initial values of the hidden variables must have exactly the same number of components as there are hidden variables"));
}

template<unsigned int spacedim>
double
ScalarFunctional<spacedim, spacedim>::get_maximum_step(	const Vector<double>& 			/*e_omega*/,
														const vector<Vector<double>>&	/*e_omega_ref_sets*/,
														const Vector<double>& 			/*delta_e_omega*/,
														const Vector<double>& 			/*hidden_vars*/,
														const Point<spacedim>& 			/*x*/)
const
{
	return DBL_MAX;
}

template<unsigned int spacedim>
void
ScalarFunctional<spacedim, spacedim>::compare_derivatives_with_numerical_derivatives(	Vector<double>&					e_omega,
																						const vector<Vector<double>>&	e_omega_ref_sets,
																						Vector<double>&					hidden_vars,
																						const Point<spacedim>&			x,
																						const string					detailed_printout_file,
																						const double					epsilon)
const
{
	cout << "START CHECKING DERIVATIVES\n";

	double h_omega, h_omega_perturbed;
	Vector<double> h_omega_1(e_omega.size()), h_omega_1_perturbed(e_omega.size());
	FullMatrix<double> h_omega_2(e_omega.size()), h_omega_2_perturbed(e_omega.size());

	//compute derivatives with implementation
	Vector<double> hidden_vars_copy=hidden_vars;	//important: use a copy of hidden_vars, because these are overwritten
	bool error_code=get_h_omega(e_omega,
								e_omega_ref_sets,
								hidden_vars_copy,
								x,
								h_omega,
								h_omega_1,
								h_omega_2,
								make_tuple(true, true, true) );
	(void)error_code;	//silence unused parameter compiler warnings
	Assert(error_code == false, ExcMessage("Error in scalar functional computation!"));

	//compute first numerical derivative
	Vector<double> h_omega_1_numerical(h_omega_1.size());
	for(unsigned int e_omega_n=0; e_omega_n<e_omega.size(); ++e_omega_n)
	{
		e_omega[e_omega_n]+=epsilon;
		hidden_vars_copy=hidden_vars;
		bool error_code=get_h_omega(e_omega,
									e_omega_ref_sets,
									hidden_vars_copy,
									x,
									h_omega_perturbed,
									h_omega_1_perturbed,
									h_omega_2_perturbed,
									make_tuple(true, false, false) );
		(void)error_code;	//silence unused parameter compiler warnings
		Assert(error_code == false, ExcMessage("Error in scalar functional computation!"));
		h_omega_1_numerical[e_omega_n]=(h_omega_perturbed - h_omega)/epsilon;
		e_omega[e_omega_n]-=epsilon;
	}

	//compute second numerical derivative
	FullMatrix<double> h_omega_2_numerical(h_omega_2.m());
	for(unsigned int e_omega_n=0; e_omega_n<e_omega.size(); ++e_omega_n)
	{
		e_omega[e_omega_n]+=epsilon;
		hidden_vars_copy=hidden_vars;
		bool error_code=get_h_omega(e_omega,
									e_omega_ref_sets,
									hidden_vars_copy,
									x,
									h_omega_perturbed,
									h_omega_1_perturbed,
									h_omega_2_perturbed,
									make_tuple(false, true, false) );
		(void)error_code;	//silence unused parameter compiler warnings
		Assert(error_code==0, ExcMessage("Error in scalar functional computation!"));
		for(unsigned int h_omega_2_m=0; h_omega_2_m<h_omega_2_numerical.m(); ++h_omega_2_m)
			h_omega_2_numerical(h_omega_2_m,e_omega_n)=(h_omega_1_perturbed[h_omega_2_m]-h_omega_1[h_omega_2_m])/epsilon;
		e_omega[e_omega_n]-=epsilon;
	}

	//check first derivative
	double h_omega_1_norm=h_omega_1.linfty_norm();
	double max_dev_1=0.0;
	for(unsigned int h_omega_1_n=0; h_omega_1_n<h_omega_1.size(); ++h_omega_1_n)
	{
		double relative_error=(h_omega_1_numerical[h_omega_1_n]-h_omega_1[h_omega_1_n])/h_omega_1_norm;
		if(fabs(relative_error)>max_dev_1)
			max_dev_1=fabs(relative_error);
	}

	//check second derivative
	double h_omega_2_norm=h_omega_2.linfty_norm();
	double max_dev_2=0.0;
	for(unsigned int h_omega_2_m=0; h_omega_2_m<h_omega_2.m(); ++h_omega_2_m)
	{
		for(unsigned int h_omega_2_n=0; h_omega_2_n<h_omega_2.n(); ++h_omega_2_n)
		{
			double relative_error=(h_omega_2_numerical(h_omega_2_m, h_omega_2_n)-h_omega_2(h_omega_2_m,h_omega_2_n))/h_omega_2_norm;
			if(fabs(relative_error)>max_dev_2)
				max_dev_2=fabs(relative_error);
		}
	}

	//output to file if requested
	if(detailed_printout_file!="")
	{
		FILE* printout = fopen(detailed_printout_file.c_str(),"w");

		//first derivative
		fprintf(printout, "-->FIRST DERIVATIVE:\n");
		fprintf(printout, "    m  |  numerical      |  exact          |  relative deviation (in infinity norm)\n");
		for(unsigned int h_omega_1_n=0; h_omega_1_n<h_omega_1.size(); ++h_omega_1_n)
		{
			double relative_error=(h_omega_1_numerical[h_omega_1_n]-h_omega_1[h_omega_1_n])/h_omega_1_norm;
			fprintf(printout, "%5i  | %- 1.8e | %- 1.8e | %- 1.8e\n", h_omega_1_n, h_omega_1_numerical[h_omega_1_n], h_omega_1[h_omega_1_n], relative_error);
		}
		fprintf(printout, "                                             %- 1.8e MAX\n", max_dev_1);

		//second derivative
		fprintf(printout, "\n-->SECOND DERIVATIVE:\n");
		fprintf(printout, "    m,     n  |  numerical      |  exact          |  relative deviation (in infinity norm)\n");
		for(unsigned int h_omega_2_m=0; h_omega_2_m<h_omega_2.m(); ++h_omega_2_m)
		{
			for(unsigned int h_omega_2_n=0; h_omega_2_n<h_omega_2.n(); ++h_omega_2_n)
			{
				double relative_error=(h_omega_2_numerical(h_omega_2_m, h_omega_2_n)-h_omega_2(h_omega_2_m,h_omega_2_n))/h_omega_2_norm;
				fprintf(printout, "%5i  %5i  | %- 1.8e | %- 1.8e | %- 1.8e\n", h_omega_2_m, h_omega_2_n, h_omega_2_numerical(h_omega_2_m,h_omega_2_n), h_omega_2(h_omega_2_m,h_omega_2_n), relative_error);
			}
		}
		fprintf(printout, "                                                    %- 1.8e MAX\n", max_dev_2);

		fclose(printout);
	}

	printf("  Maximum relative deviation between numerical first derivative and directly computed first derivative (in infinity norm)   = % 10.8f\n", max_dev_1);
	printf("  Maximum relative deviation between numerical second derivative and directly computed second derivative (in infinity norm) = % 10.8f\n", max_dev_2);
	printf("                                                                                                         Value of integrand = % 10.8f\n", h_omega);

	cout << "FINISHED CHECKING DERIVATIVES\n";

	return;
}

template<unsigned int dim, unsigned int spacedim>
double
ScalarFunctional<dim, spacedim>::get_maximum_step(	const Vector<double>& 			/*e_sigma*/,
													const vector<Vector<double>>&	/*e_sigma_ref_sets*/,
													const Vector<double>& 			/*delta_e_sigma*/,
													const Vector<double>&			/*hidden_vars*/,
													const Point<spacedim>&			/*x*/,
													const Tensor<1, spacedim>&		/*n*/)
const
{
	return DBL_MAX;
}

template<unsigned int dim, unsigned int spacedim>
void
ScalarFunctional<dim, spacedim>::compare_derivatives_with_numerical_derivatives(Vector<double>&					e_sigma,
																				const vector<Vector<double>>&	e_sigma_ref_sets,
																				Vector<double>&					hidden_vars,
																				const Point<spacedim>&			x,
																				const Tensor<1,spacedim>&		n,
																				const string					detailed_printout_file,
																				const double					epsilon)
const
{
	cout << "START CHECKING DERIVATIVES\n";

	double h_sigma, h_sigma_perturbed;
	Vector<double> h_sigma_1(e_sigma.size()), h_sigma_1_perturbed(e_sigma.size());
	FullMatrix<double> h_sigma_2(e_sigma.size()), h_sigma_2_perturbed(e_sigma.size());

	Vector<double> hidden_vars_copy=hidden_vars;	//important: use a copy of hidden_vars, because these are overwritten
	bool error_code=get_h_sigma(e_sigma,
								e_sigma_ref_sets,
								hidden_vars_copy,
								x,
								n,
								h_sigma,
								h_sigma_1,
								h_sigma_2,
								make_tuple(true, true, true) );
	(void)error_code;	//silence unused parameter compiler warnings
	Assert(error_code == false, ExcMessage("Error in scalar functional computation!"));

	//compute first numerical derivative
	Vector<double> h_sigma_1_numerical(h_sigma_1.size());
	for(unsigned int e_sigma_n=0; e_sigma_n<e_sigma.size(); ++e_sigma_n)
	{
		e_sigma[e_sigma_n]+=epsilon;
		hidden_vars_copy=hidden_vars;
		bool error_code=get_h_sigma(e_sigma,
									e_sigma_ref_sets,
									hidden_vars_copy,
									x,
									n,
									h_sigma_perturbed,
									h_sigma_1_perturbed,
									h_sigma_2_perturbed,
									make_tuple(true, false, false) );
		(void)error_code;	//silence unused parameter compiler warnings
		Assert(error_code == false, ExcMessage("Error in scalar functional computation!"));
		h_sigma_1_numerical[e_sigma_n]=(h_sigma_perturbed - h_sigma)/epsilon;
		e_sigma[e_sigma_n]-=epsilon;
	}

	//compute second numerical derivative
	FullMatrix<double> h_sigma_2_numerical(h_sigma_2.m());
	for(unsigned int e_sigma_n=0; e_sigma_n<e_sigma.size(); ++e_sigma_n)
	{
		e_sigma[e_sigma_n]+=epsilon;
		hidden_vars_copy=hidden_vars;
		bool error_code=get_h_sigma(e_sigma,
									e_sigma_ref_sets,
									hidden_vars_copy,
									x,
									n,
									h_sigma_perturbed,
									h_sigma_1_perturbed,
									h_sigma_2_perturbed,
									make_tuple(false, true, false) );
		(void)error_code;	//silence unused parameter compiler warnings
		Assert(error_code == false, ExcMessage("Error in scalar functional computation!"));
		for(unsigned int h_sigma_2_m=0; h_sigma_2_m<h_sigma_2_numerical.m(); ++h_sigma_2_m)
			h_sigma_2_numerical(h_sigma_2_m,e_sigma_n)=(h_sigma_1_perturbed[h_sigma_2_m]-h_sigma_1[h_sigma_2_m])/epsilon;
		e_sigma[e_sigma_n]-=epsilon;
	}

	//check first derivative
	double h_sigma_1_norm=h_sigma_1.linfty_norm();
	double max_dev_1=0.0;
	for(unsigned int h_sigma_1_m=0; h_sigma_1_m<h_sigma_1.size(); ++h_sigma_1_m)
	{
		double relative_error=(h_sigma_1_numerical[h_sigma_1_m]-h_sigma_1[h_sigma_1_m])/h_sigma_1_norm;
		if(fabs(relative_error)>max_dev_1)
			max_dev_1=fabs(relative_error);
	}

	//check second derivative
	double h_sigma_2_norm=h_sigma_2.linfty_norm();
	double max_dev_2=0.0;
	for(unsigned int h_sigma_2_m=0; h_sigma_2_m<h_sigma_2.m(); ++h_sigma_2_m)
	{
		for(unsigned int h_sigma_2_n=0; h_sigma_2_n<h_sigma_2.n(); ++h_sigma_2_n)
		{
			double relative_error=(h_sigma_2_numerical(h_sigma_2_m, h_sigma_2_n)-h_sigma_2(h_sigma_2_m,h_sigma_2_n))/h_sigma_2_norm;
			if(fabs(relative_error)>max_dev_2)
				max_dev_2=fabs(relative_error);
		}
	}

	//output to file if requested
	if(detailed_printout_file!="")
	{
		FILE* printout = fopen(detailed_printout_file.c_str(),"w");

		//first derivative
		fprintf(printout, "-->FIRST DERIVATIVE:\n");
		fprintf(printout, "    m  |  numerical      |  exact          |  relative deviation (in infinity norm)\n");
		for(unsigned int h_sigma_1_m=0; h_sigma_1_m<h_sigma_1.size(); ++h_sigma_1_m)
		{
			double relative_error=(h_sigma_1_numerical[h_sigma_1_m]-h_sigma_1[h_sigma_1_m])/h_sigma_1_norm;
			fprintf(printout, "%5i  | %- 1.8e | %- 1.8e | %- 1.8e\n", h_sigma_1_m, h_sigma_1_numerical[h_sigma_1_m], h_sigma_1[h_sigma_1_m], relative_error);
		}
		fprintf(printout, "                                             %- 1.8e MAX\n", max_dev_1);

		//second derivative
		fprintf(printout, "\n-->SECOND DERIVATIVE:\n");
		fprintf(printout, "    m,     n  |  numerical      |  exact          |  relative deviation (in infinity norm)\n");
		for(unsigned int h_sigma_2_m=0; h_sigma_2_m<h_sigma_2.m(); ++h_sigma_2_m)
		{
			for(unsigned int h_sigma_2_n=0; h_sigma_2_n<h_sigma_2.n(); ++h_sigma_2_n)
			{
				double relative_error=(h_sigma_2_numerical(h_sigma_2_m, h_sigma_2_n)-h_sigma_2(h_sigma_2_m,h_sigma_2_n))/h_sigma_2_norm;
				fprintf(printout, "%5i  %5i  | %- 1.8e | %- 1.8e | %- 1.8e\n", h_sigma_2_m, h_sigma_2_n, h_sigma_2_numerical(h_sigma_2_m,h_sigma_2_n), h_sigma_2(h_sigma_2_m,h_sigma_2_n), relative_error);
			}
		}
		fprintf(printout, "                                                    %- 1.8e MAX\n", max_dev_2);

		fclose(printout);
	}

	printf("  Maximum relative deviation between numerical first derivative and directly computed first derivative (in infinity norm)   = % 10.8f\n", max_dev_1);
	printf("  Maximum relative deviation between numerical second derivative and directly computed second derivative (in infinity norm) = % 10.8f\n", max_dev_2);
	printf("                                                                                                         Value of integrand = % 10.8f\n", h_sigma);

	cout << "FINISHED CHECKING DERIVATIVES\n";

	return;
}

template<unsigned int spacedim>
ScalarFunctional<spacedim, spacedim>::~ScalarFunctional()
{
	Assert(n_subscriptions() == 0, ExcMessage("You are about to destroy a ScalarFunctional, which is currently in use! Make sure that all ScalarFunctional objects live as least as long as the objects using them!"));
}

template<unsigned int dim, unsigned int spacedim>
ScalarFunctional<dim,spacedim>::~ScalarFunctional()
{
	Assert(n_subscriptions() == 0, ExcMessage("You are about to destroy a ScalarFunctional, which is currently in use! Make sure that all ScalarFunctional objects live as least as long as the objects using them!"));
}

template class ScalarFunctional<3,3>;
template class ScalarFunctional<2,2>;
template class ScalarFunctional<2,3>;
template class ScalarFunctional<1,2>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
