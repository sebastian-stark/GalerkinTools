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

#include <math.h>

#include <galerkin_tools/scalar_functional.h>
#include <deal.II/lac/lapack_full_matrix.h>
using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

namespace Implementation
{

template<unsigned int dim, unsigned int spacedim>
	vector<DependentField<dim,spacedim>>
	get_e_combined(const vector<ScalarFunctional<dim,spacedim>*> scalar_functionals)
	{
		Assert( (scalar_functionals.size() > 0), ExcMessage("The combined scalar functional must be formed by at least one underlying scalar functional"));
		set<DependentField<dim,spacedim>> dependent_fields;
		for(const auto& sf : scalar_functionals)
		{
			if(dim < spacedim)
			{
				for(const auto& df : sf->e_sigma)
					dependent_fields.insert(df);
			}
			else
			{
				for(const auto& df : sf->e_omega)
					dependent_fields.insert(df);
			}
		}
		vector<DependentField<dim,spacedim>> ret;
		for(const auto& df : dependent_fields)
			ret.push_back(df);
		return ret;
	}

	template<unsigned int dim, unsigned int spacedim>
	unsigned int
	get_n_ref_sets(const vector<ScalarFunctional<dim, spacedim>*> scalar_functionals)
	{
		unsigned int ret = 0;
		for(const auto& sf : scalar_functionals)
			if(sf->n_ref_sets > ret)
				ret = sf->n_ref_sets;
		return ret;
	}

	template<unsigned int dim, unsigned int spacedim>
	unsigned int
	get_n_hidden(const vector<ScalarFunctional<dim, spacedim>*> scalar_functionals)
	{
		unsigned int ret = 0;
		for(const auto& sf : scalar_functionals)
			ret += sf->n_hidden;
		return ret;
	}

	template<unsigned int dim, unsigned int spacedim>
	bool
	get_system(	Vector<double>&											e,
				const vector<Vector<double>>&							e_ref_sets,
				Vector<double>& 										hidden_vars,
				const Point<spacedim>& 									x,
				const Tensor<1,spacedim>* 								n,
				double& 												h,
				Vector<double>& 										h_1,
				FullMatrix<double>& 									h_2,
				Vector<double>&											f_A,
				Vector<double>&											f_B,
				LAPACKFullMatrix<double>&								AA,
				LAPACKFullMatrix<double>&								AB,
				LAPACKFullMatrix<double>&								BA,
				LAPACKFullMatrix<double>&								BB,
				const vector<ScalarFunctional<dim,spacedim>*>&			scalar_functionals,
				const vector<std::vector<unsigned int>>&				map_dependent_fields,
				const vector<unsigned int>&								indices_nonlocal_dependent_fields,
				const vector<unsigned int>&								indices_local_dependent_fields)
	{

		const tuple<bool, bool, bool> requested_quantities = make_tuple(true, true, true);

		// copy data into scalar functional wise data structures
		const unsigned int N_sf = scalar_functionals.size();
		vector<Vector<double>> e_sf(N_sf);
		vector<vector<Vector<double>>> e_ref_sets_sf(N_sf);
		vector<Vector<double>> hidden_vars_sf(N_sf);
		vector<double> h_sf(N_sf, 0.0);
		vector<Vector<double>> h_1_sf(N_sf);
		vector<FullMatrix<double>> h_2_sf(N_sf);
		unsigned int counter_hv = 0;
		for(unsigned int sf = 0; sf < N_sf; ++sf)
		{
			const unsigned int N_df = (spacedim == dim) ? scalar_functionals[sf]->e_omega.size() : scalar_functionals[sf]->e_sigma.size();
			const unsigned int N_hv = scalar_functionals[sf]->n_hidden;
			e_sf[sf].reinit(N_df);
			e_ref_sets_sf[sf].resize(scalar_functionals[sf]->n_ref_sets);
			for(auto& ref_set : e_ref_sets_sf[sf])
				ref_set.reinit(N_df);
			hidden_vars_sf[sf].reinit(N_hv);
			h_1_sf[sf].reinit(N_df);
			h_2_sf[sf].reinit(N_df, N_df);

			for(unsigned int df = 0; df < N_df; ++df)
			{
				e_sf[sf][df] = e[map_dependent_fields[sf][df]];
				for(unsigned int ref_set = 0; ref_set < e_ref_sets_sf[sf].size(); ++ref_set)
					e_ref_sets_sf[sf][ref_set][df] = e_ref_sets[ref_set][map_dependent_fields[sf][df]];
			}

			for(unsigned int hv = 0; hv < N_hv; ++hv)
			{
				hidden_vars_sf[sf][hv] = hidden_vars[counter_hv];
				++counter_hv;
			}
		}

		// call all scalar functionals
		for(unsigned int sf = 0; sf < N_sf; ++sf)
		{
			if(spacedim == dim)
			{
				if(scalar_functionals[sf]->get_h_omega(e_sf[sf], e_ref_sets_sf[sf], hidden_vars_sf[sf], x, h_sf[sf], h_1_sf[sf], h_2_sf[sf], requested_quantities))
					return true;
			}
			else
			{
				if(scalar_functionals[sf]->get_h_sigma(e_sf[sf], e_ref_sets_sf[sf], hidden_vars_sf[sf], x, *n, h_sf[sf], h_1_sf[sf], h_2_sf[sf], requested_quantities))
					return true;
			}
		}

		// copy data into combined data structures
		counter_hv = 0;
		h = 0.0;
		h_1 = 0.0;
		h_2 = 0.0;
		for(unsigned int sf = 0; sf < N_sf; ++sf)
		{
			const unsigned int N_df = (spacedim == dim) ? scalar_functionals[sf]->e_omega.size() : scalar_functionals[sf]->e_sigma.size();
			const unsigned int N_hv = scalar_functionals[sf]->n_hidden;

			h += h_sf[sf];
			for(unsigned int df = 0; df < N_df; ++df)
				h_1[map_dependent_fields[sf][df]] += h_1_sf[sf][df];

			for(unsigned int df_1 = 0; df_1 < N_df; ++df_1)
				for(unsigned int df_2 = 0; df_2 < N_df; ++df_2)
					h_2(map_dependent_fields[sf][df_1],map_dependent_fields[sf][df_2]) += h_2_sf[sf](df_1, df_2);

			for(unsigned int hv = 0; hv < N_hv; ++hv)
			{
				hidden_vars[counter_hv] = hidden_vars_sf[sf][hv];
				++counter_hv;
			}
		}

		//copy data into LAPACK structures
		const unsigned int NA = indices_local_dependent_fields.size();
		const unsigned int NB = indices_nonlocal_dependent_fields.size();
		for(unsigned int m = 0; m < NA; ++m)
		{
			f_A[m] = -h_1[indices_local_dependent_fields[m]];
			for(unsigned int n = 0; n < NA; ++n)
				AA(m,n) = h_2(indices_local_dependent_fields[m], indices_local_dependent_fields[n]);
			for(unsigned int n = 0; n < NB; ++n)
				AB(m,n) = h_2(indices_local_dependent_fields[m], indices_nonlocal_dependent_fields[n]);
		}
		for(unsigned int m = 0; m < NB; ++m)
		{
			f_B[m] = -h_1[indices_nonlocal_dependent_fields[m]];
			for(unsigned int n = 0; n < NA; ++n)
				BA(m,n) = h_2(indices_nonlocal_dependent_fields[m], indices_local_dependent_fields[n]);
			for(unsigned int n = 0; n < NB; ++n)
				BB(m,n) = h_2(indices_nonlocal_dependent_fields[m], indices_nonlocal_dependent_fields[n]);
		}
		return false;
	}

	template<unsigned int dim, unsigned int spacedim>
	double
	get_maximum_step(	const Vector<double>&									e,
						const vector<Vector<double>>&							e_ref_sets,
						const Vector<double>&									delta_e,
						const Vector<double>& 									hidden_vars,
						const Point<spacedim>& 									x,
						const Tensor<1,spacedim>* 								n,
						const vector<ScalarFunctional<dim,spacedim>*>&			scalar_functionals,
						const vector<std::vector<unsigned int>>&				map_dependent_fields,
						const vector<unsigned int>&								/*indices_nonlocal_dependent_fields*/,
						const vector<unsigned int>&								/*indices_local_dependent_fields*/)
	{
		// copy data over into scalar functional wise data structures
		const unsigned int N_sf = scalar_functionals.size();
		vector<Vector<double>> e_sf(N_sf), delta_e_sf(N_sf);
		vector<vector<Vector<double>>> e_ref_sets_sf(N_sf);
		vector<Vector<double>> hidden_vars_sf(N_sf);
		unsigned int counter_hv = 0;
		for(unsigned int sf = 0; sf < N_sf; ++sf)
		{
			const unsigned int N_df = (spacedim == dim) ? scalar_functionals[sf]->e_omega.size() : scalar_functionals[sf]->e_sigma.size();
			const unsigned int N_hv = scalar_functionals[sf]->n_hidden;
			delta_e_sf[sf].reinit(N_df);
			e_sf[sf].reinit(N_df);
			e_ref_sets_sf[sf].resize(scalar_functionals[sf]->n_ref_sets);
			for(auto& ref_set : e_ref_sets_sf[sf])
				ref_set.reinit(N_df);
			hidden_vars_sf[sf].reinit(N_hv);

			for(unsigned int df = 0; df < N_df; ++df)
			{
				e_sf[sf][df] = e[map_dependent_fields[sf][df]];
				delta_e_sf[sf][df] = delta_e[map_dependent_fields[sf][df]];
				for(unsigned int ref_set = 0; ref_set < e_ref_sets_sf[sf].size(); ++ref_set)
					e_ref_sets_sf[sf][ref_set][df] = e_ref_sets[ref_set][map_dependent_fields[sf][df]];
			}

			for(unsigned int hv = 0; hv < N_hv; ++hv)
			{
				hidden_vars_sf[sf][hv] = hidden_vars[counter_hv];
				++counter_hv;
			}
		}

		// call all scalar functionals
		double max_step = DBL_MAX;
		for(unsigned int sf = 0; sf < N_sf; ++sf)
		{
			double max_step_;
			if(spacedim == dim)
				max_step_ =  scalar_functionals[sf]->get_maximum_step(e_sf[sf], e_ref_sets_sf[sf], delta_e_sf[sf], hidden_vars_sf[sf], x);
			else
				max_step_ =  scalar_functionals[sf]->get_maximum_step(e_sf[sf], e_ref_sets_sf[sf], delta_e_sf[sf], hidden_vars_sf[sf], x, *n);
			if( (max_step_ > 0.0) && (max_step_ < max_step) )
				max_step = max_step_;
		}
		return max_step;

	}

	template<unsigned int dim, unsigned int spacedim>
	bool
	get_h(	Vector<double>&											e,
			const vector<Vector<double>>&							e_ref_sets,
			Vector<double>& 										hidden_vars,
			const Point<spacedim>& 									x,
			const Tensor<1,spacedim>* 								n,
			double& 												h,
			Vector<double>& 										h_1,
			FullMatrix<double>& 									h_2,
			const tuple<bool, bool, bool>							/*requested_quantities*/,
			const vector<ScalarFunctional<dim,spacedim>*>&			scalar_functionals,
			const vector<std::vector<unsigned int>>&				map_dependent_fields,
			const vector<unsigned int>&								indices_nonlocal_dependent_fields,
			const vector<unsigned int>&								indices_local_dependent_fields,
			const double&											safety_distance,
			const double&											threshold_residual,
			const double&											threshold_step_size,
			const unsigned int&										max_iter,
			const unsigned int&										max_cutbacks,
			const bool&												use_line_search)
	{
		const unsigned int NA = indices_local_dependent_fields.size();
		const unsigned int NB = indices_nonlocal_dependent_fields.size();
		LAPACKFullMatrix<double> AA(NA), AB(NA, NB), BA(NB, NA), BB(NB);
		Vector<double> f_A(NA), f_B(NB), delta_A(NA), delta_e(NA + NB), scaling(NA), f_A_scaled(NA);

		// compute the system and determine solution increment for local dependent fields
		if(get_system<dim, spacedim>(e, e_ref_sets, hidden_vars, x, n, h, h_1, h_2, f_A, f_B, AA, AB, BA, BB, scalar_functionals, map_dependent_fields, indices_nonlocal_dependent_fields, indices_local_dependent_fields))
			return true;


		if(NA > 0)
		{
			// compute row scaling for calculation of residual
			for(unsigned int m = 0; m < NA; ++m)
				for(unsigned int n = 0; n < NA; ++n)
					if(fabs(AA(m,n)) > scaling[m])
						scaling[m] = fabs(AA(m,n));

			// compute current residual
			for(unsigned int m = 0; m < NA; ++m)
				f_A_scaled[m] = f_A[m] / scaling[m];
			double residual = f_A_scaled.l2_norm();

			// iterate
			unsigned int iter = 0;
			double residual_old = 0.0;
			for(;;)
			{
				// determine increment to solution
				AA.compute_lu_factorization();
				delta_A = f_A;
				AA.solve(delta_A);
				for(unsigned int m = 0; m < NA; ++m)
					delta_e[indices_local_dependent_fields[m]] = delta_A[m];
				const double maximum_step = Implementation::get_maximum_step<dim,spacedim>(e, e_ref_sets, delta_e, hidden_vars, x, n, scalar_functionals, map_dependent_fields, indices_nonlocal_dependent_fields, indices_local_dependent_fields);
				if(maximum_step < 1.0/safety_distance)
					delta_A *= safety_distance * maximum_step;

				// update solution
				for(unsigned int m = 0; m < NA; ++m)
					e[indices_local_dependent_fields[m]] += delta_A[m];

				// do line search
				unsigned int cutbacks = 0;
				for(;;)
				{
					// rebuild system
					AA.reinit(NA, NA);
					if(get_system<dim, spacedim>(e, e_ref_sets, hidden_vars, x, n, h, h_1, h_2, f_A, f_B, AA, AB, BA, BB, scalar_functionals, map_dependent_fields, indices_nonlocal_dependent_fields, indices_local_dependent_fields))
						return true;

					// recalculate residual
					for(unsigned int m = 0; m < NA; ++m)
						f_A_scaled[m] = f_A[m] / scaling[m];
					residual = f_A_scaled.l2_norm();

					// stop if no line search requested
					if(!use_line_search)
						break;

					// also don't line search in the first iteration or if the residual has decreased
					if( (iter == 0) || ((threshold_step_size > 0.0) && (delta_A.linfty_norm() < threshold_step_size)) || (residual < residual_old) )
						break;
					else
					{
						++cutbacks;
//						cout << "Local cutback" << endl;
						if(cutbacks > max_cutbacks)
						{
//							cout << "Exceeded number of cutbacks. Quit line search." << endl;
							return true;
						}

						delta_A *= -0.5;
						for(unsigned int m = 0; m < NA; ++m)
							e[indices_local_dependent_fields[m]] += delta_A[m];
						delta_A *= -1.0;
					}
				}


				residual_old = residual;
				++iter;

				// check convergence
				bool converged_by_residual = true;
				if(threshold_residual > 0.0)
				{
					if(!(residual < sqrt(NA) * threshold_residual))
						converged_by_residual = false;
				}

				bool converged_by_step_size = true;
				double step_size = 0.0;
				if(threshold_step_size > 0.0)
				{
					step_size = delta_A.linfty_norm();
					if(!(fabs(step_size) < threshold_step_size))
						converged_by_step_size = false;
				}

				if(converged_by_residual && converged_by_step_size)
				{
					break;
				}

				if(iter > max_iter)
				{
//					cout << "No local convergence, stopping with residual " << residual << endl;
					return true;
				}
			}

			// write consistent tangent, etc.
			LAPACKFullMatrix<double> BA_AA_inv_AB(NB);
			AA.compute_lu_factorization();
			AA.solve(AB);
			BA.mmult(BA_AA_inv_AB, AB);

			for(unsigned int m = 0; m < NA; ++m)
				h_1[indices_local_dependent_fields[m]] = 0.0;

			for(unsigned int m = 0; m < NA; ++m)
			{
				for(unsigned int n = 0; n < NA; ++n)
					h_2(indices_local_dependent_fields[m], indices_local_dependent_fields[n]) = (m == n) ? 1.0 : 0.0;
				for(unsigned int n = 0; n < NB; ++n)
					h_2(indices_local_dependent_fields[m], indices_nonlocal_dependent_fields[n]) = h_2(indices_nonlocal_dependent_fields[n], indices_local_dependent_fields[m]) = 0.0;
			}
			for(unsigned int m = 0; m < NB; ++m)
			{
				for(unsigned int n = 0; n < NB; ++n)
					h_2(indices_nonlocal_dependent_fields[m], indices_nonlocal_dependent_fields[n]) += -BA_AA_inv_AB(m,n);
			}

		}
		return false;
	}


}

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

template<unsigned int spacedim>
bool
ScalarFunctional<spacedim, spacedim>::get_h_sigma(	Vector<double>& 				/*e_sigma*/,
													const vector<Vector<double>>&	/*e_sigma_ref_sets*/,
													Vector<double>&					/*hidden_vars*/,
													const Point<spacedim>&			/*x*/,
													const Tensor<1, spacedim>&		/*n*/,
													double&							/*h_sigma*/,
													Vector<double>&					/*h_sigma_1*/,
													FullMatrix<double>& 			/*h_sigma_2*/,
													const tuple<bool, bool, bool> 	/*requested_quantities*/)
const
{
	Assert(false, ExcMessage("This function should never be called. This indicates a bug!"));
	return true;
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

template<unsigned int dim, unsigned int spacedim>
bool
ScalarFunctional<dim, spacedim>::get_h_omega(	Vector<double>& 				/*e_omega*/,
												const vector<Vector<double>>&	/*e_omega_ref_sets*/,
												Vector<double>&					/*hidden_vars*/,
												const Point<spacedim>&			/*x*/,
												double&							/*h_omega*/,
												Vector<double>&					/*h_omega_1*/,
												FullMatrix<double>& 			/*h_omega_2*/,
												const tuple<bool, bool, bool> 	/*requested_quantities*/)
const
{
	Assert(false, ExcMessage("This function should never be called. This indicates a bug!"));
	return true;
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
double
ScalarFunctional<spacedim, spacedim>::get_maximum_step(	const Vector<double>& 			/*e_sigma*/,
														const vector<Vector<double>>&	/*e_sigma_ref_sets*/,
														const Vector<double>& 			/*delta_e_sigma*/,
														const Vector<double>& 			/*hidden_vars*/,
														const Point<spacedim>& 			/*x*/,
														const Tensor<1,spacedim>&		/*n*/)
const
{
	Assert(false, ExcMessage("This function should never be called. This indicates a bug!"));
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
double
ScalarFunctional<dim, spacedim>::get_maximum_step(	const Vector<double>& 				/*e_omega*/,
														const vector<Vector<double>>&	/*e_omega_ref_sets*/,
														const Vector<double>& 			/*delta_e_omega*/,
														const Vector<double>& 			/*hidden_vars*/,
														const Point<spacedim>& 			/*x*/)
const
{
	Assert(false, ExcMessage("This function should never be called. This indicates a bug!"));
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

template<unsigned int spacedim>
ScalarFunctionalLocalElimination<spacedim, spacedim>::ScalarFunctionalLocalElimination(	const vector<ScalarFunctional<spacedim,spacedim>*>	scalar_functionals,
																						const string 										name)
:
ScalarFunctional<spacedim,spacedim>(Implementation::get_e_combined(scalar_functionals), scalar_functionals[0]->domain_of_integration, scalar_functionals[0]->quadrature, name, Implementation::get_n_ref_sets(scalar_functionals), Implementation::get_n_hidden(scalar_functionals)),
scalar_functionals(scalar_functionals)
{
#ifdef DEBUG
	for(const auto& sf : scalar_functionals)
	{
		Assert( (sf->domain_of_integration == this->domain_of_integration), ExcMessage("Domains of integration of the underlying scalar functionals inconsistent!"));
		Assert( (sf->quadrature == this->quadrature), ExcMessage("Quadrature schemes of the underlying scalar functionals inconsistent!"));
		Assert( (sf->initial_vals_hidden == nullptr), ExcMessage("Initial values for the hidden variables are not currently implemented for ScalarFunctionalLocalElimination!") );
	}
#endif

	for(unsigned int df = 0; df < this->e_omega.size(); ++df)
	{
		if(this->e_omega[df].get_is_locally_eliminated())
			indices_locally_eliminated_dependent_fields.push_back(df);
		else
			indices_not_eliminated_dependent_fields.push_back(df);
	}

	map<DependentField<spacedim,spacedim>, unsigned int> map_df_index;
	for(unsigned int df = 0; df < this->e_omega.size(); ++df)
		map_df_index[this->e_omega[df]] = df;

	map_dependent_fields.resize(scalar_functionals.size());
	for(unsigned int sf = 0; sf < scalar_functionals.size(); ++sf)
		for(const auto& df : scalar_functionals[sf]->e_omega)
			map_dependent_fields[sf].push_back(map_df_index[df]);
}

template<unsigned int spacedim>
bool
ScalarFunctionalLocalElimination<spacedim, spacedim>::get_h_omega(	Vector<double>& 				e_omega,
																	const vector<Vector<double>>&	e_omega_ref_sets,
																	Vector<double>& 				hidden_vars,
																	const Point<spacedim>& 			x,
																	double& 						h_omega,
																	Vector<double>& 				h_omega_1,
																	FullMatrix<double>& 			h_omega_2,
																	const tuple<bool, bool, bool> 	requested_quantities)
const
{
	Vector<double> e_old = e_omega;
	Vector<double> hidden_vars_old = hidden_vars;
	bool error = Implementation::get_h<spacedim, spacedim>(e_omega, e_omega_ref_sets, hidden_vars, x, nullptr, h_omega, h_omega_1, h_omega_2, requested_quantities, scalar_functionals, map_dependent_fields, indices_not_eliminated_dependent_fields, indices_locally_eliminated_dependent_fields, safety_distance, threshold_residual, threshold_step_size, max_iter, max_cutbacks, use_line_search);
	if(!error)
	{
		Vector<double> de = e_omega;
		for(unsigned int m = 0; m < de.size(); ++m)
			de[m] = de[m] - e_old[m];

		const double max_step = get_maximum_step(e_old, e_omega_ref_sets, de, hidden_vars_old, x);
		if(max_step < 1.0/safety_distance_step)
		{
			cout << "Error in ScalarFunctionalLocalElimination due to too quick approach of the boundary of admissibility" << endl;
			error = true;
		}
	}
	else
	{
//		cout << "Error in ScalarFunctionalLocalElimination due to non-convergence of Newton-Raphson iteration" << endl;

		/*
		cout << "e_omega" << " " << "e_omega_ref" << endl;
		for(unsigned int m = 0; m < e_omega.size(); ++m)
			printf("%- 1.16e %- 1.16e\n", e_omega[m], e_omega_ref_sets[0][m]);
		cout << endl;
		cout << "hidden vars" << endl;
		for(unsigned int m = 0; m < hidden_vars.size(); ++m)
		{
			printf("%- 1.16e\n", hidden_vars[m]);
		}
		cout << endl;
		 */
	}
	return error;
}

template<unsigned int spacedim>
double
ScalarFunctionalLocalElimination<spacedim, spacedim>::get_maximum_step(	const Vector<double>& 			e_omega,
																		const vector<Vector<double>>&	e_omega_ref_sets,
																		const Vector<double>& 			delta_e_omega,
																		const Vector<double>& 			hidden_vars,
																		const Point<spacedim>& 			x)
const
{
	return Implementation::get_maximum_step<spacedim, spacedim>(e_omega, e_omega_ref_sets, delta_e_omega, hidden_vars, x, nullptr, scalar_functionals, map_dependent_fields, indices_not_eliminated_dependent_fields, indices_locally_eliminated_dependent_fields);
}

template<unsigned int spacedim>
void
ScalarFunctionalLocalElimination<spacedim, spacedim>::set_safety_distance(const double safety_distance)
{
	Assert(safety_distance < 1.0, ExcMessage("Safety distance must be smaller than 1.0"));
	Assert(safety_distance > 0.0, ExcMessage("Safety must be larger than 0.0"));
	this->safety_distance = safety_distance;
}

template<unsigned int spacedim>
void
ScalarFunctionalLocalElimination<spacedim, spacedim>::set_safety_distance_step(const double safety_distance_step)
{
	Assert(safety_distance_step <= 1.0, ExcMessage("Safety distance for step must be <= 1.0"));
	Assert(safety_distance_step> 0.0, ExcMessage("Safety must be larger than 0.0"));
	this->safety_distance_step = safety_distance_step;
}


template<unsigned int spacedim>
void
ScalarFunctionalLocalElimination<spacedim, spacedim>::set_threshold_residual(const double threshold_residual)
{
	this->threshold_residual = threshold_residual;
}

template<unsigned int spacedim>
void
ScalarFunctionalLocalElimination<spacedim, spacedim>::set_threshold_step_size(const double threshold_step_size)
{
	this->threshold_step_size = threshold_step_size;
}

template<unsigned int spacedim>
void
ScalarFunctionalLocalElimination<spacedim, spacedim>::set_max_iter(const unsigned int max_iter)
{
	this->max_iter = max_iter;
}

template<unsigned int spacedim>
void
ScalarFunctionalLocalElimination<spacedim, spacedim>::set_max_cutbacks(const unsigned int max_cutbacks)
{
	this->max_cutbacks = max_cutbacks;
}

template<unsigned int spacedim>
void
ScalarFunctionalLocalElimination<spacedim, spacedim>::set_use_line_search(const bool use_line_search)
{
	this->use_line_search = use_line_search;
}

template<unsigned int dim, unsigned int spacedim>
ScalarFunctionalLocalElimination<dim, spacedim>::ScalarFunctionalLocalElimination(	const vector<ScalarFunctional<dim,spacedim>*>	scalar_functionals,
																					const string 									name)
:
ScalarFunctional<dim,spacedim>(Implementation::get_e_combined(scalar_functionals), scalar_functionals[0]->domain_of_integration, scalar_functionals[0]->quadrature, name, Implementation::get_n_ref_sets(scalar_functionals), Implementation::get_n_hidden(scalar_functionals)),
scalar_functionals(scalar_functionals)
{
#ifdef DEBUG
	for(const auto& sf : scalar_functionals)
	{
		Assert( (sf->domain_of_integration == this->domain_of_integration), ExcMessage("Domains of integration of the underlying scalar functionals inconsistent!"));
		Assert( (sf->quadrature == this->quadrature), ExcMessage("Quadrature schemes of the underlying scalar functionals inconsistent!"));
		Assert( (sf->initial_vals_hidden == nullptr), ExcMessage("Initial values for the hidden variables are not currently implemented for ScalarFunctionalLocalElimination!") );
	}
#endif

	for(unsigned int df = 0; df < this->e_sigma.size(); ++df)
	{
		if(this->e_sigma[df].get_is_locally_eliminated())
			indices_locally_eliminated_dependent_fields.push_back(df);
		else
			indices_not_eliminated_dependent_fields.push_back(df);
	}

	map<DependentField<dim,spacedim>, unsigned int> map_df_index;
	for(unsigned int df = 0; df < this->e_sigma.size(); ++df)
		map_df_index[this->e_sigma[df]] = df;

	map_dependent_fields.resize(scalar_functionals.size());
	for(unsigned int sf = 0; sf < scalar_functionals.size(); ++sf)
		for(const auto& df : scalar_functionals[sf]->e_sigma)
			map_dependent_fields[sf].push_back(map_df_index[df]);

}

template<unsigned int dim, unsigned int spacedim>
bool
ScalarFunctionalLocalElimination<dim, spacedim>::get_h_sigma(Vector<double>&		 		e_sigma,
															const vector<Vector<double>>&	e_sigma_ref_sets,
															Vector<double>& 				hidden_vars,
															const Point<spacedim>& 			x,
															const Tensor<1,spacedim>& 		n,
															double& 						h_sigma,
															Vector<double>& 				h_sigma_1,
															FullMatrix<double>& 			h_sigma_2,
															const tuple<bool, bool, bool>	requested_quantities)
const
{
	Vector<double> e_old = e_sigma;
	Vector<double> hidden_vars_old = hidden_vars;
	bool error = Implementation::get_h<dim, spacedim>(e_sigma, e_sigma_ref_sets, hidden_vars, x, &n, h_sigma, h_sigma_1, h_sigma_2, requested_quantities, scalar_functionals, map_dependent_fields, indices_not_eliminated_dependent_fields, indices_locally_eliminated_dependent_fields, safety_distance, threshold_residual, threshold_step_size, max_iter, max_cutbacks, use_line_search);
	if(!error)
	{
		Vector<double> de = e_sigma;
		for(unsigned int m = 0; m < de.size(); ++m)
			de[m] = de[m] - e_old[m];
		const double max_step = get_maximum_step(e_old, e_sigma_ref_sets, de, hidden_vars_old, x, n);
		if(max_step < 1.0/safety_distance_step)
		{
			cout << "Error in ScalarFunctionalLocalElimination due to too quick approach of the boundary of admissibility" << endl;
			error = true;
		}
	}
	else
	{
//		cout << "Error in ScalarFunctionalLocalElimination due to non-convergence of Newton-Raphson iteration" << endl;

		/*
		cout << "e_omega" << " " << "e_omega_ref" << endl;
		for(unsigned int m = 0; m < e_omega.size(); ++m)
			printf("%- 1.16e %- 1.16e\n", e_omega[m], e_omega_ref_sets[0][m]);
		cout << endl;
		cout << "hidden vars" << endl;
		for(unsigned int m = 0; m < hidden_vars.size(); ++m)
		{
			printf("%- 1.16e\n", hidden_vars[m]);
		}
		cout << endl;
		 */
	}
	return error;

}

template<unsigned int dim, unsigned int spacedim>
double
ScalarFunctionalLocalElimination<dim, spacedim>::get_maximum_step(	const Vector<double>& 			e_sigma,
																	const vector<Vector<double>>&	e_sigma_ref_sets,
																	const Vector<double>& 			delta_e_sigma,
																	const Vector<double>&			hidden_vars,
																	const Point<spacedim>&			x,
																	const Tensor<1, spacedim>&		n)
const
{
	return Implementation::get_maximum_step<dim, spacedim>(e_sigma, e_sigma_ref_sets, delta_e_sigma, hidden_vars, x, &n, scalar_functionals, map_dependent_fields, indices_not_eliminated_dependent_fields, indices_locally_eliminated_dependent_fields);
}

template<unsigned int dim, unsigned int spacedim>
void
ScalarFunctionalLocalElimination<dim, spacedim>::set_safety_distance(const double safety_distance)
{
	Assert(safety_distance < 1.0, ExcMessage("Safety distance must be smaller than 1.0"));
	Assert(safety_distance > 0.0, ExcMessage("Safety must be larger than 0.0"));
	this->safety_distance = safety_distance;
}

template<unsigned int dim, unsigned int spacedim>
void
ScalarFunctionalLocalElimination<dim, spacedim>::set_safety_distance_step(const double safety_distance_step)
{
	Assert(safety_distance_step <= 1.0, ExcMessage("Safety distance for step must be <= 1.0"));
	Assert(safety_distance_step> 0.0, ExcMessage("Safety must be larger than 0.0"));
	this->safety_distance_step = safety_distance_step;
}


template<unsigned int dim, unsigned int spacedim>
void
ScalarFunctionalLocalElimination<dim, spacedim>::set_threshold_residual(const double threshold_residual)
{
	this->threshold_residual = threshold_residual;
}

template<unsigned int dim, unsigned int spacedim>
void
ScalarFunctionalLocalElimination<dim, spacedim>::set_threshold_step_size(const double threshold_step_size)
{
	this->threshold_step_size = threshold_step_size;
}


template<unsigned int dim, unsigned int spacedim>
void
ScalarFunctionalLocalElimination<dim, spacedim>::set_max_iter(const unsigned int max_iter)
{
	this->max_iter = max_iter;
}

template<unsigned int dim, unsigned int spacedim>
void
ScalarFunctionalLocalElimination<dim, spacedim>::set_max_cutbacks(const unsigned int max_cutbacks)
{
	this->max_cutbacks = max_cutbacks;
}

template<unsigned int dim, unsigned int spacedim>
void
ScalarFunctionalLocalElimination<dim, spacedim>::set_use_line_search(const bool use_line_search)
{
	this->use_line_search = use_line_search;
}

template class ScalarFunctional<3,3>;
template class ScalarFunctional<2,2>;
template class ScalarFunctional<2,3>;
template class ScalarFunctional<1,2>;
template class ScalarFunctionalLocalElimination<2,2>;
template class ScalarFunctionalLocalElimination<3,3>;
template class ScalarFunctionalLocalElimination<1,2>;
template class ScalarFunctionalLocalElimination<2,3>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
