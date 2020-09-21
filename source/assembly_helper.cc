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

#include <galerkin_tools/assembly_helper.h>

#include <math.h>
#include <fstream>

#include <deal.II/numerics/data_out.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/sparsity_tools.h>

#include <galerkin_tools/ldr.h>
#include <galerkin_tools/tools.h>
#include <galerkin_tools/two_block_matrix.h>
#include <galerkin_tools/two_block_sparsity_pattern.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/block_vector.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
AssemblyHelper<spacedim>::AssemblyHelper(	const TotalPotential<spacedim>&						total_potential,
											TriangulationSystem<spacedim>&						tria_system,
											const Mapping<spacedim, spacedim>&					mapping_domain,
											const Mapping<spacedim-1, spacedim>&				mapping_interface,
											const set<const IndependentField<0, spacedim>*>&	independent_scalars):
total_potential(total_potential),
tria_system(tria_system),
mapping_domain(&mapping_domain, typeid(*this).name()),
mapping_interface(&mapping_interface, typeid(*this).name()),
dof_handler_system(tria_system),
this_proc(tria_system.get_this_proc_n_procs().first),
n_procs(tria_system.get_this_proc_n_procs().second),
pout(cout, this_proc == 0)
{

/*****************************************************************************
 * make maps between material_ids and internal (consecutive) indices 		 *
 * 																	 		 *
 * member variables initialized:	material_id_to_internal_index_domain	 *
 * 									material_ids_to_internal_index_interface *
 *****************************************************************************/

	//set storing the different material id's of domains actually occurring in the triangulation
	set<types::material_id> domain_ids;

	for(const auto& domain_cell : tria_system.get_triangulation_domain().cell_iterators_on_level(0))
		domain_ids.insert(domain_cell->material_id());

	//set up material_id_to_internal_index_domain
	for(const auto& domain_id : domain_ids)
		material_id_to_internal_index_domain.insert(make_pair(domain_id, material_id_to_internal_index_domain.size()));
	const unsigned int n_domain_portions = material_id_to_internal_index_domain.size();

	//map storing the interface subportions associated with a particular interface portion
	//(the first element of the pair is the material_id of the - side of the interface and the second
	// the material_id of the + side)
	map<
		types::material_id,
		set<pair<const types::material_id, const types::material_id>>
	> interface_ids_to_domain_neighbors;

	for(const auto& interface_cell_domain_cells : tria_system.interface_coarse_iterators())
	{
		const types::material_id material_id=interface_cell_domain_cells.interface_cell->material_id();
		//add interface sub portion type to interface_ids_to_domain_neighbors
		if(interface_cell_domain_cells.refinement_case == InterfaceRefinementCase::at_boundary)
			interface_ids_to_domain_neighbors[material_id].insert(make_pair(interface_cell_domain_cells.domain_cell_minus->material_id(), numbers::invalid_material_id));
		else
			interface_ids_to_domain_neighbors[material_id].insert(make_pair(interface_cell_domain_cells.domain_cell_minus->material_id(), interface_cell_domain_cells.domain_cell_plus->material_id()));
	}

	//set up material_ids_to_internal_index_interface
	for(const auto& interface_id_to_domain_neighbors : interface_ids_to_domain_neighbors)
		for(const auto& neighbors : interface_id_to_domain_neighbors.second)
			material_ids_to_internal_index_interface.insert(make_pair(	make_tuple(	interface_id_to_domain_neighbors.first,
																					neighbors.first,
																					neighbors.second),
																		material_ids_to_internal_index_interface.size()));
	const unsigned int n_interface_subportions = material_ids_to_internal_index_interface.size();

/************************************************************************************************
 * create sets of independent fields, maps relating scalar functionals to total potential,		*
 * indexing of scalar functionals																*
 * 																								*
 * member variables initialized:	contributions_scalar_functionals_domain_total_potential		*
 * 									contributions_scalar_functionals_interface_total_potential	*
 * 									n_scalar_functionals_nonprimitive							*
 * 									n_scalar_functionals_primitive								*
 * 									scalar_functionals_domain_nonprimitive_indices				*
 * 									scalar_functionals_interface_nonprimitive_indices			*
 * 									scalar_functionals_domain_primitive_indices					*
 * 									scalar_functionals_interface_primitive_indices				*
 ************************************************************************************************/

	//set of independent fields defined on domain
	set<const IndependentField<spacedim, spacedim>*> u_omega_set;

	//set of independent fields defined on interface
	set<const IndependentField<spacedim-1, spacedim>*> u_sigma_set;

	//set of independent scalars
	set<const IndependentField<0,spacedim>*> C_set;

	//counter for number of scalar functionals entering into at least one primitive TotalPotentialContribution
	n_scalar_functionals_primitive = 0;

	//counter for number of scalar functionals entering into at least one non-primitive TotalPotentialContribution
	n_scalar_functionals_nonprimitive = 0;

	//loop over total potential contributions
	for(unsigned int total_potential_contribution_n = 0; total_potential_contribution_n<this->total_potential.total_potential_contributions.size(); total_potential_contribution_n++)
	{
		const auto total_potential_contribution = this->total_potential.total_potential_contributions[total_potential_contribution_n];
		const unsigned int n_H_omega = total_potential_contribution->H_omega.size();
		const unsigned int n_H_sigma = total_potential_contribution->H_sigma.size();
		const unsigned int n_C = total_potential_contribution->C.size();
		const bool is_primitive      = total_potential_contribution->is_primitive;

		//loop over domain related scalar functionals of total potential contribution
		for(unsigned int H_omega_n = 0; H_omega_n < n_H_omega; ++H_omega_n)
		{
			const auto H_omega = total_potential_contribution->H_omega[H_omega_n];
			//assign non-primitive index of scalar functional if it enters the total potential non-primitively and increment n_scalar_functionals_nonprimitive
			if( (!is_primitive) && (scalar_functionals_domain_nonprimitive_indices.find(H_omega) == scalar_functionals_domain_nonprimitive_indices.end()) )
			{
				scalar_functionals_domain_nonprimitive_indices.insert(make_pair(H_omega, n_scalar_functionals_nonprimitive));
				++n_scalar_functionals_nonprimitive;
			}
			//assign primitive index of scalar functional if it enters the total potential as a summand and increment n_scalar_functionals_primitive
			else if (scalar_functionals_domain_primitive_indices.find(H_omega) == scalar_functionals_domain_primitive_indices.end())
			{
				scalar_functionals_domain_primitive_indices.insert(make_pair(H_omega, n_scalar_functionals_primitive));
				++n_scalar_functionals_primitive;
			}
			//insert into map relating scalar functionals to total potential
			contributions_scalar_functionals_domain_total_potential[H_omega].push_back(make_pair(total_potential_contribution_n, H_omega_n));
			//insert into set of independent fields defined on the domain and set of independent scalars
			for(const auto& e_omega : H_omega->e_omega)
			{
				const auto u_omega_set_ = e_omega.get_independent_fields_domain();
				u_omega_set.insert(u_omega_set_.begin(), u_omega_set_.end());
				const auto C_set_ = e_omega.get_independent_scalars();
				C_set.insert(C_set_.begin(), C_set_.end());
			}
		}

		//loop over interface related scalar functionals of total potential contribution
		for(unsigned int H_sigma_n = 0; H_sigma_n < n_H_sigma; ++H_sigma_n)
		{
			const auto& H_sigma = this->total_potential.total_potential_contributions[total_potential_contribution_n]->H_sigma[H_sigma_n];
			//assign non-primitive index of scalar functional if it enters the total potential non-primitively and increment n_scalar_functionals_nonprimitive
			if((!is_primitive) && ( scalar_functionals_interface_nonprimitive_indices.find(H_sigma) == scalar_functionals_interface_nonprimitive_indices.end()) )
			{
				scalar_functionals_interface_nonprimitive_indices.insert(make_pair(H_sigma, n_scalar_functionals_nonprimitive));
				++n_scalar_functionals_nonprimitive;
			}
			//assign primitive_index of scalar functional if it enters the total potential as a summand and increment n_scalar_functionals_primitive
			else if (scalar_functionals_interface_primitive_indices.find(H_sigma) == scalar_functionals_interface_primitive_indices.end())
			{
				scalar_functionals_interface_primitive_indices.insert(make_pair(H_sigma, n_scalar_functionals_primitive));
				++n_scalar_functionals_primitive;
			}
			//insert into map relating scalar functionals to total potential
			contributions_scalar_functionals_interface_total_potential[H_sigma].push_back(make_pair(total_potential_contribution_n, H_sigma_n + n_H_omega));
			//insert into sets of independent fields defined on the domain and the interface, respectively, as well as into set of independent scalars
			for(const auto& e_sigma : H_sigma->e_sigma)
			{
				const auto u_omega_set_ = e_sigma.get_independent_fields_neighbors();
				u_omega_set.insert(u_omega_set_.begin(), u_omega_set_.end());
				const auto u_sigma_set_ = e_sigma.get_independent_fields_interface();
				u_sigma_set.insert(u_sigma_set_.begin(), u_sigma_set_.end());
				const auto C_set_ = e_sigma.get_independent_scalars();
				C_set.insert(C_set_.begin(), C_set_.end());
			}
		}

		//insert into set of independent scalar variables
		for(unsigned int C_n = 0; C_n < n_C; ++C_n)
			C_set.insert(total_potential_contribution->C[C_n]);

	}

	//insert C's into C_set which are not appearing in total potential
	C_set.insert(independent_scalars.begin(), independent_scalars.end());

/**************************************************************************************
 * Assign global indices to independent fields and assemble independent field vectors *
 * 																					  *
 * member variables initialized:	component_names_domain							  *
 * 									component_names_interface						  *
 * 									global_component_indices_u_omega				  *
 * 									global_component_indices_u_sigma				  *
 * 									global_indices_C								  *
 * 									u_omega											  *
 * 									u_sigma											  *
 * 									C												  *
 **************************************************************************************/

	//sort independent fields by name to make sure there is a well defined ordering
	map<string, const IndependentField<spacedim,spacedim>*> u_omega_ordered_by_name;
	for(const auto& u_omega : u_omega_set)
	{
		Assert(u_omega_ordered_by_name.find(u_omega->name) == u_omega_ordered_by_name.end(), ExcMessage("Two independent fields have the same name. This is not allowed!"));
		u_omega_ordered_by_name[u_omega->name] = u_omega;
	}
	map<string, const IndependentField<spacedim-1,spacedim>*> u_sigma_ordered_by_name;
	for(const auto& u_sigma : u_sigma_set)
	{
		Assert(u_sigma_ordered_by_name.find(u_sigma->name) == u_sigma_ordered_by_name.end(), ExcMessage("Two independent fields have the same name. This is not allowed!"));
		u_sigma_ordered_by_name[u_sigma->name] = u_sigma;
	}
	map<string, const IndependentField<0,spacedim>*> C_ordered_by_name;
	for(const auto& C : C_set){
		Assert(C_ordered_by_name.find(C->name) == C_ordered_by_name.end(), ExcMessage("Two independent fields have the same name. This is not allowed!"));
		C_ordered_by_name[C->name] = C;
	}

	//global indices of domain related independent fields
	map<const IndependentField<spacedim,spacedim>*, const unsigned int> global_indices_u_omega;

	//global indices of interface related independent fields
	map<const IndependentField<spacedim-1,spacedim>*, const unsigned int> global_indices_u_sigma;


	unsigned int global_index_u_omega = 0;
	unsigned int global_component_index_u_omega = 0;
	for(auto u_omega : u_omega_ordered_by_name)
	{
		global_indices_u_omega.insert(make_pair(u_omega.second, global_index_u_omega));
		global_component_indices_u_omega.insert(make_pair(u_omega.second, global_component_index_u_omega));
		++global_index_u_omega;
        global_component_index_u_omega += (u_omega.second)->n_components;
        for(unsigned int component=0; component<(u_omega.second)->n_components; component++)
        	component_names_domain.push_back(make_pair((u_omega.second)->name, component));
        this->u_omega.push_back(u_omega.second);
	}
	unsigned int global_index_u_sigma = 0;
	unsigned int global_component_index_u_sigma = 0;
	for(auto u_sigma : u_sigma_ordered_by_name)
	{
		global_indices_u_sigma.insert(make_pair(u_sigma.second, global_index_u_sigma));
		global_component_indices_u_sigma.insert(make_pair(u_sigma.second, global_component_index_u_sigma));
		++global_index_u_sigma;
		global_component_index_u_sigma += (u_sigma.second)->n_components;
        for(unsigned int component=0; component<(u_sigma.second)->n_components; component++)
        	component_names_interface.push_back(make_pair((u_sigma.second)->name, component));
        this->u_sigma.push_back(u_sigma.second);
	}
	unsigned int global_index_C=0;
	for(auto C : C_ordered_by_name)
	{
		global_indices_C.insert(make_pair(C.second, global_index_C));
		++global_index_C;
		this->C.push_back(C.second);
	}

	pout << "Independent fields defined on domain:" << endl;
	for(const auto& u_omega : u_omega_ordered_by_name)
		pout << u_omega.first << " (global index="<< global_indices_u_omega.at(u_omega.second) <<")" << endl;
	if(u_omega_ordered_by_name.size() == 0)
		pout << "None" << endl;

	pout << endl << "Independent fields defined on interfaces:" << endl;
	for(auto u_sigma : u_sigma_ordered_by_name)
		pout << u_sigma.first << " (global index="<< global_indices_u_sigma.at(u_sigma.second) <<")" << endl;
	if(u_sigma_ordered_by_name.size() == 0)
		pout << "None" << endl;

	pout << endl << "Scalar independent variables:" << endl;
	for(auto C : C_ordered_by_name)
		pout << C.first << " (global index="<< global_indices_C.at(C.second) <<")" << endl;
	if(C_ordered_by_name.size() == 0)
		pout << "None" << endl;
	pout << endl;

/*************************************************************************
 * Finite elements on domain						   					 *
 * 													   					 *
 * member variables initialized:	material_id_to_fe_system_id_domain   *
 * 									fe_collection_domain 				 *
 *************************************************************************/

	FE_Nothing<spacedim, spacedim> fe_nothing_domain;

	//this map contains, for each domain portion, the finite elements of the u_omega together with the number of components
	map<
		types::material_id,
		pair<
			 vector< const FiniteElement<spacedim, spacedim>* >,
			 vector< unsigned int >
		>
	> fes_domain;

	unsigned int no_elements = u_omega_set.size();
	//insert a dummy element, if there are no domain related independent fields and finite elements (only to avoid errors)
	if(no_elements == 0)
		no_elements = 1;

	//we will start with FE_Nothing elements everywhere
	vector< const FiniteElement<spacedim, spacedim>* > fe_nothings_domain(no_elements, &fe_nothing_domain);
	//the number of vector components of the respective finite elements
	vector < unsigned int > n_components_domain(fe_nothings_domain.size());
	for(const auto& u_omega : u_omega_ordered_by_name)
		n_components_domain[global_indices_u_omega.at(u_omega.second)] = (u_omega.second)->n_components;
	//that's again the case that there are no domain related independent fields and finite elements
	if(u_omega_set.size() == 0)
		n_components_domain[0] = 1;
	//assemble fes_domain
	for(const auto& u_omega : u_omega_ordered_by_name)
	{
		for(const auto& non_zero_region : (u_omega.second)->non_zero_regions)
		{
			Assert(	domain_ids.find(non_zero_region) != domain_ids.end(),
					ExcMessage("An independent variable has been defined to be non-zero on a domain portion for which no corresponding mesh exists!"));
			//insert a set of FE_Nothing elements if we see this domain portion for the first time
			if(fes_domain.find(non_zero_region) == fes_domain.end())
				fes_domain.insert({non_zero_region, make_pair(fe_nothings_domain, n_components_domain)});
			//replace the FE_Nothing by the appropriate element
			(fes_domain[non_zero_region].first)[global_indices_u_omega.at(u_omega.second)] = ((u_omega.second)->fe).get();
			//if the vector element is vector valued by itself, we have to set the multiplicity explicitly to 1 because
			//in this case a single element is used to represent all components of u_omega
			if( ((u_omega.second)->fe).get()->n_components() > 1 )
				(fes_domain[non_zero_region].second)[global_indices_u_omega.at(u_omega.second)] = 1;
		}
	}
	//now assemble fe_collection_domain
	for(const auto& fes_domain_n : fes_domain)
	{
		FESystem<spacedim, spacedim> fe_system(fes_domain_n.second.first, fes_domain_n.second.second);
		material_id_to_fe_system_id_domain.insert({fes_domain_n.first, fe_collection_domain.size()});
		fe_collection_domain.push_back(fe_system);
	}

	//make sure that there is an FESystem assigned for all material id's even if all fields are zero for a certain material id
	//(just add an element to fe_collection_domain with all elements being FE_Nothing and use this for material id's with all fields zero)
	fe_collection_domain.push_back(FESystem<spacedim, spacedim>(fe_nothings_domain, n_components_domain));
	for(const auto& domain_id : domain_ids)
		if(material_id_to_fe_system_id_domain.find(domain_id) == material_id_to_fe_system_id_domain.end())
			material_id_to_fe_system_id_domain.insert({domain_id, fe_collection_domain.size()-1});
	//make sure that there is an FESystem for the environment (which is identified by numbers::invalid_material_id)
	material_id_to_fe_system_id_domain.insert({numbers::invalid_material_id, fe_collection_domain.size()-1});

/***************************************************************************
 * Finite elements on interfaces						  				   *
 * 														  				   *
 * member variables initialized:	material_id_to_fe_system_id_interface  *
 * 									fe_collection_interface 			   *
 ***************************************************************************/

	FE_Nothing<spacedim-1, spacedim> fe_nothing_interface;
	map<
		types::material_id,
		pair<
			 vector< const FiniteElement<spacedim-1, spacedim>* >,
			 vector< unsigned int >
		>
	> fes_interface;

	no_elements = u_sigma_set.size();
	//insert a dummy element, if there are no interface related independent fields and finite elements (only to avoid errors)
	if(no_elements == 0)
		no_elements = 1;

	//we will start with FE_Nothing elements everywhere
	vector< const FiniteElement<spacedim-1, spacedim>* > fe_nothings_interface(no_elements, &fe_nothing_interface);
	//the number of vector components of the respective finite elements
	vector < unsigned int > n_components_interface(fe_nothings_interface.size());
	for(const auto& u_sigma : u_sigma_ordered_by_name)
		n_components_interface[global_indices_u_sigma.at(u_sigma.second)]=(u_sigma.second)->n_components;
	//that's again the case that there are no interface related independent fields and finite elements
	if(u_sigma_set.size() == 0)
		n_components_interface[0] = 1;
	//assemble fes_interface
	for(const auto& u_sigma : u_sigma_ordered_by_name)
	{
		for(const auto& non_zero_region : (u_sigma.second)->non_zero_regions)
		{
			Assert(	interface_ids_to_domain_neighbors.find(non_zero_region) != interface_ids_to_domain_neighbors.end(),
					ExcMessage("An independent variable has been defined to be non-zero on an interface portion for which no corresponding mesh exists!"));
			if(fes_interface.find(non_zero_region) == fes_interface.end())
				fes_interface.insert({non_zero_region, make_pair(fe_nothings_interface, n_components_interface)});
			(fes_interface[non_zero_region].first)[global_indices_u_sigma.at(u_sigma.second)] = ((u_sigma.second)->fe).get();
			//if element is vector valued, multiplicity is 1
			if( ((u_sigma.second)->fe).get()->n_components()>1 )
				(fes_interface[non_zero_region].second)[global_indices_u_sigma.at(u_sigma.second)] = 1;
		}
	}
	//now assemble fe_collection_interface
	for(const auto& fes_interface_n : fes_interface)
	{
		FESystem<spacedim-1, spacedim> fe_system(fes_interface_n.second.first, fes_interface_n.second.second);
		material_id_to_fe_system_id_interface.insert({fes_interface_n.first, fe_collection_interface.size()});
		fe_collection_interface.push_back(fe_system);
	}

	//make sure that there is an FESystem assigned for all material id's even if all fields are zero for a certain material id
	//(just add an element to fe_collection_interface with all elements being FE_Nothing and use this for material id's with all fields zero)
	fe_collection_interface.push_back(FESystem<spacedim-1, spacedim>(fe_nothings_interface, n_components_interface));
	for(const auto& interface_id_to_domain_neighbors : interface_ids_to_domain_neighbors)
		if(material_id_to_fe_system_id_interface.find(interface_id_to_domain_neighbors.first) == material_id_to_fe_system_id_interface.end())
			material_id_to_fe_system_id_interface.insert({interface_id_to_domain_neighbors.first, fe_collection_interface.size()-1});

/****************************************************************************
 * assemble domain related scalar functional vectors and fe values vectors	*
 * 																			*
 * member variables initialized:	scalar_functionals_domain				*
 * 									scalar_functionals_domain_nonprimitive	*
 * 									fe_values_domain						*
 * 									fe_values_domain_reinit					*
 * 									fe_values_domain_reinit_nonprimitive	*
 ****************************************************************************/

	scalar_functionals_domain.resize(n_domain_portions);
	scalar_functionals_domain_nonprimitive.resize(n_domain_portions);
	fe_values_domain.resize(n_domain_portions);
	fe_values_domain_reinit.resize(n_domain_portions);
	fe_values_domain_reinit_nonprimitive.resize(n_domain_portions);

	for(const auto& scalar_functional_it : contributions_scalar_functionals_domain_total_potential)
	{
		const auto scalar_functional = scalar_functional_it.first;
		for(const auto& domain_portion : scalar_functional->domain_of_integration)
		{
			Assert(	domain_ids.find(domain_portion)!=domain_ids.end(),
					ExcMessage("There is a domain portion defined in a scalar functional for which no corresponding mesh exists!"));
			const unsigned int internal_index = material_id_to_internal_index_domain[domain_portion];
			shared_ptr<FEValues<spacedim,spacedim>> fe_values(nullptr);
			//check whether there is already a scalar functional on this domain portion which is associated with the same quadrature rule
			//(if yes: the corresponding FEValues object can be recycled!)
			for(const auto& fe_values_n : fe_values_domain[internal_index])
			{
				if(fe_values_n->get_quadrature() == scalar_functional->quadrature)
				{
					fe_values = fe_values_n;
					break;
				}
			}
			//no recyclable FEValues object was found->create a new one
			if(fe_values == nullptr)
			{
				Assert(material_id_to_fe_system_id_domain.find(domain_portion) != material_id_to_fe_system_id_domain.end(), ExcMessage("Internal error!"));
				fe_values = shared_ptr<FEValues<spacedim,spacedim>>(new FEValues<spacedim, spacedim>(	*(this->mapping_domain),
																										fe_collection_domain[material_id_to_fe_system_id_domain[domain_portion]],
																										scalar_functional->quadrature,
																										update_quadrature_points|
																										update_JxW_values|
																										update_values|
																										update_gradients) );
				fe_values_domain_reinit[internal_index].insert(fe_values);
			}
			scalar_functionals_domain[internal_index].push_back(scalar_functional);
			//if the scalar functional enters non-primitively into the total potential
			if( scalar_functionals_domain_nonprimitive_indices.find(scalar_functional) != scalar_functionals_domain_nonprimitive_indices.end() )
			{
				scalar_functionals_domain_nonprimitive[internal_index].push_back(scalar_functionals_domain[internal_index].size()-1);
				fe_values_domain_reinit_nonprimitive[internal_index].insert(fe_values);
			}
			fe_values_domain[internal_index].push_back(fe_values);
		}
	}

/******************************************************************************
 * assemble interface related scalar functional vector,						  *
 * assemble interface related fe values object vectors and fe values vectors  *
 * 																			  *
 * member variables initialized:	scalar_functionals_interface			  *
 * 									scalar_functionals_interface_nonprimitive *
 * 									fe_values_interface						  *
 * 									fe_values_interface_reinit				  *
 * 									fe_values_interface_reinit_nonprimitive	  *
 ******************************************************************************/

	scalar_functionals_interface.resize(n_interface_subportions);
	scalar_functionals_interface_nonprimitive.resize(n_interface_subportions);
	fe_values_interface.resize(n_interface_subportions);
	fe_values_interface_reinit.resize(n_interface_subportions);
	fe_values_interface_reinit_nonprimitive.resize(n_interface_subportions);

	for(const auto& scalar_functional_it : contributions_scalar_functionals_interface_total_potential)
	{
		const auto scalar_functional = scalar_functional_it.first;
		for(const auto& interface_portion : scalar_functional->domain_of_integration)
		{
			Assert(	interface_ids_to_domain_neighbors.find(interface_portion) != interface_ids_to_domain_neighbors.end(),
					ExcMessage("There is an interface defined in a scalar functional for which no corresponding mesh exists!"));
			for(const auto& interface_subportion : interface_ids_to_domain_neighbors[interface_portion])
			{
				const unsigned int internal_index = material_ids_to_internal_index_interface[make_tuple(interface_portion, interface_subportion.first, interface_subportion.second)];
				shared_ptr<FEValuesInterface<spacedim>> fe_values(nullptr);
				//check whether there is already a scalar functional on this interface sub portion which is associated with the same quadrature rule
				//(if yes: the corresponding FEValues and FEFaceValues objects can be recycled!)
				for(const auto& fe_values_n : fe_values_interface[internal_index])
				{
					if(fe_values_n->get_quadrature() == scalar_functional->quadrature)
					{
						fe_values = fe_values_n;
						break;
					}
				}
				//no recyclable FEValues objects were found->create new ones
				if(fe_values == nullptr)
				{
					fe_values = shared_ptr<FEValuesInterface<spacedim>>(new FEValuesInterface<spacedim>(*(this->mapping_interface),
																										*(this->mapping_domain),
																										fe_collection_interface[material_id_to_fe_system_id_interface[interface_portion]],
																										fe_collection_domain[material_id_to_fe_system_id_domain[interface_subportion.first]],
																										fe_collection_domain[material_id_to_fe_system_id_domain[interface_subportion.second]],
																										scalar_functional->quadrature,
																										update_quadrature_points|
																										update_JxW_values|
																										update_values|
																										update_gradients,
																										update_quadrature_points|
																										update_JxW_values|
																										update_values|
																										update_gradients|
																										//note: normal vectors are always evaluated on
																										//      the FE of the - side (normal direction
																										//      is not consistently defined on surface
																										//      mesh)
																										update_normal_vectors,
																										update_quadrature_points|
																										update_JxW_values|
																										update_values|
																										update_gradients));
					fe_values_interface_reinit[internal_index].insert(fe_values);
				}
				scalar_functionals_interface[internal_index].push_back(scalar_functional);
				//if the scalar functional enters non-primitively into the total potential
				if( scalar_functionals_interface_nonprimitive_indices.find(scalar_functional) != scalar_functionals_interface_nonprimitive_indices.end() )
				{
					scalar_functionals_interface_nonprimitive[internal_index].push_back(scalar_functionals_interface[internal_index].size()-1);
					fe_values_interface_reinit_nonprimitive[internal_index].insert(fe_values);
				}
				fe_values_interface[internal_index].push_back(fe_values);
			}
		}
	}

/****************************************************************************************
 * generate the data structures which relate the global components to shape functions	*
 *																						*
 * member variables initialized:	components_to_shapefuns_domain						*
 * 									components_to_shapefuns_interface					*
 * 									components_to_shapefuns_domain_facewise				*
 ****************************************************************************************/

	components_to_shapefuns_domain.resize(fe_collection_domain.size());
	for(unsigned int fe_system_n = 0; fe_system_n<fe_collection_domain.size(); ++fe_system_n)
	{
		components_to_shapefuns_domain[fe_system_n].resize(fe_collection_domain.n_components());
		for(unsigned int shapefun = 0; shapefun<fe_collection_domain[fe_system_n].n_dofs_per_cell(); ++shapefun)
		{
			const auto& nonzero_components = fe_collection_domain[fe_system_n].get_nonzero_components(shapefun);
			for(unsigned int component = 0; component<nonzero_components.size(); ++component)
				if(nonzero_components[component])
					components_to_shapefuns_domain[fe_system_n][ component ].push_back(shapefun);
		}
	}

	components_to_shapefuns_interface.resize(fe_collection_interface.size());
	for(unsigned int fe_system_n = 0; fe_system_n<fe_collection_interface.size(); ++fe_system_n)
	{
		components_to_shapefuns_interface[fe_system_n].resize(fe_collection_interface.n_components());
		for(unsigned int shapefun = 0; shapefun<fe_collection_interface[fe_system_n].n_dofs_per_cell(); ++shapefun)
		{
			const auto& nonzero_components = fe_collection_interface[fe_system_n].get_nonzero_components(shapefun);
			for(unsigned int component = 0; component<nonzero_components.size(); ++component)
				if(nonzero_components[component])
					components_to_shapefuns_interface[fe_system_n][ component ].push_back(shapefun);
		}
	}

	components_to_shapefuns_domain_facewise.resize(fe_collection_domain.size());
	for(unsigned int fe_system_n = 0; fe_system_n < fe_collection_domain.size(); ++fe_system_n)
	{
		components_to_shapefuns_domain_facewise[fe_system_n].resize(fe_collection_domain.n_components());
		for(unsigned int component = 0; component<fe_collection_domain.n_components(); ++component)
			components_to_shapefuns_domain_facewise[fe_system_n][component].resize(GeometryInfo<spacedim>::faces_per_cell);
		for(unsigned int face = 0; face<GeometryInfo<spacedim>::faces_per_cell; ++face)
		{
			for(unsigned int shapefun = 0; shapefun<fe_collection_domain[fe_system_n].dofs_per_cell; ++shapefun)
			{
				const unsigned int base_element = fe_collection_domain[fe_system_n].system_to_base_index(shapefun).first.first;
				//we check here that the finite element is defined in terms of support points, because it seems to be possible
				//that has_support_on_face() returns true even if a finite element has no support points defined
				if(fe_collection_domain[fe_system_n].base_element(base_element).has_support_points() && fe_collection_domain[fe_system_n].has_support_on_face(shapefun, face))
				{
					//we can safely do this because we have checked that the finite element is defined in terms of support points
					const unsigned int component = fe_collection_domain[fe_system_n].system_to_component_index(shapefun).first;
					components_to_shapefuns_domain_facewise[fe_system_n][ component ][face].push_back(shapefun);
				}
			}
		}
	}

/****************************************************************************************
 * converting the dependent field definitions into a format suitable 					*
 * for assembly of the finite element system                         					*
 * member variables initialized: coupled_dof_indices_scalar_functionals_domain			*
 *								 coupled_dof_indices_scalar_functionals_interface		*
 *								 coupled_dof_indices_scalar_functionals_interface_minus	*
 *								 coupled_dof_indices_scalar_functionals_interface_plus	*
 *								 coupled_C_indices_scalar_functionals_domain			*
 *								 coupled_C_indices_scalar_functionals_interface			*
 *								 a_omega												*
 *								 b_omega												*
 *								 c_omega												*
 *								 d_omega												*
 *								 a_sigma												*
 *								 b_sigma												*
 *								 a_minus												*
 *								 b_minus												*
 *								 a_plus													*
 *								 b_plus													*
 *								 c_sigma												*
 *								 d_sigma												*
 ****************************************************************************************/

	convert_dependent_fields_to_shapefunctions();

/*****************************************************************************************************
 * distribute dof's and make sure that dof's are redistributed after any change of the triangulation *
 *****************************************************************************************************/

	distribute_dofs();
	tria_listeners.push_back(this->tria_system.post_refinement.connect(boost::bind(&AssemblyHelper<spacedim>::distribute_dofs, this)));
	//this is important to make sure that distribute_dofs() is really the last call after the triangulation is changed
	//(in particular this ensures that distribute_dofs() is called AFTER the post_refinement_action() of the dof handlers)
	tria_system.push_tria_listeners_to_end();

/*******************************
 * initialize hidden variables *
 *******************************/

	initialize_hidden_variables();
}

template<unsigned int spacedim>
AssemblyHelper<spacedim>::~AssemblyHelper()
{
	//deallocate memory of hidden variables on domain
	for(const auto& domain_cell : tria_system.get_triangulation_domain().active_cell_iterators())
	{
		if(domain_cell->user_pointer() != nullptr)
		{
			const auto cell_data = static_cast<vector<vector<Vector<double>>>*>(domain_cell->user_pointer());
			delete cell_data;
			domain_cell->clear_user_pointer();
		}
	}

	//deallocate memory of hidden variables on interface
	for(const auto& interface_cell_domain_cells : tria_system.interface_active_iterators())
	{
		if(interface_cell_domain_cells.interface_cell->user_pointer() != nullptr)
		{
			const auto cell_data = static_cast<vector<vector<Vector<double>>>*>(interface_cell_domain_cells.interface_cell->user_pointer());
			delete cell_data;
			interface_cell_domain_cells.interface_cell->clear_user_pointer();
		}
	}

	for(auto &connection : tria_listeners)
		connection.disconnect();
	tria_listeners.clear();
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::convert_dependent_fields_to_shapefunctions()
{

	//inverse of material_id_to_internal_index_domain
	map<unsigned int, const types::material_id> internal_index_to_material_id_domain;
	for(const auto& material_id_to_internal_index_domain_n : material_id_to_internal_index_domain)
		internal_index_to_material_id_domain.insert(make_pair(material_id_to_internal_index_domain_n.second, material_id_to_internal_index_domain_n.first));

	//number of different domain portions
	const unsigned int n_domain_portions = material_id_to_internal_index_domain.size();

	a_omega.resize(n_domain_portions);
	b_omega.resize(n_domain_portions);
	c_omega.resize(n_domain_portions);
	d_omega.resize(n_domain_portions);
	coupled_dof_indices_scalar_functionals_domain.resize(n_domain_portions);
	coupled_C_indices_scalar_functionals_domain.resize(n_domain_portions);

	for(const auto& internal_index_material_id_n : internal_index_to_material_id_domain)
	{

		const unsigned int internal_index = internal_index_material_id_n.first;
		const unsigned int fe_system_id = material_id_to_fe_system_id_domain.at(internal_index_material_id_n.second);
		const unsigned int n_scalar_functionals = scalar_functionals_domain[internal_index].size();

		a_omega[internal_index].resize(n_scalar_functionals);
		b_omega[internal_index].resize(n_scalar_functionals);
		c_omega[internal_index].resize(n_scalar_functionals);
		d_omega[internal_index].resize(n_scalar_functionals);

		for(unsigned int scalar_functional_n = 0; scalar_functional_n < n_scalar_functionals; ++scalar_functional_n)
		{
			set<unsigned int> coupled_dofs_omega, coupled_dofs_C;
			const unsigned int n_dependent_fields = scalar_functionals_domain[internal_index][scalar_functional_n]->e_omega.size();

			a_omega[internal_index][scalar_functional_n].resize(n_dependent_fields);
			b_omega[internal_index][scalar_functional_n].resize(n_dependent_fields);
			c_omega[internal_index][scalar_functional_n].resize(n_dependent_fields);
			d_omega[internal_index][scalar_functional_n].resize(n_dependent_fields);

			for(unsigned int dependent_field_n = 0; dependent_field_n < n_dependent_fields; ++dependent_field_n)
			{
				//a_omega and b_omega
				const auto& terms_domain = scalar_functionals_domain[internal_index][scalar_functional_n]->e_omega[dependent_field_n].get_terms_domain();
				map<unsigned int, const double> a_;
				map<pair<const unsigned int, const unsigned int>, const double> b_;
				for(const auto& term : terms_domain)
				{
					const unsigned int global_component_index = global_component_indices_u_omega.at(term.independent_field) + term.component;
					if(term.n_derivatives() == 0)
						a_.insert(make_pair(global_component_index, term.coefficient));
					else if(term.n_derivatives() == 1)
						b_.insert(make_pair(make_pair(global_component_index, term.derivatives[0]), term.coefficient));
					else
						Assert(term.n_derivatives() < 2, ExcMessage("Dependent fields with second or higher derivatives are currently no allowed!"))
				}

				for(const auto& term_a_n : a_)
				{
					a_omega[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_a_n.second, term_a_n.first, vector<unsigned int>()));
					const auto& shapefuns = components_to_shapefuns_domain[fe_system_id][term_a_n.first];
					for(unsigned int shapefun_n = 0; shapefun_n < shapefuns.size(); ++shapefun_n)
					{
						get<2>(a_omega[internal_index][scalar_functional_n][dependent_field_n].back()).push_back(shapefuns[shapefun_n]);
						coupled_dofs_omega.insert(shapefuns[shapefun_n]);
					}
				}
				for(const auto& term_b_n : b_)
				{
					b_omega[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_b_n.second, term_b_n.first.first, term_b_n.first.second, vector<unsigned int>()));
					const auto& shapefuns = components_to_shapefuns_domain[fe_system_id][term_b_n.first.first];
					for(unsigned int shapefun_n = 0; shapefun_n < shapefuns.size(); ++shapefun_n)
					{
						get<3>(b_omega[internal_index][scalar_functional_n][dependent_field_n].back()).push_back(shapefuns[shapefun_n]);
						coupled_dofs_omega.insert(shapefuns[shapefun_n]);
					}
				}

				//c_omega
				const auto& terms_independent_scalars = scalar_functionals_domain[internal_index][scalar_functional_n]->e_omega[dependent_field_n].get_terms_independent_scalars();
				map<unsigned int, const double> c_;
				for(const auto& term : terms_independent_scalars)
				{
					const unsigned int global_component_index = global_indices_C.at(term.independent_field);
					c_.insert(make_pair(global_component_index, term.coefficient));
				}
				for(const auto& term_c_n : c_)
				{
					c_omega[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_c_n.second, term_c_n.first));
					coupled_dofs_C.insert(term_c_n.first);
				}

				//d_omega
				d_omega[internal_index][scalar_functional_n][dependent_field_n] = scalar_functionals_domain[internal_index][scalar_functional_n]->e_omega[dependent_field_n].get_constant();
			}

			//set up coupled_dof_indices_scalar_functionals_domain (and, for temporary use, the inverse thereof)
			coupled_dof_indices_scalar_functionals_domain[internal_index].resize(n_scalar_functionals);
			map<unsigned int, unsigned int>  dof_indices_cell_to_dof_indices_scalar_functional_domain;
			for(const auto& coupled_dof_n : coupled_dofs_omega)
			{
				coupled_dof_indices_scalar_functionals_domain[internal_index][scalar_functional_n].push_back(coupled_dof_n);
				dof_indices_cell_to_dof_indices_scalar_functional_domain[coupled_dof_n]=coupled_dof_indices_scalar_functionals_domain[internal_index][scalar_functional_n].size()-1;
			}

			//set up coupled_C_indices_scalar_functionals_domain (and, for temporary use, the inverse thereof)
			coupled_C_indices_scalar_functionals_domain[internal_index].resize(n_scalar_functionals);
			map<unsigned int, const unsigned int>  dof_indices_C_to_dof_indices_scalar_functional_domain;
			for(const auto& coupled_dof_n : coupled_dofs_C)
			{
				coupled_C_indices_scalar_functionals_domain[internal_index][scalar_functional_n].push_back(coupled_dof_n);
				dof_indices_C_to_dof_indices_scalar_functional_domain.insert(make_pair(coupled_dof_n, coupled_C_indices_scalar_functionals_domain[internal_index][scalar_functional_n].size()-1));
			}

			//convert a, b and c to scalar functional related dof indexing (initially it was set up with cell related dof indexing)
			for(unsigned int dependent_field_n = 0; dependent_field_n<n_dependent_fields; ++dependent_field_n)
			{
				for(auto& a_omega_n : a_omega[internal_index][scalar_functional_n][dependent_field_n])
					for(auto& a_omega_n_sf : get<2>(a_omega_n))
						a_omega_n_sf = dof_indices_cell_to_dof_indices_scalar_functional_domain[a_omega_n_sf];
				for(auto& b_omega_n : b_omega[internal_index][scalar_functional_n][dependent_field_n])
					for(auto& b_omega_n_sf : get<3>(b_omega_n))
						b_omega_n_sf = dof_indices_cell_to_dof_indices_scalar_functional_domain[b_omega_n_sf];
				for(auto& c_omega_n : c_omega[internal_index][scalar_functional_n][dependent_field_n])
					get<1>(c_omega_n) = dof_indices_C_to_dof_indices_scalar_functional_domain[get<1>(c_omega_n)];
			}
		}
	}

	//inverse of material_ids_to_internal_index_interface
	map<unsigned int, tuple<const types::material_id, const types::material_id, const types::material_id>> internal_index_to_material_id_interface;
	for(const auto& it : material_ids_to_internal_index_interface)
		internal_index_to_material_id_interface.insert(make_pair(it.second, it.first));

	//number of different interface subportions
	const unsigned int n_interface_sub_portions = material_ids_to_internal_index_interface.size();

	a_sigma.resize(n_interface_sub_portions);
	b_sigma.resize(n_interface_sub_portions);
	a_minus.resize(n_interface_sub_portions);
	b_minus.resize(n_interface_sub_portions);
	a_plus.resize(n_interface_sub_portions);
	b_plus.resize(n_interface_sub_portions);
	c_sigma.resize(n_interface_sub_portions);
	d_sigma.resize(n_interface_sub_portions);
	coupled_dof_indices_scalar_functionals_interface.resize(n_interface_sub_portions);
	coupled_dof_indices_scalar_functionals_interface_minus.resize(n_interface_sub_portions);
	coupled_dof_indices_scalar_functionals_interface_plus.resize(n_interface_sub_portions);
	coupled_C_indices_scalar_functionals_interface.resize(n_interface_sub_portions);

	for(const auto& internal_index_material_id_n : internal_index_to_material_id_interface)
	{

		const unsigned int internal_index = internal_index_material_id_n.first;
		const unsigned int fe_system_id_interface = material_id_to_fe_system_id_interface.at(get<0>(internal_index_material_id_n.second));
		const unsigned int fe_system_id_minus = material_id_to_fe_system_id_domain.at(get<1>(internal_index_material_id_n.second));
		const unsigned int fe_system_id_plus = material_id_to_fe_system_id_domain.at(get<2>(internal_index_material_id_n.second));
		const unsigned int n_scalar_functionals = scalar_functionals_interface[internal_index].size();

		a_sigma[internal_index].resize(n_scalar_functionals);
		b_sigma[internal_index].resize(n_scalar_functionals);
		a_minus[internal_index].resize(n_scalar_functionals);
		b_minus[internal_index].resize(n_scalar_functionals);
		a_plus[internal_index].resize(n_scalar_functionals);
		b_plus[internal_index].resize(n_scalar_functionals);
		c_sigma[internal_index].resize(n_scalar_functionals);
		d_sigma[internal_index].resize(n_scalar_functionals);

		for(unsigned int scalar_functional_n = 0; scalar_functional_n < n_scalar_functionals; ++scalar_functional_n)
		{

			set<unsigned int> coupled_dofs_minus, coupled_dofs_plus, coupled_dofs_sigma, coupled_dofs_C;
			const unsigned int n_dependent_fields = scalar_functionals_interface[internal_index][scalar_functional_n]->e_sigma.size();

			a_sigma[internal_index][scalar_functional_n].resize(n_dependent_fields);
			b_sigma[internal_index][scalar_functional_n].resize(n_dependent_fields);
			a_minus[internal_index][scalar_functional_n].resize(n_dependent_fields);
			b_minus[internal_index][scalar_functional_n].resize(n_dependent_fields);
			a_plus[internal_index][scalar_functional_n].resize(n_dependent_fields);
			b_plus[internal_index][scalar_functional_n].resize(n_dependent_fields);
			c_sigma[internal_index][scalar_functional_n].resize(n_dependent_fields);
			d_sigma[internal_index][scalar_functional_n].resize(n_dependent_fields);

			for(unsigned int dependent_field_n = 0; dependent_field_n < n_dependent_fields; ++dependent_field_n)
			{
				//a_sigma and b_sigma
				const auto& terms_interface = scalar_functionals_interface[internal_index][scalar_functional_n]->e_sigma[dependent_field_n].get_terms_interface();
				map<unsigned int, const double> a_;
				map<pair<const unsigned int, const unsigned int>, double> b_;
				for(const auto& term : terms_interface)
				{
					const unsigned int global_component_index = global_component_indices_u_sigma.at(term.independent_field) + term.component;
					if(term.n_derivatives() == 0)
						a_.insert(make_pair(global_component_index, term.coefficient));
					else if(term.n_derivatives() == 1)
						b_.insert(make_pair(make_pair(global_component_index, term.derivatives[0]), term.coefficient));
					else
						Assert(term.n_derivatives() < 2, ExcMessage("Dependent fields with second or higher derivatives are currently no allowed!"))
				}

				for(const auto& term_a_n : a_)
				{
					a_sigma[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_a_n.second, term_a_n.first, vector<unsigned int>()));
					const auto& shapefuns = components_to_shapefuns_interface[fe_system_id_interface][term_a_n.first];
					for(unsigned int shapefun_n = 0; shapefun_n < shapefuns.size(); ++shapefun_n)
					{
						get<2>(a_sigma[internal_index][scalar_functional_n][dependent_field_n].back()).push_back(shapefuns[shapefun_n]);
						coupled_dofs_sigma.insert(shapefuns[shapefun_n]);
					}
				}
				for(const auto& term_b_n : b_)
				{
					b_sigma[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_b_n.second, term_b_n.first.first, term_b_n.first.second, vector<unsigned int>()));
					const auto& shapefuns = components_to_shapefuns_interface[fe_system_id_interface][term_b_n.first.first];
					for(unsigned int shapefun_n = 0; shapefun_n < shapefuns.size(); ++shapefun_n)
					{
						get<3>(b_sigma[internal_index][scalar_functional_n][dependent_field_n].back()).push_back(shapefuns[shapefun_n]);
						coupled_dofs_sigma.insert(shapefuns[shapefun_n]);
					}
				}

				//a_minus and b_minus
				const auto& terms_minus = scalar_functionals_interface[internal_index][scalar_functional_n]->e_sigma[dependent_field_n].get_terms_neighbor(InterfaceSide::minus);
				a_.clear();
				b_.clear();
				for(const auto& term : terms_minus)
				{
					const unsigned int global_component_index = global_component_indices_u_omega.at(term.independent_field) + term.component;
					if(term.n_derivatives() == 0)
						a_.insert(make_pair(global_component_index, term.coefficient));
					else if(term.n_derivatives() == 1)
						b_.insert(make_pair(make_pair(global_component_index, term.derivatives[0]), term.coefficient));
					else
						Assert(term.n_derivatives() < 2, ExcMessage("Dependent fields with second or higher derivatives are currently no allowed!"))
				}

				for(const auto& term_a_n : a_)
				{
					a_minus[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_a_n.second, term_a_n.first, vector<unsigned int>()));
					const auto& shapefuns = components_to_shapefuns_domain[fe_system_id_minus][term_a_n.first];
					for(unsigned int shapefun_n = 0; shapefun_n < shapefuns.size(); ++shapefun_n)
					{
						get<2>(a_minus[internal_index][scalar_functional_n][dependent_field_n].back()).push_back(shapefuns[shapefun_n]);
						coupled_dofs_minus.insert(shapefuns[shapefun_n]);
					}
				}
				for(const auto& term_b_n : b_)
				{
					b_minus[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_b_n.second, term_b_n.first.first, term_b_n.first.second, vector<unsigned int>()));
					const auto& shapefuns = components_to_shapefuns_domain[fe_system_id_minus][term_b_n.first.first];
					for(unsigned int shapefun_n = 0; shapefun_n < shapefuns.size(); ++shapefun_n)
					{
						get<3>(b_minus[internal_index][scalar_functional_n][dependent_field_n].back()).push_back(shapefuns[shapefun_n]);
						coupled_dofs_minus.insert(shapefuns[shapefun_n]);
					}
				}

				//a_plus and b_plus
				const auto& terms_plus = scalar_functionals_interface[internal_index][scalar_functional_n]->e_sigma[dependent_field_n].get_terms_neighbor(InterfaceSide::plus);
				a_.clear();
				b_.clear();
				for(const auto& term : terms_plus)
				{
					const unsigned int global_component_index = global_component_indices_u_omega.at(term.independent_field) + term.component;
					if(term.n_derivatives() == 0)
						a_.insert(make_pair(global_component_index, term.coefficient));
					else if(term.n_derivatives() == 1)
						b_.insert(make_pair(make_pair(global_component_index, term.derivatives[0]), term.coefficient));
					else
						Assert(term.n_derivatives() < 2, ExcMessage("Dependent fields with second or higher derivatives are currently no allowed!"))
				}

				for(const auto& term_a_n : a_)
				{
					a_plus[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_a_n.second, term_a_n.first, vector<unsigned int>()));
					const auto& shapefuns = components_to_shapefuns_domain[fe_system_id_plus][term_a_n.first];
					for(unsigned int shapefun_n = 0; shapefun_n < shapefuns.size(); ++shapefun_n)
					{
						get<2>(a_plus[internal_index][scalar_functional_n][dependent_field_n].back()).push_back(shapefuns[shapefun_n]);
						coupled_dofs_plus.insert(shapefuns[shapefun_n]);
					}
				}
				for(const auto& term_b_n : b_)
				{
					b_plus[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_b_n.second, term_b_n.first.first, term_b_n.first.second, vector<unsigned int>()));
					const auto& shapefuns = components_to_shapefuns_domain[fe_system_id_plus][term_b_n.first.first];
					for(unsigned int shapefun_n = 0; shapefun_n < shapefuns.size(); ++shapefun_n)
					{
						get<3>(b_plus[internal_index][scalar_functional_n][dependent_field_n].back()).push_back(shapefuns[shapefun_n]);
						coupled_dofs_plus.insert(shapefuns[shapefun_n]);
					}
				}

				//c_omega
				const auto& terms_independent_scalars = scalar_functionals_interface[internal_index][scalar_functional_n]->e_sigma[dependent_field_n].get_terms_independent_scalars();
				map<unsigned int, const double> c_;
				for(const auto& term : terms_independent_scalars)
				{
					const unsigned int global_component_index = global_indices_C.at(term.independent_field);
					c_.insert(make_pair(global_component_index, term.coefficient));
				}
				for(const auto& term_c_n : c_)
				{
					c_sigma[internal_index][scalar_functional_n][dependent_field_n].push_back(make_tuple(term_c_n.second, term_c_n.first));
					coupled_dofs_C.insert(term_c_n.first);
				}

				//d_sigma
				d_sigma[internal_index][scalar_functional_n][dependent_field_n]=scalar_functionals_interface[internal_index][scalar_functional_n]->e_sigma[dependent_field_n].get_constant();
			}

			//set up coupled_dof_indices_scalar_functionals_interface, coupled_dof_indices_scalar_functionals_interface_minus, coupled_dof_indices_scalar_functionals_interface_plus
			//(and, for temporary use, the inverses thereof)
			coupled_dof_indices_scalar_functionals_interface[internal_index].resize(n_scalar_functionals);
			coupled_dof_indices_scalar_functionals_interface_minus[internal_index].resize(n_scalar_functionals);
			coupled_dof_indices_scalar_functionals_interface_plus[internal_index].resize(n_scalar_functionals);
			map<unsigned int, const unsigned int>  dof_indices_cell_to_dof_indices_scalar_functional_interface;
			map<unsigned int, const unsigned int>  dof_indices_cell_to_dof_indices_scalar_functional_minus;
			map<unsigned int, const unsigned int>  dof_indices_cell_to_dof_indices_scalar_functional_plus;
			for(const auto& coupled_dof_n : coupled_dofs_sigma)
			{
				coupled_dof_indices_scalar_functionals_interface[internal_index][scalar_functional_n].push_back(coupled_dof_n);
				dof_indices_cell_to_dof_indices_scalar_functional_interface.insert(make_pair(coupled_dof_n, coupled_dof_indices_scalar_functionals_interface[internal_index][scalar_functional_n].size()-1));
			}
			for(const auto& coupled_dof_n : coupled_dofs_minus)
			{
				coupled_dof_indices_scalar_functionals_interface_minus[internal_index][scalar_functional_n].push_back(coupled_dof_n);
				dof_indices_cell_to_dof_indices_scalar_functional_minus.insert(make_pair(coupled_dof_n, coupled_dof_indices_scalar_functionals_interface_minus[internal_index][scalar_functional_n].size()-1));
			}
			for(const auto& coupled_dof_n : coupled_dofs_plus)
			{
				coupled_dof_indices_scalar_functionals_interface_plus[internal_index][scalar_functional_n].push_back(coupled_dof_n);
				dof_indices_cell_to_dof_indices_scalar_functional_plus.insert(make_pair(coupled_dof_n, coupled_dof_indices_scalar_functionals_interface_plus[internal_index][scalar_functional_n].size()-1));
			}

			//set up coupled_C_indices_scalar_functionals_interface (and, for temporary use, the inverse thereof)
			coupled_C_indices_scalar_functionals_interface[internal_index].resize(n_scalar_functionals);
			map<unsigned int, const unsigned int>  dof_indices_C_to_dof_indices_scalar_functional_interface;
			for(const auto& coupled_dof_n : coupled_dofs_C)
			{
				coupled_C_indices_scalar_functionals_interface[internal_index][scalar_functional_n].push_back(coupled_dof_n);
				dof_indices_C_to_dof_indices_scalar_functional_interface.insert(make_pair(coupled_dof_n, coupled_C_indices_scalar_functionals_interface[internal_index][scalar_functional_n].size()-1));
			}

			//convert a, b and c to scalar functional related dof indexing (initially it was set up with cell related dof indexing)
			for(unsigned int dependent_field_n = 0; dependent_field_n < n_dependent_fields; ++dependent_field_n)
			{
				for(auto& a_sigma_n : a_sigma[internal_index][scalar_functional_n][dependent_field_n])
					for(auto& a_sigma_n_sf : get<2>(a_sigma_n))
						a_sigma_n_sf = dof_indices_cell_to_dof_indices_scalar_functional_interface[a_sigma_n_sf];
				for(auto& b_sigma_n : b_sigma[internal_index][scalar_functional_n][dependent_field_n])
					for(auto& b_sigma_n_sf : get<3>(b_sigma_n))
						b_sigma_n_sf = dof_indices_cell_to_dof_indices_scalar_functional_interface[b_sigma_n_sf];
				for(auto& a_minus_n : a_minus[internal_index][scalar_functional_n][dependent_field_n])
					for(auto& a_minus_n_sf : get<2>(a_minus_n))
						a_minus_n_sf = dof_indices_cell_to_dof_indices_scalar_functional_minus[a_minus_n_sf];
				for(auto& b_minus_n : b_minus[internal_index][scalar_functional_n][dependent_field_n])
					for(auto& b_minus_n_sf : get<3>(b_minus_n))
						b_minus_n_sf = dof_indices_cell_to_dof_indices_scalar_functional_minus[b_minus_n_sf];
				for(auto& a_plus_n : a_plus[internal_index][scalar_functional_n][dependent_field_n])
					for(auto& a_plus_n_sf : get<2>(a_plus_n))
						a_plus_n_sf = dof_indices_cell_to_dof_indices_scalar_functional_plus[a_plus_n_sf];
				for(auto& b_plus_n : b_plus[internal_index][scalar_functional_n][dependent_field_n])
					for(auto& b_plus_n_sf : get<3>(b_plus_n))
						b_plus_n_sf = dof_indices_cell_to_dof_indices_scalar_functional_plus[b_plus_n_sf];
				for(auto& c_sigma_n : c_sigma[internal_index][scalar_functional_n][dependent_field_n])
					get<1>(c_sigma_n) = dof_indices_C_to_dof_indices_scalar_functional_interface[get<1>(c_sigma_n)];
			}
		}
	}

	return;
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::initialize_hidden_variables()
const
{
	//allocate memory for hidden variables on domain and set initial values
	for(const auto& domain_cell : tria_system.get_triangulation_domain().active_cell_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			//structure of hidden variables for a cell: scalar functional->quad point->data
			auto cell_data = new vector<vector<Vector<double>>>();
			domain_cell->set_user_pointer(cell_data);
			const unsigned int internal_index = material_id_to_internal_index_domain.at(domain_cell->material_id());
			cell_data->resize(scalar_functionals_domain[internal_index].size());
			for(unsigned int scalar_functional_n = 0; scalar_functional_n < scalar_functionals_domain[internal_index].size(); ++scalar_functional_n)
			{
				const auto scalar_functional = scalar_functionals_domain[internal_index][scalar_functional_n];
				const unsigned int n_hidden = scalar_functional->n_hidden;
				const unsigned int n_q_points = fe_values_domain[internal_index][scalar_functional_n]->n_quadrature_points;
				(*cell_data)[scalar_functional_n].resize(n_q_points);
				for(unsigned int q_point = 0; q_point < n_q_points; ++q_point)
				{
					(*cell_data)[scalar_functional_n][q_point].reinit(n_hidden);
					if(scalar_functional->initial_vals_hidden != nullptr)
					{
						fe_values_domain[internal_index][scalar_functional_n]->reinit(domain_cell);
						scalar_functional->initial_vals_hidden->vector_value(	fe_values_domain[internal_index][scalar_functional_n]->quadrature_point(q_point),
																				(*cell_data)[scalar_functional_n][q_point]);
					}
				}
			}
		}
	}

	//allocate memory for hidden variables on interface and set initial values
	for(const auto& interface_cell_domain_cells : tria_system.interface_active_iterators())
	{
		if(interface_cell_domain_cells.interface_cell->is_locally_owned())
		{
			//structure of hidden variables for a cell: scalar functional->quad point->data
			auto cell_data = new vector<vector<Vector<double>>>();
			interface_cell_domain_cells.interface_cell->set_user_pointer(cell_data);
			const unsigned int internal_index = material_ids_to_internal_index_interface.at(interface_cell_domain_cells.get_material_ids());
			cell_data->resize(scalar_functionals_interface[internal_index].size());
			for(unsigned int scalar_functional_n = 0; scalar_functional_n < scalar_functionals_interface[internal_index].size(); ++scalar_functional_n)
			{
				const auto scalar_functional = scalar_functionals_interface[internal_index][scalar_functional_n];
				const unsigned int n_hidden=scalar_functionals_interface[internal_index][scalar_functional_n]->n_hidden;
				const unsigned int n_q_points=fe_values_interface[internal_index][scalar_functional_n]->n_quadrature_points;
				(*cell_data)[scalar_functional_n].resize(n_q_points);
				for(unsigned int q_point = 0; q_point < n_q_points; ++q_point)
				{
					(*cell_data)[scalar_functional_n][q_point].reinit(n_hidden);
					if(scalar_functional->initial_vals_hidden != nullptr)
					{
						fe_values_interface[internal_index][scalar_functional_n]->reinit(interface_cell_domain_cells);
						scalar_functional->initial_vals_hidden->vector_value(	fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_interface().quadrature_point(q_point),
																				(*cell_data)[scalar_functional_n][q_point]);
					}
				}
			}
		}
	}
}

template<unsigned int spacedim>
template<class VectorType>
void
AssemblyHelper<spacedim>::get_initial_fields_vector(	VectorType& 								initial_fields,
														const dealii::AffineConstraints<double>*	constraints)
const
{

	Assert( initial_fields.size() == system_size(), ExcMessage("Vector not initialized to correct size"));

	//support point in real space, unit support point of domain cell
	Point<spacedim> support_point, unit_support_point_domain;

	//vector mapping local dof indices to global dof indices
	vector<unsigned int> dof_indices_local_global;
	dof_indices_local_global.reserve(fe_collection_domain.max_dofs_per_cell());

	for(const auto& domain_cell : dof_handler_system.domain_active_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			const unsigned int fe_system_id = material_id_to_fe_system_id_domain.at(domain_cell->material_id());
			dof_indices_local_global.resize(domain_cell->get_fe().dofs_per_cell);
			domain_cell.get_dof_indices(dof_indices_local_global);
			for(const auto& u_omega_n : u_omega)
			{
				if(u_omega_n->initial_vals != nullptr)
				{
					for(unsigned int component = 0; component < u_omega_n->n_components; ++component)
					{
						const unsigned int global_component_index = global_component_indices_u_omega.at(u_omega_n) + component;
						const auto& shapefuns = components_to_shapefuns_domain[fe_system_id][global_component_index];
						for(const auto& shapefun : shapefuns)
						{
							unit_support_point_domain = domain_cell->get_fe().unit_support_point(shapefun);
							support_point = mapping_domain->transform_unit_to_real_cell(domain_cell, unit_support_point_domain);
							const double u_omega_val = u_omega_n->initial_vals->value(support_point, component);
							const unsigned int global_dof_index = dof_indices_local_global[shapefun];
							initial_fields[global_dof_index] = u_omega_val;
						}
					}
				}
			}
		}
	}

	//unit support point of interface cell
	Point<spacedim-1> unit_support_point_interface;

	//re-use the dof_indices_local_global vector
	dof_indices_local_global.reserve(fe_collection_interface.max_dofs_per_cell());

	for(const auto& interface_cell_domain_cells : dof_handler_system.interface_active_iterators())
	{
		if(interface_cell_domain_cells.interface_cell->is_locally_owned())
		{
			const auto& interface_cell = interface_cell_domain_cells.interface_cell;
			const unsigned int fe_system_id = material_id_to_fe_system_id_interface.at(interface_cell->material_id());
			interface_cell_domain_cells.get_dof_indices_local_global_interface(dof_indices_local_global);
			for(const auto& u_sigma_n : u_sigma)
			{
				if(u_sigma_n->initial_vals != nullptr)
				{
					for(unsigned int component = 0; component < u_sigma_n->n_components; ++component)
					{
						const unsigned int global_component_index = global_component_indices_u_sigma.at(u_sigma_n) + component;
						const auto& shapefuns = components_to_shapefuns_interface[fe_system_id][global_component_index];
						for(const auto& shapefun : shapefuns)
						{
							unit_support_point_interface = interface_cell->get_fe().unit_support_point(shapefun);
							support_point = mapping_interface->transform_unit_to_real_cell(interface_cell, unit_support_point_interface);
							const double u_sigma_val = u_sigma_n->initial_vals->value(support_point, component);
							const unsigned int global_dof_index = dof_indices_local_global[shapefun];
							initial_fields[global_dof_index] = u_sigma_val;
						}
					}
				}
			}
		}
	}

	//initial values for scalar independent variables
	if(this_proc == n_procs - 1)
	{
		for(const auto& C_n : C)
			initial_fields[get_global_dof_index_C(C_n)] = C_n->initial_value;
	}
	initial_fields.compress(VectorOperation::insert);

	//apply constraints to solution vector in order to get initial state which is consistent with constraints
	if(constraints != nullptr)
		constraints->distribute(initial_fields);
	initial_fields.compress(VectorOperation::insert);

	//make sure that ghosts are imported here
	if(n_procs > 1)
	{
#ifdef DEAL_II_WITH_MPI
		const auto initial_fields_ptr = dynamic_cast<dealii::LinearAlgebra::distributed::Vector<double>*>(&initial_fields);
		if(initial_fields_ptr != nullptr)
			initial_fields_ptr->update_ghost_values();
#else
		Assert(n_procs == 1, ExcMessage("Internal error: deal.II not configured with MPI, but computation run with more than one processor"));
#endif //DEAL_II_WITH_MPI
	}
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::make_dirichlet_constraints_recursion(	const typename TriangulationSystem<spacedim>::DomainCell&	domain_cell,
																const unsigned int											face,
																const std::vector<unsigned int>&							shapefuns,
																const DirichletConstraint<spacedim>&						constraint,
																AffineConstraints<double>&									constraint_matrix,
																const AffineConstraints<double>&							constraints_ignore)
const
{

	static vector<unsigned int> dof_indices_local_global;
	dof_indices_local_global.reserve(fe_collection_domain.max_dofs_per_cell());

	unsigned int global_dof_index_C = 0;
	if(constraint.independent_scalar != nullptr)
		global_dof_index_C = get_global_dof_index_C(constraint.independent_scalar);

	if( domain_cell->has_children())
	{
		for(unsigned int subface = 0; subface < domain_cell->face(face)->n_children(); ++subface)
		{
			//we don't have to worry about face orientations here, because all we want is to get all the child cells behind the face, while
			//the order in which we are traversing them does not matter.
			const unsigned int child = GeometryInfo<spacedim>::child_cell_on_face(domain_cell->refinement_case(), face, subface/*,  domain_cell->face(face)->face_orientation(),
																																	domain_cell->face(face)->face_flip(),
																																	domain_cell->face(face)->face_rotation(),
																																	domain_cell->face(face)->refinement_case()*/);
			make_dirichlet_constraints_recursion(domain_cell->child(child), face, shapefuns, constraint, constraint_matrix, constraints_ignore);
		}
		// this may happen for anisotropic refinement
		if(domain_cell->face(face)->n_children() == 0)
		{
			const unsigned int child = GeometryInfo<spacedim>::child_cell_on_face(domain_cell->refinement_case(), face, 0/*, domain_cell->face(face)->face_orientation(),
																															 domain_cell->face(face)->face_flip(),
																															 domain_cell->face(face)->face_rotation(),
																															 domain_cell->face(face)->refinement_case()*/);
			make_dirichlet_constraints_recursion(domain_cell->child(child), face, shapefuns, constraint, constraint_matrix, constraints_ignore);
		}

	}
	else if(!domain_cell->is_artificial())
	{
		//unfortunately that conversion is necessary to get the numbering of dofs right if they have been renumbered
		DomainCellDoFIterator<spacedim> domain_cell_dof(domain_cell, dof_handler_system);
		dof_indices_local_global.resize(domain_cell_dof->get_fe().dofs_per_cell);
		domain_cell_dof.get_dof_indices(dof_indices_local_global);

		static unsigned int constrained_dof_index_global;
		static Point<spacedim> support_point, unit_support_point_domain;
		for(const auto& shapefun : shapefuns)
		{
			constrained_dof_index_global = dof_indices_local_global[shapefun];

			//dof's which are already constrained cannot be constrained a second time
			if( (!constraint_matrix.is_constrained(constrained_dof_index_global)) && (!constraints_ignore.is_constrained(constrained_dof_index_global)) )
			{
				constraint_matrix.add_line(constrained_dof_index_global);
				if(constraint.constraint_inhomogeneity != nullptr)
				{
					unit_support_point_domain = domain_cell_dof->get_fe().unit_support_point(shapefun);
					support_point = mapping_domain->transform_unit_to_real_cell(domain_cell_dof, unit_support_point_domain);
					const double inhomogeneity = constraint.constraint_inhomogeneity->value(support_point);
					constraint_matrix.set_inhomogeneity(constrained_dof_index_global, inhomogeneity);
				}
				if(constraint.independent_scalar != nullptr)
				{
					double coefficient_c;
					if(constraint.coefficient_c != nullptr)
					{
						unit_support_point_domain = domain_cell_dof->get_fe().unit_support_point(shapefun);
						support_point = mapping_domain->transform_unit_to_real_cell(domain_cell_dof, unit_support_point_domain);
						coefficient_c = constraint.coefficient_c->value(support_point);
					}
					else
						coefficient_c = 1.0;
					if(coefficient_c != 0.0)
						constraint_matrix.add_entry(constrained_dof_index_global, global_dof_index_C, coefficient_c);
				}
			}
		}
	}
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::make_dirichlet_constraints(	AffineConstraints<double>&							constraint_matrix,
														const vector<const DirichletConstraint<spacedim>*>&	dirichlet_constraints,
														const AffineConstraints<double>&					constraints_ignore)
const
{
	constraint_matrix.clear();
	if(n_procs > 1)
		constraint_matrix.reinit(get_locally_relevant_indices());

	//first generate data structures more suitable for assembly
	map<types::material_id, vector< pair< const DirichletConstraint<spacedim>*, const unsigned int > >> constraints_domain_internal_;
	for(const auto& constraint : dirichlet_constraints)
		for(const auto& material_id : constraint->domain_of_constraint)
			constraints_domain_internal_[material_id].push_back(make_pair(constraint, global_component_indices_u_omega.at(constraint->independent_field) + (constraint->component)));

	for(const auto& interface_cell_domain_cell : tria_system.interface_coarse_iterators())
	{
		const auto& constraints_domain_internal_mat = constraints_domain_internal_[interface_cell_domain_cell.interface_cell->material_id()];

		if( constraints_domain_internal_mat.size() == 0 )
			continue;

		//loop over the domain related constraints
		for(const auto& constraint_n: constraints_domain_internal_mat)
		{
			types::material_id material_id_domain;
			unsigned int face;
			if (constraint_n.first->side == InterfaceSide::minus)
			{
				material_id_domain = interface_cell_domain_cell.domain_cell_minus->material_id();
				face = interface_cell_domain_cell.face_minus;
			}
			else{
				Assert(	interface_cell_domain_cell.refinement_case != InterfaceRefinementCase::at_boundary,
						ExcMessage("You have tried to apply a boundary condition at the + side of a boundary, which is not possible!"))
				material_id_domain = interface_cell_domain_cell.domain_cell_plus->material_id();
				face = interface_cell_domain_cell.face_plus;
			}

			const unsigned int fe_system_id = material_id_to_fe_system_id_domain.at(material_id_domain);
			const auto& shapefuns = components_to_shapefuns_domain_facewise[fe_system_id][constraint_n.second][face];

			if (constraint_n.first->side == InterfaceSide::minus)
				make_dirichlet_constraints_recursion(interface_cell_domain_cell.domain_cell_minus, face, shapefuns, *(constraint_n.first), constraint_matrix, constraints_ignore);
			else
				make_dirichlet_constraints_recursion(interface_cell_domain_cell.domain_cell_plus, face, shapefuns, *(constraint_n.first), constraint_matrix, constraints_ignore);
		}
	}
}

template<unsigned int spacedim>
template<class SparsityPatternType>
void
AssemblyHelper<spacedim>::generate_sparsity_pattern_by_simulation(	SparsityPatternType&				dsp_K,
																	const AffineConstraints<double>&	constraints)
const
{


	//determine system size
	const unsigned int n_dofs_domain = dof_handler_system.n_dofs_domain();
	const unsigned int n_dofs_interface = dof_handler_system.n_dofs_interface();
	const unsigned int n_C = dof_handler_system.n_dofs_additional();
	const unsigned int n_system = system_size();
	const unsigned int index_first_nonlinear_scalar_functional = n_dofs_domain + n_dofs_interface + n_C;


	//this Vector stores which entries are coupled due to the handling of non-primitive total potential contributions
	//in the extra off-diagonal block of the stretched system
	IndexSet coupled_due_to_nonprimitivity(n_system);

	//////////////////////////////////////////////////
	// first add entries related to dof's on domain //
	//////////////////////////////////////////////////

	//stores the material_id of the cell visited previously
	types::material_id material_id_previous = numbers::invalid_material_id;

	//holds the internal index of a cell
	unsigned int internal_index = 0;

	//global dof indices coupling on the cell for a particular scalar functional
	vector<unsigned int> dof_indices_global;
	dof_indices_global.reserve(fe_collection_domain.max_dofs_per_cell() + global_indices_C.size());

	//mapping between local dof indices and global ones (dof_indices_local_global[i] is the global dof index corresponding to the local dof index i)
	vector<unsigned int> dof_indices_local_global;
	dof_indices_local_global.reserve(fe_collection_domain.max_dofs_per_cell());

	//global dof indices of C's
	vector<unsigned int> dof_indices_local_global_C;
	dof_handler_system.get_dof_indices(dof_indices_local_global_C);

	//global dof indices of C's for a particular scalar functional
	vector<unsigned int> dof_indices_global_C;
	dof_indices_global_C.reserve(global_indices_C.size());

	//the actual loop over the cells
	for(const auto& domain_cell : dof_handler_system.domain_active_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			//the internal material index needs only to be updated if the cell has another material_id than the one visited before
			if(domain_cell->material_id() != material_id_previous)
			{
				material_id_previous = domain_cell->material_id();
				internal_index = material_id_to_internal_index_domain.at(material_id_previous);
			}

			//get the mapping between local and global dof indices
			dof_indices_local_global.resize(domain_cell->get_fe().dofs_per_cell);
			domain_cell.get_dof_indices(dof_indices_local_global);

			//loop over scalar functionals on domain_cell
			for(unsigned int scalar_functional_n = 0; scalar_functional_n < scalar_functionals_domain[internal_index].size(); ++scalar_functional_n)
			{
				//this contains the local dof indices coupling on the domain cell for the scalar functional under consideration
				const auto& coupled_dof_indices_local = coupled_dof_indices_scalar_functionals_domain[internal_index][scalar_functional_n];

				//this contains the local C indices coupling on the domain cell for the scalar functional under consideration
				const auto& coupled_dof_indices_local_C = coupled_C_indices_scalar_functionals_domain[internal_index][scalar_functional_n];

				//map the local dof indices to the global ones
				Auxiliary::convert_local_indices_to_global_indices(coupled_dof_indices_local, dof_indices_global, dof_indices_local_global);
				Auxiliary::convert_local_indices_to_global_indices(coupled_dof_indices_local_C, dof_indices_global_C, dof_indices_local_global_C);

				//append the global dof indices of the C's
				if(coupled_dof_indices_local_C.size() > 0)
					for(const auto& dof_index_global_C : dof_indices_global_C)
						dof_indices_global.push_back(dof_index_global_C);

				//add the entries to the sparsity pattern
				if(dof_indices_global.size() > 0)
					constraints.add_entries_local_to_global(dof_indices_global, dsp_K, false);

				//if the scalar functional enters non-primitively in one or more of the contributions to the total potential
				//set the respective entries in coupled_due_to_nonprimitivity to true (the entries will be added to the sparsity pattern later)
				if(scalar_functionals_domain_nonprimitive_indices.find(scalar_functionals_domain[internal_index][scalar_functional_n]) != scalar_functionals_domain_nonprimitive_indices.end())
					for(const auto& dof_index :  dof_indices_global)
						coupled_due_to_nonprimitivity.add_index(dof_index);
			}
		}
	}

	///////////////////////////////////////////
	// now add entries related to interfaces //
	///////////////////////////////////////////

	//stores the material_id's of the cell visited previously
	tuple<types::material_id, types::material_id, types::material_id> material_ids_previous	=	make_tuple(	numbers::invalid_material_id,
																											numbers::invalid_material_id,
																											numbers::invalid_material_id);
	//global dof indices coupling on the interface cell for a particular scalar functional
	dof_indices_global.reserve(fe_collection_interface.max_dofs_per_cell());

	//global dof indices coupling on - side for a particular scalar functional
	vector<unsigned int> dof_indices_global_minus;
	dof_indices_global_minus.reserve(fe_collection_domain.max_dofs_per_cell());

	//global dof indices coupling on + side for a particular scalar functional
	vector<unsigned int> dof_indices_global_plus;
	dof_indices_global_plus.reserve(fe_collection_domain.max_dofs_per_cell());

	//combined global dof indices coupling on the interface cell and the adjacent domain cells for a particular scalar functional
	//(without any duplicate dof's -> this quantity is NOT just a concatenation of dof_indices_global, dof_indices_global_minus, dof_indices_global_plus)
	vector<unsigned int> dof_indices_global_combined;
	dof_indices_global_combined.reserve(	  fe_collection_interface.max_dofs_per_cell()
											+ 2*fe_collection_domain.max_dofs_per_cell()
											+ global_indices_C.size());

	//defined such that dof_indices_global_combined[dof_indices_interface_dof_indices_combined[i]] = dof_indices_global[i]
	vector<unsigned int> dof_indices_interface_dof_indices_combined;
	dof_indices_interface_dof_indices_combined.reserve(fe_collection_interface.max_dofs_per_cell());

	//defined such that dof_indices_global_combined[dof_indices_minus_dof_indices_combined[i]] = dof_indices_global_minus[i]
	vector<unsigned int> dof_indices_minus_dof_indices_combined;
	dof_indices_minus_dof_indices_combined.reserve(fe_collection_domain.max_dofs_per_cell());

	//defined such that dof_indices_global_combined[dof_indices_plus_dof_indices_combined[i]] = dof_indices_global_plus[i]
	vector<unsigned int> dof_indices_plus_dof_indices_combined;
	dof_indices_plus_dof_indices_combined.reserve(fe_collection_domain.max_dofs_per_cell());

	//defined such that dof_indices_global_combined[dof_indices_C_dof_indices_combined[i]] = dof_indices_global_C[i]
	vector<unsigned int> dof_indices_C_dof_indices_combined;
	dof_indices_C_dof_indices_combined.reserve(global_indices_C.size());

	//mapping between local dof indices and global ones on interface
	dof_indices_local_global.reserve(fe_collection_interface.max_dofs_per_cell());

	//mapping between local dof indices and global ones on - side
	vector<unsigned int> dof_indices_local_global_minus;
	dof_indices_local_global_minus.reserve(fe_collection_domain.max_dofs_per_cell());

	//mapping between local dof indices and global ones on + side
	vector<unsigned int> dof_indices_local_global_plus;
	dof_indices_local_global_plus.reserve(fe_collection_domain.max_dofs_per_cell());

	//the actual loop over the cells
	for(const auto& interface_cell_domain_cell : dof_handler_system.interface_active_iterators())
	{
		if(interface_cell_domain_cell.interface_cell->is_locally_owned())
		{
			//find out the material id's characterizing the interface
			const auto material_ids = interface_cell_domain_cell.get_material_ids();

			//the internal material index needs only to be updated if the cell is associated with different material_id's than the one visited before
			if(material_ids != material_ids_previous)
			{
				internal_index = material_ids_to_internal_index_interface.at(material_ids);
				material_ids_previous = material_ids;
			}

			//get the mappings between local and global dof indices
			interface_cell_domain_cell.get_dof_indices_local_global_interface(dof_indices_local_global, dof_indices_local_global_minus, dof_indices_local_global_plus);

			//loop over scalar functionals on cell
			for(unsigned int scalar_functional_n = 0; scalar_functional_n < scalar_functionals_interface[internal_index].size(); ++scalar_functional_n)
			{
				//this contains the local dof indices coupling on the interface cell for the scalar functional under consideration
				const vector<unsigned int>& coupled_dof_indices_local=coupled_dof_indices_scalar_functionals_interface[internal_index][scalar_functional_n];

				//this contains the local dof indices coupling on the - side for the scalar functional under consideration
				const vector<unsigned int>& coupled_dof_indices_local_minus = coupled_dof_indices_scalar_functionals_interface_minus[internal_index][scalar_functional_n];

				//this contains the local dof indices coupling on the + side for the scalar functional under consideration
				const vector<unsigned int>& coupled_dof_indices_local_plus = coupled_dof_indices_scalar_functionals_interface_plus[internal_index][scalar_functional_n];

				//this contains the local C indices coupling on the domain cell for the scalar functional under consideration
				const vector<unsigned int>& dof_indices_local_C = coupled_C_indices_scalar_functionals_interface[internal_index][scalar_functional_n];

				//map the local dof indices to the global ones
				Auxiliary::convert_local_indices_to_global_indices(coupled_dof_indices_local, dof_indices_global, dof_indices_local_global);
				Auxiliary::convert_local_indices_to_global_indices(coupled_dof_indices_local_minus, dof_indices_global_minus, dof_indices_local_global_minus);
				Auxiliary::convert_local_indices_to_global_indices(coupled_dof_indices_local_plus, dof_indices_global_plus, dof_indices_local_global_plus);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_C, dof_indices_global_C, dof_indices_local_global_C);

				//assemble combined dof indices
				Auxiliary::combine_dof_indices(	dof_indices_global,
												dof_indices_global_minus,
												dof_indices_global_plus,
												dof_indices_global_C,
												dof_indices_interface_dof_indices_combined,
												dof_indices_minus_dof_indices_combined,
												dof_indices_plus_dof_indices_combined,
												dof_indices_C_dof_indices_combined,
												dof_indices_global_combined);

				//add the entries to the sparsity pattern
				if(dof_indices_global_combined.size() > 0)
					constraints.add_entries_local_to_global(dof_indices_global_combined, dsp_K, false);

				//if the scalar functional enters non-primitively in one or more of the contributions to the total potential
				//set the respective entries in coupled_due_to_nonprimitivity to true (the entries will be added to the sparsity pattern later)
				if(scalar_functionals_interface_nonprimitive_indices.find(scalar_functionals_interface[internal_index][scalar_functional_n]) != scalar_functionals_interface_nonprimitive_indices.end())
					for(const auto& dof_index :  dof_indices_global_combined)
						coupled_due_to_nonprimitivity.add_index(dof_index);
			}
		}
	}


	////////////////////////////////////////////////////////////////////////////////
	// finally add entries related to non-primitive total potential contributions //
	////////////////////////////////////////////////////////////////////////////////

	//the C's enter in any case non-primitively into the total potential
	vector<unsigned int> global_dof_indices_C(n_C);
	dof_handler_system.get_dof_indices(global_dof_indices_C);
	for(const auto& dof_index : global_dof_indices_C)
		coupled_due_to_nonprimitivity.add_index(dof_index);

	if(n_C + n_scalar_functionals_nonprimitive > 0)
	{
		//vector holding global dof indices which are coupled due to handling of scalar functionals entering the total potential non-primitively
		vector<unsigned int> nonprimitively_coupled_dof_indices_global;
		nonprimitively_coupled_dof_indices_global.reserve(coupled_due_to_nonprimitivity.n_elements());
		for(const auto& index : coupled_due_to_nonprimitivity)
			nonprimitively_coupled_dof_indices_global.push_back(index);

		//these are the rows (and columns) in the system matrix which are related to the scalar functionals entering the total potential non-primitively
		vector<unsigned int> rows_to_add_to;
		for(unsigned int row = 0; row < n_scalar_functionals_nonprimitive + n_C; ++row)
		{
			rows_to_add_to.push_back(row + index_first_nonlinear_scalar_functional);
			//use this chance to add the diagonal entries which are needed for the handling of scalar functionals entering the total potential non-primitively
			dsp_K.add(row + index_first_nonlinear_scalar_functional, row + index_first_nonlinear_scalar_functional);
		}
		//add the entries
		constraints.add_entries_local_to_global(rows_to_add_to, nonprimitively_coupled_dof_indices_global, dsp_K, false);
		constraints.add_entries_local_to_global(nonprimitively_coupled_dof_indices_global, rows_to_add_to, dsp_K, false);
	}
}

template<unsigned int spacedim>
template<class SolutionVectorType, class RHSVectorType, class MatrixType>
bool
AssemblyHelper<spacedim>::assemble_system(	const SolutionVectorType&				solution,
											const vector<const SolutionVectorType*>	solution_ref_sets,
											const AffineConstraints<double>&		constraints,
											double&									potential_value,
											RHSVectorType&							f,
											MatrixType&								K,
											const tuple<bool,bool,bool>				requested_quantities)
const
{

	Assert( solution.size() == system_size(), ExcMessage("The solution vector has not the correct size!"));
	for(const auto& solution_ref : solution_ref_sets)
	{
		(void) solution_ref;	//silence unused parameter warning of compiler
		Assert( solution_ref->size() == system_size(), ExcMessage("One of the reference solution vectors has not the correct size!"));
	}
	Assert( f.size() == system_size(), ExcMessage("The vector f has not the correct size!"));
	Assert( (K.m() == system_size()) && (K.n() == system_size()), ExcMessage("The matrix K  has not the correct size!"));

	//the return value, false if everything is ok
	bool error = false;

	//requested quantities
	bool compute_potential = get<0>(requested_quantities);
	bool assemble_rhs = get<1>(requested_quantities);
	bool assemble_matrix = get<2>(requested_quantities);

	Assert( compute_potential || assemble_rhs || assemble_matrix, ExcMessage("It does not make sense to call assemble_system() with not requesting anything to be calculated!") );

	Assert( !((assemble_matrix == true) && (assemble_rhs == false)),
		ExcMessage("If the matrix is assembled, the rhs must be assembled as well!") );

	//initialize requested quantities to zero
	if(compute_potential)
		potential_value = 0.0;
	if(assemble_rhs)
	{
		f = 0.0;
		f.compress(VectorOperation::insert);
	}
	if(assemble_matrix)
	{
		K = 0.0;
		//this appears not to be necessary after K=0.0, but to be sure
		K.compress(VectorOperation::insert);
	}

	//vector with values of scalar functionals entering at least one total potential contribution primitively
	//(sorted according to scalar_functionals_domain_primitive_indices and scalar_functionals_interface_primitive_indices)
	Vector<double> primitive_scalar_functionals_values(n_scalar_functionals_primitive);

/*************************************************************************************************
 * first calculate values of scalar functionals entering total potential non-primitively, if any *
 *************************************************************************************************/

	//vector with values of scalar functionals entering the total potential non-primitively
	//(sorted according to scalar_functionals_domain_nonprimitive_indices and scalar_functionals_interface_nonprimitive_indices)
	Vector<double> nonprimitive_scalar_functionals_values(n_scalar_functionals_nonprimitive);
	if(n_scalar_functionals_nonprimitive > 0)
	{
		if(get_nonprimitive_scalar_functional_values(solution, solution_ref_sets, nonprimitive_scalar_functionals_values))
			error = true;
	}

/******************************************************************************************
 * calculate values and derivatives of total potential w.r.t. scalar functionals and C's, *
 * add non-primitive total potential contributions already here to potential value 		  *
 ******************************************************************************************/

	//value of a total potential contribution
	double Pi;

	//first derivative of a total potential contribution w.r.t. scalar functionals and C's
	Vector<double> Pi_1;

	//second derivative of a total potential contribution w.r.t. scalar functionals and C's
	FullMatrix<double> Pi_2;

	//values of scalar functionals and C's
	Vector<double> H_omega_H_sigma_C;

	//reference values of C's
	vector<Vector<double>> C_ref(solution_ref_sets.size());

	//first derivative of total potential w.r.t. scalar functionals and C's entering at least one total potential contribution primitively
	Vector<double> dPi_dH_omega_H_sigma_C_nonprimitive(n_scalar_functionals_nonprimitive + C.size());

	//first derivative of total potential part w.r.t. scalar functionals entering all total potential contributions only primitively
	//(we don't need the second derivatives because they are zero)
	Vector<double> dPi_dH_omega_H_sigma_primitive(n_scalar_functionals_primitive);
	
	//second derivative of total potential w.r.t. scalar functionals and C's entering at least one total potential contribution primitively
	FullMatrix<double> d2Pi_d2H_omega_H_sigma_C_nonprimitive(n_scalar_functionals_nonprimitive + C.size());

	//mapping between indexing of scalar functionals within total potential contribution and
	//global indexing of scalar functionals entering the total potential non-primitively)
	vector<unsigned int> nonprimitive_indices_local_global;
	nonprimitive_indices_local_global.reserve(n_scalar_functionals_nonprimitive + C.size());

	//just used to distribute the scalar functional distributions from Pi_1 and Pi_2 to dPi_dH_omega_H_sigma_C_nonprimitive and d2Pi_d2H_omega_H_sigma_C_nonprimitive
	AffineConstraints<double> dummy_constraints;

	//loop over the total potential contributions
	for(const auto& total_potential_contribution : total_potential.total_potential_contributions)
	{
		//total_potential_contribution is non-primitive (i.e., it is not just equal to a scalar functional)
		if(!total_potential_contribution->is_primitive)
		{
			const unsigned int n_H_omega = total_potential_contribution->H_omega.size();
			const unsigned int n_H_sigma = total_potential_contribution->H_sigma.size();
			const unsigned int n_C = total_potential_contribution->C.size();
			const unsigned int n_H_omega_H_sigma_C = n_H_omega + n_H_sigma + n_C;

			//initializations
			H_omega_H_sigma_C.reinit(n_H_omega_H_sigma_C);
			for(auto& C_ref_n : C_ref)
				C_ref_n.reinit(n_C);
			nonprimitive_indices_local_global.resize(n_H_omega_H_sigma_C);
			if(compute_potential)
				Pi = 0.0;
			if(assemble_rhs)
				Pi_1.reinit(H_omega_H_sigma_C.size());
			if(assemble_matrix)
				Pi_2.reinit(H_omega_H_sigma_C.size(), H_omega_H_sigma_C.size());


			//assemble nonprimitive_indices_local_global and H_omega_H_sigma_C
			for(unsigned int H_omega_n = 0; H_omega_n < n_H_omega; ++H_omega_n)
			{
				nonprimitive_indices_local_global[H_omega_n] = scalar_functionals_domain_nonprimitive_indices.at(total_potential_contribution->H_omega[H_omega_n]);
				H_omega_H_sigma_C[H_omega_n] = nonprimitive_scalar_functionals_values[ nonprimitive_indices_local_global[H_omega_n] ];
			}
			for(unsigned int H_sigma_n = 0; H_sigma_n < n_H_sigma; H_sigma_n++)
			{
				const unsigned int H_sigma_n_index = n_H_omega + H_sigma_n;
				nonprimitive_indices_local_global[H_sigma_n_index] = scalar_functionals_interface_nonprimitive_indices.at(total_potential_contribution->H_sigma[H_sigma_n]);
				H_omega_H_sigma_C[H_sigma_n_index] = nonprimitive_scalar_functionals_values[ nonprimitive_indices_local_global[H_sigma_n_index] ];
			}
			for(unsigned int C_n = 0; C_n < n_C; ++C_n)
			{
				const unsigned int C_n_index = n_H_omega + n_H_sigma + C_n;
				const unsigned int C_n_global_index = get_global_dof_index_C(total_potential_contribution->C[C_n]);
				nonprimitive_indices_local_global[C_n_index] = n_scalar_functionals_nonprimitive + global_indices_C.at(total_potential_contribution->C[C_n]);
				H_omega_H_sigma_C[C_n_index] = solution[C_n_global_index];
				for(unsigned int ref_set=0; ref_set<solution_ref_sets.size(); ++ref_set)
					C_ref[ref_set][C_n] = (*solution_ref_sets[ref_set])[C_n_global_index];
			}

			//call function get_potential_contribution to get Pi, Pi_1, Pi_2
			if(total_potential_contribution->get_potential_contribution(H_omega_H_sigma_C, C_ref, Pi, Pi_1, Pi_2, requested_quantities))
				error = true;

			//add results to potential_value, dPi_dH_omega_H_sigma_C_nonprimitive, d2Pi_d2H_omega_H_sigma_C_nonprimitive
			if(compute_potential)
				potential_value += Pi;
			if(assemble_rhs)
				dummy_constraints.distribute_local_to_global(Pi_1, nonprimitive_indices_local_global, dPi_dH_omega_H_sigma_C_nonprimitive);
			if(assemble_matrix)
				dummy_constraints.distribute_local_to_global(Pi_2, nonprimitive_indices_local_global, d2Pi_d2H_omega_H_sigma_C_nonprimitive);


		}
		//total_potential_contribution is primitive (i.e., it is a single scalar functional)
		else
		{
			if(assemble_rhs)
			{
				int nonprimitive_index, primitive_index;
				if(total_potential_contribution->H_omega.size() == 1)
				{
					const auto index_pair = get_scalar_functional_indices(total_potential_contribution->H_omega[0]);
					nonprimitive_index = index_pair.first;
					primitive_index = index_pair.second;
				}
				else
				{
					const auto index_pair = get_scalar_functional_indices(total_potential_contribution->H_sigma[0]);
					nonprimitive_index = index_pair.first;
					primitive_index = index_pair.second;
				}
				//the total_potential_contribution is primitive, but the scalar functional forming the total potential contribution
				//enters another total potential contribution non-primitively
				if(nonprimitive_index > -1)
					dPi_dH_omega_H_sigma_C_nonprimitive[nonprimitive_index] += 1.0;
				else
					dPi_dH_omega_H_sigma_primitive[primitive_index] += 1.0;
			}
		}
	}

/************************************************************************
 * calculate ldr decomposition of d2Pi_d2H_omega_H_sigma_C_nonprimitive *
 ************************************************************************/

	//todo: strictly this is required only on a single processor ... but also it doesn't hurt to do it everywhere as it's a small computation
	Vector<double> D_;
	vector<Vector<double>> L, R;
	if(assemble_matrix && (d2Pi_d2H_omega_H_sigma_C_nonprimitive.size()[0] > 0) )
		Auxiliary::compute_ldr(d2Pi_d2H_omega_H_sigma_C_nonprimitive, D_, L, R);

/*******************************************
 * do the actual assembly of the RHS and K *
 *******************************************/

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// diagonal entries related to scalar functionals and independent scalars entering the total potential non-primitively //
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//this vector holds the global indices of the rows (and columns) in the system matrix
	//which are related to the scalar functionals and independent scalars entering the total potential non-primitively,
	//the vector is also used later
	vector<unsigned int> nonlinear_scalar_functional_indices;
	const unsigned int index_first_nonlinear_scalar_functional = dof_handler_system.n_dofs();
	nonlinear_scalar_functional_indices.resize(D_.size());
	if(assemble_matrix)
	{
		for(unsigned int nonlinear_scalar_functional_index = 0; nonlinear_scalar_functional_index < D_.size(); ++nonlinear_scalar_functional_index)
		{
			nonlinear_scalar_functional_indices[nonlinear_scalar_functional_index] = nonlinear_scalar_functional_index + index_first_nonlinear_scalar_functional;
			//set the diagonal entry (which is actually either +1 or -1)
			K.set(	nonlinear_scalar_functional_indices[nonlinear_scalar_functional_index],
					nonlinear_scalar_functional_indices[nonlinear_scalar_functional_index],
					-D_[nonlinear_scalar_functional_index]);
		}
		K.compress(VectorOperation::insert);
	}

	/////////////////////////////////////
	// contributions related to domain //
	/////////////////////////////////////

	//stores the material_id of the cell visited previously
	types::material_id material_id_previous=numbers::invalid_material_id;

	//holds the internal index of a cell
	unsigned int internal_index = 0;

	//global dof indices coupling on the cell for a particular scalar functional
	vector<unsigned int> dof_indices_global;
	dof_indices_global.reserve(fe_collection_domain.max_dofs_per_cell() + global_indices_C.size());

	//mapping between local dof indices and global ones (dof_indices_local_global[i] is the global dof index corresponding to the local dof index i)
	vector<unsigned int> dof_indices_local_global;
	dof_indices_local_global.reserve(fe_collection_domain.max_dofs_per_cell());

	//global dof indices of C's
	vector<unsigned int> dof_indices_local_global_C;
	dof_handler_system.get_dof_indices(dof_indices_local_global_C);

	//global dof indices of C's for a particular scalar functional
	vector<unsigned int> dof_indices_global_C;
	dof_indices_global_C.reserve(global_indices_C.size());

	//solution vector restricted to dof indices coupling on a cell for a particular scalar functional
	Vector<double> solution_local;
	solution_local.reinit(fe_collection_domain.max_dofs_per_cell());

	//solution vector restricted to C indices coupling on a cell for a particular scalar functional
	Vector<double> solution_local_C;
	solution_local_C.reinit(global_indices_C.size());

	//reference solution vector restricted to dof indices coupling on a cell for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_n : solution_ref_sets_local)
		solution_ref_sets_local_n.reinit(fe_collection_domain.max_dofs_per_cell());

	//reference solution vector restricted to C indices coupling on a cell for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local_C(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_C_n : solution_ref_sets_local_C)
		solution_ref_sets_local_C_n.reinit(global_indices_C.size());

	//vectors for dependent variables
	Vector<double> e_omega(total_potential.max_dependent_vars);
	vector<Vector<double>> e_omega_ref_sets(solution_ref_sets.size());

	//matrix for derivatives of dependent variables w.r.t. local dof's
	FullMatrix<double> de_omega_dsol_T;

	//value of a scalar functional
	double h_omega;

	//first derivative of a scalar functional
	Vector<double> h_omega_1;

	//second derivative of a scalar functional
	FullMatrix<double> h_omega_2;

	//h_omega_2*de_omega_dsol_T
	FullMatrix<double> h_omega_2_de_omega_dsol_T;

	//rhs contribution associated with a particular scalar functional on a cell
	Vector<double> f_cell;

	//this dummy vector will be required during assembly to take the effect of constraint inhomogeneities into account
	Vector<double> rhs_cell_dummy;

	//system matrix contribution associated with a particular scalar functional on a cell
	FullMatrix<double> K_cell;

	//contribution of a particular scalar functional on a cell to the lower left block of the system matrix
	//(these are the rows related to the handling of the scalar functionals entering the total potential non-primitively)
	FullMatrix<double> V_T_L_T_cell;

	//contribution of a particular scalar functional on a cell to the upper right block of the system matrix
	//(these are the columns related to the handling of the scalar functionals entering the total potential non-primitively)
	FullMatrix<double> L_U_cell;

	//the actual loop over the cells
	for(const auto& domain_cell : dof_handler_system.domain_active_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			//the internal material index needs only to be updated if the cell has another material_id than the one visited before
			if(domain_cell->material_id() != material_id_previous)
			{
				material_id_previous = domain_cell->material_id();
				internal_index = material_id_to_internal_index_domain.at(material_id_previous);
			}

			//initialize relevant FEValues objects with cell
			initialize_fe_values_domain(domain_cell, internal_index, false);

			//get the mapping between local and global dof indices
			dof_indices_local_global.resize(domain_cell->get_fe().dofs_per_cell);
			domain_cell.get_dof_indices(dof_indices_local_global);

			//loop over scalar functionals on cell
			for(unsigned int scalar_functional_n = 0; scalar_functional_n<scalar_functionals_domain[internal_index].size(); ++scalar_functional_n)
			{

				//this contains the local dof indices coupling on the domain cell for the scalar functional under consideration
				const auto& dof_indices_local = coupled_dof_indices_scalar_functionals_domain[internal_index][scalar_functional_n];

				//this contains the local C indices coupling on the domain cell for the scalar functional under consideration
				const auto& dof_indices_local_C = coupled_C_indices_scalar_functionals_domain[internal_index][scalar_functional_n];

				//nothing to do if there are no dof's present
				if((dof_indices_local.size() + dof_indices_local_C.size()) == 0)
					continue;

				//map the local dof indices to the global ones
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local, dof_indices_global, dof_indices_local_global);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_C, dof_indices_global_C, dof_indices_local_global_C);

				//restrict solution vectors to scalar functional on cell
				solution.extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_local.begin());
				solution.extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_local_C.begin());
				for(unsigned int ref_set=0; ref_set<solution_ref_sets.size(); ++ref_set)
				{
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_ref_sets_local[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_ref_sets_local_C[ref_set].begin());
				}

				//append global dof indices of C to global dof index vector
				for(const auto& dof_index_global_C : dof_indices_global_C)
					dof_indices_global.push_back(dof_index_global_C);

				//initialize F and K to correct size
				//attention: F and K are not only assembled cell-wise, but also scalar functional wise,
				//			 therefore, the initialization cannot be done outside this loop
				if(assemble_rhs)
					f_cell.reinit(dof_indices_global.size());
				if(assemble_matrix)
					K_cell.reinit(dof_indices_global.size(), dof_indices_global.size());

				//get a pointer to the scalar functional under consideration and its nonprimitive_index and primitive_index
				const auto scalar_functional = scalar_functionals_domain[internal_index][scalar_functional_n];
				const auto index_pair = get_scalar_functional_indices(scalar_functional);
				const int nonprimitive_index = index_pair.first;
				const int primitive_index = index_pair.second;

				//initialize dependent variable vectors appropriately
				e_omega.reinit(scalar_functional->e_omega.size());
				for(auto& e_omega_ref_n : e_omega_ref_sets)
					e_omega_ref_n.reinit(e_omega.size());

				//initialize derivatives of scalar functional w.r.t. dependent variables appropriately
				if(assemble_rhs)
					h_omega_1.reinit(e_omega.size());
				if(assemble_matrix)
				{
					h_omega_2.reinit(e_omega.size(), e_omega.size());
					h_omega_2_de_omega_dsol_T.reinit(e_omega.size(), f_cell.size());
				}

				//get JxW values and quadrature points
				const auto& JxW = fe_values_domain[internal_index][scalar_functional_n]->get_JxW_values();
				const auto& q_points = fe_values_domain[internal_index][scalar_functional_n]->get_quadrature_points();

				//loop over integration points
				for(unsigned int q_point = 0; q_point < scalar_functional->quadrature.size(); ++q_point)
				{

					//compute dependent variables
					compute_e_omega(internal_index, scalar_functional_n, q_point, solution_local, solution_local_C, e_omega, de_omega_dsol_T, assemble_rhs);

					//compute reference values of dependent variables
					for(unsigned int ref_set=0; ref_set<solution_ref_sets.size(); ++ref_set)
						compute_e_omega(internal_index, scalar_functional_n, q_point, solution_ref_sets_local[ref_set], solution_ref_sets_local_C[ref_set], e_omega_ref_sets[ref_set], de_omega_dsol_T, false);

					//hidden variables
					Assert(domain_cell->user_pointer() != nullptr, ExcMessage("Internal error - it seems that the hidden variables are not propery allocated!"));
					Vector<double>& hidden_vars = (*static_cast<vector<vector<Vector<double>>>*>(domain_cell->user_pointer()))[scalar_functional_n][q_point];

					//the location of the quadrature point in real space
					const auto& x = q_points[q_point];

					//evaluate the integrand of the scalar functional and its derivatives w.r.t. the independent variables
					if(scalar_functional->get_h_omega(e_omega, e_omega_ref_sets, hidden_vars, x, h_omega, h_omega_1, h_omega_2, requested_quantities))
						error = true;

					//add contribution to the value of the scalar functional if it enters the total potential primitively
					//(values of scalar functionals entering the total potential non-primitively have been computed earlier
					// and need not be taken into account again)
					if(compute_potential && (primitive_index > -1))
						primitive_scalar_functionals_values[primitive_index] += JxW[q_point]*h_omega;

					//add contributions to f_cell and K_cell
					if(assemble_rhs)
					{
						h_omega_1 *= JxW[q_point];
						de_omega_dsol_T.vmult(f_cell, h_omega_1, true);
					}
					if(assemble_matrix)
					{
						h_omega_2 *= JxW[q_point];
						h_omega_2.mTmult(h_omega_2_de_omega_dsol_T, de_omega_dsol_T);
						de_omega_dsol_T.mmult(K_cell, h_omega_2_de_omega_dsol_T, true);
					}
				}
				scalar_functional->modify_K_cell_f_cell(domain_cell,
														K_cell,
														f_cell,
														this->coupled_dof_indices_scalar_functionals_domain[internal_index][scalar_functional_n],
														this->coupled_C_indices_scalar_functionals_domain[internal_index][scalar_functional_n]);

				//additional work to do for scalar functionals entering the total potential non-primitively
				//(lower left and upper right part of stretched system matrix)
				if((nonprimitive_index > -1) && assemble_matrix)
				{
					//compute  V_T_L_T_cell and L_U_cell
					V_T_L_T_cell.reinit(n_scalar_functionals_nonprimitive + C.size(), f_cell.size());
					L_U_cell.reinit(f_cell.size(), n_scalar_functionals_nonprimitive + C.size());
					for(unsigned int k = 0; k < V_T_L_T_cell.size()[0]; ++k)
						for(unsigned int l = 0; l < V_T_L_T_cell.size()[1]; ++l)
							V_T_L_T_cell(k,l) = R[nonprimitive_index](k) * f_cell(l);
					for(unsigned int k = 0; k < L_U_cell.size()[0]; ++k)
						for(unsigned int l = 0; l < L_U_cell.size()[1]; ++l)
							L_U_cell(k,l) = L[nonprimitive_index](l) * f_cell(k);

					//distribute result to system matrix (or respective vectors)
					constraints.distribute_local_to_global(V_T_L_T_cell, nonlinear_scalar_functional_indices, dof_indices_global, K);
					constraints.distribute_local_to_global(L_U_cell, dof_indices_global, nonlinear_scalar_functional_indices, K);
					//take the effect of constraint inhomogeneities into account (though the rhs related to V_T_L_T_cell and L_U_cell, respectively,
					//is zero, there may still be an effect of inhomogeneous constraints on the system rhs)
					rhs_cell_dummy.reinit(nonlinear_scalar_functional_indices.size());
					constraints.distribute_local_to_global(rhs_cell_dummy, nonlinear_scalar_functional_indices, dof_indices_global, f, V_T_L_T_cell);
					rhs_cell_dummy.reinit(dof_indices_global.size());
					constraints.distribute_local_to_global(rhs_cell_dummy, dof_indices_global, nonlinear_scalar_functional_indices, f, L_U_cell);
				}

				//scale rhs_cell and K_cell according to derivative of total potential contribution w.r.t. scalar functional
				if(assemble_rhs)
				{
					if(nonprimitive_index>-1)
						f_cell *= dPi_dH_omega_H_sigma_C_nonprimitive[nonprimitive_index];
					else
						f_cell *= dPi_dH_omega_H_sigma_primitive[primitive_index];
				}
				if(assemble_matrix)
				{
					if(nonprimitive_index>-1)
						K_cell *= dPi_dH_omega_H_sigma_C_nonprimitive[nonprimitive_index];
					else
						K_cell *= dPi_dH_omega_H_sigma_primitive[primitive_index];
				}

				//rhs is negative of gradient
				//(up to this point the contribution of the scalar functional on the cell to the gradient has been computed)
				if(assemble_rhs)
					f_cell *= -1.0;

				//distribute local contributions to global f and K
				if(assemble_rhs && assemble_matrix)
					constraints.distribute_local_to_global(K_cell, f_cell, dof_indices_global, K, f, false);
				else if(assemble_rhs)
					constraints.distribute_local_to_global(f_cell, dof_indices_global, f);
			}
		}
	}

	////////////////////////////////////////
	// contributions related to interface //
	////////////////////////////////////////

	//stores the material_ids of the cell visited previously
	tuple<types::material_id, types::material_id, types::material_id> material_ids_previous = make_tuple(	numbers::invalid_material_id,
																											numbers::invalid_material_id,
																											numbers::invalid_material_id);

	//global dof indices coupling on the interface cell for a particular scalar functional
	dof_indices_global.reserve(fe_collection_interface.max_dofs_per_cell());

	//global dof indices coupling on - side for a particular scalar functional
	vector<unsigned int> dof_indices_global_minus;
	dof_indices_global_minus.reserve(fe_collection_domain.max_dofs_per_cell());

	//global dof indices coupling on + side for a particular scalar functional
	vector<unsigned int> dof_indices_global_plus;
	dof_indices_global_plus.reserve(fe_collection_domain.max_dofs_per_cell());

	//mapping between local dof indices and global ones on interface
	dof_indices_local_global.reserve(fe_collection_interface.max_dofs_per_cell());

	//mapping between local dof indices and global ones on - side
	vector<unsigned int> dof_indices_local_global_minus;
	dof_indices_local_global_minus.reserve(fe_collection_domain.max_dofs_per_cell());

	//mapping between local dof indices and global ones on + side
	vector<unsigned int> dof_indices_local_global_plus;
	dof_indices_local_global_plus.reserve(fe_collection_domain.max_dofs_per_cell());

	//combined global dof indices coupling on the interface cell and the adjacent domain cells for a particular scalar functional
	//(without any duplicate dof's -> this quantity is NOT just a concatenation of dof_indices_global, dof_indices_global_minus, dof_indices_global_plus)
	vector<unsigned int> dof_indices_global_combined;
	dof_indices_global_combined.reserve(	  fe_collection_interface.max_dofs_per_cell()
											+ 2*fe_collection_domain.max_dofs_per_cell()
											+ global_indices_C.size());

	//defined such that dof_indices_global_combined[dof_indices_interface_dof_indices_combined[i]] = dof_indices_global[i]
	vector<unsigned int> dof_indices_interface_dof_indices_combined;
	dof_indices_interface_dof_indices_combined.reserve(fe_collection_interface.max_dofs_per_cell());

	//defined such that dof_indices_global_combined[dof_indices_minus_dof_indices_combined[i]] = dof_indices_global_minus[i]
	vector<unsigned int> dof_indices_minus_dof_indices_combined;
	dof_indices_minus_dof_indices_combined.reserve(fe_collection_domain.max_dofs_per_cell());

	//defined such that dof_indices_global_combined[dof_indices_plus_dof_indices_combined[i]] = dof_indices_global_plus[i]
	vector<unsigned int> dof_indices_plus_dof_indices_combined;
	dof_indices_plus_dof_indices_combined.reserve(fe_collection_domain.max_dofs_per_cell());

	//defined such that dof_indices_global_combined[dof_indices_C_dof_indices_combined[i]] = dof_indices_global_C[i]
	vector<unsigned int> dof_indices_C_dof_indices_combined;
	dof_indices_C_dof_indices_combined.reserve(global_indices_C.size());

	//solution vector restricted to dof indices coupling on an interface cell for a particular scalar functional
	solution_local.reinit(fe_collection_interface.max_dofs_per_cell());

	//reference solution vector restricted to dof indices coupling on an interface cell for a particular scalar functional
	for(auto& solution_ref_sets_local_n : solution_ref_sets_local)
		solution_ref_sets_local_n.reinit(fe_collection_interface.max_dofs_per_cell());

	//solution vector restricted to dof indices coupling on a domain cell located on the - side of the interface for a particular scalar functional
	Vector<double> solution_local_minus;
	solution_local_minus.reinit(fe_collection_domain.max_dofs_per_cell());

	//reference solution vector restricted to dof indices coupling on a domain cell located on the - side of the interface for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local_minus(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_minus_n : solution_ref_sets_local_minus)
		solution_ref_sets_local_minus_n.reinit(fe_collection_domain.max_dofs_per_cell());

	//solution vector restricted to dof indices coupling on a domain cell located on the + side of the interface for a particular scalar functional
	Vector<double> solution_local_plus;
	solution_local_plus.reinit(fe_collection_domain.max_dofs_per_cell());

	//reference solution vector restricted to dof indices coupling on a domain cell located on the + side of the interface for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local_plus(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_plus_n : solution_ref_sets_local_plus)
		solution_ref_sets_local_plus_n.reinit(fe_collection_domain.max_dofs_per_cell());

	//vectors for dependent variables
	Vector<double> e_sigma(total_potential.max_dependent_vars);
	vector<Vector<double>> e_sigma_ref_sets(solution_ref_sets.size());

	//matrix for derivatives of dependent variables w.r.t. local dof's
	FullMatrix<double> de_sigma_dsol_T;

	//value of a scalar functional
	double h_sigma;

	//first derivative of a scalar functional
	Vector<double> h_sigma_1;

	//second derivative of a scalar functional
	FullMatrix<double> h_sigma_2;

	//h_sigma_2*de_sigma_dsol_T
	FullMatrix<double> h_sigma_2_de_sigma_dsol_T;

	//the actual loop over the cells
	for(const auto& interface_cell_domain_cells : dof_handler_system.interface_active_iterators())
	{

		if(interface_cell_domain_cells.interface_cell->is_locally_owned())
		{
			//find out the material id's characterizing the interface
			const auto material_ids = interface_cell_domain_cells.get_material_ids();
			//the internal material index needs only to be updated if the cell is associated with different material_id's than the one visited before
			if(material_ids != material_ids_previous)
			{
				internal_index = material_ids_to_internal_index_interface.at(material_ids);
				material_ids_previous = material_ids;
			}

			//initialize FEValues objects with cell
			initialize_fe_values_interface(interface_cell_domain_cells, internal_index, false);

			//get the mapping between local and global dof indices
			interface_cell_domain_cells.get_dof_indices_local_global_interface(dof_indices_local_global, dof_indices_local_global_minus, dof_indices_local_global_plus);

			//loop over scalar functionals on cell
			for(unsigned int scalar_functional_n = 0; scalar_functional_n < scalar_functionals_interface[internal_index].size(); ++scalar_functional_n)
			{

				//this contains the local dof indices coupling on the interface cell for the scalar functional under consideration
				const auto& dof_indices_local = coupled_dof_indices_scalar_functionals_interface[internal_index][scalar_functional_n];

				//this contains the local dof indices coupling on the - side for the scalar functional under consideration
				const auto& dof_indices_local_minus = coupled_dof_indices_scalar_functionals_interface_minus[internal_index][scalar_functional_n];

				//this contains the local dof indices coupling on the + side for the scalar functional under consideration
				const auto& dof_indices_local_plus = coupled_dof_indices_scalar_functionals_interface_plus[internal_index][scalar_functional_n];

				//this contains the local C indices coupling on the interface cell for the scalar functional under consideration
				const auto& dof_indices_local_C = coupled_C_indices_scalar_functionals_interface[internal_index][scalar_functional_n];

				//nothing to do if there are no dof's present
				const unsigned int n_dofs_cell = dof_indices_local.size() + dof_indices_local_minus.size() + dof_indices_local_plus.size() + dof_indices_local_C.size();
				if(n_dofs_cell == 0)
					continue;

				//map the local dof indices to the global ones
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local, dof_indices_global, dof_indices_local_global);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_minus, dof_indices_global_minus, dof_indices_local_global_minus);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_plus, dof_indices_global_plus, dof_indices_local_global_plus);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_C, dof_indices_global_C, dof_indices_local_global_C);

				//assemble combined dof indices
				Auxiliary::combine_dof_indices(	dof_indices_global,
												dof_indices_global_minus,
												dof_indices_global_plus,
												dof_indices_global_C,
												dof_indices_interface_dof_indices_combined,
												dof_indices_minus_dof_indices_combined,
												dof_indices_plus_dof_indices_combined,
												dof_indices_C_dof_indices_combined,
												dof_indices_global_combined);

				//restrict solution vectors to cell
				solution.extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_local.begin());
				solution.extract_subvector_to(dof_indices_global_minus.begin(), dof_indices_global_minus.end(), solution_local_minus.begin());
				solution.extract_subvector_to(dof_indices_global_plus.begin(), dof_indices_global_plus.end(), solution_local_plus.begin());
				solution.extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_local_C.begin());
				for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
				{
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_ref_sets_local[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_minus.begin(), dof_indices_global_minus.end(), solution_ref_sets_local_minus[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_plus.begin(), dof_indices_global_plus.end(), solution_ref_sets_local_plus[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_ref_sets_local_C[ref_set].begin());
				}

				//initialize f_cell and K_cell to correct size
				//attention: f and K are not only assembled cell-wise, but also scalar functional wise,
				//			 therefore, the initialization cannot be done outside this loop
				if(assemble_rhs)
					f_cell.reinit(dof_indices_global_combined.size());
				if(assemble_matrix)
					K_cell.reinit(dof_indices_global_combined.size(), dof_indices_global_combined.size());

				//get a pointer to the scalar functional under consideration and its nonprimitive_index and primitive_index
				const auto scalar_functional = scalar_functionals_interface[internal_index][scalar_functional_n];
				const auto index_pair = get_scalar_functional_indices(scalar_functional);
				const int nonprimitive_index = index_pair.first;
				const int primitive_index = index_pair.second;

				//initialize dependent variable vectors appropriately
				e_sigma.reinit(scalar_functional->e_sigma.size());
				for(auto& e_sigma_ref_n : e_sigma_ref_sets)
					e_sigma_ref_n.reinit(e_sigma.size());

				//initialize derivatives of scalar functional w.r.t. dependent variables appropriately
				if(assemble_rhs)
					h_sigma_1.reinit(e_sigma.size());
				if(assemble_matrix)
				{
					h_sigma_2.reinit(e_sigma.size(), e_sigma.size());
					h_sigma_2_de_sigma_dsol_T.reinit(e_sigma.size(), dof_indices_global_combined.size());
				}

				//get JxW values and quadrature points
				const auto& JxW = fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_interface().get_JxW_values();
				const auto& q_points=fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_interface().get_quadrature_points();

				//get normal vectors at quadrature points
				const auto& normals = fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_domain(InterfaceSide::minus).get_normal_vectors();

				//loop over integration points
				for(unsigned int q_point = 0; q_point < scalar_functional->quadrature.size(); ++q_point)
				{
					//compute dependent variables
					compute_e_sigma(internal_index,
									scalar_functional_n,
									q_point,
									solution_local,
									solution_local_minus,
									solution_local_plus,
									solution_local_C,
									dof_indices_interface_dof_indices_combined,
									dof_indices_minus_dof_indices_combined,
									dof_indices_plus_dof_indices_combined,
									dof_indices_C_dof_indices_combined,
									dof_indices_global_combined,
									e_sigma,
									de_sigma_dsol_T,
									assemble_rhs);

					//compute reference values of dependent variables
					for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
						compute_e_sigma(internal_index,
										scalar_functional_n,
										q_point,
										solution_ref_sets_local[ref_set],
										solution_ref_sets_local_minus[ref_set],
										solution_ref_sets_local_plus[ref_set],
										solution_ref_sets_local_C[ref_set],
										dof_indices_interface_dof_indices_combined,
										dof_indices_minus_dof_indices_combined,
										dof_indices_plus_dof_indices_combined,
										dof_indices_C_dof_indices_combined,
										dof_indices_global_combined,
										e_sigma_ref_sets[ref_set],
										de_sigma_dsol_T,
										false);

					//hidden variables
					Assert(interface_cell_domain_cells.interface_cell->user_pointer() != nullptr, ExcMessage("Internal error - it seems that the hidden variables are not propery allocated!"));
					Vector<double>& hidden_vars = (*static_cast<vector<vector<Vector<double>>>*>(interface_cell_domain_cells.interface_cell->user_pointer()))[scalar_functional_n][q_point];;

					//the location of the quadrature point in real space
					const auto& x = q_points[q_point];

					//the normal vector at the quadrature point
					const auto& n = normals[q_point];

					//evaluate the integrand of the scalar functional and its derivatives w.r.t. the independent variables
					if(scalar_functional->get_h_sigma(e_sigma, e_sigma_ref_sets, hidden_vars, x, n, h_sigma, h_sigma_1, h_sigma_2, requested_quantities))
						error = true;

					//add contribution to the value of the scalar functional if it enters the total potential primitively
					//(values of scalar functionals entering the total potential non-primitively have been computed earlier
					// and need not be taken into account again)
					if(compute_potential && (primitive_index > -1))
						primitive_scalar_functionals_values[primitive_index] += JxW[q_point] * h_sigma;

					//add contributions to f_cell and K_cell
					if(assemble_rhs)
					{
						h_sigma_1 *= JxW[q_point];
						de_sigma_dsol_T.vmult(f_cell, h_sigma_1, true);
					}
					if(assemble_matrix)
					{
						h_sigma_2 *= JxW[q_point];
						h_sigma_2.mTmult(h_sigma_2_de_sigma_dsol_T, de_sigma_dsol_T);
						de_sigma_dsol_T.mmult(K_cell, h_sigma_2_de_sigma_dsol_T, true);
					}

					//check that quadrature points are aligned with each other
					Assert(	x.distance(fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_domain(InterfaceSide::minus).quadrature_point(q_point)) < 1e-8,
							ExcMessage("Internal error: Quadrature points are not aligned on interface or boundary!"));
					if(interface_cell_domain_cells.refinement_case != InterfaceRefinementCase::at_boundary)
					{
						Assert(	x.distance(fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_domain(InterfaceSide::plus).quadrature_point(q_point)) < 1e-8,
								ExcMessage("Internal error: Quadrature points are not aligned on interface or boundary!"));
					}
				}

				//additional work to do for scalar functionals entering the total potential non-primitively
				//(lower left and upper right part of stretch system matrix)
				if((nonprimitive_index > -1) && assemble_matrix)
				{
					//compute  V_T_L_T_cell and L_U_cell
					V_T_L_T_cell.reinit(n_scalar_functionals_nonprimitive + C.size(), f_cell.size());
					L_U_cell.reinit(f_cell.size(), n_scalar_functionals_nonprimitive + C.size());
					for(unsigned int k = 0; k < V_T_L_T_cell.size()[0]; ++k)
						for(unsigned int l = 0; l < V_T_L_T_cell.size()[1]; ++l)
							V_T_L_T_cell(k,l) = R[nonprimitive_index](k) * f_cell(l);
					for(unsigned int k = 0; k < L_U_cell.size()[0]; ++k)
						for(unsigned int l = 0; l < L_U_cell.size()[1]; ++l)
							L_U_cell(k,l) = L[nonprimitive_index](l) * f_cell(k);

					//distribute result to system matrix (or respective vectors)
					constraints.distribute_local_to_global(V_T_L_T_cell, nonlinear_scalar_functional_indices, dof_indices_global_combined, K);
					constraints.distribute_local_to_global(L_U_cell, dof_indices_global_combined, nonlinear_scalar_functional_indices, K);
					//take the effect of constraint inhomogeneities into account (though the rhs related to V_T_L_T_cell and L_U_cell, respectively,
					//is zero, there may still be an effect of inhomogeneous constraints on the system rhs)
					rhs_cell_dummy.reinit(nonlinear_scalar_functional_indices.size());
					constraints.distribute_local_to_global(rhs_cell_dummy, nonlinear_scalar_functional_indices, dof_indices_global_combined, f, V_T_L_T_cell);
					rhs_cell_dummy.reinit(dof_indices_global_combined.size());
					constraints.distribute_local_to_global(rhs_cell_dummy, dof_indices_global_combined, nonlinear_scalar_functional_indices, f, L_U_cell);
				}

				//scale f_cell and K_cell according to derivative of total potential contribution w.r.t. scalar functional
				if(assemble_rhs)
				{
					if(nonprimitive_index>-1)
						f_cell *= dPi_dH_omega_H_sigma_C_nonprimitive[nonprimitive_index];
					else
						f_cell *= dPi_dH_omega_H_sigma_primitive[primitive_index];
				}

				if(assemble_matrix)
				{
					if(nonprimitive_index>-1)
						K_cell *= dPi_dH_omega_H_sigma_C_nonprimitive[nonprimitive_index];
					else
						K_cell *= dPi_dH_omega_H_sigma_primitive[primitive_index];
				}

				//rhs is negative of gradient
				//(up to this point the contribution of the scalar functional on the cell to the gradient has been computed)
				if(assemble_rhs)
					f_cell *= -1.0;

				//distribute local contributions to K and f
				if(assemble_rhs && assemble_matrix)
					constraints.distribute_local_to_global(K_cell, f_cell, dof_indices_global_combined, K, f, false);
				else if(assemble_rhs)
					constraints.distribute_local_to_global(f_cell, dof_indices_global_combined, f);
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// contributions related to C's (similar to contributions entering total potential non-primitvely) //
	/////////////////////////////////////////////////////////////////////////////////////////////////////

	//do this on last processor only
	if( assemble_rhs && (this_proc == n_procs - 1) )
	{
		//loop over C's
		for(const auto& C_n : C)
		{
			//get the global dof index of the C under consideration
			unsigned int global_dof_index_C = get_global_dof_index_C(C_n);

			//index into dPi_dH_omega_H_sigma_C_nonprimitive
			const int nonprimitive_index = n_scalar_functionals_nonprimitive + global_indices_C.at(C_n);

			//contribution of C to the RHS
			f[global_dof_index_C] -= dPi_dH_omega_H_sigma_C_nonprimitive[nonprimitive_index];

			//lower left and upper right part of system matrix
			if(assemble_matrix)
			{
				//compute  V_T_L_T_cell and L_U_cell
				V_T_L_T_cell.reinit(nonlinear_scalar_functional_indices.size(), 1);
				L_U_cell.reinit(1, nonlinear_scalar_functional_indices.size());
				for(unsigned int k = 0; k < V_T_L_T_cell.size()[0]; ++k)
					V_T_L_T_cell(k,0) = R[nonprimitive_index](k);
				for(unsigned int l = 0; l < L_U_cell.size()[1]; ++l)
					L_U_cell(0,l) = L[nonprimitive_index](l);

				//distribute result to system matrix (or respective vectors)
				dof_indices_global.resize(1);
				dof_indices_global[0] = global_dof_index_C;
				constraints.distribute_local_to_global(V_T_L_T_cell, nonlinear_scalar_functional_indices, dof_indices_global, K);
				constraints.distribute_local_to_global(L_U_cell, dof_indices_global, nonlinear_scalar_functional_indices, K);
				//take the effect of constraint inhomogeneities into account (though the rhs related to V_T_L_T_cell and L_U_cell, respectively,
				//is zero, there may still be an effect of inhomogeneous constraints on the system rhs)
				rhs_cell_dummy.reinit(nonlinear_scalar_functional_indices.size());
				constraints.distribute_local_to_global(rhs_cell_dummy, nonlinear_scalar_functional_indices, dof_indices_global, f, V_T_L_T_cell);
				rhs_cell_dummy.reinit(dof_indices_global.size());
				constraints.distribute_local_to_global(rhs_cell_dummy, dof_indices_global, nonlinear_scalar_functional_indices, f, L_U_cell);

				// don't forget to set diagonal to 1 if dof is constrained
				if(constraints.is_constrained(global_dof_index_C))
					K.add(global_dof_index_C, global_dof_index_C, 1.0);
			}
		}
	}

	////////////////////////////////////
	// finalize RHS and system matrix //
	////////////////////////////////////

	if(assemble_rhs)
	{
		f.compress(VectorOperation::add);
		constraints.set_zero(f);
		f.compress(VectorOperation::insert);
	}
	if(assemble_matrix)
		K.compress(VectorOperation::add);

/***************************************************************
 * add linear total potential contributions to potential value *
 ***************************************************************/

	if(compute_potential)
	{
		//add up contributions of different processors of primitive_scalar_functionals_values
		if(n_procs > 1)
		{
#ifdef DEAL_II_WITH_MPI
			const auto tria_domain_ptr = dynamic_cast<const dealii::parallel::Triangulation<spacedim, spacedim>*>(&(tria_system.get_triangulation_domain()));
			Assert(tria_domain_ptr != nullptr, ExcMessage("Internal error!"));
			int ierr = MPI_Allreduce(MPI_IN_PLACE, primitive_scalar_functionals_values.data(), primitive_scalar_functionals_values.size(), MPI_DOUBLE, MPI_SUM, tria_domain_ptr->get_communicator());
			AssertThrowMPI(ierr);
#else
		Assert(n_procs == 1, ExcMessage("Internal error: deal.II not configured with MPI, but computation run with more than one processor"));
#endif //DEAL_II_WITH_MPI
		}

		//loop over the total potential contributions
		for(const auto total_potential_contribution : total_potential.total_potential_contributions)
		{
			//only work to do if total_potential_contribution is primitive (the non-primitive contributions have been added in the very beginning)
			if(total_potential_contribution->is_primitive)
			{
				//the scalar functional is domain related
				if(total_potential_contribution->H_omega.size() == 1)
					Pi = primitive_scalar_functionals_values[get_scalar_functional_indices(total_potential_contribution->H_omega[0]).second];
				//the scalar functional is interface related
				else
					Pi = primitive_scalar_functionals_values[get_scalar_functional_indices(total_potential_contribution->H_sigma[0]).second];
				//add result to potential_value
				potential_value += Pi;
			}
		}
	}

	if(n_procs > 1)
	{
#ifdef DEAL_II_WITH_MPI
		const auto tria_domain_ptr = dynamic_cast<const dealii::parallel::Triangulation<spacedim, spacedim>*>(&(tria_system.get_triangulation_domain()));
		Assert(tria_domain_ptr != nullptr, ExcMessage("Internal error!"));
		return Auxiliary::communicate_bool(error, tria_domain_ptr->get_communicator());
#else
		Assert(n_procs == 1, ExcMessage("Internal error: deal.II not configured with MPI, but computation run with more than one processor"));
#endif //DEAL_II_WITH_MPI
	}

	return error;

}

template<unsigned int spacedim>
template<class VectorType>
bool
AssemblyHelper<spacedim>::get_nonprimitive_scalar_functional_values(const VectorType&											solution,
																	const vector<const VectorType*>								solution_ref_sets,
																	map<const ScalarFunctional<spacedim, spacedim>*, double>& 	nonprimitive_scalar_functional_values_domain,
																	map<const ScalarFunctional<spacedim-1, spacedim>*, double>& nonprimitive_scalar_functional_values_interface)
const
{
	Vector<double> nonprimitive_scalar_functionals_values(n_scalar_functionals_nonprimitive);
	const bool error = get_nonprimitive_scalar_functional_values(solution, solution_ref_sets, nonprimitive_scalar_functionals_values);
	for(const auto& domain_scalar_functional : scalar_functionals_domain_nonprimitive_indices)
		nonprimitive_scalar_functional_values_domain.insert(make_pair(domain_scalar_functional.first, nonprimitive_scalar_functionals_values[domain_scalar_functional.second]));
	for(const auto& interface_scalar_functional : scalar_functionals_interface_nonprimitive_indices)
		nonprimitive_scalar_functional_values_interface.insert(make_pair(interface_scalar_functional.first, nonprimitive_scalar_functionals_values[interface_scalar_functional.second]));
	return error;
}

template<unsigned int spacedim>
template<class VectorType>
double
AssemblyHelper<spacedim>::get_maximum_step_length(	const VectorType&				solution,
													const vector<const VectorType*>	solution_ref_sets,
													const VectorType&				delta_solution)
const
{

	double max_step_local;
	double max_step_global = DBL_MAX;

	Assert( solution.size() == system_size(), ExcMessage("The solution vector has not the correct size!"));
	for(const auto& solution_ref : solution_ref_sets)
	{
		(void) solution_ref;	//silence unused parameter warning of compiler
		Assert( solution_ref->size() == system_size(), ExcMessage("One of the reference solution vectors has not the correct size!"));
	}
	Assert( delta_solution.size() == system_size(), ExcMessage("The solution vector has not the correct size!"));

	//////////////////
	// domain cells //
	//////////////////

	//stores the material_id of the cell visited previously
	types::material_id material_id_previous = numbers::invalid_material_id;

	//holds the internal index of a cell
	unsigned int internal_index = 0;

	//global dof indices coupling on the cell for a particular scalar functional
	vector<unsigned int> dof_indices_global;
	dof_indices_global.reserve(fe_collection_domain.max_dofs_per_cell() + global_indices_C.size());

	//mapping between local dof indices and global ones (dof_indices_local_global[i] is the global dof index corresponding to the local dof index i)
	vector<unsigned int> dof_indices_local_global;
	dof_indices_local_global.reserve(fe_collection_domain.max_dofs_per_cell());

	//global dof indices of C's
	vector<unsigned int> dof_indices_local_global_C;
	dof_handler_system.get_dof_indices(dof_indices_local_global_C);

	//global dof indices of C's for a particular scalar functional
	vector<unsigned int> dof_indices_global_C;
	dof_indices_global_C.reserve(global_indices_C.size());

	//solution vector restricted to dof indices coupling on a cell for a particular scalar functional
	Vector<double> solution_local;
	solution_local.reinit(fe_collection_domain.max_dofs_per_cell());

	//solution vector restricted to C indices coupling on a cell for a particular scalar functional
	Vector<double> solution_local_C;
	solution_local_C.reinit(global_indices_C.size());

	//reference solution vector restricted to dof indices coupling on a cell for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_n : solution_ref_sets_local)
		solution_ref_sets_local_n.reinit(fe_collection_domain.max_dofs_per_cell());

	//reference solution vector restricted to C indices coupling on a cell for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local_C(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_C_n : solution_ref_sets_local_C)
		solution_ref_sets_local_C_n.reinit(global_indices_C.size());

	//delta solution vector restricted to dof indices coupling on a cell for a particular scalar functional
	Vector<double> delta_solution_local;
	delta_solution_local.reinit(fe_collection_domain.max_dofs_per_cell());

	//delta solution vector restricted to C indices coupling on a cell for a particular scalar functional
	Vector<double> delta_solution_local_C;
	delta_solution_local_C.reinit(global_indices_C.size());


	//vectors for dependent variables
	Vector<double> e_omega(total_potential.max_dependent_vars);
	vector<Vector<double>> e_omega_ref_sets(solution_ref_sets.size());
	Vector<double> delta_e_omega(total_potential.max_dependent_vars);

	//matrix for derivatives of dependent variables w.r.t. local dof's
	FullMatrix<double> de_omega_dsol_T;

	//the actual loop over the cells
	for(const auto& domain_cell : dof_handler_system.domain_active_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			//the internal material index needs only to be updated if the cell has another material_id than the one visited before
			if(domain_cell->material_id() != material_id_previous)
			{
				material_id_previous = domain_cell->material_id();
				internal_index = material_id_to_internal_index_domain.at(material_id_previous);
			}

			//initialize relevant FEValues objects with cell
			initialize_fe_values_domain(domain_cell, internal_index, false);

			//get the mapping between local and global dof indices
			dof_indices_local_global.resize(domain_cell->get_fe().dofs_per_cell);
			domain_cell.get_dof_indices(dof_indices_local_global);

			//loop over scalar functionals on cell
			for(unsigned int scalar_functional_n = 0; scalar_functional_n<scalar_functionals_domain[internal_index].size(); ++scalar_functional_n)
			{

				//this contains the local dof indices coupling on the domain cell for the scalar functional under consideration
				const auto& dof_indices_local = coupled_dof_indices_scalar_functionals_domain[internal_index][scalar_functional_n];

				//this contains the local C indices coupling on the domain cell for the scalar functional under consideration
				const auto& dof_indices_local_C = coupled_C_indices_scalar_functionals_domain[internal_index][scalar_functional_n];

				//nothing to do if there are no dof's present
				if((dof_indices_local.size() + dof_indices_local_C.size()) == 0)
					continue;

				//map the local dof indices to the global ones
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local, dof_indices_global, dof_indices_local_global);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_C, dof_indices_global_C, dof_indices_local_global_C);

				//restrict solution vectors to scalar functional on cell
				solution.extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_local.begin());
				solution.extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_local_C.begin());
				for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
				{
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_ref_sets_local[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_ref_sets_local_C[ref_set].begin());
				}
				delta_solution.extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), delta_solution_local.begin());
				delta_solution.extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), delta_solution_local_C.begin());

				//append global dof indices of C to global dof index vector
				for(const auto& dof_index_global_C : dof_indices_global_C)
					dof_indices_global.push_back(dof_index_global_C);

				//get a pointer to the scalar functional under consideration and its nonprimitive_index and primitive_index
				const auto scalar_functional = scalar_functionals_domain[internal_index][scalar_functional_n];

				//initialize dependent variable vectors appropriately
				e_omega.reinit(scalar_functional->e_omega.size());
				for(auto& e_omega_ref_n : e_omega_ref_sets)
					e_omega_ref_n.reinit(e_omega.size());
				delta_e_omega.reinit(scalar_functional->e_omega.size());

				//get quadrature points
				const auto& q_points = fe_values_domain[internal_index][scalar_functional_n]->get_quadrature_points();

				//loop over integration points
				for(unsigned int q_point = 0; q_point < scalar_functional->quadrature.size(); ++q_point)
				{

					//compute dependent variables
					compute_e_omega(internal_index, scalar_functional_n, q_point, solution_local, solution_local_C, e_omega, de_omega_dsol_T, false);

					//compute reference values of dependent variables
					for(unsigned int ref_set=0; ref_set<solution_ref_sets.size(); ++ref_set)
						compute_e_omega(internal_index, scalar_functional_n, q_point, solution_ref_sets_local[ref_set], solution_ref_sets_local_C[ref_set], e_omega_ref_sets[ref_set], de_omega_dsol_T, false);

					//compute dependent variable increments
					compute_e_omega(internal_index, scalar_functional_n, q_point, delta_solution_local, delta_solution_local_C, delta_e_omega, de_omega_dsol_T, false, true);

					//hidden variables
					Assert(domain_cell->user_pointer() != nullptr, ExcMessage("Internal error - it seems that the hidden variables are not propery allocated!"));
					Vector<double>& hidden_vars = (*static_cast<vector<vector<Vector<double>>>*>(domain_cell->user_pointer()))[scalar_functional_n][q_point];

					//the location of the quadrature point in real space
					const auto& x = q_points[q_point];

					//evaluate the integrand of the scalar functional and its derivatives w.r.t. the independent variables
					max_step_local = scalar_functional->get_maximum_step(e_omega, e_omega_ref_sets, delta_e_omega, hidden_vars, x);
					if(max_step_local < max_step_global)
						max_step_global = max_step_local;
				}
			}
		}
	}

	////////////////////////////////////////
	// contributions related to interface //
	////////////////////////////////////////

	//stores the material_ids of the cell visited previously
	tuple<types::material_id, types::material_id, types::material_id> material_ids_previous = make_tuple(	numbers::invalid_material_id,
																											numbers::invalid_material_id,
																											numbers::invalid_material_id);

	//global dof indices coupling on the interface cell for a particular scalar functional
	dof_indices_global.reserve(fe_collection_interface.max_dofs_per_cell());

	//global dof indices coupling on - side for a particular scalar functional
	vector<unsigned int> dof_indices_global_minus;
	dof_indices_global_minus.reserve(fe_collection_domain.max_dofs_per_cell());

	//global dof indices coupling on + side for a particular scalar functional
	vector<unsigned int> dof_indices_global_plus;
	dof_indices_global_plus.reserve(fe_collection_domain.max_dofs_per_cell());

	//mapping between local dof indices and global ones on interface
	dof_indices_local_global.reserve(fe_collection_interface.max_dofs_per_cell());

	//mapping between local dof indices and global ones on - side
	vector<unsigned int> dof_indices_local_global_minus;
	dof_indices_local_global_minus.reserve(fe_collection_domain.max_dofs_per_cell());

	//mapping between local dof indices and global ones on + side
	vector<unsigned int> dof_indices_local_global_plus;
	dof_indices_local_global_plus.reserve(fe_collection_domain.max_dofs_per_cell());

	//combined global dof indices coupling on the interface cell and the adjacent domain cells for a particular scalar functional
	//(without any duplicate dof's -> this quantity is NOT just a concatenation of dof_indices_global, dof_indices_global_minus, dof_indices_global_plus)
	vector<unsigned int> dof_indices_global_combined;
	dof_indices_global_combined.reserve(	  fe_collection_interface.max_dofs_per_cell()
											+ 2*fe_collection_domain.max_dofs_per_cell()
											+ global_indices_C.size());

	//defined such that dof_indices_global_combined[dof_indices_interface_dof_indices_combined[i]] = dof_indices_global[i]
	vector<unsigned int> dof_indices_interface_dof_indices_combined;
	dof_indices_interface_dof_indices_combined.reserve(fe_collection_interface.max_dofs_per_cell());

	//defined such that dof_indices_global_combined[dof_indices_minus_dof_indices_combined[i]] = dof_indices_global_minus[i]
	vector<unsigned int> dof_indices_minus_dof_indices_combined;
	dof_indices_minus_dof_indices_combined.reserve(fe_collection_domain.max_dofs_per_cell());

	//defined such that dof_indices_global_combined[dof_indices_plus_dof_indices_combined[i]] = dof_indices_global_plus[i]
	vector<unsigned int> dof_indices_plus_dof_indices_combined;
	dof_indices_plus_dof_indices_combined.reserve(fe_collection_domain.max_dofs_per_cell());

	//defined such that dof_indices_global_combined[dof_indices_C_dof_indices_combined[i]] = dof_indices_global_C[i]
	vector<unsigned int> dof_indices_C_dof_indices_combined;
	dof_indices_C_dof_indices_combined.reserve(global_indices_C.size());

	//solution vector restricted to dof indices coupling on an interface cell for a particular scalar functional
	solution_local.reinit(fe_collection_interface.max_dofs_per_cell());

	//reference solution vector restricted to dof indices coupling on an interface cell for a particular scalar functional
	for(auto& solution_ref_sets_local_n : solution_ref_sets_local)
		solution_ref_sets_local_n.reinit(fe_collection_interface.max_dofs_per_cell());

	//delta solution vector restricted to dof indices coupling on an interface cell for a particular scalar functional
	delta_solution_local.reinit(fe_collection_interface.max_dofs_per_cell());

	//solution vector restricted to dof indices coupling on a domain cell located on the - side of the interface for a particular scalar functional
	Vector<double> solution_local_minus;
	solution_local_minus.reinit(fe_collection_domain.max_dofs_per_cell());

	//reference solution vector restricted to dof indices coupling on a domain cell located on the - side of the interface for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local_minus(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_minus_n : solution_ref_sets_local_minus)
		solution_ref_sets_local_minus_n.reinit(fe_collection_domain.max_dofs_per_cell());

	//delta solution vector restricted to dof indices coupling on a domain cell located on the - side of the interface for a particular scalar functional
	Vector<double> delta_solution_local_minus;
	delta_solution_local_minus.reinit(fe_collection_domain.max_dofs_per_cell());

	//solution vector restricted to dof indices coupling on a domain cell located on the + side of the interface for a particular scalar functional
	Vector<double> solution_local_plus;
	solution_local_plus.reinit(fe_collection_domain.max_dofs_per_cell());

	//reference solution vector restricted to dof indices coupling on a domain cell located on the + side of the interface for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local_plus(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_plus_n : solution_ref_sets_local_plus)
		solution_ref_sets_local_plus_n.reinit(fe_collection_domain.max_dofs_per_cell());

	//solution vector restricted to dof indices coupling on a domain cell located on the + side of the interface for a particular scalar functional
	Vector<double> delta_solution_local_plus;
	delta_solution_local_plus.reinit(fe_collection_domain.max_dofs_per_cell());

	//vectors for dependent variables
	Vector<double> e_sigma(total_potential.max_dependent_vars);
	vector<Vector<double>> e_sigma_ref_sets(solution_ref_sets.size());
	Vector<double> delta_e_sigma(total_potential.max_dependent_vars);

	//matrix for derivatives of dependent variables w.r.t. local dof's
	FullMatrix<double> de_sigma_dsol_T;

	//the actual loop over the cells
	for(const auto& interface_cell_domain_cells : dof_handler_system.interface_active_iterators())
	{

		if(interface_cell_domain_cells.interface_cell->is_locally_owned())
		{
			//find out the material id's characterizing the interface
			const auto material_ids = interface_cell_domain_cells.get_material_ids();
			//the internal material index needs only to be updated if the cell is associated with different material_id's than the one visited before
			if(material_ids != material_ids_previous)
			{
				internal_index = material_ids_to_internal_index_interface.at(material_ids);
				material_ids_previous = material_ids;
			}

			//initialize FEValues objects with cell
			initialize_fe_values_interface(interface_cell_domain_cells, internal_index, false);

			//get the mapping between local and global dof indices
			interface_cell_domain_cells.get_dof_indices_local_global_interface(dof_indices_local_global, dof_indices_local_global_minus, dof_indices_local_global_plus);

			//loop over scalar functionals on cell
			for(unsigned int scalar_functional_n = 0; scalar_functional_n < scalar_functionals_interface[internal_index].size(); ++scalar_functional_n)
			{

				//this contains the local dof indices coupling on the interface cell for the scalar functional under consideration
				const auto& dof_indices_local = coupled_dof_indices_scalar_functionals_interface[internal_index][scalar_functional_n];

				//this contains the local dof indices coupling on the - side for the scalar functional under consideration
				const auto& dof_indices_local_minus = coupled_dof_indices_scalar_functionals_interface_minus[internal_index][scalar_functional_n];

				//this contains the local dof indices coupling on the + side for the scalar functional under consideration
				const auto& dof_indices_local_plus = coupled_dof_indices_scalar_functionals_interface_plus[internal_index][scalar_functional_n];

				//this contains the local C indices coupling on the interface cell for the scalar functional under consideration
				const auto& dof_indices_local_C = coupled_C_indices_scalar_functionals_interface[internal_index][scalar_functional_n];

				//nothing to do if there are no dof's present
				const unsigned int n_dofs_cell = dof_indices_local.size() + dof_indices_local_minus.size() + dof_indices_local_plus.size() + dof_indices_local_C.size();
				if(n_dofs_cell == 0)
					continue;

				//map the local dof indices to the global ones
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local, dof_indices_global, dof_indices_local_global);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_minus, dof_indices_global_minus, dof_indices_local_global_minus);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_plus, dof_indices_global_plus, dof_indices_local_global_plus);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_C, dof_indices_global_C, dof_indices_local_global_C);

				//assemble combined dof indices
				Auxiliary::combine_dof_indices(	dof_indices_global,
												dof_indices_global_minus,
												dof_indices_global_plus,
												dof_indices_global_C,
												dof_indices_interface_dof_indices_combined,
												dof_indices_minus_dof_indices_combined,
												dof_indices_plus_dof_indices_combined,
												dof_indices_C_dof_indices_combined,
												dof_indices_global_combined);

				//restrict solution vectors to cell
				solution.extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_local.begin());
				solution.extract_subvector_to(dof_indices_global_minus.begin(), dof_indices_global_minus.end(), solution_local_minus.begin());
				solution.extract_subvector_to(dof_indices_global_plus.begin(), dof_indices_global_plus.end(), solution_local_plus.begin());
				solution.extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_local_C.begin());
				for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
				{
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_ref_sets_local[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_minus.begin(), dof_indices_global_minus.end(), solution_ref_sets_local_minus[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_plus.begin(), dof_indices_global_plus.end(), solution_ref_sets_local_plus[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_ref_sets_local_C[ref_set].begin());
				}
				delta_solution.extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), delta_solution_local.begin());
				delta_solution.extract_subvector_to(dof_indices_global_minus.begin(), dof_indices_global_minus.end(), delta_solution_local_minus.begin());
				delta_solution.extract_subvector_to(dof_indices_global_plus.begin(), dof_indices_global_plus.end(), delta_solution_local_plus.begin());
				delta_solution.extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), delta_solution_local_C.begin());

				//get a pointer to the scalar functional under consideration and its nonprimitive_index and primitive_index
				const auto scalar_functional = scalar_functionals_interface[internal_index][scalar_functional_n];

				//initialize dependent variable vectors appropriately
				e_sigma.reinit(scalar_functional->e_sigma.size());
				for(auto& e_sigma_ref_n : e_sigma_ref_sets)
					e_sigma_ref_n.reinit(e_sigma.size());
				delta_e_sigma.reinit(scalar_functional->e_sigma.size());

				//get quadrature points
				const auto& q_points=fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_interface().get_quadrature_points();

				//get normal vectors at quadrature points
				const auto& normals = fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_domain(InterfaceSide::minus).get_normal_vectors();

				//loop over integration points
				for(unsigned int q_point = 0; q_point < scalar_functional->quadrature.size(); ++q_point)
				{
					//compute dependent variables
					compute_e_sigma(internal_index,
									scalar_functional_n,
									q_point,
									solution_local,
									solution_local_minus,
									solution_local_plus,
									solution_local_C,
									dof_indices_interface_dof_indices_combined,
									dof_indices_minus_dof_indices_combined,
									dof_indices_plus_dof_indices_combined,
									dof_indices_C_dof_indices_combined,
									dof_indices_global_combined,
									e_sigma,
									de_sigma_dsol_T,
									false);

					//compute reference values of dependent variables
					for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
						compute_e_sigma(internal_index,
										scalar_functional_n,
										q_point,
										solution_ref_sets_local[ref_set],
										solution_ref_sets_local_minus[ref_set],
										solution_ref_sets_local_plus[ref_set],
										solution_ref_sets_local_C[ref_set],
										dof_indices_interface_dof_indices_combined,
										dof_indices_minus_dof_indices_combined,
										dof_indices_plus_dof_indices_combined,
										dof_indices_C_dof_indices_combined,
										dof_indices_global_combined,
										e_sigma_ref_sets[ref_set],
										de_sigma_dsol_T,
										false);

					//compute dependent variable increments
					compute_e_sigma(internal_index,
									scalar_functional_n,
									q_point,
									delta_solution_local,
									delta_solution_local_minus,
									delta_solution_local_plus,
									delta_solution_local_C,
									dof_indices_interface_dof_indices_combined,
									dof_indices_minus_dof_indices_combined,
									dof_indices_plus_dof_indices_combined,
									dof_indices_C_dof_indices_combined,
									dof_indices_global_combined,
									delta_e_sigma,
									de_sigma_dsol_T,
									false,
									true);

					//hidden variables
					Assert(interface_cell_domain_cells.interface_cell->user_pointer() != nullptr, ExcMessage("Internal error - it seems that the hidden variables are not propery allocated!"));
					Vector<double>& hidden_vars = (*static_cast<vector<vector<Vector<double>>>*>(interface_cell_domain_cells.interface_cell->user_pointer()))[scalar_functional_n][q_point];;

					//the location of the quadrature point in real space
					const auto& x = q_points[q_point];

					//the normal vector at the quadrature point
					const auto& n = normals[q_point];

					//evaluate the integrand of the scalar functional and its derivatives w.r.t. the independent variables
					max_step_local = scalar_functional->get_maximum_step(e_sigma, e_sigma_ref_sets, delta_e_sigma, hidden_vars, x, n);
					if(max_step_local < max_step_global)
						max_step_global = max_step_local;
				}
			}
		}
	}

	//exchange local step lengths
	if(n_procs > 1)
	{
#ifdef DEAL_II_WITH_MPI
		const auto tria_domain_ptr = dynamic_cast<const dealii::parallel::Triangulation<spacedim, spacedim>*>(&tria_system.get_triangulation_domain());
		Assert(tria_domain_ptr != nullptr, ExcMessage("Internal error!"));

		vector<double> max_step_global_vect(n_procs);
		max_step_global_vect[this_proc] = max_step_global;

		double send_value = max_step_global_vect[this_proc];
		int ierr =  MPI_Allgather(&send_value, 1, MPI_DOUBLE, max_step_global_vect.data(), 1, MPI_DOUBLE, tria_domain_ptr->get_communicator());
		AssertThrowMPI(ierr);

		max_step_global = *min_element(max_step_global_vect.begin(), max_step_global_vect.end());
#else
		Assert(n_procs == 1, ExcMessage("Internal error: deal.II not configured with MPI, but computation run with more than one processor"));
#endif //DEAL_II_WITH_MPI
		}


	return max_step_global;

}


template<unsigned int spacedim>
template<class VectorType>
void
AssemblyHelper<spacedim>::compare_derivatives_with_numerical_derivatives(	const VectorType&				solution,
																			const vector<const VectorType*> solution_ref_sets,
																			const string					detailed_printout_file,
																			const double					epsilon)
const
{
	Assert(solution.size() == system_size(), ExcMessage("The solution vector has not the correct size!"));
	for(const auto& solution_ref : solution_ref_sets)
	{
		(void) solution_ref;	//silence unused parameter warning of compiler
		Assert(solution_ref->size() == system_size(), ExcMessage("One of the reference solution vectors has not the correct size!"));
	}

	//copy solution vector and reference solution vectors
	Vector<double> solution_copy(solution.size());
	for(unsigned int n = 0; n < solution.size(); ++n)
		solution_copy[n] = solution[n];
	vector<Vector<double>> solution_ref_sets_copy(solution_ref_sets.size());
	for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
	{
		solution_ref_sets_copy[ref_set].reinit(system_size());
		for(unsigned int n = 0; n < system_size(); ++n)
			solution_ref_sets_copy[ref_set][n] = (*solution_ref_sets[ref_set])[n];
	}
	vector<const Vector<double>*> solution_ref_sets_copy_pointers(solution_ref_sets.size());
	for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
	{
		solution_ref_sets_copy_pointers[ref_set] = &(solution_ref_sets_copy[ref_set]);
	}

	pout << "START CHECKING DERIVATIVES\n";

	//first of all compute system_matrix and rhs without taking into account any constraints
	AffineConstraints<double> constraints;
	constraints.close();
	DynamicSparsityPattern dsp;
	dsp.reinit(system_size(), system_size());
	generate_sparsity_pattern_by_simulation(dsp, constraints);
	SparsityPattern sp;
	sp.copy_from(dsp);
	double potential_value = 0.0;
	Vector<double> rhs(system_size());
	SparseMatrix<double> system_matrix(sp);
	bool assembly_error = assemble_system(solution_copy, solution_ref_sets_copy_pointers, constraints, potential_value, rhs, system_matrix);
	(void)assembly_error;	//silence unused parameter warning of compiler
	Assert(assembly_error == false, ExcMessage("Error during assembly!"));

	//the number of dof's
	const unsigned int n_dofs = dof_handler_system.n_dofs();

	//the extra number of rows in the stretched system
	const unsigned int n_nonprimitive_rows = system_matrix.m()-n_dofs;

	//the derivative of the total potential w.r.t. the dof's
	Vector<double> dPi_du(n_dofs);
	for(unsigned int n = 0; n < n_dofs; n++)
		dPi_du[n] = -rhs[n];

	//the second derivative of the total potential w.r.t. the dof's
	FullMatrix<double> 	d2Pi_du2(n_dofs);
	if(n_nonprimitive_rows > 0)
	{
		//that's necessary to be able to do the following operations without
		//thinking about the sparsity pattern of the system matrix
		FullMatrix<double> system_matrix_copy(n_dofs + n_nonprimitive_rows);
		system_matrix_copy.copy_from(system_matrix);

		FullMatrix<double> system_matrix_01(n_dofs, n_nonprimitive_rows);
		for(unsigned int m = 0; m < n_dofs; m++)
			for(unsigned int n = 0; n < n_nonprimitive_rows; n++)
				system_matrix_01(m, n)=system_matrix_copy(m, n+n_dofs);

		FullMatrix<double> system_matrix_11_inv(n_nonprimitive_rows, n_nonprimitive_rows);
		for(unsigned int m = 0; m < n_nonprimitive_rows; m++)
				system_matrix_11_inv(m, m) = 1.0/system_matrix_copy(m+n_dofs, m+n_dofs);

		FullMatrix<double> system_matrix_11_inv_system_matrix_01_T(n_nonprimitive_rows, n_dofs);
		system_matrix_11_inv.mTmult(system_matrix_11_inv_system_matrix_01_T, system_matrix_01);

		system_matrix_01.mmult(d2Pi_du2, system_matrix_11_inv_system_matrix_01_T);
		for(unsigned int m = 0; m < n_dofs; m++)
			for(unsigned int n = 0; n < n_dofs; n++)
				d2Pi_du2(m, n) = system_matrix_copy(m, n) - d2Pi_du2(m, n);
	}
	else
		d2Pi_du2.copy_from(system_matrix);

	//compute first numerical derivative
	Vector<double> dPi_du_numerical(n_dofs);
	double potential_value_ref = potential_value;
	for(unsigned int m = 0; m < n_dofs; m++)
	{
		solution_copy[m] += epsilon;
		solution_copy.compress(VectorOperation::add);

		assembly_error = assemble_system(solution_copy, solution_ref_sets_copy_pointers, constraints, potential_value, rhs, system_matrix, make_tuple(true, false, false));
		Assert(assembly_error == false, ExcMessage("Error during assembly!"));

		dPi_du_numerical[m] = (potential_value - potential_value_ref)/epsilon;
		solution_copy[m] -= epsilon;
		solution_copy.compress(VectorOperation::add);
	}

	//compute second numerical derivative
	FullMatrix<double> d2Pi_du2_numerical(n_dofs);
	Vector<double> dPi_du_ref = dPi_du;
	for(unsigned int m = 0; m < n_dofs; m++)
	{
		solution_copy[m] += epsilon;
		solution_copy.compress(VectorOperation::add);

		assembly_error = assemble_system(solution_copy, solution_ref_sets_copy_pointers, constraints, potential_value, rhs, system_matrix, make_tuple(false, true, false));
		Assert( assembly_error == false, ExcMessage("Error during assembly!"));

		for(unsigned int n = 0; n < n_dofs; n++)
			d2Pi_du2_numerical(n, m) = (-rhs[n] - dPi_du_ref[n])/epsilon;
		solution_copy[m] -= epsilon;
		solution_copy.compress(VectorOperation::add);
	}

	//check first derivative
	const double dPi_du_norm = dPi_du.linfty_norm();
	double max_relative_error_1 = 0.0;
	for(unsigned int m = 0; m < n_dofs; m++)
	{
		const double relative_error = (dPi_du_numerical[m] - dPi_du[m])/dPi_du_norm;
		if(fabs(relative_error) > max_relative_error_1)
			max_relative_error_1 = fabs(relative_error);
	}

	//check second derivative
	double d2Pi_du2_norm = d2Pi_du2.linfty_norm();
	double max_relative_error_2 = 0.0;
	for(unsigned int m = 0; m < n_dofs; m++)
	{
		for(unsigned int n = 0; n < n_dofs; n++)
		{
			const double relative_error = (d2Pi_du2_numerical(m, n) - d2Pi_du2(m,n))/d2Pi_du2_norm;
			if(fabs(relative_error) > max_relative_error_2)
				max_relative_error_2 = fabs(relative_error);
		}
	}

	//output to file if requested
	if(detailed_printout_file != "")
	{
		FILE* printout = fopen(detailed_printout_file.c_str(), "w");

		//first derivative
		fprintf(printout, "-->FIRST DERIVATIVE:\n");
		fprintf(printout, "    m  |  numerical      |  exact          |  relative deviation (in infinity norm)\n");
		for(unsigned int m = 0; m < n_dofs; m++)
		{
			const double relative_error=(dPi_du_numerical[m] - dPi_du[m])/dPi_du_norm;
			fprintf(printout, "%5i  | %- 1.8e | %- 1.8e | %- 1.8e\n", m, dPi_du_numerical[m], dPi_du[m], relative_error);
		}
		fprintf(printout, "                                             %- 1.8e MAX\n", max_relative_error_1);

		//second derivative
		fprintf(printout, "\n-->SECOND DERIVATIVE:\n");
		fprintf(printout, "    m,     n  |  numerical      |  exact          |  relative deviation (in infinity norm)\n");
		for(unsigned int m = 0; m < n_dofs; m++)
		{
			for(unsigned int n = 0; n < n_dofs; n++)
			{
				const double relative_error = (d2Pi_du2_numerical(m, n) - d2Pi_du2(m, n))/d2Pi_du2_norm;
				fprintf(printout, "%5i  %5i  | %- 1.8e | %- 1.8e | %- 1.8e\n", m, n, d2Pi_du2_numerical(m, n), d2Pi_du2(m, n), relative_error);
			}
		}
		fprintf(printout, "                                                    %- 1.8e MAX\n", max_relative_error_2);

		fclose(printout);
	}

	printf("  Maximum relative deviation between numerical first derivative and directly computed first derivative (in infinity norm)   = % 10.8f\n", max_relative_error_1);
	printf("  Maximum relative deviation between numerical second derivative and directly computed second derivative (in infinity norm) = % 10.8f\n", max_relative_error_2);
	printf("                                                                                                   Value of total potential = % 10.8f\n", potential_value_ref);

	pout << "FINISHED CHECKING DERIVATIVES\n";

	return;
}

template<unsigned int spacedim>
template<class VectorType>
pair<const string, const string>
AssemblyHelper<spacedim>::write_output_independent_fields(	const VectorType&												solution,
															const string													file_name_domain,
															const string													file_name_interface,
															const unsigned int												file_index,
															const vector<SmartPointer<const DataPostprocessor<spacedim>>>&	dp_domain,
															const vector<SmartPointer<const DataPostprocessor<spacedim>>>&	dp_interface,
															const unsigned int												n_subdivisions)
const
{
	Assert( solution.size() == system_size(), ExcMessage("The solution vector has not the correct size!"));

	VectorType solution_domain, solution_interface, solution_C;
	dof_handler_system.split_vector(solution, solution_domain, solution_interface, solution_C);

	string processor_number = "";
	if(n_procs > 1)
		processor_number = "-proc" + Utilities::to_string(this_proc, 4);

	//file number and file names
	const string file_number = Utilities::to_string(file_index, 3);
	const string file_name_domain_base = file_name_domain + "-" + file_number;
	const string file_name_interface_base = file_name_interface + "-" + file_number;

	size_t pos = file_name_domain_base.find_last_of("/\\") + 1;
	const string file_name_domain_base_wo_folder= file_name_domain_base.substr(pos);
	pos = file_name_interface_base.find_last_of("/\\") + 1;
	const string file_name_interface_base_wo_folder= file_name_interface_base.substr(pos);


/*************************
 * domain related output *
 *************************/

	if(file_name_domain != "")
	{
		//DataOut object for domain
		DataOut<spacedim,hp::DoFHandler<spacedim> > data_out_domain;

		//solution names and component interpretation of domain quantities
		vector<string> solution_names_domain;
		vector<DataComponentInterpretation::DataComponentInterpretation> interpretation_domain;
		for(const auto& u_omega_n : u_omega)
		{
			for(unsigned int component = 0; component < u_omega_n->n_components; ++component)
			{
				solution_names_domain.push_back(u_omega_n->name + "_" + Utilities::to_string(component));
				if(u_omega_n->n_components == spacedim)
					interpretation_domain.push_back(DataComponentInterpretation::component_is_part_of_vector);
				else
					interpretation_domain.push_back(DataComponentInterpretation::component_is_scalar);
			}
		}

		//write output of domain related fields
		data_out_domain.attach_dof_handler(dof_handler_system.get_dof_handler_domain());
		if(solution_domain.size() > 0)
		{
			data_out_domain.add_data_vector(dof_handler_system.get_dof_handler_domain(), solution_domain, solution_names_domain, interpretation_domain);
			for(const auto& dp_sptr : dp_domain)
				data_out_domain.add_data_vector(dof_handler_system.get_dof_handler_domain(), solution_domain, *dp_sptr);
		}
		data_out_domain.build_patches(*mapping_domain, n_subdivisions, DataOut<spacedim,hp::DoFHandler<spacedim> >::curved_inner_cells);
		ofstream output_domain((file_name_domain_base + processor_number + ".vtu").c_str());
		data_out_domain.write_vtu(output_domain);

		//if there is more than one processor, also write *.pvtu record on processor 0
		if( (n_procs > 1) && (this_proc == 0 ))
		{
			vector<std::string> file_names;
			for(unsigned int proc = 0; proc < n_procs; ++proc)
				file_names.push_back(file_name_domain_base_wo_folder + "-proc" + Utilities::to_string(proc, 4) + ".vtu");
			ofstream pvtu_output(file_name_domain_base + ".pvtu");
			data_out_domain.write_pvtu_record(pvtu_output, file_names);
      }
	}

/****************************
 * interface related output *
 ****************************/

	if(file_name_interface != "")
	{
		//DataOut object for interface
		DataOut<spacedim-1, hp::DoFHandler<spacedim-1, spacedim> > data_out_interface;

		//solution names and component interpretation of interface quantities
		vector<string> solution_names_interface;
		vector<DataComponentInterpretation::DataComponentInterpretation> interpretation_interface;
		for(const auto& u_sigma_n : u_sigma)
		{
			for(unsigned int component = 0; component < u_sigma_n->n_components; ++component)
			{
				solution_names_interface.push_back(u_sigma_n->name + "_" + Utilities::to_string(component));
				if(u_sigma_n->n_components == spacedim)
					interpretation_interface.push_back(DataComponentInterpretation::component_is_part_of_vector);
				else
					interpretation_interface.push_back(DataComponentInterpretation::component_is_scalar);
			}
		}

		//write output of interface related fields
		data_out_interface.attach_dof_handler(dof_handler_system.get_dof_handler_interface());
		if(solution_interface.size() > 0)
		{
			data_out_interface.add_data_vector(dof_handler_system.get_dof_handler_interface(), solution_interface, solution_names_interface, interpretation_interface);
			for(const auto& dp_sptr : dp_interface)
				data_out_interface.add_data_vector(dof_handler_system.get_dof_handler_interface(), solution_interface, *dp_sptr);
		}
		data_out_interface.build_patches(*mapping_interface, n_subdivisions, DataOut<spacedim-1, hp::DoFHandler<spacedim-1, spacedim> >::curved_inner_cells);
		ofstream output_interface((file_name_interface_base + processor_number + ".vtu").c_str());
		data_out_interface.write_vtu(output_interface);

		//if there is more than one processor, also write *.pvtu record on processor 0
		if( (n_procs > 1) && (this_proc == 0 ))
		{
			vector<std::string> file_names;
			for(unsigned int proc = 0; proc < n_procs; ++proc)
				file_names.push_back(file_name_interface_base_wo_folder + "-proc" + Utilities::to_string(proc, 4) + ".vtu");
			ofstream pvtu_output(file_name_interface_base + ".pvtu");
			data_out_interface.write_pvtu_record(pvtu_output, file_names);
      }

	}

	if(n_procs == 1)
		return make_pair(file_name_domain_base + ".vtu", file_name_interface_base + ".vtu");
	else
		return make_pair(file_name_domain_base + ".pvtu", file_name_interface_base + ".pvtu");
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::print_assembly_helper_definition(const bool detailed_printout_shapefuns)
const
{
	//inverse of material_id_to_internal_index_domain
	map<unsigned int, const types::material_id> map_internal_index_to_material_id_domain;
	for(const auto& material_id_to_internal_index_n : material_id_to_internal_index_domain)
		map_internal_index_to_material_id_domain.insert(make_pair(material_id_to_internal_index_n.second, material_id_to_internal_index_n.first));

	//scalar functionals on domain
	pout << "******************************************************" << endl
		 << "****** scalar functionals DEFINED ON THE DOMAIN ******" << endl
		 << "******************************************************" << endl << endl;
	for(unsigned int internal_index = 0; internal_index < scalar_functionals_domain.size(); ++internal_index)
	{
		//print domain portion
		pout	<< "  *** DOMAIN PORTION #"
				<< map_internal_index_to_material_id_domain[internal_index]
				<< " ***"
				<< endl;

		//print finite element in use on the domain portion
		const unsigned int fe_system_id_domain = material_id_to_fe_system_id_domain.at(map_internal_index_to_material_id_domain[internal_index]);
		pout	<< endl
				<< "    -> FE = "
				<< fe_collection_domain[fe_system_id_domain].get_name()
				<< endl
				<< endl;

		//print scalar functional related information on domain portion
		if(scalar_functionals_domain[internal_index].size() == 0)
			pout << endl;
		for(unsigned int scalar_functional_n = 0; scalar_functional_n < scalar_functionals_domain[internal_index].size(); ++scalar_functional_n)
		{
			//name of scalar functional
			pout	<< "    -> scalar functional "
					<< scalar_functionals_domain[internal_index][scalar_functional_n]->name
					<< endl
					<< endl;

			//relation of scalar functional to total potential
			pout	<< "       Relation to total potential : ";
			for(unsigned int contribution = 0; contribution<contributions_scalar_functionals_domain_total_potential.at(scalar_functionals_domain[internal_index][scalar_functional_n]).size(); ++contribution)
				pout	<< "("
						<< contributions_scalar_functionals_domain_total_potential.at(scalar_functionals_domain[internal_index][scalar_functional_n])[contribution].first
						<< "," << contributions_scalar_functionals_domain_total_potential.at(scalar_functionals_domain[internal_index][scalar_functional_n])[contribution].second
						<< ")"
						<< ", "
						<< endl
						<< endl;

			//quadrature scheme used for scalar functional
			pout	<< "       Quadrature = "
					<< fe_values_domain[internal_index][scalar_functional_n]->get_quadrature().size()
					<< endl;
			Assert(	fe_values_domain[internal_index][scalar_functional_n]->get_fe().get_name() == fe_collection_domain[fe_system_id_domain].get_name(),
					ExcMessage("Internal error!"));
			pout << endl;

			//information about dependent fields associated with scalar functional
			for(unsigned int e_omega_n = 0; e_omega_n < scalar_functionals_domain[internal_index][scalar_functional_n]->e_omega.size(); e_omega_n++)
			{
				pout	<< "       "
						<< scalar_functionals_domain[internal_index][scalar_functional_n]->e_omega[e_omega_n].name;
				//detailed printout how dependent fields are related to shape functions
				if(detailed_printout_shapefuns)
				{
					pout	<< " = "
							<< endl;
					for(const auto& a_omega_term : a_omega[internal_index][scalar_functional_n][e_omega_n])
					{
						const double a = get<0>(a_omega_term);
						const unsigned int component = get<1>(a_omega_term);
						const auto& shapefuns = get<2>(a_omega_term);
						const auto& name = component_names_domain[component];
						for(const auto& shapefun : shapefuns)
						{
							const unsigned int shapefun_in_element = fe_collection_domain[fe_system_id_domain].system_to_base_index(shapefun).second;
							pout	<< "              "
									<< a
									<< " * "
									<< name.first
									<< "_"
									<< name.second
									<< "("
									<< shapefun_in_element
									<< ") + "
									<< endl;
						}
					}
					for(const auto& b_omega_term : b_omega[internal_index][scalar_functional_n][e_omega_n])
					{
						const double b = get<0>(b_omega_term);
						const unsigned int component = get<1>(b_omega_term);
						const unsigned int derivative = get<2>(b_omega_term);
						const auto& shapefuns = get<3>(b_omega_term);
						const auto& name=component_names_domain[component];
						for(const auto& shapefun : shapefuns)
						{
							const unsigned int shapefun_in_element = fe_collection_domain[fe_system_id_domain].system_to_base_index(shapefun).second;
							pout	<< "              "
									<< b
									<< " * "
									<< name.first
									<< "_"
									<< name.second
									<< ","
									<< derivative
									<< "("
									<< shapefun_in_element
									<< ") + "
									<< endl;
						}
					}
				}
				pout << endl;
			}
		}
	}

	//inverse of material_ids_to_internal_index_interface
	map<unsigned int, tuple<const types::material_id, const types::material_id, const types::material_id>> map_internal_index_material_ids_interface;
	for(const auto& material_ids_to_internal_index_n : material_ids_to_internal_index_interface)
		map_internal_index_material_ids_interface.insert(make_pair(material_ids_to_internal_index_n.second, material_ids_to_internal_index_n.first));

	//scalar functionals on interface
	pout << "*********************************************************" << endl;
	pout << "****** scalar functionals DEFINED ON THE INTERFACE ******" << endl;
	pout << "*********************************************************" << endl << endl;
	for(unsigned int internal_index = 0; internal_index < scalar_functionals_interface.size(); ++internal_index)
	{
		//print domain (sub)portion
		pout	<< "  *** INTERFACE PORTION #("
				<< get<0>(map_internal_index_material_ids_interface[internal_index])
				<< ","
				<< get<1>(map_internal_index_material_ids_interface[internal_index])
				<< ","
				<< get<2>(map_internal_index_material_ids_interface[internal_index])
				<< ")"
				<< " ***" << endl;

		//print finite element in use on the interface (sub)portion
		const unsigned int fe_system_id_interface=material_id_to_fe_system_id_interface.at(get<0>(map_internal_index_material_ids_interface[internal_index]));
		const unsigned int fe_system_id_domain_minus=material_id_to_fe_system_id_domain.at(get<1>(map_internal_index_material_ids_interface[internal_index]));
		unsigned int fe_system_id_domain_plus;
		if(material_id_to_fe_system_id_domain.find(get<2>(map_internal_index_material_ids_interface[internal_index])) == material_id_to_fe_system_id_domain.end())
			fe_system_id_domain_plus = fe_collection_domain.size()-1;
		else
			fe_system_id_domain_plus = material_id_to_fe_system_id_domain.at(get<2>(map_internal_index_material_ids_interface[internal_index]));
		pout	<< endl
				<< "    -> FE    = "
				<< fe_collection_interface[fe_system_id_interface].get_name()
				<< endl
				<< "    -> FE(-) = "
				<< fe_collection_domain[fe_system_id_domain_minus].get_name()
				<< endl
				<< "    -> FE(+) = "
				<< fe_collection_domain[fe_system_id_domain_plus].get_name()
				<< endl
				<< endl;

		//print scalar functional related information on domain portion
		if(scalar_functionals_interface[internal_index].size()==0)
			pout << endl;
		for(unsigned int scalar_functional_n = 0; scalar_functional_n < scalar_functionals_interface[internal_index].size(); ++scalar_functional_n)
		{
			pout	<< "    -> scalar functional "
					<< scalar_functionals_interface[internal_index][scalar_functional_n]->name
					<< " <-"
					<< endl
					<< endl;

			//relation of scalar functional to total potential
			pout	<< "       Relation to total potential : ";
			for(unsigned int contribution = 0; contribution<contributions_scalar_functionals_interface_total_potential.at(scalar_functionals_interface[internal_index][scalar_functional_n]).size(); ++contribution)
				pout	<< "("
				<< contributions_scalar_functionals_interface_total_potential.at(scalar_functionals_interface[internal_index][scalar_functional_n])[contribution].first
				<< ","
				<< contributions_scalar_functionals_interface_total_potential.at(scalar_functionals_interface[internal_index][scalar_functional_n])[contribution].second
				<< ")"
				<< ", "
				<< endl
				<< endl;

			//quadrature scheme used for scalar functional
			pout	<<"       Quadrature = "
					<< fe_values_interface[internal_index][scalar_functional_n]->get_quadrature().size()
					<< endl
					<< endl;

			//information about dependent fields associated with scalar functional
			for(unsigned int e_sigma_n = 0; e_sigma_n < scalar_functionals_interface[internal_index][scalar_functional_n]->e_sigma.size(); ++e_sigma_n)
			{
				pout	<< "       "
						<< scalar_functionals_interface[internal_index][scalar_functional_n]->e_sigma[e_sigma_n].name;
				//detailed printout how dependent fields are related to shape functions
				if(detailed_printout_shapefuns)
				{
					pout	<< " = "
							<< endl;
					for(const auto& a_sigma_term : a_sigma[internal_index][scalar_functional_n][e_sigma_n])
					{
						const double a = get<0>(a_sigma_term);
						const unsigned int component = get<1>(a_sigma_term);
						const auto& shapefuns = get<2>(a_sigma_term);
						const auto& name = component_names_interface[component];
						for(const auto& shapefun : shapefuns)
						{
							const unsigned int shapefun_in_element = fe_collection_interface[fe_system_id_interface].system_to_base_index(shapefun).second;
							pout	<< "              "
									<< a
									<< " * "
									<< name.first
									<< "_"
									<< name.second
									<< "("
									<< shapefun_in_element
									<< ") + "
									<< endl;
						}
					}
					for(const auto& b_sigma_term : b_sigma[internal_index][scalar_functional_n][e_sigma_n])
					{
						const double b = get<0>(b_sigma_term);
						const unsigned int component = get<1>(b_sigma_term);
						const unsigned int derivative = get<2>(b_sigma_term);
						const auto& shapefuns = get<3>(b_sigma_term);
						const auto& name = component_names_interface[component];
						for(const auto& shapefun : shapefuns)
						{
							const unsigned int shapefun_in_element = fe_collection_interface[fe_system_id_interface].system_to_base_index(shapefun).second;
							pout	<< "              "
									<< b << " * "
									<< name.first
									<< "_"
									<< name.second
									<< ","
									<< derivative
									<< "("
									<< shapefun_in_element
									<< ") + "
									<< endl;
						}
					}
					for(const auto& a_minus_term : a_minus[internal_index][scalar_functional_n][e_sigma_n])
					{
						const double a = get<0>(a_minus_term);
						const unsigned int component = get<1>(a_minus_term);
						const auto& shapefuns = get<2>(a_minus_term);
						const auto& name = component_names_domain[component];
						for(const auto& shapefun : shapefuns)
						{
							const unsigned int shapefun_in_element = fe_collection_domain[fe_system_id_domain_minus].system_to_base_index(shapefun).second;
							pout	<< "              "
									<< a
									<< " * "
									<< name.first
									<< "^-"
									<< "_"
									<< name.second
									<< "("
									<< shapefun_in_element
									<< ") + "
									<< endl;
						}
					}
					for(const auto& b_minus_term : b_minus[internal_index][scalar_functional_n][e_sigma_n])
					{
						const double b = get<0>(b_minus_term);
						const unsigned int component = get<1>(b_minus_term);
						const unsigned int derivative = get<2>(b_minus_term);
						const auto& shapefuns = get<3>(b_minus_term);
						const auto& name = component_names_domain[component];
						for(const auto& shapefun : shapefuns)
						{
							const unsigned int shapefun_in_element = fe_collection_domain[fe_system_id_domain_minus].system_to_base_index(shapefun).second;
							pout	<< "              "
									<< b
									<< " * "
									<< name.first
									<< "^-"
									<< "_"
									<< name.second
									<< ","
									<< derivative
									<< "("
									<< shapefun_in_element
									<< ") + "
									<< endl;
						}
					}
					for(const auto& a_plus_term : a_plus[internal_index][scalar_functional_n][e_sigma_n])
					{
						const double a = get<0>(a_plus_term);
						const unsigned int component = get<1>(a_plus_term);
						const auto& shapefuns = get<2>(a_plus_term);
						const auto& name = component_names_domain[component];
						for(const auto& shapefun : shapefuns)
						{
							const unsigned int shapefun_in_element = fe_collection_domain[fe_system_id_domain_plus].system_to_base_index(shapefun).second;
							pout	<< "              "
									<< a
									<< " * "
									<< name.first
									<< "^+"
									<< "_"
									<< name.second
									<< "("
									<< shapefun_in_element
									<< ") + "
									<< endl;
						}
					}
					for(const auto& b_plus_term : b_plus[internal_index][scalar_functional_n][e_sigma_n])
					{
						const double b = get<0>(b_plus_term);
						const unsigned int component = get<1>(b_plus_term);
						const unsigned int derivative = get<2>(b_plus_term);
						const auto& shapefuns = get<3>(b_plus_term);
						const auto& name = component_names_domain[component];
						for(const auto& shapefun : shapefuns)
						{
							const unsigned int shapefun_in_element = fe_collection_domain[fe_system_id_domain_plus].system_to_base_index(shapefun).second;
							pout	<< "              "
									<< b
									<< " * "
									<< name.first
									<< "^+"
									<< "_"
									<< name.second
									<< ","
									<< derivative
									<< "("
									<< shapefun_in_element
									<< ") + "
									<< endl;
						}
					}
				}
				pout << endl;
			}
		}
	}

}

template<unsigned int spacedim>
template<class VectorType>
std::pair<const double, const double>
AssemblyHelper<spacedim>::compute_distance_to_other_solution(	const VectorType&				solution,
																const VectorType&				other_solution,
																const AssemblyHelper<spacedim>& other_assembly_helper,
																const Quadrature<spacedim>		quadrature_domain,
																const Quadrature<spacedim-1>	quadrature_interface,
																const VectorTools::NormType		norm_type,
																const ComponentMask				component_mask_domain,
																const ComponentMask				component_mask_interface,
																const double					/*exponent*/,
																const Vector<double>			scaling_domain,
																const Vector<double>			scaling_interface)
const
{

	Assert(solution.size() == system_size(), ExcMessage("The solution vector has not the correct size!"));
	Assert(other_solution.size() == other_assembly_helper.system_size(), ExcMessage("The other solution vector has not the correct size!"));
	Assert(	(norm_type == VectorTools::NormType::L2_norm) ||
			(norm_type == VectorTools::NormType::Linfty_norm),
			ExcMessage("Up to now only the L2 norm and the Linfty norm are implemented!"));

	double norm_domain = 0.0;
	double norm_interface = 0.0;

	//get dof handlers and current solution of other AssemblyHelper
	const auto& dof_handler_system_other_assembly_helper = other_assembly_helper.get_dof_handler_system();
	const auto& other_dof_handler_domain = dof_handler_system_other_assembly_helper.get_dof_handler_domain();
	const auto& other_dof_handler_interface = dof_handler_system_other_assembly_helper.get_dof_handler_interface();

	//sizes
	const unsigned int n_dofs_domain = dof_handler_system.n_dofs_domain();
	const unsigned int n_dofs_interface = dof_handler_system.n_dofs_interface();
	const unsigned int other_n_dofs_domain = other_dof_handler_domain.n_dofs();
	const unsigned int other_n_dofs_interface = other_dof_handler_interface.n_dofs();
	const unsigned int n_C = C.size();

	//vectors
	Vector<double> other_solution_domain(other_n_dofs_domain);
	Vector<double> other_solution_interface(other_n_dofs_interface);
	Vector<double> e_domain(n_dofs_domain);
	Vector<double> e_interface(n_dofs_interface);
	Vector<double> e_C(n_C);

	//transfer data into vectors
	for(unsigned int m = 0; m < other_n_dofs_domain; ++m)
		other_solution_domain[m] = other_solution[m];
	for(unsigned int m = 0; m < other_n_dofs_interface; ++m)
		other_solution_interface[m] = other_solution[m + other_n_dofs_domain];

	//interpolate data to mesh of this AssemblyHelper
	VectorTools::interpolate_to_different_mesh(other_dof_handler_domain, other_solution_domain, dof_handler_system.get_dof_handler_domain(), e_domain);
	VectorTools::interpolate_to_different_mesh(other_dof_handler_interface, other_solution_interface, dof_handler_system.get_dof_handler_interface(), e_interface);

	//@todo: Think about hanging node constraints here

	//compute difference vector
	for(unsigned int m = 0; m < n_dofs_domain; ++m)
		e_domain[m] -= solution[m];
	for(unsigned int m = 0; m < n_dofs_interface; ++m)
		e_interface[m] -= solution[m + n_dofs_domain];
	for(unsigned int m = 0; m < n_C; ++m)
		e_C[m] = other_solution[m + other_n_dofs_domain + other_n_dofs_interface] - solution[m + n_dofs_domain + n_dofs_interface];

	//compute domain part of the norm
	hp::MappingCollection<spacedim, spacedim> mapping_collection_domain(*mapping_domain);
	hp::QCollection<spacedim> quadrature_collection_domain(quadrature_domain);
	hp::FEValues<spacedim, spacedim> fe_values_domain(	mapping_collection_domain,
														fe_collection_domain,
														quadrature_collection_domain,
														update_quadrature_points|
														update_JxW_values|
														update_values);

	const unsigned int n_components_domain = fe_collection_domain.n_components();
	Vector<double> scaling_domain_(n_components_domain);
	if(scaling_domain.size() == 0)
	{
		for(auto& el : scaling_domain_)
			el = 1.0;
	}
	else
	{
		Assert(scaling_domain.size() == n_components_domain, ExcMessage("Scaling vector does not match the number of solution components"));
		scaling_domain_ = scaling_domain;
	}
	const unsigned int n_q_points_domain = quadrature_domain.size();
	vector<Vector<double>> function_values_domain;
	function_values_domain.resize(n_q_points_domain);
	for(auto& function_values_domain_n : function_values_domain)
		function_values_domain_n.reinit(n_components_domain);
	for(const auto& cell : dof_handler_system.domain_active_iterators())
	{
		fe_values_domain.reinit(cell);
		const FEValues<spacedim, spacedim>& fe_values_domain_present = fe_values_domain.get_present_fe_values();
		fe_values_domain_present.get_function_values(e_domain, function_values_domain);
		for(unsigned int q = 0; q < n_q_points_domain; ++q)
		{
			if(norm_type == VectorTools::NormType::L2_norm)
			{
				double norm_contribution = 0.0;
				for(unsigned int component = 0; component < function_values_domain[q].size(); ++component)
					if( (component_mask_domain.size() == 0) || (component_mask_domain[component]) )
						norm_contribution += function_values_domain[q][component] * function_values_domain[q][component] * scaling_domain_[component] * scaling_domain_[component];
				norm_domain += fe_values_domain_present.JxW(q) * norm_contribution;
			}
			else if(norm_type == VectorTools::NormType::Linfty_norm)
			{
				for(unsigned int component = 0; component < function_values_domain[q].size(); ++component)
					if( (component_mask_domain.size() == 0) || (component_mask_domain[component]) )
						if(fabs(function_values_domain[q][component] * scaling_domain_[component]) > norm_domain)
							norm_domain = fabs(function_values_domain[q][component] * scaling_domain_[component]);
			}
		}
	}

	//compute interface part of the norm
	hp::MappingCollection<spacedim-1, spacedim> mapping_collection_interface(*mapping_interface);
	hp::QCollection<spacedim-1> quadrature_collection_interface(quadrature_interface);
	hp::FEValues<spacedim-1, spacedim> fe_values_interface(	mapping_collection_interface,
															fe_collection_interface,
															quadrature_collection_interface,
															update_quadrature_points|
															update_JxW_values|
															update_values);

	const unsigned int n_components_interface = fe_collection_interface.n_components();
	Vector<double> scaling_interface_(n_components_interface);
	if(scaling_interface.size() == 0)
	{
		for(auto& el : scaling_interface_)
			el = 1.0;
	}
	else
	{
		Assert(scaling_interface.size() == n_components_interface, ExcMessage("Scaling vector does not match the number of solution components"));
		scaling_interface_ = scaling_interface;
	}
	const unsigned int n_q_points_interface = quadrature_interface.size();
	vector<Vector<double>> function_values_interface;
	function_values_interface.resize(n_q_points_interface);
	for(auto& function_values_interface_n : function_values_interface)
		function_values_interface_n.reinit(n_components_interface);
	for(const auto& interface_cell_domain_cell : dof_handler_system.interface_active_iterators())
	{
		fe_values_interface.reinit(interface_cell_domain_cell.interface_cell);
		const FEValues<spacedim-1, spacedim>& fe_values_interface_present = fe_values_interface.get_present_fe_values();
		fe_values_interface_present.get_function_values(e_interface, function_values_interface);
		for(unsigned int q = 0; q < n_q_points_interface; q++)
		{
			if(norm_type == VectorTools::NormType::L2_norm)
			{
				double norm_contribution = 0.0;
				for(unsigned int component = 0; component < function_values_interface[q].size(); ++component)
					if( (component_mask_interface.size() == 0) || (component_mask_interface[component]) )
						norm_contribution += function_values_interface[q][component] * function_values_interface[q][component] * scaling_interface_[component] * scaling_interface_[component];
				norm_interface += fe_values_interface_present.JxW(q) * norm_contribution;
			}
			else if(norm_type == VectorTools::NormType::Linfty_norm)
			{
				for(unsigned int component = 0; component < function_values_interface[q].size(); ++component)
					if( (component_mask_interface.size() == 0) || (component_mask_interface[component]) )
						if(fabs(function_values_interface[q][component] * scaling_interface_[component]) > norm_interface)
							norm_interface = fabs(function_values_interface[q][component] * scaling_interface_[component]);
			}
		}
	}

	if(norm_type == VectorTools::NormType::L2_norm)
	{
		norm_domain = sqrt(norm_domain);
		norm_interface = sqrt(norm_interface);
	}

	return make_pair(norm_domain, norm_interface);
}

template<unsigned int spacedim>
template<class VectorType>
pair<const double, const double>
AssemblyHelper<spacedim>::compute_distance_to_exact_solution(	const VectorType&				solution,
																const Function<spacedim>&		exact_solution_domain,
																const Function<spacedim>&		exact_solution_interface,
																const Quadrature<spacedim>		quadrature_domain,
																const Quadrature<spacedim-1>	quadrature_interface,
																const VectorTools::NormType		norm_type,
																const ComponentMask				component_mask_domain,
																const ComponentMask				component_mask_interface,
																const double					exponent)
const
{

	Assert(solution.size() == system_size(), ExcMessage("The solution vector has not the correct size!"));

	double norm_domain = 0.0;
	double norm_interface = 0.0;

	//sizes
	const unsigned int n_dofs_domain = dof_handler_system.n_dofs_domain();
	const unsigned int n_dofs_interface = dof_handler_system.n_dofs_interface();
	const unsigned int n_C = C.size();

	//aolurion vector parts
	Vector<double> solution_domain(n_dofs_domain);
	Vector<double> solution_interface(n_dofs_interface);
	Vector<double> solution_C(n_C);

	//transfer data into vectors
	for(unsigned int m = 0; m < n_dofs_domain; ++m)
		solution_domain[m] = solution[m];
	for(unsigned int m = 0; m < n_dofs_interface; ++m)
		solution_interface[m] = solution[m + n_dofs_domain];;
	for(unsigned int m = 0; m < n_C; ++m)
		solution_C[m] = solution[m + n_dofs_domain + n_dofs_interface];

	//compute norms
	if(u_omega.size() > 0)
	{
		Vector<double> cell_differences(tria_system.get_triangulation_domain().n_active_cells());
		const hp::MappingCollection<spacedim> mapping_collection_domain(*mapping_domain);
		const hp::QCollection<spacedim> quadrature_collection_domain(quadrature_domain);
		if(component_mask_domain.represents_the_all_selected_mask())
		{
			const Function<spacedim>* weight_fun = nullptr;
			VectorTools::integrate_difference(mapping_collection_domain, dof_handler_system.get_dof_handler_domain(), solution_domain, exact_solution_domain, cell_differences, quadrature_collection_domain, norm_type, weight_fun, exponent);
		}
		else
		{
			Assert(component_mask_domain.size() == component_names_domain.size(), ExcMessage("The component mask has not the correct size!"));
			const unsigned int first_component = component_mask_domain.first_selected_component();
			const unsigned int n_components = component_mask_domain.n_selected_components();
			Assert(n_components > 0, ExcMessage("There is no component selected by the component mask!"));
			const unsigned int last_component = first_component + n_components - 1;
			for(unsigned int component = first_component; component <= last_component; ++component)
			{
				if(component_mask_domain[component] != true)
				{
					Assert(false, ExcMessage("Currently there is only support for component masks where the selected components form a contiguous range!"));
				}
			}
			const ComponentSelectFunction<spacedim, double> weight_fun(pair<unsigned int, unsigned int>(first_component, last_component+1), component_mask_domain.size());
			VectorTools::integrate_difference(mapping_collection_domain, dof_handler_system.get_dof_handler_domain(), solution_domain, exact_solution_domain, cell_differences, quadrature_collection_domain, norm_type, &weight_fun, exponent);
		}
		norm_domain = VectorTools::compute_global_error(tria_system.get_triangulation_domain(), cell_differences, norm_type, exponent);
	}

	if(u_sigma.size()>0)
	{
		Vector<double> cell_differences(tria_system.get_triangulation_interface().n_active_cells());
		const hp::MappingCollection<spacedim-1, spacedim> mapping_collection_interface(*mapping_interface);
		const hp::QCollection<spacedim-1> quadrature_collection_interface(quadrature_interface);
		if(component_mask_interface.represents_the_all_selected_mask())
		{
			const Function<spacedim>* weight_fun = nullptr;
			VectorTools::integrate_difference(mapping_collection_interface, dof_handler_system.get_dof_handler_interface(), solution_interface, exact_solution_interface, cell_differences, quadrature_collection_interface, norm_type, weight_fun, exponent);
		}
		else
		{
			Assert(component_mask_interface.size() == component_names_interface.size(), ExcMessage("The component mask has not the correct size!"));
			const unsigned int first_component = component_mask_interface.first_selected_component();
			const unsigned int n_components = component_mask_interface.n_selected_components();
			Assert(n_components > 0, ExcMessage("There is no component selected by the component mask!"));
			const unsigned int last_component = first_component + n_components - 1;
			for(unsigned int component = first_component; component <= last_component; ++component)
			{
				if(component_mask_interface[component] != true)
				{
					Assert(false, ExcMessage("Currently there is only support for component masks where the selected components form a contiguous range!"));
				}
			}
			const ComponentSelectFunction<spacedim, double> weight_fun(pair<unsigned int, unsigned int>(first_component, last_component+1), component_mask_interface.size());
			VectorTools::integrate_difference(mapping_collection_interface, dof_handler_system.get_dof_handler_interface(), solution_interface, exact_solution_interface, cell_differences, quadrature_collection_interface, norm_type, &weight_fun, exponent);
		}
		norm_interface = VectorTools::compute_global_error(tria_system.get_triangulation_interface(), cell_differences, norm_type, exponent);
	}

	return make_pair(norm_domain, norm_interface);
}

template<unsigned int spacedim>
const TriangulationSystem<spacedim>&
AssemblyHelper<spacedim>::get_triangulation_system()
const
{
	return tria_system;
}

template<unsigned int spacedim>
TriangulationSystem<spacedim>&
AssemblyHelper<spacedim>::get_triangulation_system()
{
	return tria_system;
}

template<unsigned int spacedim>
const DoFHandlerSystem<spacedim>&
AssemblyHelper<spacedim>::get_dof_handler_system()
const
{
	return dof_handler_system;
}

template<unsigned int spacedim>
DoFHandlerSystem<spacedim>&
AssemblyHelper<spacedim>::get_dof_handler_system()
{
	return dof_handler_system;
}

template<unsigned int spacedim>
map<const IndependentField<spacedim,spacedim>*, const unsigned int>
AssemblyHelper<spacedim>::get_u_omega_global_component_indices()
const
{
	map<const IndependentField<spacedim,spacedim>*, const unsigned int> return_;
	for(const auto& component_index : global_component_indices_u_omega)
		return_.emplace(component_index.first, component_index.second);
	return return_;
}

template<unsigned int spacedim>
map<const IndependentField<spacedim-1,spacedim>*, const unsigned int>
AssemblyHelper<spacedim>::get_u_sigma_global_component_indices()
const
{
	map<const IndependentField<spacedim-1,spacedim>*, const unsigned int> return_;
	for(const auto& component_index : global_component_indices_u_sigma)
		return_.emplace(component_index.first, component_index.second);
	return return_;
}

template<unsigned int spacedim>
unsigned int
AssemblyHelper<spacedim>::get_u_omega_global_component_index(const IndependentField<spacedim, spacedim>& u_omega)
const
{
	return global_component_indices_u_omega.at(&u_omega);
}


template<unsigned int spacedim>
unsigned int
AssemblyHelper<spacedim>::get_u_sigma_global_component_index(const IndependentField<spacedim-1, spacedim>& u_sigma)
const
{
	return global_component_indices_u_sigma.at(&u_sigma);
}

template<unsigned int spacedim>
unsigned int
AssemblyHelper<spacedim>::system_size()
const
{
	return dof_handler_system.n_dofs() + C.size() + n_scalar_functionals_nonprimitive;
}

template<unsigned int spacedim>
unsigned int
AssemblyHelper<spacedim>::get_n_stretched_rows()
const
{
	return n_scalar_functionals_nonprimitive + C.size();
}

template<unsigned int spacedim>
unsigned int
AssemblyHelper<spacedim>::get_n_C()
const
{
	return C.size();
}

template<unsigned int spacedim>
unsigned int
AssemblyHelper<spacedim>::get_global_dof_index_C(const IndependentField<0,spacedim>* independent_scalar)
const
{
	return dof_handler_system.get_dof_index(global_indices_C.at(independent_scalar));
}

template<unsigned int spacedim>
const IndexSet
AssemblyHelper<spacedim>::get_locally_relevant_indices()
const
{
	const unsigned int n_system = system_size();
	IndexSet locally_relevant_indices(n_system);
	locally_relevant_indices.add_indices(dof_handler_system.get_locally_relevant_dofs());
	locally_relevant_indices.add_range(n_system - get_n_stretched_rows(), n_system);
	return locally_relevant_indices;
}

template<unsigned int spacedim>
const IndexSet
AssemblyHelper<spacedim>::get_locally_owned_indices()
const
{
	const unsigned int n_system = system_size();
	IndexSet locally_owned_indices(n_system);
	locally_owned_indices.add_indices(dof_handler_system.get_locally_owned_dofs());
	if( this_proc == (n_procs - 1) )
	locally_owned_indices.add_range(n_system - get_n_stretched_rows(), n_system);
	return locally_owned_indices;
}

template<unsigned int spacedim>
const vector<IndexSet>
AssemblyHelper<spacedim>::get_locally_relevant_indices_blocks()
const
{
	const unsigned int n_system = system_size();
	const unsigned int first_block_size = dof_handler_system.n_dofs_domain() + dof_handler_system.n_dofs_interface();
	auto locally_relevant_indices = get_locally_relevant_indices();
	return {locally_relevant_indices.get_view(0, first_block_size), locally_relevant_indices.get_view(first_block_size, n_system)};
}

template<unsigned int spacedim>
unsigned int
AssemblyHelper<spacedim>::get_dof_index_at_point_omega(	const IndependentField<spacedim, spacedim>* u_omega,
														const unsigned int							component,
														const Point<spacedim>						p)
const
{


	//support point in real space, unit support point of domain cell
	Point<spacedim> support_point, unit_support_point_domain;

	for(const auto& domain_cell : dof_handler_system.domain_active_iterators())
	{
		if(!domain_cell->is_artificial())
		{
			for(const auto& u_omega_n : this->u_omega)
			{
				if(u_omega_n == u_omega)
				{
					const unsigned int fe_system_id = material_id_to_fe_system_id_domain.at(domain_cell->material_id());

					const unsigned int global_component_index = global_component_indices_u_omega.at(u_omega_n) + component;
					const auto& shapefuns = components_to_shapefuns_domain[fe_system_id][global_component_index];
					for(const auto& shapefun : shapefuns)
					{
						unit_support_point_domain = domain_cell->get_fe().unit_support_point(shapefun);
						support_point = mapping_domain->transform_unit_to_real_cell(domain_cell, unit_support_point_domain);
						if(support_point.distance(p) < 1e-12)
						{
							vector<unsigned int> dof_indices_local_global(domain_cell->get_fe().dofs_per_cell);
							domain_cell.get_dof_indices(dof_indices_local_global);
							return dof_indices_local_global[shapefun];
						}
					}
				}
			}
		}
	}

	return numbers::invalid_dof_index;

}

template<unsigned int spacedim>
unsigned int
AssemblyHelper<spacedim>::get_dof_index_at_point_sigma(	const IndependentField<spacedim-1, spacedim>*	u_sigma,
														const unsigned int								component,
														const Point<spacedim>							p)
const
{


	//support point in real space, unit support point of domain cell
	Point<spacedim> support_point;
	Point<spacedim-1> unit_support_point_interface;

	for(const auto& interface_cell_domain_cells : dof_handler_system.interface_active_iterators())
	{
		const auto& interface_cell = interface_cell_domain_cells.interface_cell;
		if(!interface_cell->is_artificial())
		{
			for(const auto& u_sigma_n : this->u_sigma)
			{
				if(u_sigma_n == u_sigma)
				{
					const unsigned int fe_system_id = material_id_to_fe_system_id_interface.at(interface_cell->material_id());

					const unsigned int global_component_index = global_component_indices_u_sigma.at(u_sigma_n) + component;
					const auto& shapefuns = components_to_shapefuns_interface[fe_system_id][global_component_index];
					for(const auto& shapefun : shapefuns)
					{
						unit_support_point_interface = interface_cell->get_fe().unit_support_point(shapefun);
						support_point = mapping_interface->transform_unit_to_real_cell(interface_cell, unit_support_point_interface);
						if(support_point.distance(p) < 1e-12)
						{
							vector<unsigned int> dof_indices_local_global(interface_cell->get_fe().dofs_per_cell);
							interface_cell.get_dof_indices(dof_indices_local_global);
							return dof_indices_local_global[shapefun];
						}
					}
				}
			}
		}
	}

	return numbers::invalid_dof_index;

}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::print_dof_information(const unsigned int dof_index)
const
{
	//support point in real space, unit support point of domain cell
	Point<spacedim> support_point, unit_support_point_domain;
	vector<unsigned int> dof_indices_local_global;

	for(const auto& domain_cell : dof_handler_system.domain_active_iterators())
	{
		if(!domain_cell->is_artificial())
		{
			dof_indices_local_global.resize(domain_cell->get_fe().dofs_per_cell);
			domain_cell.get_dof_indices(dof_indices_local_global);
			for(const auto& u_omega_n : this->u_omega)
			{
				const unsigned int fe_system_id = material_id_to_fe_system_id_domain.at(domain_cell->material_id());

				for(unsigned int component = 0; component < u_omega_n->n_components; ++component)
				{
					const unsigned int global_component_index = global_component_indices_u_omega.at(u_omega_n) + component;
					const auto& shapefuns = components_to_shapefuns_domain[fe_system_id][global_component_index];
					for(const auto& shapefun : shapefuns)
					{
						if(dof_indices_local_global[shapefun] == dof_index)
						{
							pout << dof_indices_local_global[shapefun] << " belongs to independent field " << u_omega_n->name << endl;
							pout << dof_indices_local_global[shapefun] << " belongs to component " << component << endl;
							pout << dof_indices_local_global[shapefun] << " belongs to cell " << domain_cell->id() << endl;
							if(domain_cell->get_fe().has_support_points())
							{
								unit_support_point_domain = domain_cell->get_fe().unit_support_point(shapefun);
								support_point = mapping_domain->transform_unit_to_real_cell(domain_cell, unit_support_point_domain);
								pout << dof_indices_local_global[shapefun] << " has support at " << support_point << endl;
							}
							return;
						}
					}
				}
			}
		}
	}
	pout << "Did not find the requested index!" << endl;

}

template<unsigned int spacedim>
const vector<IndexSet>
AssemblyHelper<spacedim>::get_locally_owned_indices_blocks()
const
{
	const unsigned int n_system = system_size();
	const unsigned int first_block_size = dof_handler_system.n_dofs_domain() + dof_handler_system.n_dofs_interface();
	auto locally_owned_indices = get_locally_owned_indices();
	return {locally_owned_indices.get_view(0, first_block_size), locally_owned_indices.get_view(first_block_size, n_system)};
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::distribute_dofs()
{
	for(const auto& cell : dof_handler_system.domain_active_iterators())
	{
		if(cell->is_locally_owned())
		{
			Assert(material_id_to_fe_system_id_domain.find(cell->material_id()) != material_id_to_fe_system_id_domain.end(), ExcMessage("Internal error!"));
			cell->set_active_fe_index(material_id_to_fe_system_id_domain[cell->material_id()]);
		}
	}
	for(const auto& interface_cell_domain_cell : dof_handler_system.interface_active_iterators())
	{
		if(interface_cell_domain_cell.interface_cell->is_locally_owned())
		{
			Assert(material_id_to_fe_system_id_interface.find(interface_cell_domain_cell.interface_cell->material_id()) != material_id_to_fe_system_id_interface.end(), ExcMessage("Internal error!"));
			interface_cell_domain_cell.interface_cell->set_active_fe_index(material_id_to_fe_system_id_interface[interface_cell_domain_cell.interface_cell->material_id()]);
		}
	}
	dof_handler_system.distribute_dofs(fe_collection_domain, fe_collection_interface, C.size());
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::initialize_fe_values_domain(	const typename hp::DoFHandler<spacedim, spacedim>::active_cell_iterator&	cell,
														const unsigned int															internal_index,
														const bool																	nonprimitive)
const
{
	if(nonprimitive)
		for(auto fe_values_n : fe_values_domain_reinit_nonprimitive[internal_index])
			fe_values_n->reinit(cell);
	else
		for(auto fe_values_n : fe_values_domain_reinit[internal_index])
			fe_values_n->reinit(cell);
	return;
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::initialize_fe_values_interface(	const InterfaceCellDomainCellsDoF<spacedim>&	interface_cell_domain_cells,
															const unsigned int								internal_index,
															const bool										nonprimitive)
const
{
	if(nonprimitive)
		for(auto fe_values_n : fe_values_interface_reinit_nonprimitive[internal_index])
			fe_values_n->reinit(interface_cell_domain_cells);
	else
		for(auto fe_values_n : fe_values_interface_reinit[internal_index])
			fe_values_n->reinit(interface_cell_domain_cells);
	return;
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::compute_e_omega(	const unsigned int		internal_index,
											const unsigned int		scalar_functional_index,
											const unsigned int		q_point,
											const Vector<double>&	solution_u_omega,
											const Vector<double>&	solution_C,
											Vector<double>&			e_omega,
											FullMatrix<double>&		de_omega_dsol_T,
											const bool				compute_derivative,
											const bool				ignore_constants)
const
{
	//@todo This function is one of the main bottlenecks during assembly -> think about optimization!

	e_omega = 0.0;
	const auto& coupled_dof_indices_scalar_functionals_domain = this->coupled_dof_indices_scalar_functionals_domain[internal_index][scalar_functional_index];
	const auto& coupled_C_indices_scalar_functionals_domain = this->coupled_C_indices_scalar_functionals_domain[internal_index][scalar_functional_index];

	const unsigned int dof_indices_cell_domain = coupled_dof_indices_scalar_functionals_domain.size();
	if(compute_derivative)
	{
		const unsigned int n_dofs = dof_indices_cell_domain + coupled_C_indices_scalar_functionals_domain.size();
		if(n_dofs>0)
			de_omega_dsol_T.reinit(n_dofs, e_omega.size());
		else
			de_omega_dsol_T.reinit(0,0);
	}

	const FEValues<spacedim,spacedim>* fe_values_domain = this->fe_values_domain[internal_index][scalar_functional_index].get();
	double shapefun_derivative;
	for(unsigned int e_omega_n = 0; e_omega_n < e_omega.size(); ++e_omega_n)
	{
		for(const auto& a_omega_n : a_omega[internal_index][scalar_functional_index][e_omega_n])
		{
			const double a = get<0>(a_omega_n);
			const unsigned int component = get<1>(a_omega_n);
			const auto& shapefuns = get<2>(a_omega_n);
			for(const auto& shapefun : shapefuns)
			{
				shapefun_derivative = a * fe_values_domain->shape_value_component(coupled_dof_indices_scalar_functionals_domain[shapefun], q_point, component);
				e_omega[e_omega_n] += shapefun_derivative * solution_u_omega[shapefun];
				if(compute_derivative)
					de_omega_dsol_T(shapefun, e_omega_n) += shapefun_derivative;
			}
		}

		for(const auto& b_omega_n : b_omega[internal_index][scalar_functional_index][e_omega_n])
		{
			const double b = get<0>(b_omega_n);
			const unsigned int component = get<1>(b_omega_n);
			const unsigned int derivative = get<2>(b_omega_n);
			const auto& shapefuns = get<3>(b_omega_n);
			for(const auto& shapefun : shapefuns)
			{
				shapefun_derivative = b * fe_values_domain->shape_grad_component(coupled_dof_indices_scalar_functionals_domain[shapefun], q_point, component)[derivative];
				e_omega[e_omega_n] += shapefun_derivative * solution_u_omega[shapefun];
				if(compute_derivative)
					de_omega_dsol_T(shapefun, e_omega_n) += shapefun_derivative;
			}
		}

		for(const auto& c_omega_n : this->c_omega[internal_index][scalar_functional_index][e_omega_n])
		{
			e_omega[e_omega_n] += get<0>(c_omega_n) * solution_C[get<1>(c_omega_n)];
			if(compute_derivative)
				de_omega_dsol_T( get<1>(c_omega_n) + dof_indices_cell_domain, e_omega_n ) += get<0>(c_omega_n);
		}

		if(!ignore_constants)
			e_omega[e_omega_n] += this->d_omega[internal_index][scalar_functional_index][e_omega_n];

	}
}

template<unsigned int spacedim>
void
AssemblyHelper<spacedim>::compute_e_sigma(	const unsigned int			internal_index,
											const unsigned int			scalar_functional_index,
											const unsigned int			q_point,
											const Vector<double>&		solution_u_sigma,
											const Vector<double>&		solution_u_omega_minus,
											const Vector<double>&		solution_u_omega_plus,
											const Vector<double>&		solution_C,
											const vector<unsigned int>&	dof_indices_interface_dof_indices_combined,
											const vector<unsigned int>& dof_indices_minus_dof_indices_combined,
											const vector<unsigned int>& dof_indices_plus_dof_indices_combined,
											const vector<unsigned int>& dof_indices_C_dof_indices_combined,
											const vector<unsigned int>& dof_indices_global_combined,
											Vector<double>&				e_sigma,
											FullMatrix<double>&			de_sigma_dsol_T,
											const bool					compute_derivative,
											const bool					ignore_constants)
const
{
	//Attention: This function is one of the main bottlenecks during assembly -> think about optimization!

	e_sigma=0.0;
	const auto& coupled_dof_indices_scalar_functionals_interface = this->coupled_dof_indices_scalar_functionals_interface[internal_index][scalar_functional_index];
	const auto& coupled_dof_indices_scalar_functionals_interface_minus = this->coupled_dof_indices_scalar_functionals_interface_minus[internal_index][scalar_functional_index];
	const auto& coupled_dof_indices_scalar_functionals_interface_plus = this->coupled_dof_indices_scalar_functionals_interface_plus[internal_index][scalar_functional_index];

	if(compute_derivative)
	{
		const unsigned int n_dofs = dof_indices_global_combined.size();
		if(n_dofs > 0)
			de_sigma_dsol_T.reinit(n_dofs, e_sigma.size());
		else
			de_sigma_dsol_T.reinit(0,0);
	}
	const FEValues<spacedim-1,spacedim>& fe_values_interface = this->fe_values_interface[internal_index][scalar_functional_index]->get_fe_values_interface();
	const FEFaceValuesBase<spacedim,spacedim>& fe_values_interface_minus = this->fe_values_interface[internal_index][scalar_functional_index]->get_fe_values_domain(InterfaceSide::minus);
	const FEFaceValuesBase<spacedim,spacedim>& fe_values_interface_plus = this->fe_values_interface[internal_index][scalar_functional_index]->get_fe_values_domain(InterfaceSide::plus);

	double shapefun_derivative;
	for(unsigned int e_sigma_n=0; e_sigma_n<e_sigma.size(); ++e_sigma_n)
	{
		for(const auto& a_sigma_n : this->a_sigma[internal_index][scalar_functional_index][e_sigma_n])
		{
			const double a = get<0>(a_sigma_n);
			const unsigned int component = get<1>(a_sigma_n);
			const auto& shapefuns = get<2>(a_sigma_n);
			for(const auto& shapefun : shapefuns)
			{
				shapefun_derivative = a * fe_values_interface.shape_value_component(coupled_dof_indices_scalar_functionals_interface[shapefun], q_point, component);
				e_sigma[e_sigma_n] += shapefun_derivative * solution_u_sigma[shapefun];
				if(compute_derivative)
					de_sigma_dsol_T(dof_indices_interface_dof_indices_combined[shapefun], e_sigma_n) += shapefun_derivative;

			}
		}

		for(const auto& b_sigma_n : this->b_sigma[internal_index][scalar_functional_index][e_sigma_n])
		{
			const double b = get<0>(b_sigma_n);
			const unsigned int component = get<1>(b_sigma_n);
			const unsigned int derivative = get<2>(b_sigma_n);
			const auto& shapefuns = get<3>(b_sigma_n);
			for(const auto& shapefun : shapefuns)
			{
				shapefun_derivative = b * fe_values_interface.shape_grad_component(coupled_dof_indices_scalar_functionals_interface[shapefun], q_point, component)[derivative];
				e_sigma[e_sigma_n] += shapefun_derivative * solution_u_sigma[shapefun];
				if(compute_derivative)
					de_sigma_dsol_T(dof_indices_interface_dof_indices_combined[shapefun], e_sigma_n) += shapefun_derivative;
			}
		}

		for(const auto& c_n : this->c_sigma[internal_index][scalar_functional_index][e_sigma_n])
		{
			e_sigma[e_sigma_n] += get<0>(c_n) * solution_C[get<1>(c_n)];
			if(compute_derivative)
				de_sigma_dsol_T( dof_indices_C_dof_indices_combined[get<1>(c_n)], e_sigma_n ) += get<0>(c_n);
		}

		if(!ignore_constants)
			e_sigma[e_sigma_n] += this->d_sigma[internal_index][scalar_functional_index][e_sigma_n];
	}

	for(unsigned int e_sigma_n = 0; e_sigma_n < e_sigma.size(); ++e_sigma_n)
	{
		for(const auto& a_minus_n : this->a_minus[internal_index][scalar_functional_index][e_sigma_n])
		{
			const double a = get<0>(a_minus_n);
			const unsigned int component = get<1>(a_minus_n);
			const auto& shapefuns = get<2>(a_minus_n);
			for(const auto& shapefun : shapefuns)
			{
				shapefun_derivative = a * fe_values_interface_minus.shape_value_component(coupled_dof_indices_scalar_functionals_interface_minus[shapefun], q_point, component);
				e_sigma[e_sigma_n] += shapefun_derivative * solution_u_omega_minus[shapefun];
				if(compute_derivative)
					de_sigma_dsol_T(dof_indices_minus_dof_indices_combined[shapefun], e_sigma_n) += shapefun_derivative;

			}
		}

		for(const auto& b_minus_n : this->b_minus[internal_index][scalar_functional_index][e_sigma_n])
		{
			const double b = get<0>(b_minus_n);
			const unsigned int component = get<1>(b_minus_n);
			const unsigned int derivative = get<2>(b_minus_n);
			const auto& shapefuns = get<3>(b_minus_n);
			for(const auto& shapefun : shapefuns)
			{
				shapefun_derivative = b * fe_values_interface_minus.shape_grad_component(coupled_dof_indices_scalar_functionals_interface_minus[shapefun], q_point, component)[derivative];
				e_sigma[e_sigma_n] += shapefun_derivative * solution_u_omega_minus[shapefun];
				if(compute_derivative)
					de_sigma_dsol_T(dof_indices_minus_dof_indices_combined[shapefun], e_sigma_n) += shapefun_derivative;
			}
		}
	}

	for(unsigned int e_sigma_n = 0; e_sigma_n < e_sigma.size(); ++e_sigma_n)
	{
		for(const auto& a_plus_n : this->a_plus[internal_index][scalar_functional_index][e_sigma_n])
		{
			const double a = get<0>(a_plus_n);
			const unsigned int component = get<1>(a_plus_n);
			const auto& shapefuns = get<2>(a_plus_n);
			for(const auto& shapefun : shapefuns)
			{
				shapefun_derivative = a * fe_values_interface_plus.shape_value_component(coupled_dof_indices_scalar_functionals_interface_plus[shapefun], q_point, component);
				e_sigma[e_sigma_n] += shapefun_derivative * solution_u_omega_plus[shapefun];
				if(compute_derivative)
					de_sigma_dsol_T(dof_indices_plus_dof_indices_combined[shapefun], e_sigma_n) += shapefun_derivative;
			}
		}

		for(const auto& b_plus_n : this->b_plus[internal_index][scalar_functional_index][e_sigma_n])
		{
			const double b = get<0>(b_plus_n);
			const unsigned int component = get<1>(b_plus_n);
			const unsigned int derivative = get<2>(b_plus_n);
			const auto& shapefuns = get<3>(b_plus_n);
			for(const auto& shapefun : shapefuns)
			{
				shapefun_derivative = b * fe_values_interface_plus.shape_grad_component(coupled_dof_indices_scalar_functionals_interface_plus[shapefun], q_point, component)[derivative];
				e_sigma[e_sigma_n] += shapefun_derivative * solution_u_omega_plus[shapefun];
				if(compute_derivative)
					de_sigma_dsol_T(dof_indices_plus_dof_indices_combined[shapefun], e_sigma_n) += shapefun_derivative;
			}
		}
	}
}


template<unsigned int spacedim>
template<class VectorType>
bool
AssemblyHelper<spacedim>::get_nonprimitive_scalar_functional_values(const VectorType& 				solution,
																	const vector<const VectorType*> solution_ref_sets,
																	Vector<double>&					nonprimitive_scalar_functional_values)
const
{

	//the return value, false if everything is ok
	bool error = false;

	//initialize nonprimitive_scalar_functional_values to zero
	nonprimitive_scalar_functional_values = 0.0;

	//requested quantities (only values of scalar functionals, no derivatives)
	const auto requested_quantities = make_tuple(true, false, false);

	///////////////////////////////////////////////
	// first add contributions related to domain //
	///////////////////////////////////////////////

	//stores the material_id of the cell visited previously
	types::material_id material_id_previous=numbers::invalid_material_id;

	//holds the internal index of a cell
	unsigned int internal_index = 0;

	//global dof indices coupling on the cell for a particular scalar functional
	vector<unsigned int> dof_indices_global;
	dof_indices_global.reserve(fe_collection_domain.max_dofs_per_cell() + global_indices_C.size());

	//mapping between local dof indices and global ones (dof_indices_local_global[i] is the global dof index corresponding to the local dof index i)
	vector<unsigned int> dof_indices_local_global;
	dof_indices_local_global.reserve(fe_collection_domain.max_dofs_per_cell());

	//global dof indices of C's
	vector<unsigned int> dof_indices_local_global_C;
	dof_handler_system.get_dof_indices(dof_indices_local_global_C);

	//global dof indices of C's for a particular scalar functional
	vector<unsigned int> dof_indices_global_C;
	dof_indices_global_C.reserve(global_indices_C.size());

	//solution vector restricted to dof indices coupling on a cell for a particular scalar functional
	Vector<double> solution_local;
	solution_local.reinit(fe_collection_domain.max_dofs_per_cell());

	//solution vector restricted to C indices coupling on a cell for a particular scalar functional
	Vector<double> solution_local_C;
	solution_local_C.reinit(global_indices_C.size());

	//reference solution vector restricted to dof indices coupling on a cell for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_n : solution_ref_sets_local)
		solution_ref_sets_local_n.reinit(fe_collection_domain.max_dofs_per_cell());

	//reference solution vector restricted to C indices coupling on a cell for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local_C(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_C_n : solution_ref_sets_local_C)
		solution_ref_sets_local_C_n.reinit(global_indices_C.size());

	//vectors for dependent variables
	Vector<double> e_omega(total_potential.max_dependent_vars);
	vector<Vector<double>> e_omega_ref_sets(solution_ref_sets.size());

	//matrix for derivatives of dependent variables w.r.t. local dof's (not needed in this function - just a dummy)
	FullMatrix<double> de_omega_dsol_T;

	//value of a scalar functional
	double h_omega;

	//first derivative of a scalar functional (not needed in this function - just a dummy)
	Vector<double> h_omega_1;

	//second derivative of a scalar functional (not needed in this function - just a dummy)
	FullMatrix<double> h_omega_2;

	//the actual loop over the cells
	for(const auto& domain_cell : dof_handler_system.domain_active_iterators())
	{
		if(domain_cell->is_locally_owned())
		{
			//the internal material index needs only to be updated if the cell has another material_id than the one visited before
			if(domain_cell->material_id() != material_id_previous)
			{
				material_id_previous = domain_cell->material_id();
				internal_index = material_id_to_internal_index_domain.at(material_id_previous);
			}

			//initialize relevant FEValues objects with cell
			initialize_fe_values_domain(domain_cell, internal_index, true);

			//get the mapping between local and global dof indices
			dof_indices_local_global.resize(domain_cell->get_fe().dofs_per_cell);
			domain_cell.get_dof_indices(dof_indices_local_global);

			//loop scalar functionals on cell which enter the total potential non-primitively
			for(unsigned int scalar_functional_nonprimitive_n = 0; scalar_functional_nonprimitive_n < scalar_functionals_domain_nonprimitive[internal_index].size(); ++scalar_functional_nonprimitive_n)
			{

				//get the actual index of the scalar functional
				const unsigned int scalar_functional_n = scalar_functionals_domain_nonprimitive[internal_index][scalar_functional_nonprimitive_n];

				//this contains the local dof indices coupling on the domain cell for the scalar functional under consideration
				const auto& dof_indices_local = coupled_dof_indices_scalar_functionals_domain[internal_index][scalar_functional_n];

				//this contains the local C indices coupling on the domain cell for the scalar functional under consideration
				const auto& dof_indices_local_C = coupled_C_indices_scalar_functionals_domain[internal_index][scalar_functional_n];

				//nothing to do if there are no dof's present
				if((dof_indices_local.size() + dof_indices_local_C.size()) == 0)
					continue;

				//map the local dof indices to the global ones
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local, dof_indices_global, dof_indices_local_global);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_C, dof_indices_global_C, dof_indices_local_global_C);

				//restrict solution vectors to scalar functional on cell
				solution.extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_local.begin());
				solution.extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_local_C.begin());
				for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
				{
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_ref_sets_local[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_ref_sets_local_C[ref_set].begin());
				}

				//get a pointer to the scalar functional under consideration and its nonprimitive_index
				const auto scalar_functional = scalar_functionals_domain[internal_index][scalar_functional_n];
				const int nonprimitive_index = scalar_functionals_domain_nonprimitive_indices.at(scalar_functional);

				//initialize dependent variable vectors appropriately
				e_omega.reinit(scalar_functional->e_omega.size());
				for(auto& e_omega_ref_n : e_omega_ref_sets)
					e_omega_ref_n.reinit(e_omega.size());

				//get JxW values and quadrature points
				const auto& JxW = fe_values_domain[internal_index][scalar_functional_n]->get_JxW_values();
				const auto& q_points = fe_values_domain[internal_index][scalar_functional_n]->get_quadrature_points();

				//loop over integration points
				for(unsigned int q_point = 0; q_point < scalar_functional->quadrature.size(); ++q_point)
				{
					//compute dependent variables
					compute_e_omega(internal_index, scalar_functional_n, q_point, solution_local, solution_local_C, e_omega, de_omega_dsol_T, false);

					//compute reference values of dependent variables
					for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
						compute_e_omega(internal_index, scalar_functional_n, q_point, solution_ref_sets_local[ref_set], solution_ref_sets_local_C[ref_set], e_omega_ref_sets[ref_set], de_omega_dsol_T, false);

					//hidden variables
					Assert(domain_cell->user_pointer() != nullptr, ExcMessage("Internal error - it seems that the hidden variables are not properly allocated!"));
					Vector<double>& hidden_vars = (*static_cast<vector<vector<Vector<double>>>*>(domain_cell->user_pointer()))[scalar_functional_n][q_point];

					//the location of the quadrature point in real space
					const auto& x = q_points[q_point];

					//evaluate the integrand of the scalar functional
					if(scalar_functional->get_h_omega(e_omega, e_omega_ref_sets, hidden_vars, x, h_omega, h_omega_1, h_omega_2, requested_quantities))
						error = true;

					//add contribution to the value of the scalar functional
					nonprimitive_scalar_functional_values[nonprimitive_index] += h_omega * JxW[q_point];
				}
			}
		}
	}

	////////////////////////////////////////////////
	// now add contributions related to interface //
	////////////////////////////////////////////////

	//stores the material_ids of the cell visited previously
	tuple<types::material_id, types::material_id, types::material_id> material_ids_previous = make_tuple(	numbers::invalid_material_id,
																											numbers::invalid_material_id,
																											numbers::invalid_material_id);
	//global dof indices coupling on the interface cell for a particular scalar functional
	dof_indices_global.reserve(fe_collection_interface.max_dofs_per_cell());

	//global dof indices coupling on - side for a particular scalar functional
	vector<unsigned int> dof_indices_global_minus;
	dof_indices_global_minus.reserve(fe_collection_domain.max_dofs_per_cell());

	//global dof indices coupling on + side for a particular scalar functional
	vector<unsigned int> dof_indices_global_plus;
	dof_indices_global_plus.reserve(fe_collection_domain.max_dofs_per_cell());

	//mapping between local dof indices and global ones on interface
	dof_indices_local_global.reserve(fe_collection_interface.max_dofs_per_cell());

	//mapping between local dof indices and global ones on - side
	vector<unsigned int> dof_indices_local_global_minus;
	dof_indices_local_global_minus.reserve(fe_collection_domain.max_dofs_per_cell());

	//mapping between local dof indices and global ones on + side
	vector<unsigned int> dof_indices_local_global_plus;
	dof_indices_local_global_plus.reserve(fe_collection_domain.max_dofs_per_cell());

	//combined global dof indices coupling on the interface cell and the adjacent domain cells for a particular scalar functional
	//(without any duplicate dof's -> this quantity is NOT just a concatenation of dof_indices_global, dof_indices_global_minus, dof_indices_global_plus)
	//just a dummy vector in this function
	vector<unsigned int> dof_indices_global_combined;

	//defined such that dof_indices_global_combined[dof_indices_interface_dof_indices_combined[i]] = dof_indices_global[i]
	//just a dummy vector in this function
	vector<unsigned int> dof_indices_interface_dof_indices_combined;

	//defined such that dof_indices_global_combined[dof_indices_minus_dof_indices_combined[i]] = dof_indices_global_minus[i]
	//just a dummy vector in this function
	vector<unsigned int> dof_indices_minus_dof_indices_combined;

	//defined such that dof_indices_global_combined[dof_indices_plus_dof_indices_combined[i]] = dof_indices_global_plus[i]
	//just a dummy vector in this function
	vector<unsigned int> dof_indices_plus_dof_indices_combined;

	//defined such that dof_indices_global_combined[dof_indices_C_dof_indices_combined[i]] = dof_indices_global_C[i]
	//just a dummy vector in this function
	vector<unsigned int> dof_indices_C_dof_indices_combined;

	//solution vector restricted to dof indices coupling on an interface cell for a particular scalar functional
	solution_local.reinit(fe_collection_interface.max_dofs_per_cell());

	//reference solution vector restricted to dof indices coupling on an interface cell for a particular scalar functional
	for(auto& solution_ref_sets_local_n : solution_ref_sets_local)
		solution_ref_sets_local_n.reinit(fe_collection_interface.max_dofs_per_cell());

	//solution vector restricted to dof indices coupling on a domain cell located on the - side of the interface for a particular scalar functional
	Vector<double> solution_local_minus;
	solution_local_minus.reinit(fe_collection_domain.max_dofs_per_cell());

	//reference solution vector restricted to dof indices coupling on a domain cell located on the - side of the interface for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local_minus(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_minus_n : solution_ref_sets_local_minus)
		solution_ref_sets_local_minus_n.reinit(fe_collection_domain.max_dofs_per_cell());

	//solution vector restricted to dof indices coupling on a domain cell located on the + side of the interface for a particular scalar functional
	Vector<double> solution_local_plus;
	solution_local_plus.reinit(fe_collection_domain.max_dofs_per_cell());

	//reference solution vector restricted to dof indices coupling on a domain cell located on the + side of the interface for a particular scalar functional
	vector<Vector<double>> solution_ref_sets_local_plus(solution_ref_sets.size());
	for(auto& solution_ref_sets_local_plus_n : solution_ref_sets_local_plus)
		solution_ref_sets_local_plus_n.reinit(fe_collection_domain.max_dofs_per_cell());

	//vectors for dependent variables
	Vector<double> e_sigma(total_potential.max_dependent_vars);
	vector<Vector<double>> e_sigma_ref_sets(solution_ref_sets.size());

	//matrix for derivatives of dependent variables w.r.t. local dof's (not needed in this function - just a dummy)
	FullMatrix<double> de_sigma_dsol_T;

	//value of a scalar functional
	double h_sigma;

	//first derivative of a scalar functional (not needed in this function - just a dummy)
	Vector<double> h_sigma_1;

	//second derivative of a scalar functional (not needed in this function - just a dummy)
	FullMatrix<double> h_sigma_2;

	//the actual loop over the cells
	for(const auto& interface_cell_domain_cells : dof_handler_system.interface_active_iterators())
	{
		if(interface_cell_domain_cells.interface_cell->is_locally_owned())
		{
			//find out the material id's characterizing the interface
			const auto material_ids = interface_cell_domain_cells.get_material_ids();

			//the internal material index needs only to be updated if the cell is associated with different material_id's than the one visited before
			if(material_ids != material_ids_previous)
			{
				internal_index = material_ids_to_internal_index_interface.at(material_ids);
				material_ids_previous = material_ids;
			}

			//initialize FEValues objects with cell
			initialize_fe_values_interface(interface_cell_domain_cells, internal_index, true);

			//get the mapping between local and global dof indices
			interface_cell_domain_cells.get_dof_indices_local_global_interface(dof_indices_local_global, dof_indices_local_global_minus, dof_indices_local_global_plus);

			//loop scalar functionals on cell which enter the total potential non-primitively
			for(unsigned int scalar_functional_nonprimitive_n = 0; scalar_functional_nonprimitive_n < scalar_functionals_interface_nonprimitive[internal_index].size(); ++scalar_functional_nonprimitive_n)
			{
				//get the actual index of the scalar functional
				const unsigned int scalar_functional_n = scalar_functionals_interface_nonprimitive[internal_index][scalar_functional_nonprimitive_n];

				//this contains the local dof indices coupling on the interface cell for the scalar functional under consideration
				const auto& dof_indices_local=coupled_dof_indices_scalar_functionals_interface[internal_index][scalar_functional_n];

				//this contains the local dof indices coupling on the - side for the scalar functional under consideration
				const auto& dof_indices_local_minus=coupled_dof_indices_scalar_functionals_interface_minus[internal_index][scalar_functional_n];

				//this contains the local dof indices coupling on the + side for the scalar functional under consideration
				const auto& dof_indices_local_plus=coupled_dof_indices_scalar_functionals_interface_plus[internal_index][scalar_functional_n];

				//this contains the local C indices coupling on the interface cell for the scalar functional under consideration
				const auto& dof_indices_local_C=coupled_C_indices_scalar_functionals_interface[internal_index][scalar_functional_n];

				//nothing to do if there are no dof's present
				unsigned int n_dofs_cell = dof_indices_local.size() + dof_indices_local_minus.size() + dof_indices_local_plus.size() + dof_indices_local_C.size();
				if(n_dofs_cell == 0)
					continue;

				//map the local dof indices to the global ones
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local, dof_indices_global, dof_indices_local_global);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_minus, dof_indices_global_minus, dof_indices_local_global_minus);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_plus, dof_indices_global_plus, dof_indices_local_global_plus);
				Auxiliary::convert_local_indices_to_global_indices(dof_indices_local_C, dof_indices_global_C, dof_indices_local_global_C);

				//restrict solution vectors to scalar functional on cell
				solution.extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_local.begin());
				solution.extract_subvector_to(dof_indices_global_minus.begin(), dof_indices_global_minus.end(), solution_local_minus.begin());
				solution.extract_subvector_to(dof_indices_global_plus.begin(), dof_indices_global_plus.end(), solution_local_plus.begin());
				solution.extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_local_C.begin());
				for(unsigned int ref_set=0; ref_set<solution_ref_sets.size(); ++ref_set)
				{
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global.begin(), dof_indices_global.end(), solution_ref_sets_local[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_minus.begin(), dof_indices_global_minus.end(), solution_ref_sets_local_minus[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_plus.begin(), dof_indices_global_plus.end(), solution_ref_sets_local_plus[ref_set].begin());
					solution_ref_sets[ref_set]->extract_subvector_to(dof_indices_global_C.begin(), dof_indices_global_C.end(), solution_ref_sets_local_C[ref_set].begin());
				}

				//get a pointer to the scalar functional under consideration and its nonprimitive_index
				const auto scalar_functional = scalar_functionals_interface[internal_index][scalar_functional_n];
				const int nonprimitive_index = scalar_functionals_interface_nonprimitive_indices.at(scalar_functional);

				//initialize dependent variable vectors appropriately
				e_sigma.reinit(scalar_functional->e_sigma.size());
				for(auto& e_sigma_ref_n : e_sigma_ref_sets)
					e_sigma_ref_n.reinit(e_sigma.size());

				//get JxW values and quadrature points
				const auto& JxW = fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_interface().get_JxW_values();
				const auto& q_points = fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_interface().get_quadrature_points();

				//get normal vectors at quadrature points
				const auto& normals = fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_domain(InterfaceSide::minus).get_normal_vectors();

				//loop over integration points
				for(unsigned int q_point = 0; q_point < scalar_functional->quadrature.size(); ++q_point)
				{
					//compute dependent variables
					compute_e_sigma(internal_index,
									scalar_functional_n,
									q_point,
									solution_local,
									solution_local_minus,
									solution_local_plus,
									solution_local_C,
									dof_indices_interface_dof_indices_combined,
									dof_indices_minus_dof_indices_combined,
									dof_indices_plus_dof_indices_combined,
									dof_indices_C_dof_indices_combined,
									dof_indices_global_combined,
									e_sigma,
									de_sigma_dsol_T,
									false);

					//compute reference values of dependent variables
					for(unsigned int ref_set = 0; ref_set < solution_ref_sets.size(); ++ref_set)
						compute_e_sigma(internal_index,
										scalar_functional_n,
										q_point,
										solution_ref_sets_local[ref_set],
										solution_ref_sets_local_minus[ref_set],
										solution_ref_sets_local_plus[ref_set],
										solution_ref_sets_local_C[ref_set],
										dof_indices_interface_dof_indices_combined,
										dof_indices_minus_dof_indices_combined,
										dof_indices_plus_dof_indices_combined,
										dof_indices_C_dof_indices_combined,
										dof_indices_global_combined,
										e_sigma_ref_sets[ref_set],
										de_sigma_dsol_T,
										false);

					//hidden variables
					Assert(interface_cell_domain_cells.interface_cell->user_pointer() != nullptr, ExcMessage("Internal error - it seems that the hidden variables are not propery allocated!"));
					Vector<double>& hidden_vars = (*static_cast<vector<vector<Vector<double>>>*>(interface_cell_domain_cells.interface_cell->user_pointer()))[scalar_functional_n][q_point];;

					//the location of the quadrature point in real space
					const auto& x = q_points[q_point];

					//the normal vector at the quadrature point
					const auto& n = normals[q_point];

					//evaluate the integrand of the scalar functional
					if(scalar_functional->get_h_sigma(e_sigma, e_sigma_ref_sets, hidden_vars, x, n, h_sigma, h_sigma_1, h_sigma_2, requested_quantities))
						error = true;

					//add contribution to the value of the scalar functional
					nonprimitive_scalar_functional_values[nonprimitive_index] += h_sigma * JxW[q_point];

					//check that quadrature points are aligned with each other
					//these checks can be removed in release
					Assert(	x.distance(fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_domain(InterfaceSide::minus).quadrature_point(q_point)) < 1e-8,
							ExcMessage("Internal error: Quadrature points are not aligned on interface or boundary!"));
					if(interface_cell_domain_cells.refinement_case != InterfaceRefinementCase::at_boundary)
						Assert(	x.distance(fe_values_interface[internal_index][scalar_functional_n]->get_fe_values_domain(InterfaceSide::plus).quadrature_point(q_point)) < 1e-8,
								ExcMessage("Internal error: Quadrature points are not aligned on interface or boundary!"));
				}
			}
		}
	}

	//add up contributions of different processors
	if(n_procs > 1)
	{
#ifdef DEAL_II_WITH_MPI
		const auto tria_domain_ptr = dynamic_cast<const dealii::parallel::Triangulation<spacedim, spacedim>*>(&(tria_system.get_triangulation_domain()));
		Assert(tria_domain_ptr != nullptr, ExcMessage("Internal error!"));
		int ierr = MPI_Allreduce(MPI_IN_PLACE, nonprimitive_scalar_functional_values.data(), nonprimitive_scalar_functional_values.size(), MPI_DOUBLE, MPI_SUM, tria_domain_ptr->get_communicator());
		AssertThrowMPI(ierr);

		return Auxiliary::communicate_bool(error, tria_domain_ptr->get_communicator());
#else
		Assert(n_procs == 1, ExcMessage("Internal error: deal.II not configured with MPI, but computation run with more than one processor"));
#endif //DEAL_II_WITH_MPI
	}

	return error;
}

template<unsigned int spacedim>
pair<const int, const int>
AssemblyHelper<spacedim>::get_scalar_functional_indices(const ScalarFunctional<spacedim, spacedim>* scalar_functional)
const{

	unsigned int nonprimitive_index, primitive_index;
	if(scalar_functionals_domain_nonprimitive_indices.find(scalar_functional) != scalar_functionals_domain_nonprimitive_indices.end())
		nonprimitive_index = scalar_functionals_domain_nonprimitive_indices.at(scalar_functional);
	else
		nonprimitive_index = -1;

	if(scalar_functionals_domain_primitive_indices.find(scalar_functional) != scalar_functionals_domain_primitive_indices.end())
		primitive_index = scalar_functionals_domain_primitive_indices.at(scalar_functional);
	else
		primitive_index = -1;

	return make_pair(nonprimitive_index, primitive_index);
}

template<unsigned int spacedim>
pair<const int, const int>
AssemblyHelper<spacedim>::get_scalar_functional_indices(const ScalarFunctional<spacedim-1, spacedim>* scalar_functional)
const{

	unsigned int nonprimitive_index, primitive_index;
	if(scalar_functionals_interface_nonprimitive_indices.find(scalar_functional) != scalar_functionals_interface_nonprimitive_indices.end())
		nonprimitive_index = scalar_functionals_interface_nonprimitive_indices.at(scalar_functional);
	else
		nonprimitive_index = -1;

	if(scalar_functionals_interface_primitive_indices.find(scalar_functional) != scalar_functionals_interface_primitive_indices.end())
		primitive_index = scalar_functionals_interface_primitive_indices.at(scalar_functional);
	else
		primitive_index = -1;

	return make_pair(nonprimitive_index, primitive_index);
}

//instantiations

//AssemblyHelper
template class AssemblyHelper<2>;
template class AssemblyHelper<3>;

//get_nonprimitive_scalar_functional_values
template
bool
AssemblyHelper<2>::get_nonprimitive_scalar_functional_values<Vector<double>>(	const Vector<double>&,
																				const vector<const Vector<double>*>,
																				Vector<double>&)
const;

template
bool
AssemblyHelper<3>::get_nonprimitive_scalar_functional_values<Vector<double>>(	const Vector<double>&,
																				const vector<const Vector<double>*>,
																				Vector<double>&)
const;

#ifdef DEAL_II_WITH_MPI
template
bool
AssemblyHelper<2>::get_nonprimitive_scalar_functional_values<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																											const vector<const LinearAlgebra::distributed::Vector<double>*>,
																											Vector<double>&)
const;

template
bool
AssemblyHelper<3>::get_nonprimitive_scalar_functional_values<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																											const vector<const LinearAlgebra::distributed::Vector<double>*>,
																											Vector<double>&)
const;
#endif // DEAL_II_WITH_MPI

//public get_nonprimitive_scalar_functional_values
template
bool
AssemblyHelper<2>::get_nonprimitive_scalar_functional_values<Vector<double>>(	const Vector<double>&,
																				const vector<const Vector<double>*>,
																				map<const ScalarFunctional<2, 2>*, double>&,
																				map<const ScalarFunctional<1, 2>*, double>&)
const;

template
bool
AssemblyHelper<3>::get_nonprimitive_scalar_functional_values<Vector<double>>(	const Vector<double>&,
																				const vector<const Vector<double>*>,
																				map<const ScalarFunctional<3, 3>*, double>&,
																				map<const ScalarFunctional<2, 3>*, double>&)
const;

#ifdef DEAL_II_WITH_MPI
template
bool
AssemblyHelper<2>::get_nonprimitive_scalar_functional_values<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																											const vector<const LinearAlgebra::distributed::Vector<double>*>,
																											map<const ScalarFunctional<2, 2>*, double>&,
																											map<const ScalarFunctional<1, 2>*, double>&)
const;

template
bool
AssemblyHelper<3>::get_nonprimitive_scalar_functional_values<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																											const vector<const LinearAlgebra::distributed::Vector<double>*>,
																											map<const ScalarFunctional<3, 3>*, double>&,
																											map<const ScalarFunctional<2, 3>*, double>&)
const;
#endif // DEAL_II_WITH_MPI

//public get_maximum_step_length
template
double
AssemblyHelper<2>::get_maximum_step_length<Vector<double>>(	const Vector<double>&,
															const vector<const Vector<double>*>,
															const Vector<double>&)
const;

template
double
AssemblyHelper<3>::get_maximum_step_length<Vector<double>>(	const Vector<double>&,
															const vector<const Vector<double>*>,
															const Vector<double>&)
const;

#ifdef DEAL_II_WITH_MPI
template
double
AssemblyHelper<2>::get_maximum_step_length<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																						const vector<const LinearAlgebra::distributed::Vector<double>*>,
																						const LinearAlgebra::distributed::Vector<double>&)
const;

template
double
AssemblyHelper<3>::get_maximum_step_length<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																						const vector<const LinearAlgebra::distributed::Vector<double>*>,
																						const LinearAlgebra::distributed::Vector<double>&)
const;
#endif // DEAL_II_WITH_MPI

//get_initial_fields_vector
template
void
AssemblyHelper<2>::get_initial_fields_vector<Vector<double>>(	Vector<double>&,
																const AffineConstraints<double>*)
const;

template
void
AssemblyHelper<3>::get_initial_fields_vector<Vector<double>>(	Vector<double>&,
																const AffineConstraints<double>*)
const;

#ifdef DEAL_II_WITH_MPI
template
void
AssemblyHelper<2>::get_initial_fields_vector<LinearAlgebra::distributed::Vector<double>>(	LinearAlgebra::distributed::Vector<double>&,
																							const AffineConstraints<double>*)
const;

template
void
AssemblyHelper<3>::get_initial_fields_vector<LinearAlgebra::distributed::Vector<double>>(	LinearAlgebra::distributed::Vector<double>&,
																							const AffineConstraints<double>*)
const;
#endif // DEAL_II_WITH_MPI

//generate_sparsity_pattern_by_simulation

template
void
AssemblyHelper<2>::generate_sparsity_pattern_by_simulation<DynamicSparsityPattern>(	DynamicSparsityPattern&,
																					const AffineConstraints<double>&)
const;

template
void
AssemblyHelper<3>::generate_sparsity_pattern_by_simulation<DynamicSparsityPattern>(	DynamicSparsityPattern&,
																					const AffineConstraints<double>&)
const;

template
void
AssemblyHelper<2>::generate_sparsity_pattern_by_simulation<TwoBlockSparsityPattern>(TwoBlockSparsityPattern&,
																					const AffineConstraints<double>&)
const;

template
void
AssemblyHelper<3>::generate_sparsity_pattern_by_simulation<TwoBlockSparsityPattern>(TwoBlockSparsityPattern&,
																					const AffineConstraints<double>&)
const;

//assemble_system

template
bool
AssemblyHelper<2>::assemble_system<Vector<double>, Vector<double>, SparseMatrix<double>>(	const Vector<double>&,
																							const vector<const Vector<double>*>,
																							const AffineConstraints<double>&,
																							double&,
																							Vector<double>&,
																							SparseMatrix<double>&,
																							const tuple<bool,bool,bool>)
const;

template
bool
AssemblyHelper<3>::assemble_system<Vector<double>, Vector<double>, SparseMatrix<double>>(	const Vector<double>&,
																							const vector<const Vector<double>*>,
																							const AffineConstraints<double>&,
																							double&,
																							Vector<double>&,
																							SparseMatrix<double>&,
																							const tuple<bool,bool,bool>)
const;

template
bool
AssemblyHelper<2>::assemble_system<Vector<double>, BlockVector<double>, TwoBlockMatrix<SparseMatrix<double>>>(	const Vector<double>&,
																												const vector<const Vector<double>*>,
																												const AffineConstraints<double>&,
																												double&,
																												BlockVector<double>&,
																												TwoBlockMatrix<SparseMatrix<double>>&,
																												const tuple<bool,bool,bool>)
const;

template
bool
AssemblyHelper<3>::assemble_system<Vector<double>, BlockVector<double>, TwoBlockMatrix<SparseMatrix<double>>>(	const Vector<double>&,
																												const vector<const Vector<double>*>,
																												const AffineConstraints<double>&,
																												double&,
																												BlockVector<double>&,
																												TwoBlockMatrix<SparseMatrix<double>>&,
																												const tuple<bool,bool,bool>)
const;

#ifdef DEAL_II_WITH_MPI
#ifdef DEAL_II_WITH_PETSC
template
bool
AssemblyHelper<2>::assemble_system<LinearAlgebra::distributed::Vector<double>, PETScWrappers::MPI::BlockVector, dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>>(	const LinearAlgebra::distributed::Vector<double>&,
																																																	const vector<const LinearAlgebra::distributed::Vector<double>*>,
																																																	const AffineConstraints<double>&,
																																																	double&,
																																																	PETScWrappers::MPI::BlockVector&,
																																																	dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>&,
																																																	const tuple<bool,bool,bool>)
const;

template
bool
AssemblyHelper<3>::assemble_system<LinearAlgebra::distributed::Vector<double>, PETScWrappers::MPI::BlockVector, dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>>(	const LinearAlgebra::distributed::Vector<double>&,
																																																	const vector<const LinearAlgebra::distributed::Vector<double>*>,
																																																	const AffineConstraints<double>&,
																																																	double&,
																																																	PETScWrappers::MPI::BlockVector&,
																																																	dealii::GalerkinTools::parallel::TwoBlockMatrix<PETScWrappers::MPI::SparseMatrix>&,
																																																	const tuple<bool,bool,bool>)
const;
#endif // DEAL_II_WITH_PETSC
#endif // DEAL_II_WITH_MPI

//write_output_independent_fields
template
std::pair<const string, const string>
AssemblyHelper<2>::write_output_independent_fields<Vector<double>>(	const Vector<double>&,
																	const string,
																	const string,
																	const unsigned int,
																	const vector<SmartPointer<const DataPostprocessor<2>>>&,
																	const vector<SmartPointer<const DataPostprocessor<2>>>&,
																	const unsigned int)
const;

template
std::pair<const string, const string>
AssemblyHelper<3>::write_output_independent_fields<Vector<double>>(	const Vector<double>&,
																	const string,
																	const string,
																	const unsigned int,
																	const vector<SmartPointer<const DataPostprocessor<3>>>&,
																	const vector<SmartPointer<const DataPostprocessor<3>>>&,
																	const unsigned int)
const;

#ifdef DEAL_II_WITH_MPI
template
std::pair<const string, const string>
AssemblyHelper<2>::write_output_independent_fields<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																								const string,
																								const string,
																								const unsigned int,
																								const vector<SmartPointer<const DataPostprocessor<2>>>&,
																								const vector<SmartPointer<const DataPostprocessor<2>>>&,
																								const unsigned int)
const;

template
std::pair<const string, const string>
AssemblyHelper<3>::write_output_independent_fields<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																								const string,
																								const string,
																								const unsigned int,
																								const vector<SmartPointer<const DataPostprocessor<3>>>&,
																								const vector<SmartPointer<const DataPostprocessor<3>>>&,
																								const unsigned int)
const;
#endif // DEAL_II_WITH_MPI

template
void
AssemblyHelper<2>::compare_derivatives_with_numerical_derivatives<Vector<double>>(	const Vector<double>&,
																					const vector<const Vector<double>*>,
																					const string,
																					const double)
const;

template
void
AssemblyHelper<3>::compare_derivatives_with_numerical_derivatives<Vector<double>>(	const Vector<double>&,
																					const vector<const Vector<double>*>,
																					const string,
																					const double)
const;

#ifdef DEAL_II_WITH_MPI
#ifdef DEAL_II_WITH_PETSC

template
void
AssemblyHelper<2>::compare_derivatives_with_numerical_derivatives<PETScWrappers::MPI::Vector>(	const PETScWrappers::MPI::Vector&,
																								const vector<const PETScWrappers::MPI::Vector*>,
																								const string,
																								const double)
const;

template
void
AssemblyHelper<3>::compare_derivatives_with_numerical_derivatives<PETScWrappers::MPI::Vector>(	const PETScWrappers::MPI::Vector&,
																								const vector<const PETScWrappers::MPI::Vector*>,
																								const string,
																								const double)
const;

#endif // DEAL_II_WITH_PETSC
#endif // DEAL_II_WITH_MPI

template
std::pair<const double, const double>
AssemblyHelper<2>::compute_distance_to_other_solution<Vector<double>>(	const Vector<double>&,
																		const Vector<double>&,
																		const AssemblyHelper<2>&,
																		const Quadrature<2>,
																		const Quadrature<1>,
																		const VectorTools::NormType,
																		const ComponentMask,
																		const ComponentMask,
																		const double,
																		const Vector<double>,
																		const Vector<double>)
const;

template
std::pair<const double, const double>
AssemblyHelper<3>::compute_distance_to_other_solution<Vector<double>>(	const Vector<double>&,
																		const Vector<double>&,
																		const AssemblyHelper<3>&,
																		const Quadrature<3>,
																		const Quadrature<2>,
																		const VectorTools::NormType,
																		const ComponentMask,
																		const ComponentMask,
																		const double,
																		const Vector<double>,
																		const Vector<double>)
const;

#ifdef DEAL_II_WITH_MPI

template
std::pair<const double, const double>
AssemblyHelper<2>::compute_distance_to_other_solution<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																									const LinearAlgebra::distributed::Vector<double>&,
																									const AssemblyHelper<2>&,
																									const Quadrature<2>,
																									const Quadrature<1>,
																									const VectorTools::NormType,
																									const ComponentMask,
																									const ComponentMask,
																									const double,
																									const Vector<double>,
																									const Vector<double>)
const;

template
std::pair<const double, const double>
AssemblyHelper<3>::compute_distance_to_other_solution<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																									const LinearAlgebra::distributed::Vector<double>&,
																									const AssemblyHelper<3>&,
																									const Quadrature<3>,
																									const Quadrature<2>,
																									const VectorTools::NormType,
																									const ComponentMask,
																									const ComponentMask,
																									const double,
																									const Vector<double>,
																									const Vector<double>)
const;

#endif // DEAL_II_WITH_MPI

template
std::pair<const double, const double>
AssemblyHelper<2>::compute_distance_to_exact_solution<Vector<double>>(	const Vector<double>&,
																		const Function<2>&,
																		const Function<2>&,
																		const Quadrature<2>,
																		const Quadrature<1>,
																		const VectorTools::NormType,
																		const ComponentMask,
																		const ComponentMask,
																		const double)
const;

template
std::pair<const double, const double>
AssemblyHelper<3>::compute_distance_to_exact_solution<Vector<double>>(	const Vector<double>&,
																		const Function<3>&,
																		const Function<3>&,
																		const Quadrature<3>,
																		const Quadrature<2>,
																		const VectorTools::NormType,
																		const ComponentMask,
																		const ComponentMask,
																		const double)
const;

#ifdef DEAL_II_WITH_MPI

template
std::pair<const double, const double>
AssemblyHelper<2>::compute_distance_to_exact_solution<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																									const Function<2>&,
																									const Function<2>&,
																									const Quadrature<2>,
																									const Quadrature<1>,
																									const VectorTools::NormType,
																									const ComponentMask,
																									const ComponentMask,
																									const double)
const;

template
std::pair<const double, const double>
AssemblyHelper<3>::compute_distance_to_exact_solution<LinearAlgebra::distributed::Vector<double>>(	const LinearAlgebra::distributed::Vector<double>&,
																									const Function<3>&,
																									const Function<3>&,
																									const Quadrature<3>,
																									const Quadrature<2>,
																									const VectorTools::NormType,
																									const ComponentMask,
																									const ComponentMask,
																									const double)
const;

#endif // DEAL_II_WITH_MPI

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
