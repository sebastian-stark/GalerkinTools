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

#ifndef GALERKINTOOLS_DOFHANDLERSYSTEM_H_
#define GALERKINTOOLS_DOFHANDLERSYSTEM_H_

#include <deal.II/hp/dof_handler.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <galerkin_tools/config.h>
#include <galerkin_tools/triangulation_system.h>
#include <galerkin_tools/dof_renumbering.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
class DoFHandlerSystem;

/**
 * This is an active cell iterator (for interface cells) pretty much like the active_cell_iterator of deal.II. However, the difference
 * is that it provides with the function InterfaceCellDoFIterator::get_dof_indices(), which returns the
 * dof indices of the cells in the global ordering of the DoFHandlerSystem instead of the dof ordering of the
 * deal.II dof handler. By using the function interface_cell_dof_iterator.get_dof_indices() you will get
 * the indices in the global ordering of the DoFHandlerSystem. In contrast, if you use interface_cell_dof_iterator->get_dof_indices(),
 * you get the indices in the ordering of the deal.II dof handler.
 *
 */
template<unsigned int spacedim>
class InterfaceCellDoFIterator : public hp::DoFHandler<spacedim-1, spacedim>::active_cell_iterator
{
private:

	/**
	 * The underlying DoFHandlerSystem
	 */
	const DoFHandlerSystem<spacedim>&
	dof_handler_system;

public:

	/**
	 * Constructor
	 *
	 * @param[in]	interface_cell		The interface cell of the triangulation
	 *
	 * @param[in]	dof_handler_system	The underlying DoFHandlerSystem
	 */
	InterfaceCellDoFIterator(	const TriaIterator<CellAccessor<spacedim-1, spacedim>>&	interface_cell,
								const DoFHandlerSystem<spacedim>& 						dof_handler_system);

	/**
	 * Constructor for the end() element
	 *
	 * @param[in]	dof_handler_system	The underlying DoFHandlerSystem
	 */
	InterfaceCellDoFIterator(const DoFHandlerSystem<spacedim>& dof_handler_system);

	/**
	 * %Function returning the dof indices in the global ordering of the dof handler system InterfaceCellDoFIterator::dof_handler_system
	 *
	 * @param[out]	dof_indices		The returned dof indices
	 */
	void
	get_dof_indices(std::vector<types::global_dof_index >& dof_indices)
	const;
};

/**
 * This is an active cell iterator (for domain cells) pretty much like the active_cell_iterator of deal.II. However, the difference
 * is that it provides the function DomainCellDoFIterator::get_dof_indices(), which returns the
 * dof indices of the cells in the global ordering of the DoFHandlerSystem instead of the dof ordering of the
 * deal.II dof handler. By using the function domain_cell_dof_iterator.get_dof_indices() you will get
 * the indices in the global ordering of the DoFHandlerSystem. In contrast, if you use domain_cell_dof_iterator->get_dof_indices(),
 * you get the indices in the ordering of the deal.II dof handler.
 *
 */
template<unsigned int spacedim>
class DomainCellDoFIterator : public hp::DoFHandler<spacedim, spacedim>::active_cell_iterator
{
private:

	/**
	 * The underlying DoFHandlerSystem
	 */
	const DoFHandlerSystem<spacedim>&
	dof_handler_system;

public:

	/**
	 * Constructor
	 *
	 * @param[in]	domain_cell			The domain cell of the triangulation
	 *
	 * @param[in]	dof_handler_system	The underlying DoFHandlerSystem
	 */
	DomainCellDoFIterator(	const TriaIterator<CellAccessor<spacedim, spacedim>>& 	domain_cell,
							const DoFHandlerSystem<spacedim>& 						dof_handler_system);

	/**
	 * Constructor for the end() element
	 *
	 * @param[in]	dof_handler_system	The underlying DoFHandlerSystem
	 */
	DomainCellDoFIterator(const DoFHandlerSystem<spacedim>& dof_handler_system);

	/**
	 * %Function returning the dof indices in the global ordering of the dof handler system DomainCellDoFIterator::dof_handler_system
	 *
	 * @param[out]	dof_indices		The returned dof indices
	 */
	void
	get_dof_indices(std::vector<types::global_dof_index >& dof_indices)
	const;
};


/**
 * This class is the equivalent to InterfaceCellDomainCells with
 * information about the degrees of freedom associated with the interface
 * cell and the neighboring domain cells
 *
 * The InterfaceCellDomainCellsDoF class inherits from Subscriptor in order to be
 * able to check that InterfaceCellDomainCellsDoF objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam	spacedim	The spatial dimension into which the interface is embedded.
 */
template<unsigned int spacedim>
class InterfaceCellDomainCellsDoF : public Subscriptor
{

public:

	/**
	 * A convenience typedef for a deal.II TriaIterator referring to an interface cell
	 */
	typedef TriaIterator<CellAccessor<spacedim-1, spacedim>>
	InterfaceCell;

	/**
	 * A convenience typedef for a deal.II cell iterator referring to an interface cell (with dof information)
	 */
	typedef typename hp::DoFHandler<spacedim-1, spacedim>::cell_iterator
	InterfaceCellDoF;

	/**
	 * A convenience typedef for a deal.II TriaIterator referring to a domain cell
	 */
	typedef TriaIterator<CellAccessor<spacedim, spacedim>>
	DomainCell;

	/**
	 * A convenience typedef for a deal.II cell iterator referring to a domain cell (with dof information)
	 */
	typedef typename hp::DoFHandler<spacedim, spacedim>::cell_iterator
	DomainCellDoF;

	/**
	 * A deal.II iterator with dof information to the interface cell
	 */
	const InterfaceCellDoFIterator<spacedim> interface_cell;

	/**
	 * A deal.II iterator with dof information to the domain cell adjacent to the interface cell
	 * on the minus side.
	 *
	 *
	 * InterfaceCellDomainCellsDoF::domain_cell_minus is refined to the same degree as
	 * InterfaceCellDomainCellsDoF::interface_cell if
	 * InterfaceCellDomainCellsDoF::refinement_case == @p minus_is_finer or if
	 * InterfaceCellDomainCellsDoF::refinement_case == @p equally_fine or if
	 * InterfaceCellDomainCellsDoF::refinement_case == @p at_boundary.
	 * Otherwise, it is once less refined than InterfaceCellDomainCellsDoF::interface_cell.
	 */
	const DomainCellDoFIterator<spacedim> domain_cell_minus;

	/**
	 * The face number of the face of the domain cell on the minus side
	 * lying on the interface
	 */
	const unsigned int face_minus;

	/**
	 * A deal.II iterator with dof information to the domain cell adjacent to the interface cell
	 * on the plus side.
	 *
	 * InterfaceCellDomainCellsDoF::domain_cell_plus is refined to the same degree as
	 * InterfaceCellDomainCellsDoF::interface_cell if
	 * InterfaceCellDomainCellsDoF::refinement_case == @p plus_is_finer or if
	 * InterfaceCellDomainCellsDoF::refinement_case == @p equally_fine;
	 * and it is once less refined than InterfaceCellDomainCellsDoF::interface_cell if
	 * InterfaceCellDomainCellsDoF::refinement_case == @p minus_is_finer.
	 * Otherwise, the interface is actually a boundary, in which case
	 * InterfaceCellDomainCellsDoF::domain_cell_plus is set to the same cell as
	 * InterfaceCellDomainCellsDoF::domain_cell_minus.
	 */
	const DomainCellDoFIterator<spacedim> domain_cell_plus;

	/**
	 * The face number of the face of the domain cell on the plus side
	 * lying on the interface
	 */
	const unsigned int face_plus;

	/**
	 * The refinement situation at the interface cell. See ::InterfaceRefinementCase
	 * for further documentation.
	 */
	const InterfaceRefinementCase refinement_case;

	/**
	 * This member applies only to the case, where the domain mesh on one side of the
	 * interface cell is refined differently than the domain mesh on the other side (i.e., where
	 * InterfaceCellDomainCellsDoF::refinement_case == @p minus_is_finer or
	 * InterfaceCellDomainCellsDoF::refinement_case == @p plus_is_finer).
	 * It represents the subface of the domain cell face on the coarser side of the interface
	 * coinciding with the interface cell (and, by assumption on the properties of the meshes,
	 * with the domain cell face on the finer side of the interface as well).
	 *
	 * See also the documentation of the CellAccessor<dim, spacedim>::neighbor_of_coarser_neighbor()
	 * method.
	 */
	const unsigned int subface;

	/**
	 * This constructor constructs an InterfaceCellDomainCellsDoF object from a corresponding
	 * InterfaceCellDomainCells object and the DoFHandlerSystem.
	 *
	 * @param[in]	interface_cell_domain_cells	InterfaceCellDomainCells object
	 * @param[in]	dof_handler_system			the DoFHandlerSystem
	 */
	InterfaceCellDomainCellsDoF(const InterfaceCellDomainCells<spacedim>&	interface_cell_domain_cells,
								const DoFHandlerSystem<spacedim>&			dof_handler_system);

	/**
	 * The destructor of InterfaceCellDomainCellsDoF essentially checks before destruction that the
	 * InterfaceCellDomainCellsDoF object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	~InterfaceCellDomainCellsDoF();


	/**
	 * @return	This method returns a tuple with the types::material_id%s associated with
	 * 			InterfaceCellDomainCellsDoF::interface_cell, InterfaceCellDomainCellsDoF::domain_cell_minus,
	 * 			InterfaceCellDomainCellsDoF::domain_cell_plus (in this order).
	 *
	 * 			If the interface cell is actually at the boundary, the third member of the tuple will
	 * 			be returned as numbers::invalid_material_id
	 */
	std::tuple<const types::material_id, const types::material_id, const types::material_id>
	get_material_ids()
	const;

	/**
	 * Method to get the maps between local and global dof indices related to the interface cell
	 * InterfaceCellDomainCellsDoF::interface_cell and the neighboring domain cells
	 * InterfaceCellDomainCellsDoF::domain_cell_minus and InterfaceCellDomainCellsDoF::domain_cell_plus.
	 *
	 * @param[out]	dof_indices_local_global		mapping between local dof indices and global ones on InterfaceCellDomainCellsDoF::interface_cell
	 *
	 * @param[out]	dof_indices_local_global_minus	mapping between local dof indices and global ones on InterfaceCellDomainCellsDoF::domain_cell_minus
	 *
	 * @param[out]	dof_indices_local_global_plus	mapping between local dof indices and global ones on InterfaceCellDomainCellsDoF::domain_cell_plus
	 */
	void
	get_dof_indices_local_global_interface(	std::vector<unsigned int>& dof_indices_local_global,
											std::vector<unsigned int>& dof_indices_local_global_minus,
											std::vector<unsigned int>& dof_indices_local_global_plus)
	const;

	/**
	 * Method to get the map between local and global dof indices related to the interface cell
	 * InterfaceCellDomainCellsDoF::interface_cell.
	 *
	 * @param[out]	dof_indices_local_global		mapping between local dof indices and global ones on InterfaceCellDomainCellsDoF::interface_cell
	 */
	void
	get_dof_indices_local_global_interface(	std::vector<unsigned int>& dof_indices_local_global)
	const;
};

/**
 * Class defining a dof handler system consisting of a dof handler
 * for domain related dofs (i.e., those dofs which are related to the
 * independent fields \f$u^\Omega_\epsilon\f$) and a dof handler
 * for interface related dofs (i.e., those dofs which are related to the
 * independent fields \f$u^\Sigma_\eta\f$). In addition, also a number of
 * dofs which are not related to a mesh is allowed.
 *
 * The standard ordering (i.e. without explicit renumbering)
 * used in this library is that the domain related dofs are followed by the
 * interface related dofs, which are in turn followed by the additional dofs not related to a mesh.
 *
 * Conceptually, all independent fields \f$u^\Omega_\epsilon\f$ and \f$u^\Sigma_\eta\f$
 * are defined on the entire domain \f$\Omega\f$ and the entire interface \f$\Sigma\f$, respectively.
 * However, certain independent fields may be extended by zero to certain portions of the domain and interface,
 * respectively (just because they are not needed there according to the model). In order to be able to
 * conveniently deal with this situation, hp::DoFHandler%s are used internally. This allows to use
 * FE_Nothing finite elements everywhere, where an independent field is zero (see also step-46 in the deal.II tutorial).
 *
 * @tparam	spacedim	spatial dimension
 *
 */
template<unsigned int spacedim>
class DoFHandlerSystem
{

public:

	/**
	 * A convenience typedef for a deal.II TriaIterator referring to an interface cell
	 */
	typedef TriaIterator<CellAccessor<spacedim-1, spacedim>>
	InterfaceCell;

	/**
	 * A convenience typedef for a deal.II cell iterator referring to an interface cell (with dof information)
	 */
	typedef typename hp::DoFHandler<spacedim-1, spacedim>::cell_iterator
	InterfaceCellDoF;

	/**
	 * A convenience typedef for a deal.II TriaIterator referring to a domain cell
	 */
	typedef TriaIterator<CellAccessor<spacedim, spacedim>>
	DomainCell;

	/**
	 * A convenience typedef for a deal.II cell iterator referring to a domain cell (with dof information)
	 */
	typedef typename hp::DoFHandler<spacedim, spacedim>::cell_iterator
	DomainCellDoF;

private:

	/**
	 * The TriangulationSystem underlying the DoFHandlerSystem
	 */
	const SmartPointer<const TriangulationSystem<spacedim>>
	tria_system;

	/**
	 * hp::DoFHandler for domain related dofs
	 */
	std::shared_ptr<hp::DoFHandler<spacedim, spacedim>>
	dof_handler_domain;

	/**
	 * hp::DoFHandler for interface related dofs
	 */
	std::shared_ptr<hp::DoFHandler<spacedim-1, spacedim>>
	dof_handler_interface;

	/**
	 * A vector with the dof indices of the additional dofs which are not related to a mesh
	 * These numbers are consecutive and start after the finite element related dofs.
	 */
	std::vector<unsigned int>
	dof_indices_C;

	/**
	 * A vector with the locally owned dof indices of the additional dofs which are not related to a mesh
	 */
	std::vector<unsigned int>
	locally_owned_dof_indices_C;

	/**
	 * This vector contains the interface cells and corresponding domain cells of the active mesh (with dof information).
	 * This member is the equivalent to TriangulationSystem::active_interface_cell_domain_cells
	 */
	std::vector<InterfaceCellDomainCellsDoF<spacedim>>
	active_interface_cell_domain_cells;

	/**
	 * A DoFRenumbering object defining the renumbering of dofs.
	 * If this is a @p nullptr, no renumbering will be applied and the dof numbering of the domain and interface related dof handlers is used with the
	 * convention that domain related dofs come first, followed by interface related dofs (which are shifted by the number
	 * of dofs associated with the domain related domain handler to avoid duplicate dofs),
	 * and finally followed by the dofs not related to any mesh. This ordering is also used
	 * as the basis for renumbering by the DoFRenumbering object.
	 *
	 * Note that the dof numbering can be influenced in two ways: (1) Renumbering of the dofs of the domain and interface related
	 * domain handlers; (2) this DoFRenumbering object. The reason for including the second possibility is that the dof
	 * handlers by themselves can only be renumbered separately, meaning that e.g. for the domain related dof handler the dof indices
	 * can only be in the range [0, number of domain related dofs). This makes it impossible to achieve a contiguous indexing of the union
	 * of domain and interface related dofs on each processor in parallel computations. However, the latter can be achieved by the
	 * DoFRenumbering object provided here.
	 */
	SmartPointer<const DoFRenumbering>
	dof_renumbering = nullptr;

	/**
	 * The index set with the locally owned dofs
	 */
	IndexSet
	locally_owned_dofs;

	/**
	 * The index set with the locally relevant dofs
	 */
	IndexSet
	locally_relevant_dofs;

	/**
	 * A vector with the numbers of locally owned dofs per processor
	 */
	std::vector<unsigned int> n_dofs_per_processor;

	/**
	 * The index set with the locally owned dofs in the standard numbering (i.e., without application of renumbering)
	 */
	IndexSet
	locally_owned_dofs_standard_numbering;

	/**
	 * The index set with the locally relevant dofs  (i.e., without application of renumbering)
	 */
	IndexSet
	locally_relevant_dofs_standard_numbering;

	/**
	 * A list of connections with which this object connects to the triangulations to get information about when the triangulations change.
	 */
	std::vector<boost::signals2::connection>
	tria_listeners;

	/**
	 * This function updates DoFHandlerSystem::active_interface_cell_domain_cells.
	 * It is called automatically after any refinement of the TriangulationSystem
	 * DoFHandlerSystem::tria_system.
	 */
	void
	update_interface_domain_relation();

	/**
	 * %Function setting DoFHandlerSystem::locally_owned_dofs
	 */
	void
	set_locally_owned_dofs();

	/**
	 * %Function setting DoFHandlerSystem::locally_relevant_dofs
	 */
	void
	set_locally_relevant_dofs();

	/**
	 * %Function setting DoFHandlerSystem::locally_owned_dofs_standard_numbering
	 */
	void
	set_locally_owned_dofs_standard_numbering();

	/**
	 * %Function setting DoFHandlerSystem::locally_relevant_dofs_standard_numbering
	 */
	void
	set_locally_relevant_dofs_standard_numbering();

#ifdef DEAL_II_WITH_MPI
	/**
	 * This function allows to renumber a distributed vector, which is renumbered according to DoFHandlerSystem::dof_renumbering,
	 * back to the standard numbering of this DoFHandlerSystem
	 *
	 * The result of this renumbering can be sliced to the interval of elements [window_begin, window_end). The latter
	 * operation creates a "view" of the vector through the window [@p window_begin, @p window_end). I.e., the resulting vector
	 * contains @p window_end - @p window_begin elements. This is used to extract different parts of the vector (domain related, etc.)
	 *
	 * The function will not change ownership of vector elements.
	 *
	 * Also, this function will only work if the locally owned range of @p in_vect and @p out_vect is contiguous!
	 *
	 * @param[in]	in_vect								The vector to be rearranged
	 *
	 * @param[out]	out_vect							The rearranged and sliced vector.
	 *
	 * @param[in]	window_begin						The first (global) element to be included in the return vector
	 *
	 * @param[in]	window_end							The past the last (global) element to be included in the return vector
	 * 													(if larger than the vector size, all elements starting from @p window_begin will be included)
	 */
	void
	split_vector_implementation(const LinearAlgebra::distributed::Vector<double>&	in_vect,
								LinearAlgebra::distributed::Vector<double>&			out_vect,
								const unsigned int									window_begin,
								const unsigned int									window_end)
	const;
#endif // DEAL_II_WITH_MPI

	/**
	 * This function allows to renumber a sequential vector, which is renumbered according to DoFHandlerSystem::dof_renumbering,
	 * back to the standard numbering of this DoFHandlerSystem
	 *
	 * The result of this renumbering can be sliced to the interval of elements [window_begin, window_end). The latter
	 * operation creates a "view" of the vector through the window [@p window_begin, @p window_end). I.e., the resulting vector
	 * contains @p window_end - @p window_begin elements. This is used to extract different parts of the vector (domain related, etc.)
	 *
	 * @param[in]	in_vect								The vector to be rearranged
	 *
	 * @param[out]	out_vect							The rearranged and sliced vector.
	 *
	 * @param[in]	window_begin						The first (global) element to be included in the return vector
	 *
	 * @param[in]	window_end							The past the last (global) element to be included in the return vector
	 * 													(if larger than the vector size, all elements starting from @p window_begin will be included)
	 */
	void
	split_vector_implementation(const Vector<double>&	in_vect,
								Vector<double>&			out_vect,
								const unsigned int		window_begin,
								const unsigned int		window_end)
	const;

	/**
	 * This is necssary to get access to the current renumbering scheme
	 */
	template<unsigned int> friend class InterfaceCellDoFIterator;

	/**
	 * This is necssary to get access to the current renumbering scheme
	 */
	template<unsigned int> friend class DomainCellDoFIterator;

public:

	/**
	 * Construct a DoFHandlerSystem from a TriangulationSystem
	 *
	 * @param[in]	tria_system		The TriangulationSystem underlying the DoFHandlerSystem
	 */
	DoFHandlerSystem(const TriangulationSystem<spacedim>& tria_system);

	/**
	 * Destructor
	 */
	~DoFHandlerSystem();

	/**
	 * This function actually distributes the dofs and must be called after any
	 * refinement of the triangulation. Before you call this function, make sure
	 * that all cells in the mesh have assigned the correct @p active_fe_index!
	 *
	 * @param[in]	fe_collection_domain	the fe collection to be used for the distribution of dofs
	 * 										on the domain
	 *
	 * @param[in]	fe_collection_interface	the fe collection to be used for the distribution of dofs
	 * 										on the interface
	 *
	 * @param[in]	n_additional_dofs		A number of additional dofs to be included in the DoFHandlerSystem, which
	 * 										are not related to a mesh
	 */
	void
	distribute_dofs(const hp::FECollection<spacedim, spacedim>&		fe_collection_domain,
					const hp::FECollection<spacedim-1, spacedim>&	fe_collection_interface,
					const unsigned int								n_additional_dofs = 0);

	/**
	 * This function sets the FECollection objects to be used. However, it does not distribute any dofs
	 *
	 * @param[in]	fe_collection_domain	the fe collection to be used for the distribution of dofs
	 * 										on the domain
	 *
	 * @param[in]	fe_collection_interface	the fe collection to be used for the distribution of dofs
	 * 										on the interface
	 */
	void
	set_fe(	const hp::FECollection<spacedim, spacedim>&		fe_collection_domain,
			const hp::FECollection<spacedim-1, spacedim>&	fe_collection_interface);

	/**
	 * @return	An iterator (with dof information) to the first active interface cell. Note that this is not
	 * 			a deal.II iterator. Rather, it is an iterator to the first element
	 * 			in DoFHandlerSystem::active_interface_cell_domain_cells.
	 * 			However, through the returned iterator, deal.II iterators to the
	 * 			interface cells and adjacent domain cells can be obtained.
	 */
	typename std::vector< InterfaceCellDomainCellsDoF<spacedim> >::iterator
	interface_begin_active();

	/**
	 * @return	An iterator (with dof information) to past the last active interface cell. Note that this is not
	 * 			a deal.II iterator. Rather, it is an iterator to the past the end element
	 * 			in DoFHandlerSystem::active_interface_cell_domain_cells.
	 */
	typename std::vector< InterfaceCellDomainCellsDoF<spacedim> >::iterator
	interface_end_active();

	/**
	 * @return	This returns a const reference to DoFHandlerSystem::active_interface_cell_domain_cells,
	 * 			which is essentially be meant for range based loops instead of using the functions
	 * 			DoFHandlerSystem::interface_begin_active() and
	 * 			DoFHandlerSystem::interface_end_active() (hence the name of the function).
	 */
	const typename std::vector< InterfaceCellDomainCellsDoF<spacedim> >&
	interface_active_iterators()
	const;

	/**
	 * @return	An iterator (with dof information) to the first active domain cell. Note that this is an extension of
	 * 			a deal.II iterator. If you ask for the dof indices by @p iterator_name.get_dof_indices(), you will get
	 * 			the dof indices according to the global ordering of this DoFHandlerSystem. If you ask for iterator_name->get_dof_indices()
	 * 			you will get the dof indices according to the deal.II ordering in the domain related dof handler.
	 */
	DomainCellDoFIterator<spacedim>
	domain_begin_active()
	const;

	/**
	 * @return	An iterator (with dof information) to past the last active domain cell.
	 */
	DomainCellDoFIterator<spacedim>
	domain_end_active()
	const;

	/**
	 * @return	This returns the iterator range between DoFHandlerSystem::domain_begin_active() and DoFHandlerSystem::domain_end_active()
	 */
	IteratorRange<DomainCellDoFIterator<spacedim>>
	domain_active_iterators()
	const;

	/**
	 * @return	Number of dofs on domain (i.e., the number of dofs currently associated with
	 * 			DoFHandlerSystem::dof_handler_domain)
	 */
	unsigned int
	n_dofs_domain()
	const;

	/**
	 * @return	Number of dofs on interface (i.e., the number of dofs currently associated with
	 * 			DoFHandlerSystem::dof_handler_interface)
	 */
	unsigned int
	n_dofs_interface()
	const;

	/**
	 * @return	Number of dofs not related to a mesh
	 */
	unsigned int
	n_dofs_additional()
	const;

	/**
	 * @return	Total number of dofs
	 */
	unsigned int
	n_dofs()
	const;

	/**
	 * %Function returning the global dof indices of the dofs not related to a mesh
	 *
	 * @param[out]	dof_indices		dof indices of the dofs not related to a mesh
	 */
	void
	get_dof_indices(std::vector<unsigned int>& dof_indices)
	const;

	/**
	 * %Function returning the global dof index of the n-th dof not related to the mesh
	 *
	 * @param[out]	dof_index		dof for which the global dof index is to be returned
	 */
	unsigned int
	get_dof_index(const unsigned int& dof_index)
	const;

	/**
	 * @return	domain related dof handler
	 */
	const hp::DoFHandler<spacedim, spacedim>&
	get_dof_handler_domain()
	const;

	/**
	 * @return	interface related dof handler
	 */
	const hp::DoFHandler<spacedim-1, spacedim>&
	get_dof_handler_interface()
	const;

	/**
	 * @return	domain related dof handler (warning: this is a non-const reference, which is mainly intended for dof renumbering)
	 */
	hp::DoFHandler<spacedim, spacedim>&
	get_dof_handler_domain();

	/**
	 * @return	interface related dof handler (warning: this is a non-const reference, which is mainly intended for dof renumbering)
	 */
	hp::DoFHandler<spacedim-1, spacedim>&
	get_dof_handler_interface();

	/**
	 * @return	The set of locally owned dofs (possibly renumbered according to DoFHandlerSystem::dof_renumbering)
	 */
	const IndexSet&
	get_locally_owned_dofs()
	const;

	/**
	 * @return					The set of locally relevant dofs (possibly renumbered according to DoFHandlerSystem::dof_renumbering)
	 */
	const IndexSet&
	get_locally_relevant_dofs()
	const;

	/**
	 * Attach a DoFRenumbering object to the DoFHandlerSystem
	 *
	 * @param[in]	dof_renumbering		The new DoFHandlerSystem::dof_renumbering to be used
	 *
	 * @warning		The renumbering scheme must not change the range of the dof indices associated with the independent scalars (the numbering of the independent scalars
	 * 				within this range can however be changed)!
	 */
	void
	attach_dof_renumbering(const DoFRenumbering& dof_renumbering);

	/**
	 * Method generating the hanging node constraints (both for the domain related as well as the interface related dofs).
	 *
	 * @param[out]	constraint_matrix	constraint matrix with the hanging node constraints
	 *
	 * @warning		The constraint matrix will be cleared before the constraints are written.
	 * 				This means that if you have other constraints, it will be necessary to merge the constraint matrices.
	 * 				The @p local_lines property will be set to what DoFHandlerSystem::get_locally_relevant_dofs() of the object AssemblyHelper::dof_handler_system
	 * 				returns.
	 */
	void
	make_hanging_node_constraints(AffineConstraints<double>& constraint_matrix)
	const;

	/**
	 * This function renumbers a distributed vector, which is renumbered according to DoFHandlerSystem::dof_renumbering,
	 * back to the standard numbering of this DoFHandlerSystem and splits it into the domain part, the interface part,
	 * and the additional dof part.
	 *
	 * @param[in]	in_vect				the original global vector
	 *
	 * @param[out]	out_vect_domain		the resulting domain related vector
	 *
	 * @param[out]	out_vect_interface	the resulting interface related vector
	 *
	 * @param[out]	out_vect_C			the resulting additional dof related vector
	 */
	template<class VectorType>
	void
	split_vector(	const VectorType&	in_vect,
					VectorType&			out_vect_domain,
					VectorType&			out_vect_interface,
					VectorType&			out_vect_C)
	const;

	/**
	 * @return					The number of dofs per processor (the n-th element of the vector corresponds to the n-th processor)
	 */
	const std::vector<unsigned int>&
	get_n_dofs_per_processor()
	const;

	/**
	 * @param[in]	component			The component for which a global dof index is to be returned
	 *
	 * @return							A single global dof index associated with @p component. If no dof index is found, numbers::invalid_dof_index is returned.
	 *
	 * @note	This function can be used to constrain a single dof of an interface related field (e.g. to fix a scalar potential)
	 *
	 * @note	This works only if the finite element corresponding to @p component is primitive.
	 *
	 * @note	This searches in the locally owned dofs. In parallel, this will select one of the dofs and return this dof on all processors for which this
	 * 			dof is locally relevant (on all other processors numbers::invalid_dof_index is returned).
	 */
	unsigned int
	get_single_dof_index_component_interface(const unsigned int component)
	const;

	/**
	 * @param[in]	component			The component for which a global dof index is to be returned
	 *
	 * @return							A single global dof index associated with @p component. If no dof index is found, numbers::invalid_dof_index is returned.
	 *
	 * @note	This function can be used to constrain a single dof of a domain related field (e.g. to fix a scalar potential)
	 *
	 * @note	This works only if the finite element corresponding to @p component is primitive.
	 *
	 * @note	This searches in the locally owned dofs. In parallel, this will select one of the dofs and return this dof on all processors for which this
	 * 			dof is locally relevant (on all other processors numbers::invalid_dof_index is returned).
	 */
	unsigned int
	get_single_dof_index_component_domain(const unsigned int component)
	const;

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_DOFHANDLERSYSTEM_H_ */
