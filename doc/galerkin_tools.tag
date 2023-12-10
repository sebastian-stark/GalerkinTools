<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.1">
  <compound kind="file">
    <name>mainpage.dox</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/doc/</path>
    <filename>mainpage_8dox.html</filename>
  </compound>
  <compound kind="file">
    <name>assembly_helper.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>assembly__helper_8h.html</filename>
    <class kind="class">FunctionCell</class>
    <class kind="class">AssemblyHelper</class>
  </compound>
  <compound kind="file">
    <name>dependent_field.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>dependent__field_8h.html</filename>
    <class kind="class">DependentFieldTerm</class>
    <class kind="struct">DependentFieldComparatorWithoutCoefficient</class>
    <class kind="class">DependentField</class>
    <class kind="class">DependentField&lt; spacedim, spacedim &gt;</class>
  </compound>
  <compound kind="file">
    <name>dirichlet_constraint.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>dirichlet__constraint_8h.html</filename>
    <class kind="class">DirichletConstraint</class>
    <class kind="class">PointConstraint</class>
    <class kind="class">PointConstraint&lt; spacedim, spacedim &gt;</class>
  </compound>
  <compound kind="file">
    <name>dof_handler_system.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>dof__handler__system_8h.html</filename>
    <class kind="class">InterfaceCellDoFIterator</class>
    <class kind="class">DomainCellDoFIterator</class>
    <class kind="class">InterfaceCellDomainCellsDoF</class>
    <class kind="class">DoFHandlerSystem</class>
  </compound>
  <compound kind="file">
    <name>dof_renumbering.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>dof__renumbering_8h.html</filename>
    <class kind="class">DoFRenumbering</class>
    <class kind="class">DoFRenumberingOffset</class>
  </compound>
  <compound kind="file">
    <name>fe_dgq_anisotropic.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>fe__dgq__anisotropic_8h.html</filename>
    <class kind="class">FE_DGQAnisotropic</class>
  </compound>
  <compound kind="file">
    <name>fe_values_interface.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>fe__values__interface_8h.html</filename>
    <class kind="class">FEValuesInterface</class>
  </compound>
  <compound kind="file">
    <name>independent_field.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>independent__field_8h.html</filename>
    <class kind="class">IndependentField</class>
    <class kind="class">IndependentField&lt; 0, spacedim &gt;</class>
  </compound>
  <compound kind="file">
    <name>ldr.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>ldr_8h.html</filename>
    <namespace>Auxiliary</namespace>
    <member kind="function">
      <type>int</type>
      <name>compute_ldr</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>abe796c1529d11eed08fd05bca82f3002</anchor>
      <arglist>(FullMatrix&lt; double &gt; &amp;C, Vector&lt; double &gt; &amp;D, std::vector&lt; Vector&lt; double &gt;&gt; &amp;L, std::vector&lt; Vector&lt; double &gt;&gt; &amp;R)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>linear_material.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>linear__material_8h.html</filename>
    <class kind="class">LinearMaterialDomain</class>
    <class kind="class">LinearMaterialInterface</class>
  </compound>
  <compound kind="file">
    <name>scalar_functional.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>scalar__functional_8h.html</filename>
    <class kind="class">ScalarFunctional</class>
    <class kind="class">ScalarFunctional&lt; spacedim, spacedim &gt;</class>
    <class kind="class">ScalarFunctionalLocalElimination</class>
    <class kind="class">ScalarFunctionalLocalElimination&lt; spacedim, spacedim &gt;</class>
  </compound>
  <compound kind="file">
    <name>solver_wrapper.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>solver__wrapper_8h.html</filename>
    <class kind="class">SolverWrapper</class>
    <class kind="class">SolverWrapperUMFPACK</class>
    <class kind="class">BlockSolverWrapperUMFPACK</class>
    <class kind="class">SolverWrapperPETSc</class>
    <class kind="class">SolverWrapperPETScIterative</class>
    <class kind="class">BlockSolverWrapperPARDISO</class>
    <class kind="class">SolverWrapperMUMPS</class>
    <class kind="class">BlockSolverWrapperMUMPS</class>
    <class kind="class">BlockSolverWrapperUMFPACK2</class>
    <class kind="class">BlockSolverWrapperMA57</class>
    <member kind="function">
      <type>void</type>
      <name>pardisoinit</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>a470d3a8fdfc8e7517388a077fbb33fb0</anchor>
      <arglist>(void *, int *, int *, int *, double *, int *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pardiso</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>a8a77a0624251c8179bebace24f4c574e</anchor>
      <arglist>(void *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *, double *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pardiso_chkmatrix</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>ac68225cc41dd97996b1dd1fe7b84f3bb</anchor>
      <arglist>(int *, int *, double *, int *, int *, int *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pardiso_chkvec</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>ae2305586e44a2ea695f472069d2bee74</anchor>
      <arglist>(int *, int *, double *, int *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pardiso_printstats</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>a6813d2a3b274bbf613ef95d965062494</anchor>
      <arglist>(int *, int *, double *, int *, int *, int *, double *, int *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pardiso_residual</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>a558fe24a52b31ce8c2f40c7d5b41a9a1</anchor>
      <arglist>(int *, int *, double *, int *, int *, double *, double *, double *, double *, double *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ma57id_</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>a60ea38bdc152b7e893bded06017dd737</anchor>
      <arglist>(double *, int *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ma57ad_</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>a60b7996ef7e34acd652fd3a618d80994</anchor>
      <arglist>(int *, int *, int *, int *, int *, int *, int *, int *, int *, double *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ma57bd_</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>ae642f60d486e9ab3309a69f874297917</anchor>
      <arglist>(int *, int *, double *, double *, int *, int *, int *, int *, int *, int *, int *, double *, int *, double *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ma57cd_</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>a1ee773d49137399358320bfa351f64f3</anchor>
      <arglist>(int *, int *, double *, int *, int *, int *, int *, double *, int *, double *, int *, int *, int *, int *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ma57dd_</name>
      <anchorfile>solver__wrapper_8h.html</anchorfile>
      <anchor>ab3ada0a527461c49cdb0755e8435400f</anchor>
      <arglist>(int *, int *, int *, double *, int *, int *, double *, int *, int *, int *, double *, double *, double *, double *, int *, int *, double *, int *, double *)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>tools.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>tools_8h.html</filename>
    <namespace>Auxiliary</namespace>
    <member kind="function">
      <type>void</type>
      <name>convert_local_indices_to_global_indices</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>ad3c42d209f0ba8e9d4ce305060634bf1</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;dof_indices_local, std::vector&lt; unsigned int &gt; &amp;dof_indices_global, const std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>combine_dof_indices</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a1d90ebc8738df3d8c70b540034137019</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_interface, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_minus, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_plus, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_C, std::vector&lt; unsigned int &gt; &amp;dof_indices_interface_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_minus_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_plus_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_C_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_global_combined)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_to_index_set</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aca5fe5966aef4ba2293fb44095bfd86c</anchor>
      <arglist>(const DoFRenumbering &amp;dof_renumbering, const IndexSet &amp;in, IndexSet &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>renumber_constraints</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aa6148bcbaf5e3003717b4dd2a4da15b3</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const DoFRenumbering &amp;dof_renumbering=DoFRenumbering(), const bool close=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_dof_renumbering_contiguous</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a4261fc726ff965166f3f53242918acea</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system, DoFRenumberingOffset &amp;dof_renumbering_offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_map_dofs</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a072f85e6d745ae3c532bb0612f4bd3ce</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dhs_1, const DoFHandlerSystem&lt; spacedim &gt; &amp;dhs_2, std::vector&lt; unsigned int &gt; &amp;map_dofs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_vector</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a59c01a6511bffd7442693e86cd194ef1</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;in, Vector&lt; double &gt; &amp;out_0, Vector&lt; double &gt; &amp;out_1, const unsigned int size_1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_matrix</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>af746a1d08b1135a3684fc990f7b1384d</anchor>
      <arglist>(const FullMatrix&lt; double &gt; &amp;in, FullMatrix&lt; double &gt; &amp;out_00, FullMatrix&lt; double &gt; &amp;out_01, FullMatrix&lt; double &gt; &amp;out_10, FullMatrix&lt; double &gt; &amp;out_11, const unsigned int size_1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>communicate_bool</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>ad0ff28386be7b54b3487ae36e5a074fa</anchor>
      <arglist>(const bool local_bool, const MPI_Comm &amp;mpi_communicator)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>get_n_locally_owned_dofs_per_processor</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aaa9fb90b2f2308f766b8130578ab9c39</anchor>
      <arglist>(const DoFHandler&lt; dim, spacedim &gt; &amp;dof_handler)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>total_potential.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>total__potential_8h.html</filename>
    <class kind="class">TotalPotential</class>
  </compound>
  <compound kind="file">
    <name>total_potential_contribution.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>total__potential__contribution_8h.html</filename>
    <class kind="class">TotalPotentialContribution</class>
  </compound>
  <compound kind="file">
    <name>triangulation_system.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>triangulation__system_8h.html</filename>
    <class kind="class">InterfaceCellDomainCells</class>
    <class kind="class">TriangulationSystem</class>
    <class kind="class">parallel::Triangulation</class>
    <class kind="class">parallel::TriangulationSystem</class>
    <namespace>parallel</namespace>
    <member kind="enumeration">
      <type></type>
      <name>InterfaceRefinementCase</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>minus_is_finer</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7a725ee72f9a8d77b879deec24e24a1d7b</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>plus_is_finer</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7aee243748f7ee1b44506873ec0b8ae191</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>equally_fine</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7a37477769ebb40217b8b2dcb406f8b564</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>at_boundary</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7ad453630f99a3f3cec0aaea2d520ac3cf</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>InterfaceSide</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a44f3c00e36c1d6e3c389ae693c09b435</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>minus</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a44f3c00e36c1d6e3c389ae693c09b435aba34d25d1a0e7f60052d7418d1ccf58b</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>plus</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a44f3c00e36c1d6e3c389ae693c09b435a772311abb4290d65c05c49d38c4da6fb</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>two_block_matrix.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>two__block__matrix_8h.html</filename>
    <class kind="class">TwoBlockMatrix</class>
    <class kind="class">parallel::TwoBlockMatrix</class>
    <namespace>parallel</namespace>
  </compound>
  <compound kind="file">
    <name>two_block_sparsity_pattern.h</name>
    <path>/home/sst/code/GalerkinTools/GalerkinTools/include/galerkin_tools/</path>
    <filename>two__block__sparsity__pattern_8h.html</filename>
    <class kind="class">TwoBlockSparsityPattern</class>
  </compound>
  <compound kind="class">
    <name>AssemblyHelper</name>
    <filename>class_assembly_helper.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>AssemblyHelper</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae164448dcc5e9e8d2849c354212d6df5</anchor>
      <arglist>(const TotalPotential&lt; spacedim &gt; &amp;total_potential, TriangulationSystem&lt; spacedim &gt; &amp;tria_system, const Mapping&lt; spacedim, spacedim &gt; &amp;mapping_domain, const Mapping&lt; spacedim-1, spacedim &gt; &amp;mapping_interface, const std::set&lt; const IndependentField&lt; 0, spacedim &gt; * &gt; &amp;independent_scalars=std::set&lt; const IndependentField&lt; 0, spacedim &gt; * &gt;())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AssemblyHelper</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac82eca7b04aedf772028c6ff77245e9b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_initial_fields_vector</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4043c993c8902ad5d3045a13ed42e6f8</anchor>
      <arglist>(VectorType &amp;initial_fields, const AffineConstraints&lt; double &gt; *constraints=nullptr) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>make_dirichlet_constraints</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a40e7eb10c5dbd5358597f38291b90d85</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const std::vector&lt; const DirichletConstraint&lt; spacedim &gt; * &gt; &amp;dirichlet_constraints, const AffineConstraints&lt; double &gt; &amp;constraints_ignore=AffineConstraints&lt; double &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>make_dirichlet_constraints</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5b2d6681755428d4abfe718e54e3c322</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const std::vector&lt; const DirichletConstraint&lt; spacedim &gt; * &gt; &amp;dirichlet_constraints, const std::vector&lt; const PointConstraint&lt; spacedim, spacedim &gt; * &gt; &amp;point_constraints_omega, const std::vector&lt; const PointConstraint&lt; spacedim-1, spacedim &gt; * &gt; &amp;point_constraints_sigma, const std::vector&lt; const PointConstraint&lt; 0, spacedim &gt; * &gt; &amp;point_constraints_C, const AffineConstraints&lt; double &gt; &amp;constraints_ignore=AffineConstraints&lt; double &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>generate_sparsity_pattern_by_simulation</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a60a2aa2aa08149682feca02e458232d4</anchor>
      <arglist>(SparsityPatternType &amp;dsp_K, const AffineConstraints&lt; double &gt; &amp;constraints) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>assemble_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5146802e9d7d4bc5df7fe35909da3e44</anchor>
      <arglist>(const SolutionVectorType &amp;solution, const std::vector&lt; const SolutionVectorType * &gt; solution_ref_sets, const AffineConstraints&lt; double &gt; &amp;constraints, double &amp;potential_value, RHSVectorType &amp;f, MatrixType &amp;K, const std::tuple&lt; bool, bool, bool &gt; requested_quantities=std::make_tuple(true, true, true), std::map&lt; unsigned int, double &gt; *local_solution=nullptr) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_nonprimitive_scalar_functional_values</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>acf6ef2dced66e223684e5df97182f428</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, std::map&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; *, double &gt; &amp;nonprimitive_scalar_functional_values_domain, std::map&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; *, double &gt; &amp;nonprimitive_scalar_functional_values_interface) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>call_scalar_functionals</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a60af8e484b7663e5fe5595a6e0e07693</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; &amp;solution_ref_sets, const std::set&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; * &gt; &amp;scalar_functionals_domain_to_call, const std::set&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; * &gt; &amp;scalar_functionals_interface_to_call, const bool call_all_functionals=false) const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_maximum_step_length</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a310594206df2622027fdc48e84600bf7</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, const VectorType &amp;delta_solution) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compare_derivatives_with_numerical_derivatives</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac54f45f37a38426db1b5f85eccc7b3e9</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, const std::string detailed_printout_file=&quot;&quot;, const double epsilon=1e-8) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>distribute_dofs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af1d2938f5c1a40ace48bc86a0a91eb44</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;renumbering_domain, const std::vector&lt; unsigned int &gt; &amp;renumbering_interface)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const std::string, const std::string &gt;</type>
      <name>write_output_independent_fields</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6120d66724f518dcfcdc30a89df01c23</anchor>
      <arglist>(const VectorType &amp;solution, const std::string file_name_domain, const std::string file_name_interface, const unsigned int file_index=0, const std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt; &amp;dp_domain=std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt;(), const std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt; &amp;dp_interface=std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt;(), const unsigned int n_subdivisions=1) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_assembly_helper_definition</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6a6f8ff7c1a8910d84beb7761b5c821b</anchor>
      <arglist>(const bool detailed_printout_shapefuns=true) const</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const double, const double &gt;</type>
      <name>compute_distance_to_other_solution</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a94bb821b6258eab0bb3a9046b6d9158a</anchor>
      <arglist>(const VectorType &amp;solution, const VectorType &amp;other_solution, const AssemblyHelper&lt; spacedim &gt; &amp;other_assembly_helper, const Quadrature&lt; spacedim &gt; quadrature_domain, const Quadrature&lt; spacedim-1 &gt; quadrature_interface, const VectorTools::NormType norm_type=VectorTools::NormType::L2_norm, const ComponentMask component_mask_domain=ComponentMask(), const ComponentMask component_mask_interface=ComponentMask(), const double exponent=2.0, const Vector&lt; double &gt; scaling_domain=dealii::Vector&lt; double &gt;(), const Vector&lt; double &gt; scaling_interface=dealii::Vector&lt; double &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const double, const double &gt;</type>
      <name>compute_distance_to_exact_solution</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0e96514d9d023949eb07d95b5a2214c4</anchor>
      <arglist>(const VectorType &amp;solution, const Function&lt; spacedim &gt; &amp;exact_solution_domain, const Function&lt; spacedim &gt; &amp;exact_solution_interface, const Quadrature&lt; spacedim &gt; quadrature_domain, const Quadrature&lt; spacedim-1 &gt; quadrature_interface, const VectorTools::NormType norm_type=VectorTools::NormType::L2_norm, const ComponentMask component_mask_domain=ComponentMask(), const ComponentMask component_mask_interface=ComponentMask(), const double exponent=2.0) const</arglist>
    </member>
    <member kind="function">
      <type>const TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>get_triangulation_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a42cc83a6b33fe48b04fa8f4c9907cbb8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>get_triangulation_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a16f9d21a79922d4879e37916b414f7d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>get_dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a095623df46217c89ee8e786f6e8a3034</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>get_dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a04523eef6062ced8c88d4c093b65df3d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; const IndependentField&lt; spacedim, spacedim &gt; *, const unsigned int &gt;</type>
      <name>get_u_omega_global_component_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a39d6fed5b90cee2e2e972e294ececffb</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; const IndependentField&lt; spacedim-1, spacedim &gt; *, const unsigned int &gt;</type>
      <name>get_u_sigma_global_component_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>afe7df3baf877b83b7c98b7f389fa2926</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_u_omega_global_component_index</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af65861aae14b724a631feaf56c82ae9a</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; &amp;u_omega) const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_u_sigma_global_component_index</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4ad9cb7482bfc3ce527ad7639a8d5843</anchor>
      <arglist>(const IndependentField&lt; spacedim-1, spacedim &gt; &amp;u_sigma) const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>system_size</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae821c8ae9c8fa6f8e85d20ecd5ad7e39</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_n_stretched_rows</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a781fcbb9a157621c8db25d8ef46aca13</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_n_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9b9603ede43f9abae845caf60e52d4a1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_global_dof_index_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3b6e5ff3a45411c2b8b42777fa94ec40</anchor>
      <arglist>(const IndependentField&lt; 0, spacedim &gt; *independent_scalar) const</arglist>
    </member>
    <member kind="function">
      <type>const IndexSet</type>
      <name>get_locally_relevant_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a30584e0ed1b2564e9b66ce9cecae40c7</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const IndexSet</type>
      <name>get_locally_owned_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8ecc73fbc0e71716805c97498a83833a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; IndexSet &gt;</type>
      <name>get_locally_relevant_indices_blocks</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5f47fad7f7f2a83a7f54cd38825c703a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; IndexSet &gt;</type>
      <name>get_locally_owned_indices_blocks</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af1a81fd16e7692501189d9c0bd96bc2d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index_at_point_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad798d22d994c4c18a05503ae9db4b408</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; *u_omega, const unsigned int component, const Point&lt; spacedim &gt; p, const std::set&lt; unsigned int &gt; ignore_dofs=std::set&lt; unsigned int &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index_at_point_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a256dc4f0fc5c62f130d75734ea94249f</anchor>
      <arglist>(const IndependentField&lt; spacedim-1, spacedim &gt; *u_sigma, const unsigned int component, const Point&lt; spacedim &gt; p, const std::set&lt; unsigned int &gt; ignore_dofs=std::set&lt; unsigned int &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_dof_information</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab0dabb84cc4a0497dbeb73a9eec3071d</anchor>
      <arglist>(const unsigned int dof_index) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_quadrature_point_alignment_tol</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a60553b183a382aff3072ea7b6f5b86dd</anchor>
      <arglist>(const double quadrature_point_alignment_tol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_cylindrical_symmetry</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5744ed966169526d69508706b7b8f2d2</anchor>
      <arglist>(const bool cylindrical_symmetry)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>convert_dependent_fields_to_shapefunctions</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aad21ced11a2d90c804827854c18f7f89</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_hidden_variables</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1c187fbb8171d6a1ff2ff6344cb454ed</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>distribute_dofs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a953859cc4cd745a1b51fdec5418be682</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_fe_values_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a95639ce05d10d6aa37ea6ae5962753e2</anchor>
      <arglist>(const typename DoFHandler&lt; spacedim, spacedim &gt;::active_cell_iterator &amp;cell, const unsigned int internal_index, const bool nonprimitive=false) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_fe_values_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a223bbebe8a1f0aa09a53ee19257b927e</anchor>
      <arglist>(const InterfaceCellDomainCellsDoF&lt; spacedim &gt; &amp;interface_cell_domain_cells, const unsigned int internal_index, const bool nonprimitive=false) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>compute_e_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5e29275d9ec1c479375707fbc670bc8c</anchor>
      <arglist>(const unsigned int internal_index, const unsigned int scalar_functional_index, const unsigned int q_point, const Vector&lt; double &gt; &amp;solution_u_omega, const Vector&lt; double &gt; &amp;solution_C, Vector&lt; double &gt; &amp;e_omega, FullMatrix&lt; double &gt; &amp;de_omega_dsol_T, const bool compute_derivative=true, const bool ignore_constants=false) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>compute_e_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1906ba80994e8bb136ccf466fd7611d5</anchor>
      <arglist>(const unsigned int internal_index, const unsigned int scalar_functional_index, const unsigned int q_point, const Vector&lt; double &gt; &amp;solution_u_sigma, const Vector&lt; double &gt; &amp;solution_u_omega_minus, const Vector&lt; double &gt; &amp;solution_u_omega_plus, const Vector&lt; double &gt; &amp;solution_C, const std::vector&lt; unsigned int &gt; &amp;dof_indices_interface_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_minus_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_plus_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_C_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_combined, Vector&lt; double &gt; &amp;e_sigma, FullMatrix&lt; double &gt; &amp;de_sigma_dsol_T, const bool compute_derivative=true, const bool ignore_constants=false) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>get_nonprimitive_scalar_functional_values</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>abe88d8ccfd69bbfcfc56551c5c7d67e9</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, Vector&lt; double &gt; &amp;nonprimitive_scalar_functional_values) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>std::pair&lt; const int, const int &gt;</type>
      <name>get_scalar_functional_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a2778924ff66c8ad8695f0cd3da5ced9f</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim, spacedim &gt; *scalar_functional) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>std::pair&lt; const int, const int &gt;</type>
      <name>get_scalar_functional_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3e3f6a06344be172d9419b26bb085073</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim-1, spacedim &gt; *scalar_functional) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>get_dof_indices_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae339631f070dbe766d84697cd9229134</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;global_dof_indices_C) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>make_dirichlet_constraints_recursion</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8a6d7a10dfa2b88ef6f2a974124b0ad5</anchor>
      <arglist>(const typename TriangulationSystem&lt; spacedim &gt;::DomainCell &amp;domain_cell, const unsigned int face, const std::vector&lt; unsigned int &gt; &amp;shapefuns, const DirichletConstraint&lt; spacedim &gt; &amp;constraint, AffineConstraints&lt; double &gt; &amp;constraint_matrix, const AffineConstraints&lt; double &gt; &amp;constraints_ignore) const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const TotalPotential&lt; spacedim &gt;</type>
      <name>total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a748eed9d73b7437a4bf2dcd73108790b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>tria_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>add08a8a7bb9c9325fcc7d92bfce525d4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const Mapping&lt; spacedim, spacedim &gt; &gt;</type>
      <name>mapping_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3fbb49461000dea8f64266f830709fad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const Mapping&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>mapping_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a055fde6217c18e62cd80188d0130c201</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, const unsigned int &gt;</type>
      <name>material_id_to_internal_index_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a10b3acf64bccc169ee14dc2505ce4b46</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; std::tuple&lt; const types::material_id, const types::material_id, const types::material_id &gt;, const unsigned int &gt;</type>
      <name>material_ids_to_internal_index_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a43b82de0ede96d03b9f7fd8740d81668</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>hp::FECollection&lt; spacedim, spacedim &gt;</type>
      <name>fe_collection_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af3803b0aad9853e6bf018c70be41e791</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>hp::FECollection&lt; spacedim-1, spacedim &gt;</type>
      <name>fe_collection_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8b4d224a9ecd2e926a8860829874d2a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, unsigned int &gt;</type>
      <name>material_id_to_fe_system_id_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3045f80801fc31920efd161a268aae8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, unsigned int &gt;</type>
      <name>material_id_to_fe_system_id_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fea54137e3c1c5a514e39c9b2ad7926</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DoFHandlerSystem&lt; spacedim &gt;</type>
      <name>dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a885e660c749e91a35e3279643ebcd87f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt; &gt;</type>
      <name>u_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1952a054a839a7a683ca108013e7d976</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; spacedim-1, spacedim &gt; &gt; &gt;</type>
      <name>u_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a696fe649b3503561235aa1ccbf2ddeef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt; &gt;</type>
      <name>C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa5234a46be82cfe7d92678169d38f326</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_component_indices_u_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6dae4b6ae7934eaec1ad7baff258ce6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_component_indices_u_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a992a53a1fcac8a393ca53fb8d504bdfe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_indices_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9a8f0e8ea8c67ce9429c16a2017cafdc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa6fa619e4c2582e95950e878cd06628e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt; &gt; &gt;</type>
      <name>scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a29aa77e0e8e6b35c94966ea88840e462</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>scalar_functionals_domain_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fe78a019aec03cbeeb336d1d2874729</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>scalar_functionals_interface_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a833383aa6d157157545204143897ed9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_domain_nonprimitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>acf05fab2ddf57769a103d82a4f2d1cd3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_interface_nonprimitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0d15b3ab0c7bec9fc4f40e532f8776f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_domain_primitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4f08790a2235e48ce19f5d8d965a7874</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_interface_primitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad99c75f32cf3f18aa1d4067ad8b56ae8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt; &gt;</type>
      <name>contributions_scalar_functionals_domain_total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aceaf7ba62dfe0fa06ecb15ee8c14da34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt; &gt;</type>
      <name>contributions_scalar_functionals_interface_total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9e76874224ab4946218fdce9bdba0e03</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>n_scalar_functionals_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af7bcfc1db651535a7aefc6071a81e124</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>n_scalar_functionals_primitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af5e03e8e47a85dbc96444ef61525c454</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0bdb6e2e2f9623f3f10dfa2ebe8e234c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>abcad51e64347fc3141d2840a2835b46c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_domain_facewise</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae92560183f1d2060265f0744a84f0349</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1a26b40224e3f04e5168accc91486493</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; std::vector&lt; unsigned int &gt;, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain_local</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a08edb01ba4c8d721862fb060860137d4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; std::vector&lt; unsigned int &gt;, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain_locally_eliminated</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a52d69b234032d5194e87e9b7e5fef759</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::pair&lt; std::vector&lt; unsigned int &gt;, std::vector&lt; unsigned int &gt; &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain_not_local</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>abcef951309d147d892660da4f62532d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::pair&lt; bool, bool &gt; &gt; &gt;</type>
      <name>has_local_locally_eliminated_dofs_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a092a86e3bce3038ac36483338afca487</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::map&lt; unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt;</type>
      <name>correspondence_e_omega_locally_eliminated</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a509326b1521d7979d55dd9b811beac52</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab346e146cf91fb7a0688076551b37355</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a12299d82365553a21fef8529c8fe8a17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af07bb528fdd350e9b467b08dc44a03e7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_C_indices_scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0edd25820c92a25ae87fc240f4916804</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_C_indices_scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a311e176038ee2b7ca0719abb384ca57b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a13fef9096e5fc7b1e922ea78d7aa2c28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fbb532e798c2427af5285c2df10c9f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>c_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a75e6f76c0b12b91c5feb230251f0137f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>d_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad93b109608d4425d318434e01cb6246c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa266cc07e9670319481da52d633d2583</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af58c9a1c7093edc306070913aa1b9be2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4461be378c9be0364ca23153c367d24c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6f51f8b4dfdae385a2e2fe2dd9e66cdb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8eb513a6239d94cc0f1a6a01a037c572</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab09dd07d3ec596525be83f8f5859bef7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>c_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab3d4370661ba726010ee687ef0e98140</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>d_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a48d7d677120eb1c84b4983f470246e02</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; bool &gt; &gt; &gt;</type>
      <name>e_sigma_local</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a16aa77ea3867740a280fd50f5578adfc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a904a24f53b66e1c1ef89f1bb7989eb32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae1e2643696005415d78421882ca80e8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain_reinit</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a28b58551c6afc68c7beaaa2604bc6e92</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface_reinit</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac9789c5a00867744dd906b85580e3091</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain_reinit_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3c592ef0a148753891cc3e03fd08324c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface_reinit_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>afaa20027ee539ca8d9c40c317127e471</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::pair&lt; std::string, unsigned int &gt; &gt;</type>
      <name>component_names_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af5e29448f133863a1859be8bfbb300c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::pair&lt; std::string, unsigned int &gt; &gt;</type>
      <name>component_names_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a7ae6ae2ec356cbb7b830d968315d280c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; boost::signals2::connection &gt;</type>
      <name>tria_listeners</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a228cec028ab5126d25c3ebf0e12a17a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>this_proc</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a2aad83ae1bfe5338794cf9b50848469a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>n_procs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a87945d87baf37637673fd124b3803fd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>ConditionalOStream</type>
      <name>pout</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a717eb6ebc7c62fe00063edcf264f3ecc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>quadrature_point_alignment_tol</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa37920e596dca3985e6d28b9d4e3d882</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>cylindrical_symmetry</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa2548bfa4a097088e3759ce2a1319aa8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const TotalPotential&lt; spacedim &gt;</type>
      <name>total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a748eed9d73b7437a4bf2dcd73108790b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>tria_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>add08a8a7bb9c9325fcc7d92bfce525d4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const Mapping&lt; spacedim, spacedim &gt; &gt;</type>
      <name>mapping_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3fbb49461000dea8f64266f830709fad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const Mapping&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>mapping_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a055fde6217c18e62cd80188d0130c201</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, const unsigned int &gt;</type>
      <name>material_id_to_internal_index_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a10b3acf64bccc169ee14dc2505ce4b46</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; std::tuple&lt; const types::material_id, const types::material_id, const types::material_id &gt;, const unsigned int &gt;</type>
      <name>material_ids_to_internal_index_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a43b82de0ede96d03b9f7fd8740d81668</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>hp::FECollection&lt; spacedim, spacedim &gt;</type>
      <name>fe_collection_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af3803b0aad9853e6bf018c70be41e791</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>hp::FECollection&lt; spacedim-1, spacedim &gt;</type>
      <name>fe_collection_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8b4d224a9ecd2e926a8860829874d2a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, unsigned int &gt;</type>
      <name>material_id_to_fe_system_id_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3045f80801fc31920efd161a268aae8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, unsigned int &gt;</type>
      <name>material_id_to_fe_system_id_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fea54137e3c1c5a514e39c9b2ad7926</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DoFHandlerSystem&lt; spacedim &gt;</type>
      <name>dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a885e660c749e91a35e3279643ebcd87f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt; &gt;</type>
      <name>u_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1952a054a839a7a683ca108013e7d976</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; spacedim-1, spacedim &gt; &gt; &gt;</type>
      <name>u_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a696fe649b3503561235aa1ccbf2ddeef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt; &gt;</type>
      <name>C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa5234a46be82cfe7d92678169d38f326</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_component_indices_u_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6dae4b6ae7934eaec1ad7baff258ce6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_component_indices_u_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a992a53a1fcac8a393ca53fb8d504bdfe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_indices_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9a8f0e8ea8c67ce9429c16a2017cafdc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa6fa619e4c2582e95950e878cd06628e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt; &gt; &gt;</type>
      <name>scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a29aa77e0e8e6b35c94966ea88840e462</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>scalar_functionals_domain_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fe78a019aec03cbeeb336d1d2874729</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>scalar_functionals_interface_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a833383aa6d157157545204143897ed9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_domain_nonprimitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>acf05fab2ddf57769a103d82a4f2d1cd3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_interface_nonprimitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0d15b3ab0c7bec9fc4f40e532f8776f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_domain_primitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4f08790a2235e48ce19f5d8d965a7874</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_interface_primitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad99c75f32cf3f18aa1d4067ad8b56ae8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt; &gt;</type>
      <name>contributions_scalar_functionals_domain_total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aceaf7ba62dfe0fa06ecb15ee8c14da34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt; &gt;</type>
      <name>contributions_scalar_functionals_interface_total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9e76874224ab4946218fdce9bdba0e03</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>n_scalar_functionals_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af7bcfc1db651535a7aefc6071a81e124</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>n_scalar_functionals_primitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af5e03e8e47a85dbc96444ef61525c454</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0bdb6e2e2f9623f3f10dfa2ebe8e234c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>abcad51e64347fc3141d2840a2835b46c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_domain_facewise</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae92560183f1d2060265f0744a84f0349</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1a26b40224e3f04e5168accc91486493</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; std::vector&lt; unsigned int &gt;, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain_local</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a08edb01ba4c8d721862fb060860137d4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::pair&lt; std::vector&lt; unsigned int &gt;, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain_locally_eliminated</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a52d69b234032d5194e87e9b7e5fef759</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::pair&lt; std::vector&lt; unsigned int &gt;, std::vector&lt; unsigned int &gt; &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain_not_local</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>abcef951309d147d892660da4f62532d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::pair&lt; bool, bool &gt; &gt; &gt;</type>
      <name>has_local_locally_eliminated_dofs_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a092a86e3bce3038ac36483338afca487</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::map&lt; unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt;</type>
      <name>correspondence_e_omega_locally_eliminated</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a509326b1521d7979d55dd9b811beac52</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab346e146cf91fb7a0688076551b37355</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a12299d82365553a21fef8529c8fe8a17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af07bb528fdd350e9b467b08dc44a03e7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_C_indices_scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0edd25820c92a25ae87fc240f4916804</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_C_indices_scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a311e176038ee2b7ca0719abb384ca57b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a13fef9096e5fc7b1e922ea78d7aa2c28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fbb532e798c2427af5285c2df10c9f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>c_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a75e6f76c0b12b91c5feb230251f0137f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>d_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad93b109608d4425d318434e01cb6246c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa266cc07e9670319481da52d633d2583</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af58c9a1c7093edc306070913aa1b9be2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4461be378c9be0364ca23153c367d24c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6f51f8b4dfdae385a2e2fe2dd9e66cdb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8eb513a6239d94cc0f1a6a01a037c572</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab09dd07d3ec596525be83f8f5859bef7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>c_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab3d4370661ba726010ee687ef0e98140</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>d_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a48d7d677120eb1c84b4983f470246e02</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; bool &gt; &gt; &gt;</type>
      <name>e_sigma_local</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a16aa77ea3867740a280fd50f5578adfc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a904a24f53b66e1c1ef89f1bb7989eb32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae1e2643696005415d78421882ca80e8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain_reinit</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a28b58551c6afc68c7beaaa2604bc6e92</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface_reinit</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac9789c5a00867744dd906b85580e3091</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain_reinit_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3c592ef0a148753891cc3e03fd08324c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface_reinit_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>afaa20027ee539ca8d9c40c317127e471</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::pair&lt; std::string, unsigned int &gt; &gt;</type>
      <name>component_names_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af5e29448f133863a1859be8bfbb300c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::pair&lt; std::string, unsigned int &gt; &gt;</type>
      <name>component_names_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a7ae6ae2ec356cbb7b830d968315d280c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; boost::signals2::connection &gt;</type>
      <name>tria_listeners</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a228cec028ab5126d25c3ebf0e12a17a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>this_proc</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a2aad83ae1bfe5338794cf9b50848469a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>n_procs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a87945d87baf37637673fd124b3803fd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>ConditionalOStream</type>
      <name>pout</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a717eb6ebc7c62fe00063edcf264f3ecc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>quadrature_point_alignment_tol</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa37920e596dca3985e6d28b9d4e3d882</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>cylindrical_symmetry</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa2548bfa4a097088e3759ce2a1319aa8</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>convert_dependent_fields_to_shapefunctions</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aad21ced11a2d90c804827854c18f7f89</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_hidden_variables</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1c187fbb8171d6a1ff2ff6344cb454ed</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>distribute_dofs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a953859cc4cd745a1b51fdec5418be682</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_fe_values_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a95639ce05d10d6aa37ea6ae5962753e2</anchor>
      <arglist>(const typename DoFHandler&lt; spacedim, spacedim &gt;::active_cell_iterator &amp;cell, const unsigned int internal_index, const bool nonprimitive=false) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_fe_values_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a223bbebe8a1f0aa09a53ee19257b927e</anchor>
      <arglist>(const InterfaceCellDomainCellsDoF&lt; spacedim &gt; &amp;interface_cell_domain_cells, const unsigned int internal_index, const bool nonprimitive=false) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>compute_e_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5e29275d9ec1c479375707fbc670bc8c</anchor>
      <arglist>(const unsigned int internal_index, const unsigned int scalar_functional_index, const unsigned int q_point, const Vector&lt; double &gt; &amp;solution_u_omega, const Vector&lt; double &gt; &amp;solution_C, Vector&lt; double &gt; &amp;e_omega, FullMatrix&lt; double &gt; &amp;de_omega_dsol_T, const bool compute_derivative=true, const bool ignore_constants=false) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>compute_e_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1906ba80994e8bb136ccf466fd7611d5</anchor>
      <arglist>(const unsigned int internal_index, const unsigned int scalar_functional_index, const unsigned int q_point, const Vector&lt; double &gt; &amp;solution_u_sigma, const Vector&lt; double &gt; &amp;solution_u_omega_minus, const Vector&lt; double &gt; &amp;solution_u_omega_plus, const Vector&lt; double &gt; &amp;solution_C, const std::vector&lt; unsigned int &gt; &amp;dof_indices_interface_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_minus_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_plus_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_C_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_combined, Vector&lt; double &gt; &amp;e_sigma, FullMatrix&lt; double &gt; &amp;de_sigma_dsol_T, const bool compute_derivative=true, const bool ignore_constants=false) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>get_nonprimitive_scalar_functional_values</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>abe88d8ccfd69bbfcfc56551c5c7d67e9</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, Vector&lt; double &gt; &amp;nonprimitive_scalar_functional_values) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>std::pair&lt; const int, const int &gt;</type>
      <name>get_scalar_functional_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a2778924ff66c8ad8695f0cd3da5ced9f</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim, spacedim &gt; *scalar_functional) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>std::pair&lt; const int, const int &gt;</type>
      <name>get_scalar_functional_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3e3f6a06344be172d9419b26bb085073</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim-1, spacedim &gt; *scalar_functional) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>get_dof_indices_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae339631f070dbe766d84697cd9229134</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;global_dof_indices_C) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>make_dirichlet_constraints_recursion</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8a6d7a10dfa2b88ef6f2a974124b0ad5</anchor>
      <arglist>(const typename TriangulationSystem&lt; spacedim &gt;::DomainCell &amp;domain_cell, const unsigned int face, const std::vector&lt; unsigned int &gt; &amp;shapefuns, const DirichletConstraint&lt; spacedim &gt; &amp;constraint, AffineConstraints&lt; double &gt; &amp;constraint_matrix, const AffineConstraints&lt; double &gt; &amp;constraints_ignore) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_initial_fields_vector</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4043c993c8902ad5d3045a13ed42e6f8</anchor>
      <arglist>(VectorType &amp;initial_fields, const AffineConstraints&lt; double &gt; *constraints=nullptr) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>make_dirichlet_constraints</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a40e7eb10c5dbd5358597f38291b90d85</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const std::vector&lt; const DirichletConstraint&lt; spacedim &gt; * &gt; &amp;dirichlet_constraints, const AffineConstraints&lt; double &gt; &amp;constraints_ignore=AffineConstraints&lt; double &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>make_dirichlet_constraints</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5b2d6681755428d4abfe718e54e3c322</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const std::vector&lt; const DirichletConstraint&lt; spacedim &gt; * &gt; &amp;dirichlet_constraints, const std::vector&lt; const PointConstraint&lt; spacedim, spacedim &gt; * &gt; &amp;point_constraints_omega, const std::vector&lt; const PointConstraint&lt; spacedim-1, spacedim &gt; * &gt; &amp;point_constraints_sigma, const std::vector&lt; const PointConstraint&lt; 0, spacedim &gt; * &gt; &amp;point_constraints_C, const AffineConstraints&lt; double &gt; &amp;constraints_ignore=AffineConstraints&lt; double &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>generate_sparsity_pattern_by_simulation</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a60a2aa2aa08149682feca02e458232d4</anchor>
      <arglist>(SparsityPatternType &amp;dsp_K, const AffineConstraints&lt; double &gt; &amp;constraints) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>assemble_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5146802e9d7d4bc5df7fe35909da3e44</anchor>
      <arglist>(const SolutionVectorType &amp;solution, const std::vector&lt; const SolutionVectorType * &gt; solution_ref_sets, const AffineConstraints&lt; double &gt; &amp;constraints, double &amp;potential_value, RHSVectorType &amp;f, MatrixType &amp;K, const std::tuple&lt; bool, bool, bool &gt; requested_quantities=std::make_tuple(true, true, true), std::map&lt; unsigned int, double &gt; *local_solution=nullptr) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_nonprimitive_scalar_functional_values</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>acf6ef2dced66e223684e5df97182f428</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, std::map&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; *, double &gt; &amp;nonprimitive_scalar_functional_values_domain, std::map&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; *, double &gt; &amp;nonprimitive_scalar_functional_values_interface) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>call_scalar_functionals</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a60af8e484b7663e5fe5595a6e0e07693</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; &amp;solution_ref_sets, const std::set&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; * &gt; &amp;scalar_functionals_domain_to_call, const std::set&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; * &gt; &amp;scalar_functionals_interface_to_call, const bool call_all_functionals=false) const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_maximum_step_length</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a310594206df2622027fdc48e84600bf7</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, const VectorType &amp;delta_solution) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compare_derivatives_with_numerical_derivatives</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac54f45f37a38426db1b5f85eccc7b3e9</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, const std::string detailed_printout_file=&quot;&quot;, const double epsilon=1e-8) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>distribute_dofs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af1d2938f5c1a40ace48bc86a0a91eb44</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;renumbering_domain, const std::vector&lt; unsigned int &gt; &amp;renumbering_interface)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const std::string, const std::string &gt;</type>
      <name>write_output_independent_fields</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6120d66724f518dcfcdc30a89df01c23</anchor>
      <arglist>(const VectorType &amp;solution, const std::string file_name_domain, const std::string file_name_interface, const unsigned int file_index=0, const std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt; &amp;dp_domain=std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt;(), const std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt; &amp;dp_interface=std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt;(), const unsigned int n_subdivisions=1) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_assembly_helper_definition</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6a6f8ff7c1a8910d84beb7761b5c821b</anchor>
      <arglist>(const bool detailed_printout_shapefuns=true) const</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const double, const double &gt;</type>
      <name>compute_distance_to_other_solution</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a94bb821b6258eab0bb3a9046b6d9158a</anchor>
      <arglist>(const VectorType &amp;solution, const VectorType &amp;other_solution, const AssemblyHelper&lt; spacedim &gt; &amp;other_assembly_helper, const Quadrature&lt; spacedim &gt; quadrature_domain, const Quadrature&lt; spacedim-1 &gt; quadrature_interface, const VectorTools::NormType norm_type=VectorTools::NormType::L2_norm, const ComponentMask component_mask_domain=ComponentMask(), const ComponentMask component_mask_interface=ComponentMask(), const double exponent=2.0, const Vector&lt; double &gt; scaling_domain=dealii::Vector&lt; double &gt;(), const Vector&lt; double &gt; scaling_interface=dealii::Vector&lt; double &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const double, const double &gt;</type>
      <name>compute_distance_to_exact_solution</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0e96514d9d023949eb07d95b5a2214c4</anchor>
      <arglist>(const VectorType &amp;solution, const Function&lt; spacedim &gt; &amp;exact_solution_domain, const Function&lt; spacedim &gt; &amp;exact_solution_interface, const Quadrature&lt; spacedim &gt; quadrature_domain, const Quadrature&lt; spacedim-1 &gt; quadrature_interface, const VectorTools::NormType norm_type=VectorTools::NormType::L2_norm, const ComponentMask component_mask_domain=ComponentMask(), const ComponentMask component_mask_interface=ComponentMask(), const double exponent=2.0) const</arglist>
    </member>
    <member kind="function">
      <type>const TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>get_triangulation_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a42cc83a6b33fe48b04fa8f4c9907cbb8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>get_triangulation_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a16f9d21a79922d4879e37916b414f7d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>get_dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a095623df46217c89ee8e786f6e8a3034</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>get_dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a04523eef6062ced8c88d4c093b65df3d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; const IndependentField&lt; spacedim, spacedim &gt; *, const unsigned int &gt;</type>
      <name>get_u_omega_global_component_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a39d6fed5b90cee2e2e972e294ececffb</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; const IndependentField&lt; spacedim-1, spacedim &gt; *, const unsigned int &gt;</type>
      <name>get_u_sigma_global_component_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>afe7df3baf877b83b7c98b7f389fa2926</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_u_omega_global_component_index</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af65861aae14b724a631feaf56c82ae9a</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; &amp;u_omega) const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_u_sigma_global_component_index</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4ad9cb7482bfc3ce527ad7639a8d5843</anchor>
      <arglist>(const IndependentField&lt; spacedim-1, spacedim &gt; &amp;u_sigma) const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>system_size</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae821c8ae9c8fa6f8e85d20ecd5ad7e39</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_n_stretched_rows</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a781fcbb9a157621c8db25d8ef46aca13</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_n_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9b9603ede43f9abae845caf60e52d4a1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_global_dof_index_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3b6e5ff3a45411c2b8b42777fa94ec40</anchor>
      <arglist>(const IndependentField&lt; 0, spacedim &gt; *independent_scalar) const</arglist>
    </member>
    <member kind="function">
      <type>const IndexSet</type>
      <name>get_locally_relevant_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a30584e0ed1b2564e9b66ce9cecae40c7</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const IndexSet</type>
      <name>get_locally_owned_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8ecc73fbc0e71716805c97498a83833a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; IndexSet &gt;</type>
      <name>get_locally_relevant_indices_blocks</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5f47fad7f7f2a83a7f54cd38825c703a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; IndexSet &gt;</type>
      <name>get_locally_owned_indices_blocks</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af1a81fd16e7692501189d9c0bd96bc2d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index_at_point_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad798d22d994c4c18a05503ae9db4b408</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; *u_omega, const unsigned int component, const Point&lt; spacedim &gt; p, const std::set&lt; unsigned int &gt; ignore_dofs=std::set&lt; unsigned int &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index_at_point_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a256dc4f0fc5c62f130d75734ea94249f</anchor>
      <arglist>(const IndependentField&lt; spacedim-1, spacedim &gt; *u_sigma, const unsigned int component, const Point&lt; spacedim &gt; p, const std::set&lt; unsigned int &gt; ignore_dofs=std::set&lt; unsigned int &gt;()) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_dof_information</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab0dabb84cc4a0497dbeb73a9eec3071d</anchor>
      <arglist>(const unsigned int dof_index) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_quadrature_point_alignment_tol</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a60553b183a382aff3072ea7b6f5b86dd</anchor>
      <arglist>(const double quadrature_point_alignment_tol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_cylindrical_symmetry</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5744ed966169526d69508706b7b8f2d2</anchor>
      <arglist>(const bool cylindrical_symmetry)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BlockSolverWrapperMA57</name>
    <filename>class_block_solver_wrapper_m_a57.html</filename>
    <base>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::BlockVector&lt; double &gt;, TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt;, TwoBlockSparsityPattern &gt;</base>
    <member kind="function">
      <type></type>
      <name>~BlockSolverWrapperMA57</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>ab21743d58c41e77d7d4689a855c57886</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a33aa6e47fb423580d0cf396e25ec4d88</anchor>
      <arglist>(const TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt;&gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::BlockVector&lt; double &gt; &amp;f_stretched, const bool symmetric=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DeclException2</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a2bdc8ae152e8af80c7e7df0ee5569415</anchor>
      <arglist>(ExcMA57Error, std::string, int,&lt;&lt; &quot;MA57 routine &quot;&lt;&lt; arg1&lt;&lt; &quot; returned error status &quot;&lt;&lt; arg2&lt;&lt; &quot;.&quot;)</arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>analyze</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a7803078dafce5d8f42a998441da3633a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>print_level</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a9a764155211658a62882671edc818d3b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>ordering_method</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a46bcec2dd946ca4ca05a149d80cde85c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>use_iterative_refinement</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>afb7094d6e048412f94d7ec7fdabf64f5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>ignore_zeros</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>af07ee5941bf8018886133978ac1d9ad2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_matrix</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>adcf3745eb9e48cdeec740bb90d56928e</anchor>
      <arglist>(const SparseMatrix&lt; double &gt; &amp;matrix)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>analyze_matrix</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>aea47ef899c6123ce5434012e119358c9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>factorize_matrix</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a982634cc81fac311f02bd49b6c88c4b3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>vmult</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>ab4fbc2f03282226b7024f657e4ca5cb8</anchor>
      <arglist>(Vector&lt; double &gt; &amp;x, const Vector&lt; double &gt; &amp;f)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>CNTL</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>ad5fb9927a34b0943f0a2b82d7ebef1c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>ICNTL</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>acb7ffdcbfbd8786630359790139d08ad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>N</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a33c87fcd2bf1ffadc5e98b82270f1f44</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>NE</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a26814446ef196f3328196785954cf450</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>IRN</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>ab72348911e0a660e186214e1f93a0057</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>JCN</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>ad032db448791365dbe27f0b3eebc7791</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>A</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>abad6583801f31d2735e179450e5bb119</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>KEEP</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a946283faa9099c589f44483c02cb011f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>LKEEP</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a1a9a0911396fe053e434e4814eda5175</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>IWORK</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a051b078fdf533ec59d3439cdbdb5dcb4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>WORK</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>ae973be23640c5e51ab7d60b3d49379a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>LWORK</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a848e0d174d9b7476b1a19124855be8ad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>INFO</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>af6102b412f012d6526ac25ef9cb87a78</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>RINFO</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a15d12be227fda3cc4580a7f64d48001e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>FACT</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a2c3ae4d295d83a1ca5a92690eff572be</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>LFACT</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a8acd2b50bf6374917b6aa8de1cbc3069</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>IFACT</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a3229abaff499ee2203afea8305ae9b47</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>LIFACT</name>
      <anchorfile>class_block_solver_wrapper_m_a57.html</anchorfile>
      <anchor>a3ff8ca3369a05fb896a69a276af397ce</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BlockSolverWrapperMUMPS</name>
    <filename>class_block_solver_wrapper_m_u_m_p_s.html</filename>
    <base>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::BlockVector&lt; double &gt;, TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt;, TwoBlockSparsityPattern &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initialize_matrix</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>abb4b34a355239e02e36ac9f6b1928c7e</anchor>
      <arglist>(const SparseMatrix&lt; double &gt; &amp;matrix)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initialize_matrix</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>aadddfe6082fe3d3321bc29ef31c1640a</anchor>
      <arglist>(std::vector&lt; int &gt; &amp;irn, std::vector&lt; int &gt; &amp;jcn, std::vector&lt; double &gt; &amp;A, unsigned int n)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>analyze_matrix</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a9eb5cbca0cf547d391332c379c3acea9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>factorize_matrix</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a39b279f0e0f1c0ecb30702763719d433</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>vmult</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a5e1a3f5918cb70bf9c73391c35b8cf9a</anchor>
      <arglist>(Vector&lt; double &gt; &amp;x, const Vector&lt; double &gt; &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BlockSolverWrapperMUMPS</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>acce2c39c5e2f289db6b85343a8ab2bdb</anchor>
      <arglist>(int sym=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BlockSolverWrapperMUMPS</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a2d2be6d9c412fc9fdc8a76c116c67104</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ae8fdb5cb19687c9e09d4328da66b7d97</anchor>
      <arglist>(const TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt;&gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::BlockVector&lt; double &gt; &amp;f_stretched, const bool symmetric=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DeclException2</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a9568838820205928af55c2a03aecd263</anchor>
      <arglist>(MUMPSError, std::string, int,&lt;&lt; &quot;MUMPS routine &quot;&lt;&lt; arg1&lt;&lt; &quot; returned error status &quot;&lt;&lt; arg2&lt;&lt; &quot;.&quot;)</arglist>
    </member>
    <member kind="variable">
      <type>DMUMPS_STRUC_C</type>
      <name>id</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ad5e82030212e5207e08d9395f015a27c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>icntl</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ad83b47666c9b0979ac729cd402005703</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>cntl</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a51302d5ca5d87e7209c7d8a3913191a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>info</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ad890f0e1ce9499cf91ae9a8d6b89b0b7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>infog</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a59f1340cb4c1f1d5baefb9b81aef6857</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>analyze</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a8be8b4fb9d6a6ebf7184327b5121f1e5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>modify_on_negative_pivot</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ab22891ea06bef2419ef6190aa99e7bf3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>beta</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a04279ba9065c9f59bd7a30c92c7804ac</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>increase_tau</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a120c310cd94f46b4186c56b37170cc25</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>irn</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ab2c3ba49a45e5570cc498e5b91ae9050</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>jcn</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a7573f600b4965a19b556651512d3b26b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>d</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>af481ade4d2832ade039399a08959f91f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>A</name>
      <anchorfile>class_block_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a1cc5cb669452e8c75033caa3a2a869fc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BlockSolverWrapperPARDISO</name>
    <filename>class_block_solver_wrapper_p_a_r_d_i_s_o.html</filename>
    <base>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::BlockVector&lt; double &gt;, TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt;, TwoBlockSparsityPattern &gt;</base>
    <member kind="function">
      <type></type>
      <name>~BlockSolverWrapperPARDISO</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a22f5bc2196681e319b323549a457178a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a73135d179e0618be550b9089c3c70689</anchor>
      <arglist>(const TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt;&gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::BlockVector&lt; double &gt; &amp;f_stretched, const bool symmetric=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DeclException2</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a8bb0292df53f83055f387c49f060499f</anchor>
      <arglist>(ExcPARDISOError, std::string, int,&lt;&lt; &quot;PARDISO routine &quot;&lt;&lt; arg1&lt;&lt; &quot; returned error status &quot;&lt;&lt; arg2&lt;&lt; &quot;.&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DeclException1</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>add4fbb653224be1348463d9bff990a9f</anchor>
      <arglist>(ExcPARDISORes, double,&lt;&lt; &quot;PARDISO residual too large: res = &quot;&lt;&lt; arg1&lt;&lt; &quot;.&quot;)</arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>analyze</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a3df0818f957b743717ce4bfd463ff8dd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>print_level</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>aa87e38dcf8e35852fbb47b4792b2540d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>matrix_type</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>adc6d3f7f878c578aff518cb907e157ff</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>ordering_method</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a3bf4e0ea8332a4e6c71eac84f385663c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>apply_scaling</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>aba68acd666b7c5c41af207eb77db201c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>iterative_scaling_tolerance</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>ae11f723031b3f4f4e9038964a2e010bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>pivoting_method</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a2924ded1d130cdd30e865ff664d8dcf7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>n_iterative_refinements</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>afdb2039acda86d2def4e0eee33635dab</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>pivot_perturbation</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a5fa1782457ce43278374cc1d5395b1d6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>user_perm</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a6146a36c55b7c7dcb0645917f92f518b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>res_max</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a32db83d42ee6bafd25ff6b01cd5452a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>use_defaults</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>ad8bf45c0568b751a2773a57a5b00ba5d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; int &gt; *</type>
      <name>perm</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a13bbf006d1fdf50f06535d9d7ab7434f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_matrix</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a73460ba651207f3b8e86699f58db896f</anchor>
      <arglist>(const SparseMatrix&lt; double &gt; &amp;matrix)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>analyze_matrix</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a53ca2b0fda79f7c29ab46b32469d058a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>factorize_matrix</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>accbe586b87cfd3d453522fee240ff37b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>vmult</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a1dfbc3e43736a1225caa0ab91a30b80e</anchor>
      <arglist>(Vector&lt; double &gt; &amp;x, const Vector&lt; double &gt; &amp;f)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>get_matrix_type</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a3ec780cd305a0fe8326a6acb5c9283e9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>scale_system</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a39b2c196237df4d28c3b806eb6d672a5</anchor>
      <arglist>(TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt;&gt; *K_stretched, dealii::BlockVector&lt; double &gt; *f_stretched)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>scale_solution</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>adfd5663c48d37c37005935176012b334</anchor>
      <arglist>(dealii::Vector&lt; double &gt; &amp;solution, TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt;&gt; *K_stretched, dealii::BlockVector&lt; double &gt; *f_stretched) const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>initialized</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a4bb7e16d07967bda228764b271fcc9d9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>iparm</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a61c49eda443075c200c9583f982ca447</anchor>
      <arglist>[64]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>dparm</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>aa158c35f9486dc043bf3485da3dec8d4</anchor>
      <arglist>[64]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>void *</type>
      <name>pt</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a4d1d5d9882ed7b4019b7412b0bbd125c</anchor>
      <arglist>[64]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>N</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a0218d17deabfb2a7b8905b90ccb6ef46</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>Ap</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>ad284b0a41b7e8f6ab2f58308dfad7c13</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>Ai</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a20249115ffb23d1c0d35c22fa76a1257</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>Ax</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>ad32493c27639c20286aea67e702f3442</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>dealii::Vector&lt; double &gt;</type>
      <name>C_inv</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>adbd0f34a5205de9e8e0c5eec7540d093</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>dealii::Vector&lt; double &gt;</type>
      <name>R_inv</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a68470095b86444392c91a1c9d6a6fa86</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>maxfct</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a3d6d894f09dcdeae1295cf5072b79451</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>nrhs</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a84ae552a3b77bb009366379148119ef8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>mnum</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>ae7c0c78a1167861fcce71bc86bddba19</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>msglvl</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a547f7dbd1fa0da9ce546a707d4e19912</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Vector&lt; double &gt;</type>
      <name>res</name>
      <anchorfile>class_block_solver_wrapper_p_a_r_d_i_s_o.html</anchorfile>
      <anchor>a7d2b9809d10bde8d8440a607522e383c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BlockSolverWrapperUMFPACK</name>
    <filename>class_block_solver_wrapper_u_m_f_p_a_c_k.html</filename>
    <base>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::BlockVector&lt; double &gt;, TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt;, TwoBlockSparsityPattern &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k.html</anchorfile>
      <anchor>a0d27521861712a3b36cc15845e093df2</anchor>
      <arglist>(const TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt;&gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::BlockVector&lt; double &gt; &amp;f_stretched, const bool symmetric=false)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BlockSolverWrapperUMFPACK2</name>
    <filename>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</filename>
    <base>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::BlockVector&lt; double &gt;, TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt;, TwoBlockSparsityPattern &gt;</base>
    <member kind="function">
      <type></type>
      <name>~BlockSolverWrapperUMFPACK2</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a8f298590602e05af072f1ce4178199f2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a651cb9299e8f368ba7294367260f9759</anchor>
      <arglist>(const TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt;&gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::BlockVector&lt; double &gt; &amp;f_stretched, const bool symmetric=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DeclException2</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>aa454139937c8f01d8e893a7e09d3b1a2</anchor>
      <arglist>(ExcUMFPACKError, std::string, int,&lt;&lt; &quot;UMFPACK routine &quot;&lt;&lt; arg1&lt;&lt; &quot; returned error status &quot;&lt;&lt; arg2&lt;&lt; &quot;.&quot;)</arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>analyze</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a5869a4301c4e7a3ae64aa8998131d8c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>print_level</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>af78261735f66012f531c3a71c9640d2c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_matrix</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a8df429ef59e02922971ed01dd4780cee</anchor>
      <arglist>(const SparseMatrix&lt; double &gt; &amp;matrix)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>analyze_matrix</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>afc9d5243e71aa81bd391e36976ef5623</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>factorize_matrix</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a1fb479f8089c3b71b2632f0267e6e7ff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>vmult</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>ac308e6181b48e7e35f533d90a89773f3</anchor>
      <arglist>(Vector&lt; double &gt; &amp;x, const Vector&lt; double &gt; &amp;f)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>control</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a295fca0b28991ae43c7bdfb72b5cdcb2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>info</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a8dfbe8f5e1f1cf26460ba584b4eb02b6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>N</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a85d397826a36330bf3602d510acdd8d4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>void *</type>
      <name>symbolic_decomposition</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a399387c1717404d92ed721fc31767f55</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>void *</type>
      <name>numeric_decomposition</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>ac3a439162324f36f5162439c075d02eb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SuiteSparse_long &gt;</type>
      <name>Ap</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a18e53f152f97e51caab8dcf3efd5eaf1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SuiteSparse_long &gt;</type>
      <name>Ai</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a088fc386567b26c63b8d1a25a9319c75</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>Ax</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k2.html</anchorfile>
      <anchor>a57a4f7722dbbdd1ec2efa416f1cbb47c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DependentField</name>
    <filename>class_dependent_field.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>DependentField</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a5a47b8f5b0e51b29bf6f3ca6e18fda9c</anchor>
      <arglist>(const std::string name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a0726cab14197769ad07bd8f071e395ea</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim, spacedim &gt; &amp;independent_field, const unsigned int component=0, const std::vector&lt; unsigned int &gt; derivatives=std::vector&lt; unsigned int &gt;())</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a852ed477f6ee71ad85a6690d4a5e56fe</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim+1, spacedim &gt; &amp;independent_field, const unsigned int component, const std::vector&lt; unsigned int &gt; derivatives, const InterfaceSide side)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>acb8979c891a45ed7326eb47df41c0129</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim+1, spacedim &gt; &amp;independent_field, const unsigned int component, const InterfaceSide side)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a1aaaf072334d4e06640383ae81766a39</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim, spacedim &gt; &amp;independent_field, const unsigned int component, const unsigned int derivative)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a43de305294f94c38c88fa4768b3f0155</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim+1, spacedim &gt; &amp;independent_field, const unsigned int component, const unsigned int derivative, const InterfaceSide side)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>adf9e4863ca954a9fa70c57aef68849f8</anchor>
      <arglist>(double coefficient, const IndependentField&lt; 0, spacedim &gt; &amp;independent_field)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a72196bff18c2c28a1e811877e739f4ee</anchor>
      <arglist>(double constant)</arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; dim, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; dim, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_interface</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a01793fd56b83754bec50e9d4e0c94cae</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; dim+1, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; dim+1, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_neighbor</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>aa88f72409a02de7505ff5ac5d694ed36</anchor>
      <arglist>(const InterfaceSide side) const</arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; 0, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; 0, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_independent_scalars</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a89fed0663944d7f392ddd4b75bcb60a8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_constant</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a213c5de7c4726eb823590de0510f1ff1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; dim, spacedim &gt; * &gt;</type>
      <name>get_independent_fields_interface</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a33760e356c5704614ed8cf7c532b4de6</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; dim+1, spacedim &gt; * &gt;</type>
      <name>get_independent_fields_neighbors</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>ac12fccee9603f2cf36990af16d8da604</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; 0, spacedim &gt; * &gt;</type>
      <name>get_independent_scalars</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a15eb8101d33f2f09ada4cdf2cde73ad7</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>adf0ad81e2d7b609b7dfde60a77f92a20</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_is_local</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>aabc31c47b3f6cb69da28758d9672d8d8</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_is_locally_eliminated</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>aa16f1559eca5060b0877dc3ec00cc552</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator&lt;</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a1c56c4e18d320e9c402f3adde546c681</anchor>
      <arglist>(const DependentField &amp;dependent_field_2) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a808b18c7b0fc7c05204aed47b1509ce4</anchor>
      <arglist>(const DependentField &amp;dependent_field_2) const</arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a698b255f131279d1edb97b3c525153a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; dim, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; dim, spacedim &gt; &gt;</type>
      <name>terms_interface</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a463f490ef397b9b65112c94932afcd5c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; dim+1, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; dim+1, spacedim &gt; &gt;</type>
      <name>terms_neighbor_plus</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a3f7812cd53686afce899438ae999999a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; dim+1, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; dim+1, spacedim &gt; &gt;</type>
      <name>terms_neighbor_minus</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a406b6490dcd17d80737500184131a3e7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; 0, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; 0, spacedim &gt; &gt;</type>
      <name>terms_independent_scalars</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a642147ff57d5412554508f293182c649</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>constant</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a32b37c78e04a16b6b606442f156c8ca9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>is_local</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>aacdabff2601b1f0efcdce7a4413b47bc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>is_locally_eliminated</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a9c41e749b67a39b68a37ad9d893fb14d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DependentField&lt; spacedim, spacedim &gt;</name>
    <filename>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>DependentField</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a3b2090f3b7c66d84dbb41ec5e677ecfe</anchor>
      <arglist>(const std::string name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a9e0b93f36ef8bc10ee9028576db57390</anchor>
      <arglist>(double coefficient, const IndependentField&lt; spacedim, spacedim &gt; &amp;independent_field, const unsigned int component=0, const std::vector&lt; unsigned int &gt; derivatives=std::vector&lt; unsigned int &gt;())</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a79f3b134f21759492662b2f7173fff8e</anchor>
      <arglist>(double coefficient, const IndependentField&lt; spacedim, spacedim &gt; &amp;independent_field, const unsigned int component, const unsigned int derivative)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a25b4485644ccc8e2aeafddceadb33950</anchor>
      <arglist>(double coefficient, const IndependentField&lt; 0, spacedim &gt; &amp;independent_field)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a24fe77c93937f9563b3dac9e4520fd9b</anchor>
      <arglist>(double constant)</arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; spacedim, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; spacedim, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_domain</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a953203ce2a0a0e5f9463184f43b19229</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; 0, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; 0, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_independent_scalars</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a93ca298c0f912330b2535e886d6ccdf0</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_constant</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ae872945005ea9f660edf2a14edea7702</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; spacedim, spacedim &gt; * &gt;</type>
      <name>get_independent_fields_domain</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ade98e0e1e1a8ebc3f8ad1ef7cf3268fa</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; 0, spacedim &gt; * &gt;</type>
      <name>get_independent_scalars</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a9e441dd97513d30309aebff87f65c3af</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ad08ab17b8ea5d0bc07a0c6755d43ab10</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_is_local</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>abf60a16e3b158c7c67c5e04d8ad41b5c</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_is_locally_eliminated</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a52de8c04f07a31cdd16c30344e953d11</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator&lt;</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a2897cd7975ebb9fa67a53c491fc8a7d8</anchor>
      <arglist>(const DependentField&lt; spacedim, spacedim &gt; &amp;dependent_field_2) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a83626491df9b2bbd6f7d2f3587cbf7e0</anchor>
      <arglist>(const DependentField&lt; spacedim, spacedim &gt; &amp;dependent_field_2) const</arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a99a47f4c10f0e472dbdc2b0945712e23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; spacedim, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; spacedim, spacedim &gt; &gt;</type>
      <name>terms_domain</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a2fc7920fe455597e911f68c56b8d4dbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; 0, spacedim &gt;, DependentFieldComparatorWithoutCoefficient&lt; 0, spacedim &gt; &gt;</type>
      <name>terms_independent_scalars</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a3a6a89ad2fe38606f09573e6e1450879</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>constant</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a24efc9c0928896be871908f050c406fc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>is_local</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ae43bf948e8f45545dd7e744689b47d5b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>is_locally_eliminated</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ae7674c26813012938b30d350d37b7722</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>DependentFieldComparatorWithoutCoefficient</name>
    <filename>struct_dependent_field_comparator_without_coefficient.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type>bool</type>
      <name>operator()</name>
      <anchorfile>struct_dependent_field_comparator_without_coefficient.html</anchorfile>
      <anchor>a0392f0fc6339ac9d803b21a13ff5cf45</anchor>
      <arglist>(const DependentFieldTerm&lt; dim, spacedim &gt; &amp;dependent_field_1, const DependentFieldTerm&lt; dim, spacedim &gt; &amp;dependent_field_2) const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DependentFieldTerm</name>
    <filename>class_dependent_field_term.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>DependentFieldTerm</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a149a09f602b0baadba9033d6a8da448c</anchor>
      <arglist>(const double coefficient, const IndependentField&lt; dim, spacedim &gt; &amp;independent_field, const unsigned int component=0, const std::vector&lt; unsigned int &gt; derivatives=std::vector&lt; unsigned int &gt;())</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>first_derivative</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a7f7b45969b939cfead095d2681d2f612</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_derivatives</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>ae901f99afbdd336f73cf550e62c93ae1</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator&lt;</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a204178ccb13a1280b96cf2ad8cf24695</anchor>
      <arglist>(const DependentFieldTerm &amp;dependent_field_2) const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a017b2119451d74e4dbbe203503d498c4</anchor>
      <arglist>(const DependentFieldTerm &amp;dependent_field_2) const</arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>coefficient</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a59a9183a32ac55fb728f3797b68a9f8f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; dim, spacedim &gt; &gt;</type>
      <name>independent_field</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a89d1c3fea36e6fe105232097a321e095</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>component</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>ac6f3ac40d4ee2c8b9f9bbdfa34079b74</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; unsigned int &gt;</type>
      <name>derivatives</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>af09c5452c3e8e71e9ee99db304b90135</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DirichletConstraint</name>
    <filename>class_dirichlet_constraint.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>DirichletConstraint</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>a14e75a7f8bcc8a18daa04a5c1164fe87</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; &amp;independent_field, const unsigned int component, const InterfaceSide side, const std::set&lt; types::material_id &gt; domain_of_constraint, const Function&lt; spacedim &gt; *const constraint_inhomogeneity=nullptr, const IndependentField&lt; 0, spacedim &gt; *independent_scalar=nullptr, const Function&lt; spacedim &gt; *const coefficient_c=nullptr)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~DirichletConstraint</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>ac6b6b7caf4694000486abff377fd1440</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_constraint_is_active</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>abdfdde6ef50639152e2a82d54b130042</anchor>
      <arglist>(const bool constraint_is_active)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_constraint_is_active</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>ad1741529befff8fca4af6ee82224abb3</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_time</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>ab8c509170d5d36f778415b88bceee99d</anchor>
      <arglist>(const double time) const</arglist>
    </member>
    <member kind="variable">
      <type>const std::set&lt; types::material_id &gt;</type>
      <name>domain_of_constraint</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>a258b6ff11b206f966bb03943bb11f469</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt;</type>
      <name>independent_field</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>abbd7945a973ed93d1d773307393ffde3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>component</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>a7e3c4d0e0906af1c81b88e05e41bdafc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceSide</type>
      <name>side</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>ae049d107664d3bf23d287ec77545b6f3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>constraint_inhomogeneity</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>af22d2bca23999bdb4f2b67e8982d29c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt;</type>
      <name>independent_scalar</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>a8793f0d41a9c6e88638d7cf5d52a6fdd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>coefficient_c</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>adea2dc6126a633b297eb7dfce832b9ef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>constraint_is_active</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>a1ae0b53768be76747f31b5bc584d9ed6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DoFHandlerSystem</name>
    <filename>class_do_f_handler_system.html</filename>
    <templarg>spacedim</templarg>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>InterfaceCell</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a64bfb9dab5f61c2583896a7d086cb1f2</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DoFHandler&lt; spacedim-1, spacedim &gt;::cell_iterator</type>
      <name>InterfaceCellDoF</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ab74b04bb37380033195cfb081f505214</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>DomainCell</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a84df955a5d8da7e1439a10f2a64811c0</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DoFHandler&lt; spacedim, spacedim &gt;::cell_iterator</type>
      <name>DomainCellDoF</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a24490ae81fcdaf7b94e1ef0cae46d096</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DoFHandlerSystem</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>acb31af36eec0abac1d6e764e3d08c7fe</anchor>
      <arglist>(const TriangulationSystem&lt; spacedim &gt; &amp;tria_system)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~DoFHandlerSystem</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a5a2f0389f4f4216b9841d9c518cb107e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>distribute_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a1a3f87726392121bf4b88e9adf09cbc7</anchor>
      <arglist>(const hp::FECollection&lt; spacedim, spacedim &gt; &amp;fe_collection_domain, const hp::FECollection&lt; spacedim-1, spacedim &gt; &amp;fe_collection_interface, const unsigned int n_additional_dofs=0, const std::vector&lt; unsigned int &gt; &amp;renumbering_domain=std::vector&lt; unsigned int &gt;(), const std::vector&lt; unsigned int &gt; &amp;renumbering_interface=std::vector&lt; unsigned int &gt;())</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCellsDoF&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_begin_active</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a88d93ac03defd1b296238de4eb122841</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCellsDoF&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_end_active</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aecbc22a86943c1bc921a13ec7e2333f2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; InterfaceCellDomainCellsDoF&lt; spacedim &gt; &gt; &amp;</type>
      <name>interface_active_iterators</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a80010f3d2b4f0f524cc4ed79e1f93762</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>DomainCellDoFIterator&lt; spacedim &gt;</type>
      <name>domain_begin_active</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>af005d3cdba9149ba44e8d0545f5d9135</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>DomainCellDoFIterator&lt; spacedim &gt;</type>
      <name>domain_end_active</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>abdc3e30591dfeffb00356468f5a0eb10</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>IteratorRange&lt; DomainCellDoFIterator&lt; spacedim &gt; &gt;</type>
      <name>domain_active_iterators</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ab1ce368741781eb042e11ac3d3b7e651</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_dofs_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a649873cd660a9888cd4ff29c35b1df2a</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_dofs_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ab0157b6d29707c6657d21610f2bfe05d</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_dofs_additional</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a3391c5a179b3d2952e7de0eb88435709</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aa0613fb931568c1f5efde3c4f6b97c9f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>add1750163fd42999c1892dce94114a43</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices) const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ac7c2592ff045c5a7cf4c2ee1de0e2c73</anchor>
      <arglist>(const unsigned int &amp;dof_index) const</arglist>
    </member>
    <member kind="function">
      <type>const DoFHandler&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_dof_handler_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a213bc9448c093924d8e44c2523492327</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const DoFHandler&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_dof_handler_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a47e09ae23ad1f653f0234afe845e4eec</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>DoFHandler&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_dof_handler_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aaa5597e66f5a62935b83f9c862e3ae9d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DoFHandler&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_dof_handler_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a4d0745d916a6ab3beb6857337b4ffc20</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const IndexSet &amp;</type>
      <name>get_locally_owned_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aa07844f26a5f98734535b29790bd9dca</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const IndexSet &amp;</type>
      <name>get_locally_relevant_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>af5e18e9200d2137f01c35a19e0890739</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>attach_dof_renumbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ab31e887efe0c3703e3a777287c250ed0</anchor>
      <arglist>(const DoFRenumbering &amp;dof_renumbering)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>make_hanging_node_constraints</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ac4c0dda51291e1a0d5ddaaa23c07475a</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_vector</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a106000629dee3e614dc87b2678d5a4e0</anchor>
      <arglist>(const VectorType &amp;in_vect, VectorType &amp;out_vect_domain, VectorType &amp;out_vect_interface, VectorType &amp;out_vect_C) const</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; unsigned int &gt; &amp;</type>
      <name>get_n_dofs_per_processor</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a91bdc19b8c8226f827aed689ad654ffe</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_single_dof_index_component_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a0df45d966fef705935ef71943f3e373d</anchor>
      <arglist>(const unsigned int component) const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_single_dof_index_component_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a93c1f41edb851a7322dd183089d11e81</anchor>
      <arglist>(const unsigned int component) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices_component_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a0d49f02ccbdd681028ce66fea9341fe1</anchor>
      <arglist>(const unsigned int component, std::set&lt; unsigned int &gt; &amp;indices) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices_component_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aa4bb80864110f7917ae86525e3c0d799</anchor>
      <arglist>(const unsigned int component, std::set&lt; unsigned int &gt; &amp;indices) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_map_dof_index_component_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ac2617866e4d92dfb309a0212e2bc116d</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;components) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_map_dof_index_component_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a99e00b1efebbfee99db0478c3eac7787</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;components) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_map_dof_index_support_point_index_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a0af24410554802a908ce2800550007b8</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;support_points, std::map&lt; unsigned int, Point&lt; spacedim &gt;&gt; &amp;map_support_point_index_support_point_location, const dealii::Mapping&lt; spacedim-1, spacedim &gt; &amp;mapping) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_map_dof_index_support_point_index_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a298db5cd2cd488abbacb72c66fc9703e</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;support_points, std::map&lt; unsigned int, Point&lt; spacedim &gt;&gt; &amp;map_support_point_index_support_point_location, const dealii::Mapping&lt; spacedim &gt; &amp;mapping) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>update_interface_domain_relation</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a237b112398864f03ef97873f276960f5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>set_locally_owned_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ab24852414baf7b391b67d21644bd66f8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>set_locally_relevant_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>affa92a0bb132688295b748e0c5436d2e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>set_locally_owned_dofs_standard_numbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a56c9cd6d17075c1b7a6d24d5d28357b2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>set_locally_relevant_dofs_standard_numbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a51f57afbe35de6bdd39a979911eaa838</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>split_vector_implementation</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a52da99405f6258795d414a2348097380</anchor>
      <arglist>(const LinearAlgebra::distributed::Vector&lt; double &gt; &amp;in_vect, LinearAlgebra::distributed::Vector&lt; double &gt; &amp;out_vect, const unsigned int window_begin, const unsigned int window_end) const</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>split_vector_implementation</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a6ef88586f9bc65676fad4a5281a489c0</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;in_vect, Vector&lt; double &gt; &amp;out_vect, const unsigned int window_begin, const unsigned int window_end) const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const TriangulationSystem&lt; spacedim &gt; &gt;</type>
      <name>tria_system</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a06d93193cb47591db138cd8f41953796</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::shared_ptr&lt; DoFHandler&lt; spacedim, spacedim &gt; &gt;</type>
      <name>dof_handler_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a342ca54f8d2447244d380db5025dab03</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::shared_ptr&lt; DoFHandler&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>dof_handler_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aa76badc239ecf29fa6021994d3411b6c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>dof_indices_C</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a15c27ca30c5905b44691779c755cd69c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>locally_owned_dof_indices_C</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a7b2e77b2c718b0b13d861da7c0530b28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; InterfaceCellDomainCellsDoF&lt; spacedim &gt; &gt;</type>
      <name>active_interface_cell_domain_cells</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>af0119b14377300f7f0457a18bb7dcd67</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SmartPointer&lt; const DoFRenumbering &gt;</type>
      <name>dof_renumbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aef6159c606a24ac7daadcd3fe082b3b6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>IndexSet</type>
      <name>locally_owned_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ad72a701a3581187eec846c831ba384c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>IndexSet</type>
      <name>locally_relevant_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a6b2de0e80cf0e67e62391cc84f8a7be5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>n_dofs_per_processor</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>acfe853db91a67d1f78e08d62d91bab1b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>IndexSet</type>
      <name>locally_owned_dofs_standard_numbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a257ddb680d9f8276f7b802884dd8e577</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>IndexSet</type>
      <name>locally_relevant_dofs_standard_numbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a0e6744c40aada76d695f2c564f4c5040</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; boost::signals2::connection &gt;</type>
      <name>tria_listeners</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a286776a935dacb3d1c79a90f4ca7e5d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>InterfaceCellDoFIterator</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a266cb42b2fa4af49d4df0d938f962749</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>DomainCellDoFIterator</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a22fa60ad60906aacbbe21d3b5704ebfc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DoFRenumbering</name>
    <filename>class_do_f_renumbering.html</filename>
    <base>Subscriptor</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>convert_dof_indices</name>
      <anchorfile>class_do_f_renumbering.html</anchorfile>
      <anchor>a3dd9e72a14b3078a9f823222eacd162c</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt;</type>
      <name>convert_range</name>
      <anchorfile>class_do_f_renumbering.html</anchorfile>
      <anchor>a39be83786966bdce332a902104993c74</anchor>
      <arglist>(const unsigned int range_begin, const unsigned int range_end) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~DoFRenumbering</name>
      <anchorfile>class_do_f_renumbering.html</anchorfile>
      <anchor>a0443d7cb614fa161e4201ff7a6bdf541</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DoFRenumberingOffset</name>
    <filename>class_do_f_renumbering_offset.html</filename>
    <base>DoFRenumbering</base>
    <member kind="function">
      <type>void</type>
      <name>add_range</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>a39333b4c3597c2b4a967f7af3cec0054</anchor>
      <arglist>(const unsigned int dof_start, const unsigned int dof_end, const int offset)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>convert_dof_indices</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>a11744a010a6502ae0aa9bac11c7ff3f0</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt;</type>
      <name>convert_range</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>adef238432ea7b8100c3d4d0e008a9f72</anchor>
      <arglist>(const unsigned int range_begin, const unsigned int range_end) const</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; std::tuple&lt; const unsigned int, const unsigned int, const int &gt; &gt; &amp;</type>
      <name>get_dof_offsets</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>a672fcd1ce3fa89174aceb105b00e7814</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clear</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>ac7156c58d0aa5dcb57fa0fe3ac437453</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>a0781d3aa455908ba49e61eafaca5f50b</anchor>
      <arglist>(std::ostream &amp;out) const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::tuple&lt; const unsigned int, const unsigned int, const int &gt; &gt;</type>
      <name>dof_offsets</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>a5df6c4b70b1394c3670ced634146c9a9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DomainCellDoFIterator</name>
    <filename>class_domain_cell_do_f_iterator.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>DomainCellDoFIterator</name>
      <anchorfile>class_domain_cell_do_f_iterator.html</anchorfile>
      <anchor>a36416883d4c91dd0849520d3e205dba1</anchor>
      <arglist>(const TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt;&gt; &amp;domain_cell, const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DomainCellDoFIterator</name>
      <anchorfile>class_domain_cell_do_f_iterator.html</anchorfile>
      <anchor>a4e987bc73bb2c8ec6974cc6427505091</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices</name>
      <anchorfile>class_domain_cell_do_f_iterator.html</anchorfile>
      <anchor>a0727db1344983d8f7c3ff40c87e0ba40</anchor>
      <arglist>(std::vector&lt; types::global_dof_index &gt; &amp;dof_indices) const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>dof_handler_system</name>
      <anchorfile>class_domain_cell_do_f_iterator.html</anchorfile>
      <anchor>a7431a0505f1e5b0b796c7f58cf9a060c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>FE_DGQAnisotropic</name>
    <filename>class_f_e___d_g_q_anisotropic.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <base>FE_Poly&lt; dim, dim &gt;</base>
    <member kind="function">
      <type></type>
      <name>FE_DGQAnisotropic</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>a1504042c3096ba9e4b958ec699667d6c</anchor>
      <arglist>(const Quadrature&lt; 1 &gt; &amp;points_x, const Quadrature&lt; 1 &gt; &amp;points_y, const Quadrature&lt; 1 &gt; &amp;points_z)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FE_DGQAnisotropic</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>a3e6e9cd0f18fe42ec99aebaff5364cce</anchor>
      <arglist>(const Quadrature&lt; 1 &gt; &amp;points_x, const Quadrature&lt; 1 &gt; &amp;points_y)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>convert_generalized_support_point_values_to_dof_values</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>a16f2216f14a001d5e16704c41ffb7f09</anchor>
      <arglist>(const std::vector&lt; Vector&lt; double &gt;&gt; &amp;support_point_values, std::vector&lt; double &gt; &amp;nodal_values) const override</arglist>
    </member>
    <member kind="function">
      <type>std::unique_ptr&lt; FiniteElement&lt; dim, spacedim &gt; &gt;</type>
      <name>clone</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>aba12ce2e8eec0366c87b26d62b3d8d89</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>get_name</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>af639aedcc90606523198e7037c142e2e</anchor>
      <arglist>() const override</arglist>
    </member>
    <member kind="function">
      <type>FiniteElementDomination::Domination</type>
      <name>compare_for_domination</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>aa1ecc8619e7841cff95c9ad226bfe373</anchor>
      <arglist>(const FiniteElement&lt; dim, spacedim &gt; &amp;fe_other, const unsigned int codim) const override</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::pair&lt; unsigned int, unsigned int &gt; &gt;</type>
      <name>hp_vertex_dof_identities</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>a14f2f28e2672694817fae015e37164db</anchor>
      <arglist>(const FiniteElement&lt; dim, spacedim &gt; &amp;fe_other) const override</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::pair&lt; unsigned int, unsigned int &gt; &gt;</type>
      <name>hp_line_dof_identities</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>a531c3309fd09eb4a5699a930ca7c3bc2</anchor>
      <arglist>(const FiniteElement&lt; dim, spacedim &gt; &amp;fe_other) const override</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::pair&lt; unsigned int, unsigned int &gt; &gt;</type>
      <name>hp_quad_dof_identities</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>ad9aa1f1d252a4db26754c8da9fd7c9e0</anchor>
      <arglist>(const FiniteElement&lt; dim, spacedim &gt; &amp;fe_other, const unsigned int face_no=0) const override</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>degree_x</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>acdf735a7ef2198784ac1771a47d22b1b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>degree_y</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>ad48f8d9428ac7719985eea4fde6c9b54</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>degree_z</name>
      <anchorfile>class_f_e___d_g_q_anisotropic.html</anchorfile>
      <anchor>aa585c966670d06d56d5b15acc3e76b2f</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>FEValuesInterface</name>
    <filename>class_f_e_values_interface.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>FEValuesInterface</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>ae48863eabc0754e1cc033143fcd2d443</anchor>
      <arglist>(const Mapping&lt; spacedim-1, spacedim &gt; &amp;mapping_interface, const Mapping&lt; spacedim, spacedim &gt; &amp;mapping_domain, const FiniteElement&lt; spacedim-1, spacedim &gt; &amp;fe_interface, const FiniteElement&lt; spacedim, spacedim &gt; &amp;fe_domain_minus, const FiniteElement&lt; spacedim, spacedim &gt; &amp;fe_domain_plus, const Quadrature&lt; spacedim-1 &gt; &amp;quadrature, const UpdateFlags update_flags_interface, const UpdateFlags update_flags_domain_minus, const UpdateFlags update_flags_domain_plus)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a334ed0ef179516854690e7e3fa273aa0</anchor>
      <arglist>(const InterfaceCellDomainCellsDoF&lt; spacedim &gt; &amp;interface_cell_domain_cells)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a875a17d471069ddb97cbd8e6f7aba2b8</anchor>
      <arglist>(const InterfaceCellDomainCells&lt; spacedim &gt; &amp;interface_cell_domain_cells)</arglist>
    </member>
    <member kind="function">
      <type>const Quadrature&lt; spacedim-1 &gt; &amp;</type>
      <name>get_quadrature</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>abbaae0efa290c18b4ecfabc90b151136</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const FEValues&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_fe_values_interface</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a15a3cc642c056496b79e41732f7012c4</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const FEFaceValuesBase&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_fe_values_domain</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a85aeb4a21ad21576d2d550f2d57f3a00</anchor>
      <arglist>(const InterfaceSide &amp;interface_side) const</arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_quadrature_points</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a100afc348b15432f9e241ff56ebffa8d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>InterfaceRefinementCase</type>
      <name>initialized_interface_refinement_case</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a4a14b4fa181df533a6c7e7971fad5c4e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const FEFaceValuesBase&lt; spacedim, spacedim &gt; *</type>
      <name>initialized_fe_values_domain_minus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a48f1a887c14615d665408f156d879980</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const FEFaceValuesBase&lt; spacedim, spacedim &gt; *</type>
      <name>initialized_fe_values_domain_plus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>ab1f7524d517a79cd5e4adbd43017f415</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FEValues&lt; spacedim-1, spacedim &gt;</type>
      <name>fe_values_interface</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a75f4ee97ebaefa2d1ded52c445953473</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FEFaceValues&lt; spacedim, spacedim &gt;</type>
      <name>fe_face_values_domain_minus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a1ce4354bfe2e852fdc5b7917614404a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FEFaceValues&lt; spacedim, spacedim &gt;</type>
      <name>fe_face_values_domain_plus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a6d6820b66a2694d327238c6b3c7fb30c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FESubfaceValues&lt; spacedim, spacedim &gt;</type>
      <name>fe_subface_values_domain_minus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>abe003c242c807dd162c396d812bbfe77</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FESubfaceValues&lt; spacedim, spacedim &gt;</type>
      <name>fe_subface_values_domain_plus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a184e705efd975db536299ad619585a65</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>FunctionCell</name>
    <filename>class_function_cell.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>FunctionCell</name>
      <anchorfile>class_function_cell.html</anchorfile>
      <anchor>a283cbfe9d1af3fcee24517e3b63b5bc8</anchor>
      <arglist>(const unsigned int n_components=1)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual</type>
      <name>~FunctionCell</name>
      <anchorfile>class_function_cell.html</anchorfile>
      <anchor>a4f91c2b4961ca9330c6e564014879670</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_cell</name>
      <anchorfile>class_function_cell.html</anchorfile>
      <anchor>a36b161ae96e789de9d03654497de03cd</anchor>
      <arglist>(const TriaIterator&lt; CellAccessor&lt; dim, spacedim &gt;&gt; &amp;cell)</arglist>
    </member>
    <member kind="function">
      <type>const TriaIterator&lt; CellAccessor&lt; dim, spacedim &gt; &gt; &amp;</type>
      <name>get_cell</name>
      <anchorfile>class_function_cell.html</anchorfile>
      <anchor>ab3f22f3fa0cb7c3eec847f883fb7fecd</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const TriaIterator&lt; CellAccessor&lt; dim, spacedim &gt; &gt; *</type>
      <name>cell</name>
      <anchorfile>class_function_cell.html</anchorfile>
      <anchor>a3ef39628d9141076b3418fe2bb529dfe</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IndependentField</name>
    <filename>class_independent_field.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>IndependentField</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>afb12ec32511568a9504312730542bfa1</anchor>
      <arglist>(const std::string name, const FiniteElement&lt; dim, spacedim &gt; &amp;fe, const unsigned int n_components, const std::set&lt; types::material_id &gt; non_zero_regions, const Function&lt; spacedim &gt; *const initial_vals=nullptr, const bool is_local=false, const bool is_locally_eliminated=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IndependentField</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>af953c9ecc5965d5bcbd61fab9e0b6e2b</anchor>
      <arglist>(const std::string name, const FiniteElement&lt; dim, spacedim &gt; &amp;fe, const std::set&lt; types::material_id &gt; non_zero_regions, const Function&lt; spacedim &gt; *const initial_vals=nullptr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IndependentField</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>af83ee8c1600bb079a3075244a7b39481</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>ae05f8565e4ce1a70b5b833555dc084b5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::unique_ptr&lt; const FiniteElement&lt; dim, spacedim &gt; &gt;</type>
      <name>fe</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a1c583665b7710bd3b815b03ba026b6d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_components</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a7b19ea8c30d72cf27f05669de61f30a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::set&lt; types::material_id &gt;</type>
      <name>non_zero_regions</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a4e09e114870c0b3761bc2e32916e5850</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>initial_vals</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a274c902785d2937a6065f7e09f3976c3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const bool</type>
      <name>is_local</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a8e77d8d321a259bec955a71f55ef41e5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const bool</type>
      <name>is_locally_eliminated</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a755514cb31014d9b39e1b3e677c4527d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IndependentField&lt; 0, spacedim &gt;</name>
    <filename>class_independent_field_3_010_00_01spacedim_01_4.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>IndependentField</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a8a486862088446941ea97fc471c5d44c</anchor>
      <arglist>(const std::string name, const double initial_value=0.0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IndependentField</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a653c01f2a6aeffba44542d9788433647</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a05ecdcc8310253f055fbc59abaa2bc90</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>initial_value</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a8c4c434806dadc5b885102322e27357c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_components</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a248c2226570c71914f5219e7e3052561</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>InterfaceCellDoFIterator</name>
    <filename>class_interface_cell_do_f_iterator.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>InterfaceCellDoFIterator</name>
      <anchorfile>class_interface_cell_do_f_iterator.html</anchorfile>
      <anchor>a41a74a50c801c67a5393d568ed2768b7</anchor>
      <arglist>(const TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt;&gt; &amp;interface_cell, const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterfaceCellDoFIterator</name>
      <anchorfile>class_interface_cell_do_f_iterator.html</anchorfile>
      <anchor>a0dc67a5bf4cbf83e786af1fb28b6dd7c</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices</name>
      <anchorfile>class_interface_cell_do_f_iterator.html</anchorfile>
      <anchor>af4d16fe93631559a8d4405352ebf1546</anchor>
      <arglist>(std::vector&lt; types::global_dof_index &gt; &amp;dof_indices) const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>dof_handler_system</name>
      <anchorfile>class_interface_cell_do_f_iterator.html</anchorfile>
      <anchor>a62045b61faee901392a7326b624d2b02</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>InterfaceCellDomainCells</name>
    <filename>class_interface_cell_domain_cells.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>InterfaceCell</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>ac6c8ada1cde14364575ef0191a278bf8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>DomainCell</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>aea8f8f65d0f5021da9c5b8e147b6d587</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterfaceCellDomainCells</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a016b875566033b0da13d774081eaeaf4</anchor>
      <arglist>(const DomainCell &amp;domain_cell, const unsigned int face, const InterfaceCell &amp;interface_cell, const InterfaceSide interface_side)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~InterfaceCellDomainCells</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a7dae568b2d3c80762d0e1af52c93f279</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::tuple&lt; const types::material_id, const types::material_id, const types::material_id &gt;</type>
      <name>get_material_ids</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a6823df36e862abe155a90e25d95cf76e</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceCell</type>
      <name>interface_cell</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>acabc5a62be3f7742c3c334ae1777fbd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const DomainCell</type>
      <name>domain_cell_minus</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>ae610a6a7b0ff4a421ee051cb873bbc22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>face_minus</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a51a9eaf54de991f5fa432220b09205cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const DomainCell</type>
      <name>domain_cell_plus</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a72c7faaed3a84c546d47960f1064f9df</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>face_plus</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a8780275f79c7137ef65df9b7d4039017</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceRefinementCase</type>
      <name>refinement_case</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>ab1b5469ca5c40256942ea179abbba92c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>subface</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>aba7f3e048c6985006988d716cb31aa7d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type></type>
      <name>InterfaceCellDomainCells</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a4c0fd73a697833a1257b9341982897de</anchor>
      <arglist>(const std::tuple&lt; const InterfaceCell, const DomainCell, const unsigned int, const DomainCell, const unsigned int, const InterfaceRefinementCase, const unsigned int &gt; input)</arglist>
    </member>
    <member kind="function" protection="private" static="yes">
      <type>static std::tuple&lt; const InterfaceCell, const DomainCell, const unsigned int, const DomainCell, const unsigned int, const InterfaceRefinementCase, const unsigned int &gt;</type>
      <name>convert_constructor_inputs</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a3f892a73f1f38f7412e549b9e62a8928</anchor>
      <arglist>(const DomainCell &amp;domain_cell, const unsigned int face, const InterfaceCell &amp;interface_cell, const InterfaceSide interface_side)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>InterfaceCellDomainCellsDoF</name>
    <filename>class_interface_cell_domain_cells_do_f.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>InterfaceCell</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a55783857b93a31baac5a6e6ee3dad22f</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DoFHandler&lt; spacedim-1, spacedim &gt;::cell_iterator</type>
      <name>InterfaceCellDoF</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>ac24e469237fddaf0310a5f818ba13afd</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>DomainCell</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a0d0b12a444adb340e1ffc2dd23d77902</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DoFHandler&lt; spacedim, spacedim &gt;::cell_iterator</type>
      <name>DomainCellDoF</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a78225607ebf9880311f4fdeb30c5a89d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterfaceCellDomainCellsDoF</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a79a88327d88717a2c3381fe7131b435d</anchor>
      <arglist>(const InterfaceCellDomainCells&lt; spacedim &gt; &amp;interface_cell_domain_cells, const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~InterfaceCellDomainCellsDoF</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a090e1536777928b78dff7fa7b1021a1e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::tuple&lt; const types::material_id, const types::material_id, const types::material_id &gt;</type>
      <name>get_material_ids</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a4823f0f9c791b61d6fb75699952464e5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices_local_global_interface</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>ab6ff11f2adfa7cd1fd13a2d3e9b3b130</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global, std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global_minus, std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global_plus) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices_local_global_interface</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a69d8ac580b83b2469ba1cc611bb5776b</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global) const</arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceCellDoFIterator&lt; spacedim &gt;</type>
      <name>interface_cell</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a5fe7922af91545598d89c77ba42a59da</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const DomainCellDoFIterator&lt; spacedim &gt;</type>
      <name>domain_cell_minus</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>ac8c84df4298d2d5be54dc8bf81e48c4b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>face_minus</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a01fc12497d44b06c7a325ea8f0b2ad28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const DomainCellDoFIterator&lt; spacedim &gt;</type>
      <name>domain_cell_plus</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a8fd9dbd3a1e53023ecca7917b004bab1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>face_plus</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>afcb7abea75dfc6e310493a6ff9aadcf6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceRefinementCase</type>
      <name>refinement_case</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>aebb7e5f13d079fc83f98f67bcfcc6de3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>subface</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>ac0f33cbe60afdc4f5856ea39a085ed1a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LinearMaterialDomain</name>
    <filename>class_linear_material_domain.html</filename>
    <templarg>spacedim</templarg>
    <base>ScalarFunctional&lt; spacedim, spacedim &gt;</base>
    <member kind="function">
      <type></type>
      <name>LinearMaterialDomain</name>
      <anchorfile>class_linear_material_domain.html</anchorfile>
      <anchor>a7f91a980e9475691513443cd6b2efbc1</anchor>
      <arglist>(const std::vector&lt; DependentField&lt; spacedim, spacedim &gt;&gt; e_omega, const std::set&lt; types::material_id &gt; domain_of_integration, const Quadrature&lt; spacedim &gt; quadrature, const FullMatrix&lt; double &gt; C, const Vector&lt; double &gt; y, const std::string name=&quot;LinearMaterialDomain&quot;)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_h_omega</name>
      <anchorfile>class_linear_material_domain.html</anchorfile>
      <anchor>a4d73bbb0e545ee21c5635e7dac071719</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, double &amp;h_omega, Vector&lt; double &gt; &amp;h_omega_1, FullMatrix&lt; double &gt; &amp;h_omega_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const FullMatrix&lt; double &gt;</type>
      <name>C</name>
      <anchorfile>class_linear_material_domain.html</anchorfile>
      <anchor>a8ac95fcf4f77790670b8520ba9593a7c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const Vector&lt; double &gt;</type>
      <name>y</name>
      <anchorfile>class_linear_material_domain.html</anchorfile>
      <anchor>a7ea9ab6930c0b0aa826e809ef245b0e2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LinearMaterialInterface</name>
    <filename>class_linear_material_interface.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>LinearMaterialInterface</name>
      <anchorfile>class_linear_material_interface.html</anchorfile>
      <anchor>ac10a5c148dd81edab41d22dd207c78a2</anchor>
      <arglist>(const std::vector&lt; DependentField&lt; spacedim-1, spacedim &gt;&gt; e_sigma, const std::set&lt; types::material_id &gt; domain_of_integration, const Quadrature&lt; spacedim-1 &gt; quadrature, const FullMatrix&lt; double &gt; C, const Vector&lt; double &gt; y, const std::string name=&quot;LinearMaterialInterface&quot;)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_h_sigma</name>
      <anchorfile>class_linear_material_interface.html</anchorfile>
      <anchor>af5cca6ebef9371fd07782e3ae0171d65</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n, double &amp;h_sigma, Vector&lt; double &gt; &amp;h_sigma_1, FullMatrix&lt; double &gt; &amp;h_sigma_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const FullMatrix&lt; double &gt;</type>
      <name>C</name>
      <anchorfile>class_linear_material_interface.html</anchorfile>
      <anchor>ae0b1ef211453e5fd3e1416ff37c315db</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const Vector&lt; double &gt;</type>
      <name>y</name>
      <anchorfile>class_linear_material_interface.html</anchorfile>
      <anchor>a3864513d7662e4d1c91e606467befd59</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>PointConstraint</name>
    <filename>class_point_constraint.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>PointConstraint</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>a3820c38236423c0bf32129258b8f24f8</anchor>
      <arglist>(const IndependentField&lt; dim, spacedim &gt; &amp;independent_field, const unsigned int component, const Point&lt; spacedim &gt; X, const Function&lt; spacedim &gt; *const constraint_inhomogeneity=nullptr, const IndependentField&lt; 0, spacedim &gt; *independent_scalar=nullptr, const Function&lt; spacedim &gt; *const coefficient_c=nullptr, const bool ignore_point=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~PointConstraint</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>afcaa6aaeb90f8128f121ac3797f77f03</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_constraint_is_active</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>afeb601e598e753c669abcd86d84ac3d3</anchor>
      <arglist>(const bool constraint_is_active)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_X</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>a748280c69c257904cad201e837a9b5d4</anchor>
      <arglist>(const Point&lt; spacedim &gt; X)</arglist>
    </member>
    <member kind="function">
      <type>Point&lt; spacedim &gt;</type>
      <name>get_X</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>a181f4e7b3042afed35db300290d46aa5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_constraint_is_active</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>a19659b5f07f5f3013ecb59eecefd56ea</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_time</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>a2178075182226dcf265e5b956190ed15</anchor>
      <arglist>(const double time) const</arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; dim, spacedim &gt; &gt;</type>
      <name>independent_field</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>a8d27c23d2d6e1e64f8b2a36f7c950efb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>component</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>af15fcbe78e2f8233ae1ff4b52b9d9de5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>constraint_inhomogeneity</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>ace3fa5f9542ee73d502abd78f48da032</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt;</type>
      <name>independent_scalar</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>a4f6c27593a39a2b1064020c256120f0a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>coefficient_c</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>a5b1a621c07e8ece928f44e2dabd01d25</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const bool</type>
      <name>ignore_point</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>a5c2605fac61aab48c65d76137c172a78</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>constraint_is_active</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>ae09c78d247df993db4fced0460889161</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Point&lt; spacedim &gt;</type>
      <name>X</name>
      <anchorfile>class_point_constraint.html</anchorfile>
      <anchor>ac021367793862b52580cc292bb7b0ea7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>PointConstraint&lt; spacedim, spacedim &gt;</name>
    <filename>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>PointConstraint</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a13052c609467fb9bfda281e37699b196</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; &amp;independent_field, const unsigned int component, const Point&lt; spacedim &gt; X, const Function&lt; spacedim &gt; *const constraint_inhomogeneity=nullptr, const IndependentField&lt; 0, spacedim &gt; *independent_scalar=nullptr, const Function&lt; spacedim &gt; *const coefficient_c=nullptr, const bool ignore_point=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~PointConstraint</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ade399e3c4a1d033505ddce38f7f5d1ce</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_constraint_is_active</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a3e06652b150bbcb759eb20b1b1a805e0</anchor>
      <arglist>(const bool constraint_is_active)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_X</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a292f383bdf4563d562dda619df76aa11</anchor>
      <arglist>(const Point&lt; spacedim &gt; X)</arglist>
    </member>
    <member kind="function">
      <type>Point&lt; spacedim &gt;</type>
      <name>get_X</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a4264d07fbc938311ea069b5a5a3b4165</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_constraint_is_active</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a9c6cde4645e498de1c74fd2c4b04df06</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_time</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a725ee8b397b6f5c2b54795e4fb5173f2</anchor>
      <arglist>(const double time) const</arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt;</type>
      <name>independent_field</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a8b124e7568331aef96870450a0113bab</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>component</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a92118ca233dadb0ff3495317cee3fe72</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>constraint_inhomogeneity</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>af6bbff9cc333d7aa678180f3b5f5d569</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt;</type>
      <name>independent_scalar</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a7406c5d95454bf68acc2f4745c1b4536</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>coefficient_c</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a15e0c3b340dca72427147c35ca40f3f3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const bool</type>
      <name>ignore_point</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ad6e6b3c5bb689c75b6aca0ef87df1b12</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>constraint_is_active</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a8f876badee193506b88dbb0959ac6776</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Point&lt; spacedim &gt;</type>
      <name>X</name>
      <anchorfile>class_point_constraint_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a8d6c7dd2a1c97869ea347ab18d57aefe</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ScalarFunctional</name>
    <filename>class_scalar_functional.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>ScalarFunctional</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>ac0618ada3b80400784ee8553a66aade9</anchor>
      <arglist>(const std::vector&lt; DependentField&lt; dim, spacedim &gt;&gt; e_sigma, const std::set&lt; types::material_id &gt; domain_of_integration, const Quadrature&lt; dim &gt; quadrature, const std::string name, const unsigned int n_ref_sets=0, const unsigned int n_hidden=0, const Function&lt; spacedim &gt; *const initial_vals_hidden=nullptr)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>get_h_sigma</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>adb295fb739a743d5a1273025eb8dae72</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n, double &amp;h_sigma, Vector&lt; double &gt; &amp;h_sigma_1, FullMatrix&lt; double &gt; &amp;h_sigma_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const =0</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_h_omega</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a0f5c776b9084f350e1c0defa69c81047</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, double &amp;h_omega, Vector&lt; double &gt; &amp;h_omega_1, FullMatrix&lt; double &gt; &amp;h_omega_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>get_maximum_step</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>ae4f5fbd69cabfda73cc2b30ae0263ca5</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, const Vector&lt; double &gt; &amp;delta_e_sigma, const Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n) const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_maximum_step</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a71310146c6cc3180e25d90f7b7329c44</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, const Vector&lt; double &gt; &amp;delta_e_omega, const Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compare_derivatives_with_numerical_derivatives</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a0c0a2095a33e657a7a984e706ed7f968</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n, const std::string detailed_printout_file=&quot;&quot;, const double epsilon=1e-8) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>set_q_point_id</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a88573ddc3d2467301a5ca30590d82303</anchor>
      <arglist>(std::string q_point_id) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ScalarFunctional</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>aed36f35e6d2c9f9a93ee2749d01f2a51</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::set&lt; types::material_id &gt;</type>
      <name>domain_of_integration</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>ae3b6dd6934e1cd55fcc55cf344179407</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>initial_vals_hidden</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a602d0bc2c945822c6b756fc63183ae2b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; DependentField&lt; dim, spacedim &gt; &gt;</type>
      <name>e_sigma</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a86662b03a63219227993a2c6c07aefc1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; DependentField&lt; dim, spacedim &gt; &gt;</type>
      <name>e_omega</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>aff467b5588c5a8a28b805ae886205f9c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const Quadrature&lt; dim &gt;</type>
      <name>quadrature</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>adea9ff214aeb2a1d8c3712a9d2433883</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_hidden</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a8b1617930242870f22eef5e306cb717f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a4d184688053b3443d10e228e4a8eba60</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_ref_sets</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a7e12423f4b29e9e0aaa0f7f9c2d1c0eb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>q_point_id</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a1609d09f1f8c2084e4c0f8b5432dc86c</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>AssemblyHelper</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>af4019c2e39cc934d646aaa35c3c52773</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ScalarFunctional&lt; spacedim, spacedim &gt;</name>
    <filename>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>ScalarFunctional</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a4c4903b402908f26b8c8b433a7cbeb76</anchor>
      <arglist>(const std::vector&lt; DependentField&lt; spacedim, spacedim &gt;&gt; e_omega, const std::set&lt; types::material_id &gt; domain_of_integration, const Quadrature&lt; spacedim &gt; quadrature, const std::string name, const unsigned int n_ref_sets=0, const unsigned int n_hidden=0, const Function&lt; spacedim &gt; *const initial_vals_hidden=nullptr)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>get_h_omega</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a8665cb5b5a57ad22217e4c112845d43b</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, double &amp;h_omega, Vector&lt; double &gt; &amp;h_omega_1, FullMatrix&lt; double &gt; &amp;h_omega_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const =0</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_h_sigma</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>af947eb73ba87510ec7667a7ae246bbfa</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n, double &amp;h_sigma, Vector&lt; double &gt; &amp;h_sigma_1, FullMatrix&lt; double &gt; &amp;h_sigma_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>get_maximum_step</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>add1852ebe7ad8b1178063ff725748856</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, const Vector&lt; double &gt; &amp;delta_e_omega, const Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x) const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_maximum_step</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a43097fde6921f2dcf85f173540eb1414</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, const Vector&lt; double &gt; &amp;delta_e_sigma, const Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compare_derivatives_with_numerical_derivatives</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ab1cbbd84088b3dae549d152e049240da</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const std::string detailed_printout_file=&quot;&quot;, const double epsilon=1e-8) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>set_q_point_id</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a8c1c9dc1488f810fc7dd0a4233e0b3d4</anchor>
      <arglist>(std::string q_point_id) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ScalarFunctional</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a42b56519b3e4338d248c1e29e29831a6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::set&lt; types::material_id &gt;</type>
      <name>domain_of_integration</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>aa192395f822a64f60df43bf9d36c2f3a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>initial_vals_hidden</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ae3282d5182360e0030e4cc5e02fbe2eb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; DependentField&lt; spacedim, spacedim &gt; &gt;</type>
      <name>e_omega</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>adfed9b70b743ba245a39c3e63b951f96</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; DependentField&lt; spacedim, spacedim &gt; &gt;</type>
      <name>e_sigma</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>afe519345583b2a9d1dcf3a29ea2172f1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const Quadrature&lt; spacedim &gt;</type>
      <name>quadrature</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ab83ee3ae077b211137824b006098382e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_hidden</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a7df6711471715f907bc9911449c5c825</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a195248af3821548af3000872e9e6d00e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_ref_sets</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>acee2c3c289e5b2b680996facc2f79e78</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>q_point_id</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a6704f42ad921bc042b5be661d5737557</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>AssemblyHelper</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>af4019c2e39cc934d646aaa35c3c52773</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ScalarFunctionalLocalElimination</name>
    <filename>class_scalar_functional_local_elimination.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <base>ScalarFunctional</base>
    <member kind="function">
      <type></type>
      <name>ScalarFunctionalLocalElimination</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>ad24025d5ac8d43d2e87479dbd9d70fa3</anchor>
      <arglist>(const std::vector&lt; ScalarFunctional&lt; dim, spacedim &gt; * &gt; scalar_functionals, const std::string name=&quot;ScalarFunctionalLocalElimination&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>get_h_sigma</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a458ebe9d4709117e66fc1e03e220fbf7</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n, double &amp;h_sigma, Vector&lt; double &gt; &amp;h_sigma_1, FullMatrix&lt; double &gt; &amp;h_sigma_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>get_maximum_step</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>abd08e15318c711bd4fdb2a66d378460f</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, const Vector&lt; double &gt; &amp;delta_e_sigma, const Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_safety_distance</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>ae75a23326905e96e1787ef10fd8af930</anchor>
      <arglist>(const double safety_distance)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_safety_distance_step</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a7de6ba97eb27336b77cf1383e79ea3ad</anchor>
      <arglist>(const double safety_step)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_threshold_residual</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a0165ca6758e3818dbd279a109ed65d3c</anchor>
      <arglist>(const double threshold_residual)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_threshold_step_size</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a2c516efc75deb2bfc649b62f1c465760</anchor>
      <arglist>(const double threshold_step_size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_max_iter</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a5b45946eed402b1865271accc58a9dc5</anchor>
      <arglist>(const unsigned int max_iter)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_max_cutbacks</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>ad4e02a446f53fe5ee185e43d99defa51</anchor>
      <arglist>(const unsigned int max_cutbacks)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_use_line_search</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a1a0fb740473442c270f6c76b146f6c61</anchor>
      <arglist>(const bool use_line_search)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_q_point_id</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>ac1bb9770d8d5b60aca7cd71ef86415b1</anchor>
      <arglist>(std::string q_point_id) const override final</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const std::vector&lt; ScalarFunctional&lt; dim, spacedim &gt; * &gt;</type>
      <name>scalar_functionals</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a473e3bc585587839f7f7349f9a25bf28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>map_dependent_fields</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a7a7cbe368936cfb0224e7885a8959c4b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>indices_not_eliminated_dependent_fields</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a1d72179de9c2d77068d090c1b483cd5a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>indices_locally_eliminated_dependent_fields</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a9a62464884767f3be49391e5ec7ec987</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>safety_distance</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a469c7e8629c32b215aba04c805617d9d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>safety_distance_step</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a75b69be34e2a90185083e0a0267781f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>threshold_residual</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a5335c4ea520dfca402e0fb3155e94406</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>threshold_step_size</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>ad9b9862ea6ad4f46f72ee4ae6dfc99b9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>max_iter</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>ad9ee2f4e07556665b70531a917a7e07d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>max_cutbacks</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>af814289c235283b9ce685c5aabb32b51</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>use_line_search</name>
      <anchorfile>class_scalar_functional_local_elimination.html</anchorfile>
      <anchor>a61f8cd4268582affd339a30ae6a5ea6e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ScalarFunctionalLocalElimination&lt; spacedim, spacedim &gt;</name>
    <filename>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</filename>
    <templarg>spacedim</templarg>
    <base>ScalarFunctional&lt; spacedim, spacedim &gt;</base>
    <member kind="function">
      <type></type>
      <name>ScalarFunctionalLocalElimination</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a03d4f058a313e7b904df1cc132b522eb</anchor>
      <arglist>(const std::vector&lt; ScalarFunctional&lt; spacedim, spacedim &gt; * &gt; scalar_functionals, const std::string name=&quot;ScalarFunctionalLocalElimination&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>get_h_omega</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>aa19331788617e6ff84a436fb2479ce03</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, double &amp;h_omega, Vector&lt; double &gt; &amp;h_omega_1, FullMatrix&lt; double &gt; &amp;h_omega_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>get_maximum_step</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a510939a9d4c429a5b779b1dbc8adf26f</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, const Vector&lt; double &gt; &amp;delta_e_omega, const Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_safety_distance</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a53352820c3f0493793196010232ef96a</anchor>
      <arglist>(const double safety_distance)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_safety_distance_step</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>aba38dd2c4de8cd593a00e8811c4be10d</anchor>
      <arglist>(const double safety_step)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_threshold_residual</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ac82db99b2b08cd32f74fadd57ef18e19</anchor>
      <arglist>(const double threshold_residual)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_threshold_step_size</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a2c0c51e69a894bd3196ac1da7f25ca33</anchor>
      <arglist>(const double threshold_step_size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_max_iter</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a7db53d103eae56be4c342d59b241d539</anchor>
      <arglist>(const unsigned int max_iter)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_max_cutbacks</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a9669b2078915bd408a42612e0b68c008</anchor>
      <arglist>(const unsigned int max_cutbacks)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_use_line_search</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ae835e407763fa55661c3405b77be1ef0</anchor>
      <arglist>(const bool use_line_search)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_q_point_id</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a3b922dbb073d2c4649e833ef52626385</anchor>
      <arglist>(std::string q_point_id) const override final</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const std::vector&lt; ScalarFunctional&lt; spacedim, spacedim &gt; * &gt;</type>
      <name>scalar_functionals</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>afff9951bc0ac7a5dfc55c01ffb36177a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>map_dependent_fields</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>aedf208eb036ffafb16871cdfd3d197b1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>indices_not_eliminated_dependent_fields</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ab7e4b409e89577ec660f7c0f83de7096</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>indices_locally_eliminated_dependent_fields</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>aa99a6e0152120b63d6052c4b597c91f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>safety_distance</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a55ab0ce1d909f96e614b84c3ee00bcb3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>safety_distance_step</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a7bb7215106da6851ec4d17aeba85230c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>threshold_residual</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a2a268af4bf9867e51446fbd33c488fe9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>threshold_step_size</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>aa3761ab1e9111468e31135fd954cb964</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>max_iter</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a329bf04de49be40085b0b70bf4d5a039</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>max_cutbacks</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a230b126730b702170f6b7edd1d23c735</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>use_line_search</name>
      <anchorfile>class_scalar_functional_local_elimination_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ac3301f2f6f6e71675ab4d21487484946</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapper</name>
    <filename>class_solver_wrapper.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <templarg></templarg>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a509594953f388e594bfdde5b927ece35</anchor>
      <arglist>(const MatrixType &amp;K_stretched, SolutionVectorType &amp;solution, const RHSVectorType &amp;f_stretched, const bool symmetric)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~SolverWrapper</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a738e65e298fb0b13327e853fc9683b6d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapper&lt; dealii::LinearAlgebra::distributed::Vector&lt; double &gt;, dealii::PETScWrappers::MPI::BlockVector, parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt;, TwoBlockSparsityPattern &gt;</name>
    <filename>class_solver_wrapper.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a509594953f388e594bfdde5b927ece35</anchor>
      <arglist>(const parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt; &amp;K_stretched, dealii::LinearAlgebra::distributed::Vector&lt; double &gt; &amp;solution, const dealii::PETScWrappers::MPI::BlockVector &amp;f_stretched, const bool symmetric)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~SolverWrapper</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a738e65e298fb0b13327e853fc9683b6d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::BlockVector&lt; double &gt;, TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt;, TwoBlockSparsityPattern &gt;</name>
    <filename>class_solver_wrapper.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a509594953f388e594bfdde5b927ece35</anchor>
      <arglist>(const TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::BlockVector&lt; double &gt; &amp;f_stretched, const bool symmetric)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~SolverWrapper</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a738e65e298fb0b13327e853fc9683b6d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::Vector&lt; double &gt;, dealii::SparseMatrix&lt; double &gt;, SparsityPattern &gt;</name>
    <filename>class_solver_wrapper.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a509594953f388e594bfdde5b927ece35</anchor>
      <arglist>(const dealii::SparseMatrix&lt; double &gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::Vector&lt; double &gt; &amp;f_stretched, const bool symmetric)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~SolverWrapper</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a738e65e298fb0b13327e853fc9683b6d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapperMUMPS</name>
    <filename>class_solver_wrapper_m_u_m_p_s.html</filename>
    <base>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::Vector&lt; double &gt;, dealii::SparseMatrix&lt; double &gt;, SparsityPattern &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initialize_matrix</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>aac0b071a4afe90f575f48283f321e9aa</anchor>
      <arglist>(const SparseMatrix&lt; double &gt; &amp;matrix)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initialize_matrix</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ac0e555858181fc6111d6302137f0a235</anchor>
      <arglist>(std::vector&lt; int &gt; &amp;irn, std::vector&lt; int &gt; &amp;jcn, std::vector&lt; double &gt; &amp;A, unsigned int n)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>analyze_matrix</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a0be6cc138caf0deb9b1066071edfbb80</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>factorize_matrix</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ab9ecf7772664524d9de206744441ba80</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>vmult</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>acdb8d522bd9816425d86e813d603f875</anchor>
      <arglist>(Vector&lt; double &gt; &amp;x, const Vector&lt; double &gt; &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SolverWrapperMUMPS</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a394f287e0761996554eed430366c8019</anchor>
      <arglist>(int sym=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SolverWrapperMUMPS</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a23ef40e6c3316b9dedaea66893e3bb29</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ab76c71fc7d8c8f0d7c796bc78e5bdb37</anchor>
      <arglist>(const dealii::SparseMatrix&lt; double &gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::Vector&lt; double &gt; &amp;f_stretched, const bool symmetric=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DeclException2</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a878271784082d91184b16780bab52906</anchor>
      <arglist>(MUMPSError, std::string, int,&lt;&lt; &quot;MUMPS routine &quot;&lt;&lt; arg1&lt;&lt; &quot; returned error status &quot;&lt;&lt; arg2&lt;&lt; &quot;.&quot;)</arglist>
    </member>
    <member kind="variable">
      <type>DMUMPS_STRUC_C</type>
      <name>id</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a919a0be0fa6173351786e0abe565359c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>icntl</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a246f1e553cda5741eede8b97cb658388</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>cntl</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a2d2a35a782b91b320f93e0962baee788</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>info</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>aeddd611d2a3287054fc3926ba7c700b8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>infog</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>ac66bb7798313f6da855617807981e529</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>analyze</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a1bd59a48f0ed4485dd4d67b810bfd701</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>modify_on_negative_pivot</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>aac357825a24e2c32f550abe5b3cfea8c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>beta</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>accf7454076e059b9ea1500e305c98cdb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>increase_tau</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a3882696d9c08572c742978f6ac5de65c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>irn</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>aa1de27d7f99bd632d6a13b6d587bfe48</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>jcn</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a697d2583ca34837591a09a896dd62060</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; int &gt;</type>
      <name>d</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a3a9bacf003d721ff3200688063818548</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; double &gt;</type>
      <name>A</name>
      <anchorfile>class_solver_wrapper_m_u_m_p_s.html</anchorfile>
      <anchor>a67a4f5d3f503d51d0f52ad9c14fe578a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapperPETSc</name>
    <filename>class_solver_wrapper_p_e_t_sc.html</filename>
    <base>SolverWrapper&lt; dealii::LinearAlgebra::distributed::Vector&lt; double &gt;, dealii::PETScWrappers::MPI::BlockVector, parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt;, TwoBlockSparsityPattern &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper_p_e_t_sc.html</anchorfile>
      <anchor>a0b4a570ebbfba3dd2fcf116b3f32948b</anchor>
      <arglist>(const parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt; &amp;K_stretched, dealii::LinearAlgebra::distributed::Vector&lt; double &gt; &amp;solution, const dealii::PETScWrappers::MPI::BlockVector &amp;f_stretched, const bool symmetric=false)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapperPETScIterative</name>
    <filename>class_solver_wrapper_p_e_t_sc_iterative.html</filename>
    <base>SolverWrapper&lt; dealii::LinearAlgebra::distributed::Vector&lt; double &gt;, dealii::PETScWrappers::MPI::BlockVector, parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt;, TwoBlockSparsityPattern &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper_p_e_t_sc_iterative.html</anchorfile>
      <anchor>ade08d6ee2dfb5a6900f99edd8d12fd86</anchor>
      <arglist>(const parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt; &amp;K_stretched, dealii::LinearAlgebra::distributed::Vector&lt; double &gt; &amp;solution, const dealii::PETScWrappers::MPI::BlockVector &amp;f_stretched, const bool symmetric=false)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapperUMFPACK</name>
    <filename>class_solver_wrapper_u_m_f_p_a_c_k.html</filename>
    <base>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::Vector&lt; double &gt;, dealii::SparseMatrix&lt; double &gt;, SparsityPattern &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper_u_m_f_p_a_c_k.html</anchorfile>
      <anchor>a68f5e2055270d172a4ffdf4a8a4f536e</anchor>
      <arglist>(const dealii::SparseMatrix&lt; double &gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::Vector&lt; double &gt; &amp;f_stretched, const bool symmetric=false)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TotalPotential</name>
    <filename>class_total_potential.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type>void</type>
      <name>add_total_potential_contribution</name>
      <anchorfile>class_total_potential.html</anchorfile>
      <anchor>ace8b22eef7789f5c7eb78606ba1f110c</anchor>
      <arglist>(const TotalPotentialContribution&lt; spacedim &gt; &amp;total_potential_contribution)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const TotalPotentialContribution&lt; spacedim &gt; &gt; &gt;</type>
      <name>total_potential_contributions</name>
      <anchorfile>class_total_potential.html</anchorfile>
      <anchor>a5a14ce0e2fabf8116566aa67fb11db35</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>max_dependent_vars</name>
      <anchorfile>class_total_potential.html</anchorfile>
      <anchor>a800f9366116679fd0f7d3173a3bfc539</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>AssemblyHelper</name>
      <anchorfile>class_total_potential.html</anchorfile>
      <anchor>af4019c2e39cc934d646aaa35c3c52773</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TotalPotentialContribution</name>
    <filename>class_total_potential_contribution.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>get_potential_contribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a515786e58c1fda7baf168be8f4a13720</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;H_omega_H_sigma_C, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;C_ref_sets, double &amp;Pi, Vector&lt; double &gt; &amp;Pi_1, FullMatrix&lt; double &gt; &amp;Pi_2, const std::tuple&lt; bool, bool, bool &gt; &amp;requested_quantities) const</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TotalPotentialContribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a1932d6d7c269344542fc8e9221de4448</anchor>
      <arglist>(const std::vector&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; * &gt; &amp;H_omega, const std::vector&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; * &gt; &amp;H_sigma, const std::vector&lt; const IndependentField&lt; 0, spacedim &gt; * &gt; &amp;C)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TotalPotentialContribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>ad64c449f4cae79fb8462f4679c50d606</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim, spacedim &gt; &amp;H_omega)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TotalPotentialContribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a4625ad0c798469b9da9bcb62cc0d00f6</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim-1, spacedim &gt; &amp;H_sigma)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TotalPotentialContribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>af9ecccfc930f6826be4a903194541d5f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const bool</type>
      <name>is_primitive</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a45bfb25a7693c949c26e223cf4a1a1e7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; * &gt;</type>
      <name>H_omega</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a15191539345978a3d0c7293bd7ecaa91</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; * &gt;</type>
      <name>H_sigma</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>aac404e3a8493d9170541e34bd96673d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; const IndependentField&lt; 0, spacedim &gt; * &gt;</type>
      <name>C</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>adea8f8f88243adec43df300e8c8d4593</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>parallel::Triangulation</name>
    <filename>classparallel_1_1_triangulation.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>Triangulation</name>
      <anchorfile>classparallel_1_1_triangulation.html</anchorfile>
      <anchor>abd6ed252f25fe0342e85e7a00aaf9fe3</anchor>
      <arglist>(MPI_Comm mpi_communicator)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>finalize_subdomain_assignment</name>
      <anchorfile>classparallel_1_1_triangulation.html</anchorfile>
      <anchor>ae058708cfbddedd508afce57b6e3b284</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_multilevel_hierarchy_constructed</name>
      <anchorfile>classparallel_1_1_triangulation.html</anchorfile>
      <anchor>abe47f077f6c88940ed25ee571e4586b5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>update_cell_relations</name>
      <anchorfile>classparallel_1_1_triangulation.html</anchorfile>
      <anchor>a9fe60e2814d4d9a5993db43daaf7f64b</anchor>
      <arglist>() override</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>save</name>
      <anchorfile>classparallel_1_1_triangulation.html</anchorfile>
      <anchor>ad722763b534f5ec0b8a768c07ab6f6d9</anchor>
      <arglist>(const std::string &amp;filename) const override</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>load</name>
      <anchorfile>classparallel_1_1_triangulation.html</anchorfile>
      <anchor>ab7745a3a3f606c6526f60038303b1453</anchor>
      <arglist>(const std::string &amp;filename) override</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>load</name>
      <anchorfile>classparallel_1_1_triangulation.html</anchorfile>
      <anchor>a843b296ec1cff052b49efb02ff02681f</anchor>
      <arglist>(const std::string &amp;filename, const bool autopartition) override</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>parallel::TriangulationSystem</name>
    <filename>classparallel_1_1_triangulation_system.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>TriangulationSystem</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a312a7db9393869aa3c9b84b1e316342b</anchor>
      <arglist>(dealii::parallel::distributed::Triangulation&lt; spacedim, spacedim &gt; &amp;tria_domain, const bool fix_vertex_positions=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TriangulationSystem</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>ad12682d17e85ee169fd53cdbe0536f16</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const dealii::parallel::distributed::Triangulation&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_triangulation_domain</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a8dfd7da98adf853fd6ec27cbfccec90b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const dealii::GalerkinTools::parallel::Triangulation&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_triangulation_interface</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>acaac838c580c71c0152b1eef7c084ebb</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual dealii::parallel::distributed::Triangulation&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_triangulation_domain</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>abadd4ec597f43a8992f6477690c9f6c4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual dealii::GalerkinTools::parallel::Triangulation&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_triangulation_interface</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>ad5c8bebc3b72cae7e1ec086790b2f76d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>close</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a5a7ef73bec1fbe698480e1bd7c6cd8a7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_meshes_per_processor_as_vtu</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a6220dec5dd926987eedf084273062f9c</anchor>
      <arglist>(const std::string file_name_domain, const std::string file_name_interface) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::pair&lt; const unsigned int, const unsigned int &gt;</type>
      <name>get_this_proc_n_procs</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>afbde3803f2243308fa33593d6b5f17b4</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual void</type>
      <name>pre_refinement_domain</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a79e9789e83e12900c85cf8de0644271f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual void</type>
      <name>post_refinement_domain</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a951181f2ad877283d458fa19db42efb2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>update_interface_subdomain_ids</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>ae214de5a3053ae320a1eb6d879f03b43</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>MPI_Comm</type>
      <name>mpi_communicator</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>aff7cdcf04d5a4fb633d714130da893b0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TriangulationSystem</name>
    <filename>class_triangulation_system.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>InterfaceCell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a1d62a56e335cf19f4a4ec7932fbaef09</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>DomainCell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>af53de5ec80a16d9cb167660b2832b240</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaActiveIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>ActiveInterfaceCell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a4c3f97884b5478ebd6771433c5273ac5</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaActiveIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>ActiveDomainCell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a6aa358f6803facdc7e7ee72698a5f4cb</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TriangulationSystem</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a2d41b2b77cc2835f81a8584ca24dccad</anchor>
      <arglist>(Triangulation&lt; spacedim, spacedim &gt; &amp;tria_domain, const bool fix_vertex_positions=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TriangulationSystem</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a9b0dc550fe1566ab1ba6668a868d0a70</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Triangulation&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_triangulation_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a70637dc46ee4a36433c3ec3c5dc544c5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Triangulation&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_triangulation_interface</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>af6c6bc5c73c44896c976bf8e118fa67c</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Triangulation&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_triangulation_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a4c4e1d0589f962d92532586d808c7da2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Triangulation&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_triangulation_interface</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a73ef25cdc88395ef6d1ebcf6792cb97f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>close</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>ad3605f3f59fbd55942288026107c4e6d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_interface_cell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a004ced64c3160c7589b6e7bf12f58019</anchor>
      <arglist>(const DomainCell &amp;domain_cell, const unsigned int face, const types::material_id material_id)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_interface_cells</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a2bb42253dc27cab092c7c76c66b83585</anchor>
      <arglist>(std::vector&lt; std::tuple&lt; const DomainCell, const unsigned int, const types::material_id &gt; &gt; interface_cells)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_interface_manifold</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a54803f193d6f237b2039ba6f5bca09b9</anchor>
      <arglist>(const types::manifold_id manifold_id, const Manifold&lt; spacedim-1, spacedim &gt; &amp;manifold)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_interface_manifolds</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>afb5564d01d3914b15b552f95648d231b</anchor>
      <arglist>(const std::map&lt; types::manifold_id, std::reference_wrapper&lt; const Manifold&lt; spacedim-1, spacedim &gt;&gt;&gt; manifolds)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_begin_active</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>abbbd02a7e813a0604d072480bda93b1a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_end_active</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>abcab5ca315d720879468cf92a820c29a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt; &amp;</type>
      <name>interface_active_iterators</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>ac46573e9f30edb21ca6f39770ceceaee</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_begin_coarse</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aaf68847ea61f015fce319bb9810833f4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_end_coarse</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a1595839f52f5f5b9d7af31836f9c1f0f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt; &amp;</type>
      <name>interface_coarse_iterators</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>afa664696e5d01f0c4b408251311f13b2</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_triangulations_vtk</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>ab8700003a38ead6d2b7cadc896b9c94c</anchor>
      <arglist>(const std::string file_name_domain, const std::string file_name_interface) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>refine_global</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a8aa29c77fd798ddd8a5e20806ea21558</anchor>
      <arglist>(const unsigned int times=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>execute_coarsening_and_refinement</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a62be2563cc8d810a71941e15490f9840</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>check_active_interface_cell_domain_cells_consistency</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a6b03827603784e9ce86ddd5484fa52be</anchor>
      <arglist>(const double tol=1e-12) const</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::pair&lt; const unsigned int, const unsigned int &gt;</type>
      <name>get_this_proc_n_procs</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a1d468c78ec2c3c57b66bf2c307201eea</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable">
      <type>boost::signals2::signal&lt; void()&gt;</type>
      <name>pre_refinement</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>ad3b16b86e7cbe4800fb5a40114cbc5fb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>boost::signals2::signal&lt; void()&gt;</type>
      <name>post_refinement</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a8813b73cab76c8b77c47c00c39233f76</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>generate_tria_interface_from_tria_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aee020d825bdf108d69f8e5b71b7ad665</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const SmartPointer&lt; Triangulation&lt; spacedim, spacedim &gt; &gt;</type>
      <name>tria_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a68bafddc70652cb7c64701c74e86279f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::unique_ptr&lt; Triangulation&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>tria_interface</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>af145fb3f6a3b3ef9ed92111f9fca5b7d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>closed</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aef5ae937cc00b7954357bb8ae88b1f73</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>generate_active_interface_cells_domain_cells</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a4214cf789cf39910305db9fe554249a0</anchor>
      <arglist>(const bool no_assert=false)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>generate_active_interface_cells_domain_cells_recursion</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a449457dd44d9c9d9aaa5754532e8537e</anchor>
      <arglist>(const DomainCell &amp;domain_cell, const unsigned int &amp;face, const InterfaceCell &amp;interface_cell, const bool no_assert)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>check_material_ids_recursion</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a0dfaa2e780dbb0f750b1a85ed9428740</anchor>
      <arglist>(const DomainCell &amp;domain_cell) const</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual void</type>
      <name>pre_refinement_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>ae1862e6da3157dc8d539fdc0439e9f48</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual void</type>
      <name>post_refinement_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a8435489384095f687363d200ccfce628</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>push_tria_listeners_to_end</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aca3dc94f1006a603f56ac12c730ae242</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; std::pair&lt; const DomainCell, const unsigned int &gt;, types::material_id &gt;</type>
      <name>coarse_domain_faces_material_ids</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a9a14a8a68bda679330a169413f0e2dfb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; dealii::types::manifold_id, const dealii::Manifold&lt; spacedim-1, spacedim &gt; * &gt;</type>
      <name>interface_manifolds</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a130322b74b66f27d8611ce92ad241c6d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;</type>
      <name>coarse_interface_cell_domain_cells</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a5b909c472fe2508da44e89b017aa146c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;</type>
      <name>active_interface_cell_domain_cells</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a516c7a253cefbc5e208714538b21424d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; boost::signals2::connection &gt;</type>
      <name>tria_listeners</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a713f97bee5570de1f5571d86a41e83f6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const bool</type>
      <name>fix_vertex_positions</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a74bfeabee5174d77a8477ba9d2de4813</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>AssemblyHelper</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>af4019c2e39cc934d646aaa35c3c52773</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>parallel::TwoBlockMatrix</name>
    <filename>classparallel_1_1_two_block_matrix.html</filename>
    <templarg></templarg>
    <member kind="function">
      <type></type>
      <name>TwoBlockMatrix</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>a93990616301e671f9a967fbc49029674</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TwoBlockMatrix</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>afa9b3f1c3d74a211c334c8f6e59b52f6</anchor>
      <arglist>(const TwoBlockSparsityPattern &amp;sp, const IndexSet &amp;locally_owned_indices, const MPI_Comm mpi_communicator=MPI_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>a0e112695034c69aa8a430efd9e7c37ef</anchor>
      <arglist>(const TwoBlockSparsityPattern &amp;sp, const IndexSet &amp;locally_owned_indices, const MPI_Comm mpi_communicator=MPI_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type>dealii::GalerkinTools::parallel::TwoBlockMatrix&lt; MatrixType &gt; &amp;</type>
      <name>operator=</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>a8e2ba56ebc9ee5ad1fe7b1d8604e9669</anchor>
      <arglist>(const double value)</arglist>
    </member>
    <member kind="function">
      <type>const MPI_Comm &amp;</type>
      <name>get_communicator</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>add0cc2ee49d139e1dbc5554fd62dd5f5</anchor>
      <arglist>() const</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TwoBlockMatrix</name>
    <filename>class_two_block_matrix.html</filename>
    <templarg></templarg>
    <member kind="typedef">
      <type>double</type>
      <name>value_type</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a9f50fb1f98df3ae0048c3c5220aa265f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TwoBlockMatrix</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>aa74bb5722415e61241b5f38c9df29d39</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TwoBlockMatrix</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a861cdc13c0fc6a65ec7d99e4d2349a92</anchor>
      <arglist>(const TwoBlockSparsityPattern &amp;sp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TwoBlockMatrix</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a0674b188ebc7bde3ee8f5eb9a94d6265</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a5eaa76e9909c600aaa256aa77b577391</anchor>
      <arglist>(const TwoBlockSparsityPattern &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>m</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>afedb959cb6812e0870abb96a36762b4f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>af05649a5859d458e3961a7aa9bafe489</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a9b5d922689d1f60cbc59d2e424df2dc4</anchor>
      <arglist>(const unsigned int i, const unsigned int j, const double value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a7373229281e1c3712681c1ac742dffb7</anchor>
      <arglist>(const unsigned int row, const unsigned int n_cols, const unsigned int *col_indices, const double *values, const bool elide_zero_values=false, const bool col_indices_are_sorted=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a5457fd3b6228bed7a731bab2569a17a0</anchor>
      <arglist>(const unsigned int i, const unsigned int j, const double value)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a5da32a6ea4b80806f03ebfa462f67169</anchor>
      <arglist>(const unsigned int i, const unsigned int j) const</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>el</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ab2a9f1280987b31fd2010c9fdaad4891</anchor>
      <arglist>(const unsigned int i, const unsigned int j) const</arglist>
    </member>
    <member kind="function">
      <type>TwoBlockMatrix&lt; MatrixType &gt; &amp;</type>
      <name>operator=</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ac008679a935fb6577931058ac0a6be36</anchor>
      <arglist>(const double value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compress</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a8de26ebbb410811059d162fb164cedce</anchor>
      <arglist>(const VectorOperation::values operation)</arglist>
    </member>
    <member kind="function">
      <type>const MatrixType &amp;</type>
      <name>get_A</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a7b410c6140d2c2b11b2cd05a1f8c732b</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const MatrixType &amp;</type>
      <name>get_B</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>aec990a6ce6f20f155cb6653a235ded95</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const MatrixType &amp;</type>
      <name>get_C</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a081c82ffbd9e70aff577088cbf501913</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const MatrixType &amp;</type>
      <name>get_D</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a4a653266b5592ebc020b954a2e8927bb</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>MatrixType &amp;</type>
      <name>get_A</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ab30f0e3b8f803456c3a2f1ce4f753081</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MatrixType &amp;</type>
      <name>get_B</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ae41ca5bbc0495c8ffb1604da488a33b7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MatrixType &amp;</type>
      <name>get_C</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>aeaa4102bb137700d7bbe0262ee16e8ee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MatrixType &amp;</type>
      <name>get_D</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a3725800849e26e2192c0ac6520777dab</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_block_0_size</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ab698897508d282162ee05e54ab59bb5f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_block_1_size</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a8d2bf6b03f753b56f9a632d3cda926b7</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>MatrixType</type>
      <name>A</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a91b02f00cf56b398e1624227e6efc8c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>MatrixType</type>
      <name>B</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ae8edfe1cd2298741586dcc46a4d09ebf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>MatrixType</type>
      <name>C</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a847505a77448124928540be0367d8747</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>MatrixType</type>
      <name>D</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a72c598f36d3e8cce12821a09920ded2d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SparsityPattern</type>
      <name>sp_A</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>af40216911c197c73687f7d3a1ff4e68e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SparsityPattern</type>
      <name>sp_B</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a854c7e77623a7c0012aa0ceeb60c147a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SparsityPattern</type>
      <name>sp_C</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a45d6e534ad775577c86e9ac9a116014e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SparsityPattern</type>
      <name>sp_D</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>adb1833352096be79279134b218659665</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned int</type>
      <name>block_0_size</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ae862db333ed43f1650d4926f237cfc43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned int</type>
      <name>block_1_size</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a60ddea1d26313aac35cfb2ae9f192cb2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned int</type>
      <name>total_dimension</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a7d63b4c4c1e499c9f9913444dc5fd6b6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TwoBlockSparsityPattern</name>
    <filename>class_two_block_sparsity_pattern.html</filename>
    <member kind="function">
      <type></type>
      <name>TwoBlockSparsityPattern</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a5bcc54c1c6b0d769807644f1c87fa136</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TwoBlockSparsityPattern</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ab7240dc72bf5fdb32fc192f70e6c68b0</anchor>
      <arglist>(const IndexSet &amp;locally_relevant_indices, const unsigned int block_0_size)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TwoBlockSparsityPattern</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>abe3cbf4d3ee10b66234c5f8b566a938f</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ae18398b187a0dba1570abfec1a4e082d</anchor>
      <arglist>(const IndexSet &amp;locally_relevant_indices, const unsigned int block_0_size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a8c9a56a20123e8eaea65f9eb6efed250</anchor>
      <arglist>(const AssemblyHelper&lt; spacedim &gt; &amp;assembly_helper)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_rows</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ab2b3683a888bedc713117771f4e052fb</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_cols</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a800e2b70b07117cdf330c41a5e9532ca</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ab3787657efc568999fbf24cb5208d46b</anchor>
      <arglist>(const unsigned int i, const unsigned int j)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_entries</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a97c69542e0e22e5653c8b4f7b4275635</anchor>
      <arglist>(const unsigned int row, ForwardIterator begin, ForwardIterator end, const bool indices_are_sorted=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>finalize</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a051b557d53ca8367228b7829c55d4af5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const SparsityPattern &amp;</type>
      <name>get_sp_A</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ac4e6c53ad297cfe1b9e0f41550d7345c</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const SparsityPattern &amp;</type>
      <name>get_sp_B</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a7826ec9b684468fa09a8b60d31c6188f</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const SparsityPattern &amp;</type>
      <name>get_sp_C</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a5dac1b5718f7d001fe3e04f8744601fd</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>const SparsityPattern &amp;</type>
      <name>get_sp_D</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ad23448c9eb58843e1c0cc3403664b6f5</anchor>
      <arglist>() const</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>exists</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a642806909b31d0cb75ca4831563c8269</anchor>
      <arglist>(const unsigned int i, const unsigned int j) const</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>distribute</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a61ef6e3f3325cffcfc1d25efc41886b0</anchor>
      <arglist>(const IndexSet &amp;locally_owned_indices, const MPI_Comm mpi_communicator)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DynamicSparsityPattern</type>
      <name>dsp_A</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ac3bafc713fec0ec2232fc084f880c50c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DynamicSparsityPattern</type>
      <name>dsp_B</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ab1b2d6dc5d92f6353b80dd7a8d32ffa4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DynamicSparsityPattern</type>
      <name>dsp_C</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a444a6465b3ba1355fb675128c6e9800a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DynamicSparsityPattern</type>
      <name>dsp_D</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a853d4000754f009a66853954d8c540bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SparsityPattern</type>
      <name>sp_A</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ad09712fc57a8cf5a2cd2cb9fc238a7cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SparsityPattern</type>
      <name>sp_B</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a95e817472f023d7237c752a18a24ac73</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SparsityPattern</type>
      <name>sp_C</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>aee8da566a426045fb8101dba91295a27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SparsityPattern</type>
      <name>sp_D</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a7a820cba98001e45ba09551a306b14d7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>block_0_size</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a81789900ac1a9f5bfd7f85ad730eb03b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>block_1_size</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ab215e3553e07b5f5872b1e01228968f2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>total_dimension</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>aa560457e8b8bf7a06acf6b19987c5ee0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>finalized</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a7c930cd0dd5e1d6fb1550d6c3939c7aa</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>Auxiliary</name>
    <filename>namespace_auxiliary.html</filename>
    <member kind="function">
      <type>int</type>
      <name>compute_ldr</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>abe796c1529d11eed08fd05bca82f3002</anchor>
      <arglist>(FullMatrix&lt; double &gt; &amp;C, Vector&lt; double &gt; &amp;D, std::vector&lt; Vector&lt; double &gt;&gt; &amp;L, std::vector&lt; Vector&lt; double &gt;&gt; &amp;R)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>convert_local_indices_to_global_indices</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>ad3c42d209f0ba8e9d4ce305060634bf1</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;dof_indices_local, std::vector&lt; unsigned int &gt; &amp;dof_indices_global, const std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>combine_dof_indices</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a1d90ebc8738df3d8c70b540034137019</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_interface, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_minus, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_plus, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_C, std::vector&lt; unsigned int &gt; &amp;dof_indices_interface_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_minus_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_plus_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_C_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_global_combined)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_to_index_set</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aca5fe5966aef4ba2293fb44095bfd86c</anchor>
      <arglist>(const DoFRenumbering &amp;dof_renumbering, const IndexSet &amp;in, IndexSet &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>renumber_constraints</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aa6148bcbaf5e3003717b4dd2a4da15b3</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const DoFRenumbering &amp;dof_renumbering=DoFRenumbering(), const bool close=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_dof_renumbering_contiguous</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a4261fc726ff965166f3f53242918acea</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system, DoFRenumberingOffset &amp;dof_renumbering_offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_map_dofs</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a072f85e6d745ae3c532bb0612f4bd3ce</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dhs_1, const DoFHandlerSystem&lt; spacedim &gt; &amp;dhs_2, std::vector&lt; unsigned int &gt; &amp;map_dofs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_vector</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a59c01a6511bffd7442693e86cd194ef1</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;in, Vector&lt; double &gt; &amp;out_0, Vector&lt; double &gt; &amp;out_1, const unsigned int size_1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_matrix</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>af746a1d08b1135a3684fc990f7b1384d</anchor>
      <arglist>(const FullMatrix&lt; double &gt; &amp;in, FullMatrix&lt; double &gt; &amp;out_00, FullMatrix&lt; double &gt; &amp;out_01, FullMatrix&lt; double &gt; &amp;out_10, FullMatrix&lt; double &gt; &amp;out_11, const unsigned int size_1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>communicate_bool</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>ad0ff28386be7b54b3487ae36e5a074fa</anchor>
      <arglist>(const bool local_bool, const MPI_Comm &amp;mpi_communicator)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>get_n_locally_owned_dofs_per_processor</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aaa9fb90b2f2308f766b8130578ab9c39</anchor>
      <arglist>(const DoFHandler&lt; dim, spacedim &gt; &amp;dof_handler)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>parallel</name>
    <filename>namespaceparallel.html</filename>
    <class kind="class">parallel::Triangulation</class>
    <class kind="class">parallel::TriangulationSystem</class>
    <class kind="class">parallel::TwoBlockMatrix</class>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>GalerkinTools library</title>
    <filename>index.html</filename>
  </compound>
</tagfile>
