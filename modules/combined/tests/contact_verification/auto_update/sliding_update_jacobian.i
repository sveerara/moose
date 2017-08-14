[Mesh]
  file = sliding_block_30000.e
  displacements = 'disp_x disp_y'
  automatic_patch_update =true
  priority_queue=true
[]

[Problem]
  update_jacobian_preallocation=true
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Kernels]
  [./TensorMechanics]
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e5
    poissons_ratio = 0.3
  [../]
  [./strain]
    type = ComputeIncrementalSmallStrain
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
  [../]
[]

[Contact]
  [./leftright]
    slave = 3
    master = 2
    system = constraint
    normalize_penalty = true
    tangential_tolerance = 1e-3
    penalty = 1e+6
    model = frictionless
    formulation = penalty
  [../]
[]

[Executioner]
   type = Transient
   solve_type = 'PJFNK'
   start_time = 0
   end_time = 0.3
   l_tol = 1e-8
   nl_rel_tol = 1e-6
   nl_abs_tol = 1e-4
   dt = 0.1
   line_search = 'none'
   petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
   petsc_options_value = 'lu superlu_dist'
   timestep_tolerance = 1e-12
[]

[BCs]
  [./fixed_1_2x]
    type = DirichletBC
    boundary = '1'
    value = 0.0
    variable = disp_x
  [../]
  [./fixed_1_2y]
    type = DirichletBC
    boundary = '1'
    value = 0.0
    variable = disp_y
  [../]
  [./sliding_1]
    type = FunctionPresetBC
    function = sliding_fn
    variable = disp_x
    boundary = '4'
  [../]
  [./normal_y]
    type = PresetBC
    variable = disp_y
    boundary = '4'
    value = -0.01
  [../]
#  [./Pressure]
#    [./normal_pressure]
#      disp_x = disp_x
#      disp_y = disp_y
#      factor = 100.0
#      boundary = 4
#    [../]
#  [../]
[]

[Functions]
  [./sliding_fn]
    type = ParsedFunction
    value = 't'
  [../]
[]

[Outputs]
  exodus = true
  print_perf_log = true
[]
