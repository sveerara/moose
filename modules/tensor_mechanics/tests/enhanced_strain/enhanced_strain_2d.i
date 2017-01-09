#This is a test of the piece-wise linear strain hardening model using the small strain formulation.
#The exact same problem was run in Abaqus with exactly the same result.

[Mesh]
  file = 2d_rect_fine.e
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
[]

[AuxKernels]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy
    execute_on = timestep_end
  [../]
 []

[NodalKernels]
  [./force1]
    type = ConstantRate
    rate = -1000.0
    boundary = 2
    variable = disp_x
  [../]
  [./force2]
    type = ConstantRate
    rate = 1000.0
    boundary = 3
    variable = disp_x
  [../]
#  [./force3]
#    type = ConstantRate
#    rate = -150.0
#    boundary = 5
#    variable = disp_y
#  [../]
[]

[BCs]
  [./x_left]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]

  [./y_left]
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
#  [./x_right_top]
#    type = DirichletBC
#    variable = disp_x
#    boundary = 2
#    value = 0.0
#  [../]
#  [./x_right_bottom]
#    type = DirichletBC
#    variable = disp_x
#    boundary = 3
#    value = 0.0
#  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1500.0
    poissons_ratio = 0.25
    block = 1
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    enhanced_strain = true
    block = 1
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 1
  [../]
[]

[Executioner]
  type = Steady

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  nl_rel_tol = 1e-4
  l_tol = 1e-4
  nl_abs_tol = 1e-4

[]

[Postprocessors]
  [./disp_1]
    type = NodalVariableValue
    nodeid = 0
    variable = disp_y
  [../]
  [./disp_2]
    type = NodalVariableValue
    nodeid = 3
    variable = disp_y
  [../]
[]

[Outputs]
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
[]
