[Mesh]
  type = GeneratedMesh
  dim = 1
[]

[Problem]
  type = FEProblem
  solve = false
[]

[Executioner]
  type = Steady
[]

[MultiApps]
  [sub]
    type = FullSolveMultiApp
    positions = '0 0 0
                 1 1 1'
    input_files = 'sub.i'
  []
[]

[Controls]
  [cmdline]
    type = MultiAppFunctionCommandLineControl
    multi_app = sub
    function = linear_test
    param_name = 'Mesh/nx'
  []
[]

[Functions]
  [./linear_test]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '1 3 5'
    axis = 0
  [../]
[]
