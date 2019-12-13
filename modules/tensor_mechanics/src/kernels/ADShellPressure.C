//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADShellPressure.h"

registerADMooseObject("TensorMechanicsApp", ADShellPressure);

/**
 * This kernel defines the residual contribution from a gravitational body force
 */
defineADValidParams(
    ADShellPressure,
    ADKernelValue,
    params.addClassDescription("Apply gravity. Value is in units of acceleration.");
    params.addParam<bool>("use_displaced_mesh", true, "Displaced mesh defaults to true");
    params.addRequiredParam<Real>(
        "value", "Value multiplied against the residual, e.g. gravitational acceleration");
    params.addRequiredParam<unsigned int>(
        "component", "Value multiplied against the residual, e.g. gravitational acceleration"););

template <ComputeStage compute_stage>
ADShellPressure<compute_stage>::ADShellPressure(const InputParameters & parameters)
  : ADKernelValue<compute_stage>(parameters),
    _value(getParam<Real>("value")),
    _component(getParam<unsigned int>("component"))
{
  printf("enter here \n");
  std::vector<const Node *> nodes(4);
  for (unsigned int i = 0; i < 4; ++i)
    nodes[i] = _current_elem->node_ptr(i);
    printf("enter here 2\n");

  ADRealVectorValue x = (*nodes[1] - *nodes[0]);
  ADRealVectorValue y = (*nodes[3] - *nodes[0]);
  _normal = x.cross(y);
  _normal /= _normal.norm();
  printf("exit here \n");
}

template <ComputeStage compute_stage>
ADReal
ADShellPressure<compute_stage>::precomputeQpResidual()
{
  return -_value * _normal(_component);
}
