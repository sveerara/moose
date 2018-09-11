//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StressDivergenceRZTensors.h"
#include "MooseMesh.h"
#include "libmesh/quadrature.h"

registerMooseObject("TensorMechanicsApp", StressDivergenceRZTensors);

template <>
InputParameters
validParams<StressDivergenceRZTensors>()
{
  InputParameters params = validParams<StressDivergence2DTensors>();
  params.addClassDescription(
      "Calculate stress divergence for an axisymmetric problem in cylinderical coordinates.");
  return params;
}

StressDivergenceRZTensors::StressDivergenceRZTensors(const InputParameters & parameters)
  : StressDivergence2DTensors(parameters), _first(!_fe_problem.mesh().hasSecondOrderElements())
{
}

void
StressDivergenceRZTensors::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_RZ)
    mooseError("The coordinate system in the Problem block must be set to RZ for axisymmetric "
               "geometries.");
}

void
StressDivergenceRZTensors::computeAverageGradientZZTest()
{
  if (_first || _legacy_volumetric_locking_correction)
  {
    _avg_grad_zz_test.resize(_test.size());
    for (_i = 0; _i < _test.size(); ++_i)
    {
      _avg_grad_zz_test[_i] = 0.0;
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
        _avg_grad_zz_test[_i] += _test[_i][_qp] / _q_point[_qp](0) * _JxW[_qp] * _coord[_qp];

      _avg_grad_zz_test[_i] /= _current_elem_volume;
    }
  }
}

void
StressDivergenceRZTensors::computeAverageGradientZZPhi()
{
  if (_first || _legacy_volumetric_locking_correction)
  {
    _avg_grad_zz_phi.resize(_phi.size());
    for (_i = 0; _i < _phi.size(); ++_i)
    {
      _avg_grad_zz_phi[_i] = 0.0;
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
        _avg_grad_zz_phi[_i] += _phi[_i][_qp] / _q_point[_qp](0) * _JxW[_qp] * _coord[_qp];

      _avg_grad_zz_phi[_i] /= _current_elem_volume;
    }
  }
}

Real
StressDivergenceRZTensors::getGradientZZTest()
{
  if (_first && !_legacy_volumetric_locking_correction)
    return _avg_grad_zz_test[_i];
  return _test[_i][_qp] / _q_point[_qp](0);
}

Real
StressDivergenceRZTensors::getGradientZZPhi()
{
  if (_first && !_legacy_volumetric_locking_correction)
    return _avg_grad_zz_phi[_j];
  return _phi[_j][_qp] / _q_point[_qp](0);
}
