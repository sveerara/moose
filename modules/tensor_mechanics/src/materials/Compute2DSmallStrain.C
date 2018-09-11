//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Compute2DSmallStrain.h"

#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<Compute2DSmallStrain>()
{
  InputParameters params = validParams<ComputeSmallStrain>();
  params.addClassDescription("Compute a small strain in a plane strain configuration.");

  MooseEnum outOfPlaneDirection("x y z", "z");
  params.addParam<MooseEnum>(
      "out_of_plane_direction", outOfPlaneDirection, "The direction of the out-of-plane strain.");
  params.addParam<bool>("legacy_volumetric_locking_correction", false, "Older version of 2D volumetric locking correction to compare results against the solid mechanics version");
  return params;
}

Compute2DSmallStrain::Compute2DSmallStrain(const InputParameters & parameters)
  : ComputeSmallStrain(parameters),
    _out_of_plane_direction(getParam<MooseEnum>("out_of_plane_direction")),
    _legacy_volumetric_locking_correction(getParam<bool>("legacy_volumetric_locking_correction")),
    _in_plane_direction(2),
    _ave_zz_strain(false)
{
  if (_out_of_plane_direction == 0)
  {
    _in_plane_direction[0] = 1;
    _in_plane_direction[1] = 2;
  }
  else if(_out_of_plane_direction == 1)
  {
    _in_plane_direction[0] = 0;
    _in_plane_direction[1] = 2;
  }
  else if(_out_of_plane_direction == 2)
  {
    _in_plane_direction[0] = 0;
    _in_plane_direction[1] = 1;
  }
}

void
Compute2DSmallStrain::initialSetup()
{
  displacementIntegrityCheck();
  for (unsigned int i = 0; i < 3; ++i)
  {
    if (_out_of_plane_direction == i)
    {
      _disp[i] = &_zero;
      _grad_disp[i] = &_grad_zero;
    }
    else
    {
      _disp[i] = &coupledValue("displacements", i);
      _grad_disp[i] = &coupledGradient("displacements", i);
    }
  }
}

void
Compute2DSmallStrain::computeProperties()
{
  Real volumetric_strain = 0.0;
  Real out_of_plane_strain = 0.0;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    _total_strain[_qp](_in_plane_direction[0], _in_plane_direction[0]) = (*_grad_disp[_in_plane_direction[0]])[_qp](_in_plane_direction[0]);
    _total_strain[_qp](_in_plane_direction[0], _in_plane_direction[1]) = ((*_grad_disp[_in_plane_direction[0]])[_qp](_in_plane_direction[1]) + (*_grad_disp[_in_plane_direction[1]])[_qp](_in_plane_direction[0])) / 2.0;
    _total_strain[_qp](_in_plane_direction[1], _in_plane_direction[0]) = _total_strain[_qp](_in_plane_direction[0], _in_plane_direction[1]); // force the symmetrical strain tensor
    _total_strain[_qp](_in_plane_direction[1], _in_plane_direction[1]) = (*_grad_disp[_in_plane_direction[1]])[_qp](_in_plane_direction[1]);
    _total_strain[_qp](_out_of_plane_direction, _out_of_plane_direction) = computeOutOfPlaneStrain();

  /*  switch (_out_of_plane_direction)
    {
      case 0: // x
      {
        _total_strain[_qp](0, 0) = computeOutOfPlaneStrain();
        _total_strain[_qp](1, 1) = (*_grad_disp[1])[_qp](1);
        _total_strain[_qp](1, 2) = ((*_grad_disp[1])[_qp](2) + (*_grad_disp[2])[_qp](1)) / 2.0;
        _total_strain[_qp](2, 1) = _total_strain[_qp](1, 2); // force the symmetrical strain tensor
        _total_strain[_qp](2, 2) = (*_grad_disp[2])[_qp](2);
        break;
      }
      case 1: // y
      {
        _total_strain[_qp](0, 0) = (*_grad_disp[0])[_qp](0);
        _total_strain[_qp](0, 2) = ((*_grad_disp[0])[_qp](2) + (*_grad_disp[2])[_qp](0)) / 2.0;
        _total_strain[_qp](1, 1) = computeOutOfPlaneStrain();
        _total_strain[_qp](2, 0) = _total_strain[_qp](0, 2); // force the symmetrical strain tensor
        _total_strain[_qp](2, 2) = (*_grad_disp[2])[_qp](2);
        break;
      }
      default: // z
      {
        _total_strain[_qp](0, 0) = (*_grad_disp[0])[_qp](0);
        _total_strain[_qp](1, 1) = (*_grad_disp[1])[_qp](1);
        _total_strain[_qp](0, 1) = ((*_grad_disp[0])[_qp](1) + (*_grad_disp[1])[_qp](0)) / 2.0;
        _total_strain[_qp](1, 0) = _total_strain[_qp](0, 1); // force the symmetrical strain tensor
        _total_strain[_qp](2, 2) = computeOutOfPlaneStrain();
        break;
      }
    } */

    if (_volumetric_locking_correction)
      volumetric_strain +=
          (_total_strain[_qp](_in_plane_direction[0], _in_plane_direction[0]) + _total_strain[_qp](_in_plane_direction[1], _in_plane_direction[1])) * _JxW[_qp] * _coord[_qp];
    else if (_legacy_volumetric_locking_correction)
      volumetric_strain += _total_strain[_qp].trace() * _JxW[_qp] * _coord[_qp];

    if (_ave_zz_strain && !_legacy_volumetric_locking_correction)
      out_of_plane_strain += _total_strain[_qp](_out_of_plane_direction, _out_of_plane_direction) * _JxW[_qp] * _coord[_qp];
  }

  if (_volumetric_locking_correction || _legacy_volumetric_locking_correction)
    volumetric_strain /= _current_elem_volume;

  if (_ave_zz_strain && !_legacy_volumetric_locking_correction)
    out_of_plane_strain /= _current_elem_volume;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    if (_volumetric_locking_correction)
    {
      const Real trace_2D = _total_strain[_qp](_in_plane_direction[0], _in_plane_direction[0]) + _total_strain[_qp](_in_plane_direction[1], _in_plane_direction[1]);
      _total_strain[_qp](_in_plane_direction[0], _in_plane_direction[0]) += (volumetric_strain - trace_2D) / 2.0;
      _total_strain[_qp](_in_plane_direction[1], _in_plane_direction[1]) += (volumetric_strain - trace_2D) / 2.0;
    }
    else if (_legacy_volumetric_locking_correction)
    {
      Real trace = _total_strain[_qp].trace();
      _total_strain[_qp](0, 0) += (volumetric_strain - trace) / 3.0;
      _total_strain[_qp](1, 1) += (volumetric_strain - trace) / 3.0;
      _total_strain[_qp](2, 2) += (volumetric_strain - trace) / 3.0;
    }

    if (_ave_zz_strain && !_legacy_volumetric_locking_correction)
      _total_strain[_qp](_out_of_plane_direction, _out_of_plane_direction) = out_of_plane_strain;

    _mechanical_strain[_qp] = _total_strain[_qp];

    // Remove the eigenstrains
    for (auto es : _eigenstrains)
      _mechanical_strain[_qp] -= (*es)[_qp];
  }
}

void
Compute2DSmallStrain::displacementIntegrityCheck()
{
  if (_out_of_plane_direction != 2 && _ndisp != 3)
    mooseError("For 2D simulations where the out-of-plane direction is x or y the number of "
               "supplied displacements must be three.");
  else if (_out_of_plane_direction == 2 && _ndisp != 2)
    mooseError("For 2D simulations where the out-of-plane direction is z the number of supplied "
               "displacements must be two.");
}
