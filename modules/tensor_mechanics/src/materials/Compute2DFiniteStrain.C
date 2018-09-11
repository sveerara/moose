//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Compute2DFiniteStrain.h"

#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<Compute2DFiniteStrain>()
{
  InputParameters params = validParams<ComputeFiniteStrain>();
  params.addClassDescription(
      "Compute a strain increment and rotation increment for finite strains in 2D geometries.");

  MooseEnum outOfPlaneDirection("x y z", "z");
  params.addParam<MooseEnum>(
      "out_of_plane_direction", outOfPlaneDirection, "The direction of the out-of-plane strain.");
  params.addParam<bool>("legacy_volumetric_locking_correction", false, "Older version of 2D volumetric locking correction to compare results against the solid mechanics version");
  return params;
}

Compute2DFiniteStrain::Compute2DFiniteStrain(const InputParameters & parameters)
  : ComputeFiniteStrain(parameters),
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
Compute2DFiniteStrain::initialSetup()
{
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

    if (_fe_problem.isTransient() && i != _out_of_plane_direction)
      _grad_disp_old[i] = &coupledGradientOld("displacements", i);
    else
      _grad_disp_old[i] = &_grad_zero;
  }
}

void
Compute2DFiniteStrain::computeProperties()
{
  RankTwoTensor ave_Fhat;
  Real ave_dfgrd_det = 0.0;
  Real ave_dfgrd_22 = 0.0;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {

    // Deformation gradient calculation for 2D problems
    RankTwoTensor A((*_grad_disp[0])[_qp],
                    (*_grad_disp[1])[_qp],
                    (*_grad_disp[2])[_qp]); // Deformation gradient
    RankTwoTensor Fbar((*_grad_disp_old[0])[_qp],
                       (*_grad_disp_old[1])[_qp],
                       (*_grad_disp_old[2])[_qp]); // Old Deformation gradient

    // Compute the displacement gradient for the out of plane direction for plane strain,
    // generalized plane strain, or axisymmetric problems

    A(_out_of_plane_direction, _out_of_plane_direction) = computeOutOfPlaneGradDisp();
    Fbar(_out_of_plane_direction, _out_of_plane_direction) = computeOutOfPlaneGradDispOld();

    // Gauss point deformation gradient
    _deformation_gradient[_qp] = A;
    _deformation_gradient[_qp].addIa(1.0);

    A -= Fbar; // very nearly A = gradU - gradUold

    Fbar.addIa(1.0); // Fbar = ( I + gradUold)

    // Incremental deformation gradient _Fhat = I + A Fbar^-1
    _Fhat[_qp] = A * Fbar.inverse();
    _Fhat[_qp].addIa(1.0);

    if (_volumetric_locking_correction)
    {
      // Calculate average _Fhat for volumetric locking correction
      ave_Fhat(_in_plane_direction[0], _in_plane_direction[0]) += _Fhat[_qp](_in_plane_direction[0], _in_plane_direction[0]) * _JxW[_qp] * _coord[_qp];
      ave_Fhat(_in_plane_direction[0], _in_plane_direction[1]) += _Fhat[_qp](_in_plane_direction[0], _in_plane_direction[1]) * _JxW[_qp] * _coord[_qp];
      ave_Fhat(_in_plane_direction[1], _in_plane_direction[0]) += _Fhat[_qp](_in_plane_direction[1], _in_plane_direction[0]) * _JxW[_qp] * _coord[_qp];
      ave_Fhat(_in_plane_direction[1], _in_plane_direction[1]) += _Fhat[_qp](_in_plane_direction[1], _in_plane_direction[1]) * _JxW[_qp] * _coord[_qp];

      // Average deformation gradient
      ave_dfgrd_det += (_deformation_gradient[_qp](_in_plane_direction[0], _in_plane_direction[0]) * _deformation_gradient[_qp](_in_plane_direction[1], _in_plane_direction[1]) -
                        _deformation_gradient[_qp](_in_plane_direction[0], _in_plane_direction[1]) * _deformation_gradient[_qp](_in_plane_direction[1], _in_plane_direction[0])) *
                       _JxW[_qp] * _coord[_qp];
    }
    else if(_legacy_volumetric_locking_correction)
    {
      // Calculate average _Fhat for volumetric locking correction
      ave_Fhat += _Fhat[_qp] * _JxW[_qp] * _coord[_qp];

      // Average deformation gradient
      ave_dfgrd_det += _deformation_gradient[_qp].det() * _JxW[_qp] * _coord[_qp];
    }
    if (_ave_zz_strain && !_legacy_volumetric_locking_correction)
    {
      ave_Fhat(_out_of_plane_direction, _out_of_plane_direction) += _Fhat[_qp](_out_of_plane_direction, _out_of_plane_direction) * _JxW[_qp] * _coord[_qp];
      ave_dfgrd_22 += _deformation_gradient[_qp](_out_of_plane_direction, _out_of_plane_direction) * _JxW[_qp] * _coord[_qp];
    }
  }
  if (_volumetric_locking_correction)
  {
    // needed for volumetric locking correction
    ave_Fhat(_in_plane_direction[0], _in_plane_direction[0]) /= _current_elem_volume;
    ave_Fhat(_in_plane_direction[0], _in_plane_direction[1]) /= _current_elem_volume;
    ave_Fhat(_in_plane_direction[1], _in_plane_direction[0]) /= _current_elem_volume;
    ave_Fhat(_in_plane_direction[1], _in_plane_direction[1]) /= _current_elem_volume;
    // average deformation gradient
    ave_dfgrd_det /= _current_elem_volume;
  }
  else if (_legacy_volumetric_locking_correction)
  {
    // needed for volumetric locking correction
    ave_Fhat /= _current_elem_volume;
    // average deformation gradient
    ave_dfgrd_det /= _current_elem_volume;
  }
  if (_ave_zz_strain && !_legacy_volumetric_locking_correction)
  {
    ave_Fhat(_out_of_plane_direction, _out_of_plane_direction) /= _current_elem_volume;
    ave_dfgrd_22 /= _current_elem_volume;
  }

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    if (_volumetric_locking_correction)
    {
      // Finalize volumetric locking correction
      const Real factor =
          std::sqrt((ave_Fhat(_in_plane_direction[0], _in_plane_direction[0]) * ave_Fhat(_in_plane_direction[1], _in_plane_direction[1]) - ave_Fhat(_in_plane_direction[0], _in_plane_direction[1]) * ave_Fhat(_in_plane_direction[1], _in_plane_direction[0])) /
                    (_Fhat[_qp](_in_plane_direction[0], _in_plane_direction[0]) * _Fhat[_qp](_in_plane_direction[1], _in_plane_direction[1]) - _Fhat[_qp](_in_plane_direction[0], _in_plane_direction[1]) * _Fhat[_qp](_in_plane_direction[1], _in_plane_direction[0])));
      _Fhat[_qp](_in_plane_direction[0], _in_plane_direction[0]) *= factor;
      _Fhat[_qp](_in_plane_direction[0], _in_plane_direction[1]) *= factor;
      _Fhat[_qp](_in_plane_direction[1], _in_plane_direction[0]) *= factor;
      _Fhat[_qp](_in_plane_direction[1], _in_plane_direction[1]) *= factor;

      // Volumetric locking correction
      const Real factor2 = std::sqrt(
          ave_dfgrd_det / (_deformation_gradient[_qp](_in_plane_direction[0], _in_plane_direction[0]) * _deformation_gradient[_qp](_in_plane_direction[1], _in_plane_direction[1]) -
                           _deformation_gradient[_qp](_in_plane_direction[0], _in_plane_direction[1]) * _deformation_gradient[_qp](_in_plane_direction[1], _in_plane_direction[0])));
      _deformation_gradient[_qp](_in_plane_direction[0], _in_plane_direction[0]) *= factor2;
      _deformation_gradient[_qp](_in_plane_direction[0], _in_plane_direction[1]) *= factor2;
      _deformation_gradient[_qp](_in_plane_direction[1], _in_plane_direction[0]) *= factor2;
      _deformation_gradient[_qp](_in_plane_direction[1], _in_plane_direction[1]) *= factor2;
    }
    else if (_legacy_volumetric_locking_correction)
    {
      // Finalize volumetric locking correction
      _Fhat[_qp] *= std::cbrt(ave_Fhat.det() / _Fhat[_qp].det());
      // Volumetric locking correction
      _deformation_gradient[_qp] *= std::cbrt(ave_dfgrd_det / _deformation_gradient[_qp].det());
    }
    if (_ave_zz_strain && !_legacy_volumetric_locking_correction)
    {
      _Fhat[_qp](_out_of_plane_direction, _out_of_plane_direction) = ave_Fhat(_out_of_plane_direction, _out_of_plane_direction);
      _deformation_gradient[_qp](_out_of_plane_direction, _out_of_plane_direction) = ave_dfgrd_22;
    }

    computeQpStrain();
  }
}

void
Compute2DFiniteStrain::displacementIntegrityCheck()
{
  if (_out_of_plane_direction != 2 && _ndisp != 3)
    mooseError("For 2D simulations where the out-of-plane direction is x or y the number of "
               "supplied displacements must be three.");
  else if (_out_of_plane_direction == 2 && _ndisp != 2)
    mooseError("For 2D simulations where the out-of-plane direction is z the number of supplied "
               "displacements must be two.");
}
