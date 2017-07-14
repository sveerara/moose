//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Compute2DIncrementalStrain.h"

#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<Compute2DIncrementalStrain>()
{
  InputParameters params = validParams<ComputeIncrementalSmallStrain>();
  params.addClassDescription("Compute strain increment for incremental strains in 2D geometries.");

  MooseEnum outOfPlaneDirection("x y z", "z");
  params.addParam<MooseEnum>(
      "out_of_plane_direction", outOfPlaneDirection, "The direction of the out-of-plane strain.");
  return params;
}

Compute2DIncrementalStrain::Compute2DIncrementalStrain(const InputParameters & parameters)
  : ComputeIncrementalSmallStrain(parameters),
    _out_of_plane_direction(getParam<MooseEnum>("out_of_plane_direction")),
    _ave_zz_strain(false)
{
}

void
Compute2DIncrementalStrain::initialSetup()
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
Compute2DIncrementalStrain::computeTotalStrainIncrement(RankTwoTensor & total_strain_increment)
{
  // Deformation gradient calculation for 2D problems
  RankTwoTensor A(
      (*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]); // Deformation gradient
  RankTwoTensor Fbar((*_grad_disp_old[0])[_qp],
                     (*_grad_disp_old[1])[_qp],
                     (*_grad_disp_old[2])[_qp]); // Old Deformation gradient

  // Compute the displacement gradient of the out of plane direction for plane strain,
  // generalized plane strain, or axisymmetric problems
  A(_out_of_plane_direction, _out_of_plane_direction) = computeOutOfPlaneGradDisp();
  Fbar(_out_of_plane_direction, _out_of_plane_direction) = computeOutOfPlaneGradDispOld();

  _deformation_gradient[_qp] = A;
  _deformation_gradient[_qp].addIa(1.0);

  A -= Fbar; // very nearly A = gradU - gradUold

  total_strain_increment = 0.5 * (A + A.transpose());
}

void
Compute2DIncrementalStrain::displacementIntegrityCheck()
{
  if (_out_of_plane_direction != 2 && _ndisp != 3)
    mooseError("For 2D simulations where the out-of-plane direction is x or y the number of "
               "supplied displacements must be three.");
  else if (_out_of_plane_direction == 2 && _ndisp != 2)
    mooseError("For 2D simulations where the out-of-plane direction is z the number of supplied "
               "displacements must be two.");
}

void
Compute2DIncrementalStrain::computeProperties()
{
  Real volumetric_strain = 0.0;
  Real out_of_plane_strain = 0.0;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    RankTwoTensor total_strain_increment;
    computeTotalStrainIncrement(total_strain_increment);

    _strain_increment[_qp] = total_strain_increment;

    if (_volumetric_locking_correction)
      volumetric_strain +=
          (total_strain_increment(0, 0) + total_strain_increment(1, 1)) * _JxW[_qp] * _coord[_qp];

    if (_ave_zz_strain)
      out_of_plane_strain += total_strain_increment(2, 2) * _JxW[_qp] * _coord[_qp];
  }
  if (_volumetric_locking_correction)
    volumetric_strain /= _current_elem_volume;

  if (_ave_zz_strain)
    out_of_plane_strain /= _current_elem_volume;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    if (_volumetric_locking_correction)
    {
      const Real trace_2D = _strain_increment[_qp](0, 0) + _strain_increment[_qp](1, 1);
      _strain_increment[_qp](0, 0) += (volumetric_strain - trace_2D) / 2.0;
      _strain_increment[_qp](1, 1) += (volumetric_strain - trace_2D) / 2.0;
    }

    if (_ave_zz_strain)
      _strain_increment[_qp](2, 2) = out_of_plane_strain;

    _total_strain[_qp] = _total_strain_old[_qp] + _strain_increment[_qp];

    // Remove the Eigen strain increment
    subtractEigenstrainIncrementFromStrain(_strain_increment[_qp]);

    // strain rate
    if (_dt > 0)
      _strain_rate[_qp] = _strain_increment[_qp] / _dt;
    else
      _strain_rate[_qp].zero();

    // Update strain in intermediate configuration: rotations are not needed
    _mechanical_strain[_qp] = _mechanical_strain_old[_qp] + _strain_increment[_qp];
  }
}
