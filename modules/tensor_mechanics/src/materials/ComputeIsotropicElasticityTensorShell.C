//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeIsotropicElasticityTensorShell.h"

registerMooseObject("TensorMechanicsApp", ComputeIsotropicElasticityTensorShell);

template <>
InputParameters
validParams<ComputeIsotropicElasticityTensorShell>()
{
  InputParameters params = validParams<ComputeElasticityTensorBase>();
  params.addClassDescription("Compute a constant isotropic elasticity tensor.");
  params.addParam<Real>("poissons_ratio", "Poisson's ratio for the material.");
  params.addParam<Real>("youngs_modulus", "Young's modulus of the material.");
  return params;
}

ComputeIsotropicElasticityTensorShell::ComputeIsotropicElasticityTensorShell(
    const InputParameters & parameters)
  : ComputeElasticityTensorBase(parameters),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _youngs_modulus(getParam<Real>("youngs_modulus"))
{
  // all tensors created by this class are always isotropic
  issueGuarantee(_elasticity_tensor_name, Guarantee::ISOTROPIC);
  if (!isParamValid("elasticity_tensor_prefactor"))
    issueGuarantee(_elasticity_tensor_name, Guarantee::CONSTANT_IN_TIME);

  _Cijkl.fillSymmetricIsotropicEandNu(_youngs_modulus, _poissons_ratio);

  // correction for plane stress
  _Cijkl(0,0,0,0) = _youngs_modulus / (1-_poissons_ratio * _poissons_ratio);
  _Cijkl(1,1,1,1) = _Cijkl(0,0,0,0);
  _Cijkl(0,0,1,1) = _Cijkl(0,0,0,0) * _poissons_ratio;
  _Cijkl(1,1,0,0) = _Cijkl(0,0,1,1);
  _Cijkl(0,0,2,2) = 0.0;
  _Cijkl(1,1,2,2) = 0.0;
  _Cijkl(2,2,2,2) = 0.0;
  _Cijkl(2,2,0,0) = 0.0;
  _Cijkl(2,2,1,1) = 0.0;
}

void
ComputeIsotropicElasticityTensorShell::computeQpElasticityTensor()
{
  // Assign elasticity tensor at a given quad point
  for (unsigned i = 0; i < _t_points.size(); ++i)
    _elasticity_tensor[_qp][i] = _Cijkl;

  // compute contravariant elasticity tensor
  for (unsigned int t = 0; t < _t_points.size(); ++t)
    for (unsinged int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        for (unsigned int k = 0; k < 3; ++k)
          for (unsigned int l = 0; l < 3; ++l)
            for (unsigned int m = 0; m < 3; ++m)
              for (unsigned int n = 0; n < 3; ++n)
                for (unsinged int o = 0; o < 3; ++o)
                  for (unsigned int p = 0; p < 3; ++p)
                    _elasticity_tensor[_qp][t][i][j][k][l] = _ge[_qp][_t](i, m) * _ge[_qp][_t](j, n) * _ge[_qp][_t](k, o) * _ge[_qp][_t](l, p) * _Cijkl(m,n,o,p);
}
