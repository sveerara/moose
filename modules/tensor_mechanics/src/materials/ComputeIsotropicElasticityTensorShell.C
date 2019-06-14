//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeIsotropicElasticityTensorShell.h"
#include "RankFourTensor.h"

registerMooseObject("TensorMechanicsApp", ComputeIsotropicElasticityTensorShell);

template <>
InputParameters
validParams<ComputeIsotropicElasticityTensorShell>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Compute a constant isotropic elasticity tensor.");
  params.addParam<Real>("poissons_ratio", "Poisson's ratio for the material.");
  params.addParam<Real>("youngs_modulus", "Young's modulus of the material.");
  return params;
}

ComputeIsotropicElasticityTensorShell::ComputeIsotropicElasticityTensorShell(
    const InputParameters & parameters)
  : Material(parameters),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _youngs_modulus(getParam<Real>("youngs_modulus")),
    _elasticity_tensor(declareProperty<std::vector<RankFourTensor>>("elasticity_tensor")),
    _ge(getMaterialProperty<std::vector<RankTwoTensor>>("ge_matrix"))
{
  _Cijkl.fillSymmetricIsotropicEandNu(_youngs_modulus, _poissons_ratio);

  // correction for plane stress
  _Cijkl(0,0,0,0) = _youngs_modulus / (1.0-_poissons_ratio * _poissons_ratio);
  _Cijkl(1,1,1,1) = _Cijkl(0,0,0,0);
  _Cijkl(0,0,1,1) = _Cijkl(0,0,0,0) * _poissons_ratio;
  _Cijkl(1,1,0,0) = _Cijkl(0,0,1,1);
  _Cijkl(0,0,2,2) = 0.0;
  _Cijkl(1,1,2,2) = 0.0;
  _Cijkl(2,2,2,2) = 0.0;
  _Cijkl(2,2,0,0) = 0.0;
  _Cijkl(2,2,1,1) = 0.0;
  printf("elasticity tensor \n");
  _Cijkl.print();
  printf("elasticity done \n");
}

void
ComputeIsotropicElasticityTensorShell::computeQpProperties()
{
  _elasticity_tensor[_qp].resize(_ge[_qp].size());
  // Assign elasticity tensor at a given quad point
  RankFourTensor a;
  for (unsigned i = 0; i < _ge[_qp].size(); ++i)
    _elasticity_tensor[_qp][i] = a;

  // compute contravariant elasticity tensor
  for (unsigned int t = 0; t < _ge[_qp].size(); ++t)
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        for (unsigned int k = 0; k < 3; ++k)
          for (unsigned int l = 0; l < 3; ++l)
            for (unsigned int m = 0; m < 3; ++m)
              for (unsigned int n = 0; n < 3; ++n)
                for (unsigned int o = 0; o < 3; ++o)
                  for (unsigned int p = 0; p < 3; ++p)
                  {
                //    if (i == m && j == n && k == o && l == p)
                //       printf(" all terms, i, j, k ,l, m, n, o, p: %u, %u, %u, %u, %u, %u, %u, %u,%e, %e, %e, %e, %e \n", i, m, j, n, k, o, l, p, _ge[_qp][t](i, m), _ge[_qp][t](j, n), _ge[_qp][t](k, o), _ge[_qp][t](l, p), _Cijkl(m,n,o,p));
                    _elasticity_tensor[_qp][t](i,j,k,l) += _ge[_qp][t](i, m) * _ge[_qp][t](j, n) * _ge[_qp][t](k, o) * _ge[_qp][t](l, p) * _Cijkl(m,n,o,p);

                  //  _elasticity_tensor[_qp][t](i,j,k,l) = _Cijkl(i,j,k,l);

                  }

/*   for (unsigned i = 0; i < _ge[_qp].size(); ++i)
   {
     printf("ge, _qp, t :%u, %u \n", _qp, i);
     _ge[_qp][i].print();
     printf("edited elasticity, _qp, t: %u, %u \n", _qp, i);
     _elasticity_tensor[_qp][i].print();
   }*/
}
