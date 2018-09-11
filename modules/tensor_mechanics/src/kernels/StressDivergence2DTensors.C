/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergence2DTensors.h"
#include "Assembly.h"
#include "ElasticityTensorTools.h"
#include "libmesh/quadrature.h"

registerMooseObject("TensorMechanicsApp", StressDivergence2DTensors);

template <>
InputParameters
validParams<StressDivergence2DTensors>()
{
  InputParameters params = validParams<StressDivergenceTensors>();
  params.addClassDescription("Calculate stress divergence for a 2D problem.");
  params.addRequiredParam<unsigned int>(
      "component",
      "An integer corresponding to the direction the variable this kernel acts in. (0 "
      "for x, 1 for y, 2 for z; note in this kernel disp_x refers to the radial "
      "displacement and disp_y refers to the axial displacement.)");
  params.addCoupledVar("out_of_plane_strain",
                       "The name of the out_of_plane_strain variable used in the "
                       "WeakPlaneStress kernel. Required only if want to provide off-diagonal "
                       "Jacobian in plane stress analysis using weak formulation.");
  MooseEnum out_of_plane_direction("x y z", "z");
  params.addParam<MooseEnum>(
      "out_of_plane_direction",
      out_of_plane_direction,
      "The direction of the out_of_plane_strain variable used in the WeakPlaneStress kernel.");
  params.addParam<bool>("legacy_volumetric_locking_correction", false, "Older version of 2D volumetric locking correction to compare results against the solid mechanics version");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

StressDivergence2DTensors::StressDivergence2DTensors(const InputParameters & parameters)
  : StressDivergenceTensors(parameters),
    _out_of_plane_strain_coupled(isCoupled("out_of_plane_strain")),
    _out_of_plane_strain_var(_out_of_plane_strain_coupled ? coupled("out_of_plane_strain") : 0),
    _out_of_plane_direction(getParam<MooseEnum>("out_of_plane_direction")),
    _legacy_volumetric_locking_correction(getParam<bool>("legacy_volumetric_locking_correction")),
    _in_plane_direction(2),
    _avg_grad_zz_test(_test.size(), 0.0),
    _avg_grad_zz_phi(_phi.size(), 0.0)
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_out_of_plane_direction != 2 && _ndisp != 3)
    mooseError("For 2D simulations where the out-of-plane direction is x or y coordinate "
               "directions the number of supplied displacements must be three.");
  else if (_out_of_plane_direction == 2 && _ndisp != _mesh.dimension())
    mooseError("The number of displacement variables supplied must match the mesh dimension");

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
StressDivergence2DTensors::computeResidual()
{
  computeAverageGradientZZTest();
  if (_legacy_volumetric_locking_correction)
    computeAverageGradientTest();

  StressDivergenceTensors::computeResidual();
}

Real
StressDivergence2DTensors::computeQpResidual()
{
  Real div = 0.0;
  if (_component == _in_plane_direction[0])
  {
    div = _grad_test[_i][_qp](_in_plane_direction[0]) * _stress[_qp](_in_plane_direction[0], _in_plane_direction[0]) + getGradientZZTest() * _stress[_qp](_out_of_plane_direction, _out_of_plane_direction) +
          +_grad_test[_i][_qp](_in_plane_direction[1]) * _stress[_qp](_in_plane_direction[0], _in_plane_direction[1]); // stress_{rz}

    // volumetric locking correction
    if (_volumetric_locking_correction)
      div += (_avg_grad_test[_i][_in_plane_direction[0]] - _grad_test[_i][_qp](_in_plane_direction[0])) *
             (_stress[_qp](_in_plane_direction[0], _in_plane_direction[0]) + _stress[_qp](_in_plane_direction[1], _in_plane_direction[1])) / 2.0;
    else if (_legacy_volumetric_locking_correction)
      div += (_avg_grad_test[_i][_in_plane_direction[0]] + _avg_grad_zz_test[_i] - _grad_test[_i][_qp](_in_plane_direction[0]) - getGradientZZTest()) * (_stress[_qp].trace()) / 3.0;
  }
  else if (_component == _in_plane_direction[1])
  {
    div = _grad_test[_i][_qp](_in_plane_direction[1]) * _stress[_qp](_in_plane_direction[1], _in_plane_direction[1]) +
          +_grad_test[_i][_qp](_in_plane_direction[0]) * _stress[_qp](_in_plane_direction[1], _in_plane_direction[0]); // stress_{zr}

    // volumetric locking correction
    if (_volumetric_locking_correction)
      div += (_avg_grad_test[_i][_in_plane_direction[1]] - _grad_test[_i][_qp](_in_plane_direction[1])) *
             (_stress[_qp](_in_plane_direction[0], _in_plane_direction[0]) + _stress[_qp](_in_plane_direction[1], _in_plane_direction[1])) / 2.0;
    else if (_legacy_volumetric_locking_correction)
      div += (_avg_grad_test[_i][_in_plane_direction[1]] - _grad_test[_i][_qp](_in_plane_direction[1])) * (_stress[_qp].trace()) / 3.0;;
  }
  else
    mooseError("Invalid component for this 2D problem.");
  return div;
}

void
StressDivergence2DTensors::computeJacobian()
{
  computeAverageGradientZZTest();
  computeAverageGradientZZPhi();

  if (_legacy_volumetric_locking_correction)
  {
    computeAverageGradientTest();
    computeAverageGradientPhi();
  }
  StressDivergenceTensors::computeJacobian();
}

Real
StressDivergence2DTensors::computeQpJacobian()
{
  return calculateJacobian(_component, _component);
}

void
StressDivergence2DTensors::computeOffDiagJacobian(MooseVariableFEBase & jvar)
{
  computeAverageGradientZZTest();
  computeAverageGradientZZPhi();

  if (_legacy_volumetric_locking_correction)
  {
    computeAverageGradientTest();
    computeAverageGradientPhi();
  }

  StressDivergenceTensors::computeOffDiagJacobian(jvar);
}

Real
StressDivergence2DTensors::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    if (jvar == _disp_var[i])
    {
      if (_out_of_plane_direction != 2)
      {
        if (i == _out_of_plane_direction)
          continue;
      }
      return calculateJacobian(_component, i);
    }
  }

  // off-diagonal Jacobian with respect to a coupled out_of_plane_strain variable
  if (_out_of_plane_strain_coupled && jvar == _out_of_plane_strain_var)
    return _Jacobian_mult[_qp](
               _component, _component, _out_of_plane_direction, _out_of_plane_direction) *
           _grad_test[_i][_qp](_component) * _phi[_j][_qp];

  if (_temp_coupled && jvar == _temp_var)
  {
    Real jac = 0.0;
    if (_component == _in_plane_direction[0])
    {
      for (unsigned k = 0; k < LIBMESH_DIM; ++k)
        for (unsigned l = 0; l < LIBMESH_DIM; ++l)
          jac -= (_grad_test[_i][_qp](_in_plane_direction[0]) * _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[0], k, l) +
                  getGradientZZTest() * _Jacobian_mult[_qp](_out_of_plane_direction, _out_of_plane_direction, k, l) +
                  _grad_test[_i][_qp](_in_plane_direction[1]) * _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[1], k, l)) *
                 (*_deigenstrain_dT)[_qp](k, l);
      return jac * _phi[_j][_qp];
    }
    else if (_component == _in_plane_direction[1])
    {
      for (unsigned k = 0; k < LIBMESH_DIM; ++k)
        for (unsigned l = 0; l < LIBMESH_DIM; ++l)
          jac -= (_grad_test[_i][_qp](_in_plane_direction[1]) * _Jacobian_mult[_qp](_in_plane_direction[1], _in_plane_direction[1], k, l) +
                  _grad_test[_i][_qp](_in_plane_direction[0]) * _Jacobian_mult[_qp](_in_plane_direction[1], _in_plane_direction[0], k, l)) *
                 (*_deigenstrain_dT)[_qp](k, l);
      return jac * _phi[_j][_qp];
    }
  }
  return 0.0;
}

Real
StressDivergence2DTensors::calculateJacobian(unsigned int ivar, unsigned int jvar)
{
  // B^T_i * C * B_j
  RealGradient test, test_z, phi, phi_z;
  Real first_term = 0.0;
  if (ivar == _in_plane_direction[0]) // Case grad_test for x, requires contributions from stress_xx, stress_xy, and stress_zz
  {
    test(_in_plane_direction[0]) = _grad_test[_i][_qp](_in_plane_direction[0]);
    test(_in_plane_direction[1]) = _grad_test[_i][_qp](_in_plane_direction[1]);
    test_z(_out_of_plane_direction) = getGradientZZTest();
  }
  else // Case grad_test for y
  {
    test(_in_plane_direction[0]) = _grad_test[_i][_qp](_in_plane_direction[0]);
    test(_in_plane_direction[1]) = _grad_test[_i][_qp](_in_plane_direction[1]);
  }

  if (jvar == _in_plane_direction[0])
  {
    phi(_in_plane_direction[0]) = _grad_phi[_j][_qp](_in_plane_direction[0]);
    phi(_in_plane_direction[1]) = _grad_phi[_j][_qp](_in_plane_direction[1]);
    phi_z(_out_of_plane_direction) = getGradientZZPhi();
  }
  else
  {
    phi(_in_plane_direction[0]) = _grad_phi[_j][_qp](_in_plane_direction[0]);
    phi(_in_plane_direction[1]) = _grad_phi[_j][_qp](_in_plane_direction[1]);
  }

  if (ivar == _in_plane_direction[0] &&
      jvar == _in_plane_direction[0]) // Case when both phi and test are functions of x and z; requires four terms
  {
    const Real first_sum = ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], ivar, jvar, test, phi); // test_x and phi_x
    const Real second_sum = ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], _out_of_plane_direction, _out_of_plane_direction, test_z, phi_z); // test_z and phi_z
    const Real mixed_sum1 = ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], ivar, _out_of_plane_direction, test, phi_z); // test_x and phi_z
    const Real mixed_sum2 = ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], _out_of_plane_direction, jvar, test_z, phi); // test_z and phi_x

    first_term = first_sum + second_sum + mixed_sum1 + mixed_sum2;
  }
  else if (ivar == _in_plane_direction[0] && jvar == _in_plane_direction[1])
  {
    const Real first_sum = ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], ivar, jvar, test, phi); // test_x and phi_y
    const Real mixed_sum2 = ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], _out_of_plane_direction, jvar, test_z, phi); // test_z and phi_y

    first_term = first_sum + mixed_sum2;
  }
  else if (ivar == _in_plane_direction[1] && jvar == _in_plane_direction[0])
  {
    const Real second_sum = ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], ivar, jvar, test, phi); // test_y and phi_x
    const Real mixed_sum1 = ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], ivar, _out_of_plane_direction, test, phi_z); // test_y and phi_z

    first_term = second_sum + mixed_sum1;
  }
  else if (ivar == _in_plane_direction[1] && jvar == _in_plane_direction[1])
    first_term = ElasticityTensorTools::elasticJacobian(
        _Jacobian_mult[_qp], ivar, jvar, test, phi); // test_y and phi_y
  else
    mooseError("Invalid component in Jacobian Calculation");

  Real val = 0.0;
  // volumetric locking correction
  // K = Bbar^T_i * C * Bbar^T_j where Bbar = B + Bvol
  // K = B^T_i * C * B_j + Bvol^T_i * C * Bvol_j + B^T_i * C * Bvol_j + Bvol^T_i * C * B_j
  if (_volumetric_locking_correction)
  {
    RealGradient new_test(3, 0.0);
    RealGradient new_phi(3, 0.0);

    new_test(_in_plane_direction[0]) = _grad_test[_i][_qp](_in_plane_direction[0]);
    new_test(_in_plane_direction[1]) = _grad_test[_i][_qp](_in_plane_direction[1]);
    new_test(_out_of_plane_direction) = getGradientZZTest();
    new_phi(_in_plane_direction[0]) = _grad_phi[_j][_qp](_in_plane_direction[0]);
    new_phi(_in_plane_direction[1]) = _grad_phi[_j][_qp](_in_plane_direction[1]);
    new_phi(_out_of_plane_direction) = getGradientZZPhi();

    // Bvol^T_i * C * Bvol_j
    Real sum = 0.0;
    for (unsigned i = 0; i < 2; ++i)
      for (unsigned j = 0; j < 2; ++j)
        sum += _Jacobian_mult[_qp](i, i, j, j);

    val += sum * (_avg_grad_test[_i][ivar] - new_test(ivar)) *
           (_avg_grad_phi[_j][jvar] - new_phi(jvar)) / 2.0;

    // B^T_i * C * Bvol_j
    RealGradient sum_2x1(3, 0.0);
    sum_2x1(_in_plane_direction[0]) = _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[0], _in_plane_direction[0], _in_plane_direction[0]) + _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[0], _in_plane_direction[1], _in_plane_direction[1]);
    sum_2x1(_in_plane_direction[1]) = _Jacobian_mult[_qp](_in_plane_direction[1], _in_plane_direction[1], _in_plane_direction[0], _in_plane_direction[0]) + _Jacobian_mult[_qp](_in_plane_direction[1], _in_plane_direction[1], _in_plane_direction[1], _in_plane_direction[1]);
    sum_2x1(_out_of_plane_direction) = _Jacobian_mult[_qp](_out_of_plane_direction, _out_of_plane_direction, _in_plane_direction[0], _in_plane_direction[0]) + _Jacobian_mult[_qp](_out_of_plane_direction, _out_of_plane_direction, _in_plane_direction[1], _in_plane_direction[1]);
    Real sum_2x1_3 = _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[1], _in_plane_direction[0], _in_plane_direction[0]) + _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[1], _in_plane_direction[1], _in_plane_direction[1]);

    if (ivar == _in_plane_direction[0] && jvar == _in_plane_direction[0])
      val += (sum_2x1(_in_plane_direction[0]) * new_test(_in_plane_direction[0]) + sum_2x1(_out_of_plane_direction) * new_test(_out_of_plane_direction) + sum_2x1_3 * new_test(_in_plane_direction[1]) / 2.0) *
             (_avg_grad_phi[_j][_in_plane_direction[0]] - new_phi(_in_plane_direction[0]));
    else if (ivar == _in_plane_direction[0] && jvar == _in_plane_direction[1])
        val += (sum_2x1(_in_plane_direction[0]) * new_test(_in_plane_direction[0]) + sum_2x1(_out_of_plane_direction) * new_test(_out_of_plane_direction) + sum_2x1_3 * new_test(_in_plane_direction[1]) / 2.0 ) * (_avg_grad_phi[_j][_in_plane_direction[1]] - new_phi(_in_plane_direction[1]));
    else if (ivar == _in_plane_direction[1] && jvar == _in_plane_direction[0])
      val += (sum_2x1(_in_plane_direction[1]) * new_test(_in_plane_direction[1]) + sum_2x1_3 * new_test(_in_plane_direction[0]) / 2.0) * (_avg_grad_phi[_j][_in_plane_direction[0]] - new_phi(_in_plane_direction[0]));
    else if (ivar == _in_plane_direction[1] && jvar == _in_plane_direction[1])
      val += (sum_2x1(_in_plane_direction[1]) * new_test(_in_plane_direction[1]) + sum_2x1_3 * new_test(_in_plane_direction[0]) / 2.0) * (_avg_grad_phi[_j][_in_plane_direction[1]] - new_phi(_in_plane_direction[1]));

    // Bvol^T_i * C * B_j
    sum_2x1(_in_plane_direction[0]) = _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[0], _in_plane_direction[0], _in_plane_direction[0]) + _Jacobian_mult[_qp](_in_plane_direction[1], _in_plane_direction[1], _in_plane_direction[0], _in_plane_direction[0]);
    sum_2x1(_in_plane_direction[1]) = _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[0], _in_plane_direction[1], _in_plane_direction[1]) + _Jacobian_mult[_qp](_in_plane_direction[1], _in_plane_direction[1], _in_plane_direction[1], _in_plane_direction[1]);
    sum_2x1(_out_of_plane_direction) = _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[0], _out_of_plane_direction, _out_of_plane_direction) + _Jacobian_mult[_qp](_in_plane_direction[1], _in_plane_direction[1], _out_of_plane_direction, _out_of_plane_direction);
    sum_2x1_3 = _Jacobian_mult[_qp](_in_plane_direction[0], _in_plane_direction[0], _in_plane_direction[0], _in_plane_direction[1]) + _Jacobian_mult[_qp](_in_plane_direction[1], _in_plane_direction[1], _in_plane_direction[0], _in_plane_direction[1]);
    if (ivar == _in_plane_direction[0] && jvar == _in_plane_direction[0])
       val += (sum_2x1(_in_plane_direction[0]) * new_phi(_in_plane_direction[0]) + sum_2x1(_out_of_plane_direction) * new_phi(_out_of_plane_direction) + sum_2x1_3 * new_phi(_in_plane_direction[1]) / 2.0) * (_avg_grad_test[_i][_in_plane_direction[0]] - new_test(_in_plane_direction[0]));
    else if (ivar == _in_plane_direction[0] && jvar == _in_plane_direction[1])
       val += (sum_2x1(_in_plane_direction[1]) * new_phi(_in_plane_direction[1]) + sum_2x1_3 * new_phi(_in_plane_direction[0]) / 2.0) * (_avg_grad_test[_i][_in_plane_direction[0]] - new_test(_in_plane_direction[0]));
    else if (ivar == _in_plane_direction[1] && jvar == _in_plane_direction[0])
       val += (sum_2x1(_in_plane_direction[0]) * new_phi(_in_plane_direction[0]) + sum_2x1(_out_of_plane_direction) * new_phi(_out_of_plane_direction) + sum_2x1_3 * new_phi(_in_plane_direction[1]) / 2.0) * (_avg_grad_test[_i][_in_plane_direction[1]] - new_test(_in_plane_direction[1]));
    else if (ivar == _in_plane_direction[1] && jvar == _in_plane_direction[1])
       val += (sum_2x1(_in_plane_direction[1]) * new_phi(_in_plane_direction[1]) + sum_2x1_3 * new_phi(_in_plane_direction[0]) / 2.0) * (_avg_grad_test[_i][_in_plane_direction[1]] - new_test(_in_plane_direction[1]));

    val /= 2.0;
  }
  else if (_legacy_volumetric_locking_correction)
  {
    RealGradient new_test(3, 0.0);
    RealGradient new_phi(3, 0.0);

    new_test(_in_plane_direction[0]) = _grad_test[_i][_qp](_in_plane_direction[0]) + getGradientZZTest();
    new_test(_in_plane_direction[1]) = _grad_test[_i][_qp](_in_plane_direction[1]);
    new_phi(_in_plane_direction[0]) = _grad_phi[_j][_qp](_in_plane_direction[0]) + getGradientZZPhi();
    new_phi(_in_plane_direction[1]) = _grad_phi[_j][_qp](_in_plane_direction[1]);

    // Bvol^T_i * C * Bvol_j
    if (ivar == _in_plane_direction[0] && jvar == _in_plane_direction[0])
      val += _Jacobian_mult[_qp].sum3x3() * (_avg_grad_test[_i][ivar] + _avg_grad_zz_test[_i] - new_test(ivar)) *
             (_avg_grad_phi[_j][jvar] + _avg_grad_zz_phi[_j] - new_phi(jvar)) / 3.0;
    else if (ivar == _in_plane_direction[0] && jvar == _in_plane_direction[1])
      val += _Jacobian_mult[_qp].sum3x3() * (_avg_grad_test[_i][ivar] + _avg_grad_zz_test[_i] - new_test(ivar)) *
             (_avg_grad_phi[_j][jvar] - new_phi(jvar)) / 3.0;
    else if (ivar == _in_plane_direction[1] && jvar == _in_plane_direction[0])
      val += _Jacobian_mult[_qp].sum3x3() * (_avg_grad_test[_i][ivar] - new_test(ivar)) *
            (_avg_grad_phi[_j][jvar] + _avg_grad_zz_phi[_j] - new_phi(jvar)) / 3.0;
    else
      val += _Jacobian_mult[_qp].sum3x3() * (_avg_grad_test[_i][ivar] - new_test(ivar)) *
            (_avg_grad_phi[_j][jvar] - new_phi(jvar)) / 3.0;

    // B^T_i * C * Bvol_j
    RealGradient sum_3x1 = _Jacobian_mult[_qp].sum3x1();
    if (ivar == _in_plane_direction[0] && jvar == _in_plane_direction[0])
      val += (sum_3x1(_in_plane_direction[0]) * test(_in_plane_direction[0]) + sum_3x1(_out_of_plane_direction) * test_z(_out_of_plane_direction)) * (_avg_grad_phi[_j][_in_plane_direction[0]] +_avg_grad_zz_phi[_j] - new_phi(_in_plane_direction[0]));
    else if (ivar == _in_plane_direction[0] && jvar == _in_plane_direction[1])
      val += (sum_3x1(_in_plane_direction[0]) * test(_in_plane_direction[0]) + sum_3x1(_out_of_plane_direction) * test_z(_out_of_plane_direction)) * (_avg_grad_phi[_j][_in_plane_direction[1]] - new_phi(_in_plane_direction[1]));
    else if (ivar == _in_plane_direction[1] && jvar == _in_plane_direction[0])
      val += sum_3x1(_in_plane_direction[1]) * test(_in_plane_direction[1]) * (_avg_grad_phi[_j][_in_plane_direction[0]] + _avg_grad_zz_phi[_j]- new_phi(_in_plane_direction[0]));
    else
      val += sum_3x1(_in_plane_direction[1]) * test(_in_plane_direction[1]) * (_avg_grad_phi[_j][_in_plane_direction[1]] - new_phi(_in_plane_direction[1]));

    // Bvol^T_i * C * B_j
    // val = trace (C * B_j) *(avg_grad_test[_i][ivar] - new_test(ivar))
    if (jvar == _in_plane_direction[0])
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        if (ivar == _in_plane_direction[0])
          val +=
             (_Jacobian_mult[_qp](i, i, _in_plane_direction[0], _in_plane_direction[0]) * phi(_in_plane_direction[0]) + _Jacobian_mult[_qp](i, i, _in_plane_direction[0], _in_plane_direction[1]) * phi(_in_plane_direction[1]) +
             _Jacobian_mult[_qp](i, i, _out_of_plane_direction, _out_of_plane_direction) * phi_z(_out_of_plane_direction)) *
             (_avg_grad_test[_i][ivar] + _avg_grad_zz_test[_i] - new_test(ivar));
        else
          val +=
            (_Jacobian_mult[_qp](i, i, _in_plane_direction[0], _in_plane_direction[0]) * phi(_in_plane_direction[0]) + _Jacobian_mult[_qp](i, i, _in_plane_direction[0], _in_plane_direction[1]) * phi(_in_plane_direction[1]) +
            _Jacobian_mult[_qp](i, i, _out_of_plane_direction, _out_of_plane_direction) * phi_z(_out_of_plane_direction)) *
            (_avg_grad_test[_i][ivar] - new_test(ivar));
      }
    }
    else if (jvar == 1)
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        if (ivar == _in_plane_direction[0])
          val +=
              (_Jacobian_mult[_qp](i, i, _in_plane_direction[0], _in_plane_direction[1]) * phi(_in_plane_direction[0]) + _Jacobian_mult[_qp](i, i, _in_plane_direction[1], _in_plane_direction[1]) * phi(_in_plane_direction[1])) *
              (_avg_grad_test[_i][ivar] + _avg_grad_zz_test[_i] - new_test(ivar));
        else
          val +=
             (_Jacobian_mult[_qp](i, i, _in_plane_direction[0], _in_plane_direction[1]) * phi(_in_plane_direction[0]) + _Jacobian_mult[_qp](i, i, _in_plane_direction[1], _in_plane_direction[1]) * phi(_in_plane_direction[1])) *
             (_avg_grad_test[_i][ivar] - new_test(ivar));
      }
    }

    val /= 3.0;
  }
  return val + first_term;
}

void
StressDivergence2DTensors::computeAverageGradientTest()
{
  // calculate volume averaged value of shape function derivative
  _avg_grad_test.resize(_test.size());
  for (_i = 0; _i < _test.size(); ++_i)
  {
    _avg_grad_test[_i].resize(2);
    _avg_grad_test[_i][_component] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _avg_grad_test[_i][_component] += _grad_test[_i][_qp](_component) * _JxW[_qp] * _coord[_qp];

    _avg_grad_test[_i][_component] /= _current_elem_volume;
  }
}

void
StressDivergence2DTensors::computeAverageGradientPhi()
{
  _avg_grad_phi.resize(_phi.size());
  for (_i = 0; _i < _phi.size(); ++_i)
  {
    _avg_grad_phi[_i].resize(3);
    for (unsigned int ii = 0; ii < 2; ++ii)
    {
      _avg_grad_phi[_i][_in_plane_direction[ii]] = 0.0;
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
        _avg_grad_phi[_i][_in_plane_direction[ii]] += _grad_phi[_i][_qp](_in_plane_direction[ii]) * _JxW[_qp] * _coord[_qp];

      _avg_grad_phi[_i][_in_plane_direction[ii]] /= _current_elem_volume;
    }
  }
}

void
StressDivergence2DTensors::computeAverageGradientZZTest()
{
  _avg_grad_zz_test.resize(_test.size());
  for (_i = 0; _i < _test.size(); ++_i)
    _avg_grad_zz_test[_i] = 0.0;
}

void
StressDivergence2DTensors::computeAverageGradientZZPhi()
{
  _avg_grad_zz_phi.resize(_phi.size());
  for (_i = 0; _i < _phi.size(); ++_i)
    _avg_grad_zz_phi[_i] = 0.0;
}

Real
StressDivergence2DTensors::getGradientZZTest()
{
  return 0.0;
}

Real
StressDivergence2DTensors::getGradientZZPhi()
{
  return 0.0;
}
