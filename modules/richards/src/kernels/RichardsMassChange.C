/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#include "RichardsMassChange.h"
#include "Material.h"

#include <iostream>


template<>
InputParameters validParams<RichardsMassChange>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addParam<bool>("use_supg", false, "True for using SUPG in this kernel, false otherwise.  This has no effect if the material does not use SUPG.");
  params.addRequiredParam<UserObjectName>("richardsVarNames_UO", "The UserObject that holds the list of Richards variable names.");
  return params;
}

RichardsMassChange::RichardsMassChange(const InputParameters & parameters) :
    TimeDerivative(parameters),
    _richards_name_UO(getUserObject<RichardsVarNames>("richardsVarNames_UO")),
    _pvar(_richards_name_UO.richards_var_num(_var.number())),

    _use_supg(getParam<bool>("use_supg")),

    _mass(getMaterialProperty<std::vector<Real> >("mass")),
    _dmass(getMaterialProperty<std::vector<std::vector<Real> > >("dmass")),
    _mass_old(getMaterialProperty<std::vector<Real> >("mass_old")),

    _tauvel_SUPG(getMaterialProperty<std::vector<RealVectorValue> >("tauvel_SUPG")),
    _dtauvel_SUPG_dgradv(getMaterialProperty<std::vector<std::vector<RealTensorValue> > >("dtauvel_SUPG_dgradv")),
    _dtauvel_SUPG_dv(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("dtauvel_SUPG_dv"))
{
}


Real
RichardsMassChange::computeQpResidual()
{
  Real test_fcn = _test[_i][_qp] ;
  if (_use_supg)
    test_fcn += _tauvel_SUPG[_qp][_pvar]*_grad_test[_i][_qp];
  return test_fcn*(_mass[_qp][_pvar] - _mass_old[_qp][_pvar])/_dt;
}

Real
RichardsMassChange::computeQpJac(unsigned int wrt_num)
{
  Real mass = _mass[_qp][_pvar];
  Real mass_old = _mass_old[_qp][_pvar];
  Real mass_prime = _phi[_j][_qp]*_dmass[_qp][_pvar][wrt_num];

  Real test_fcn = _test[_i][_qp] ;
  Real test_fcn_prime = 0;

  if (_use_supg)
  {
    test_fcn += _tauvel_SUPG[_qp][_pvar]*_grad_test[_i][_qp];
    test_fcn_prime += _grad_phi[_j][_qp]*(_dtauvel_SUPG_dgradv[_qp][_pvar][wrt_num]*_grad_test[_i][_qp]) + _phi[_j][_qp]*_dtauvel_SUPG_dv[_qp][_pvar][wrt_num]*_grad_test[_i][_qp];
  }
  return (test_fcn*mass_prime + test_fcn_prime*(mass- mass_old))/_dt;
}


Real
RichardsMassChange::computeQpJacobian()
{
  return computeQpJac(_pvar);
}

Real
RichardsMassChange::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_richards_name_UO.not_richards_var(jvar))
    return 0.0;
  unsigned int dvar = _richards_name_UO.richards_var_num(jvar);
  return computeQpJac(dvar);
}


// DEPRECATED CONSTRUCTOR
RichardsMassChange::RichardsMassChange(const std::string & deprecated_name, InputParameters parameters) :
    TimeDerivative(deprecated_name, parameters),
    _richards_name_UO(getUserObject<RichardsVarNames>("richardsVarNames_UO")),
    _pvar(_richards_name_UO.richards_var_num(_var.number())),

    _use_supg(getParam<bool>("use_supg")),

    _mass(getMaterialProperty<std::vector<Real> >("mass")),
    _dmass(getMaterialProperty<std::vector<std::vector<Real> > >("dmass")),
    _mass_old(getMaterialProperty<std::vector<Real> >("mass_old")),

    _tauvel_SUPG(getMaterialProperty<std::vector<RealVectorValue> >("tauvel_SUPG")),
    _dtauvel_SUPG_dgradv(getMaterialProperty<std::vector<std::vector<RealTensorValue> > >("dtauvel_SUPG_dgradv")),
    _dtauvel_SUPG_dv(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("dtauvel_SUPG_dv"))
{
}
