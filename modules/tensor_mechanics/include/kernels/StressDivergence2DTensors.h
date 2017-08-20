/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STRESSDIVERGENCE2DTENSORS_H
#define STRESSDIVERGENCE2DTENSORS_H

#include "StressDivergenceTensors.h"

// Forward Declarations
class StressDivergence2DTensors;

template <>
InputParameters validParams<StressDivergence2DTensors>();

/**
 * StressDivergence2DTensors is a modification of StressDivergenceTensors to
 * accommodate 2D models.  The out-of-plane (zz) response is zero.
 */
class StressDivergence2DTensors : public StressDivergenceTensors
{
public:
  StressDivergence2DTensors(const InputParameters & parameters);
  
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(MooseVariableFEBase & jvar) override;

protected:
  virtual void computeResidual() override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  virtual void computeAverageGradientTest() override;
  virtual void computeAverageGradientPhi() override;
  virtual void computeAverageGradientZZTest();
  virtual void computeAverageGradientZZPhi();
  virtual Real getGradientZZTest();
  virtual Real getGradientZZPhi();

  Real calculateJacobian(unsigned int ivar, unsigned int jvar);

  std::vector<Real> _avg_grad_zz_test;
  std::vector<Real> _avg_grad_zz_phi;
};

#endif // STRESSDIVERGENCE2DTENSORS_H
