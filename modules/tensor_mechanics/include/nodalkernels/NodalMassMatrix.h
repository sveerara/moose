/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef NodalMassMatrix_H
#define NodalMassMatrix_H

#include "NodalKernel.h"

// Forward Declarations
class NodalMassMatrix;

template <>
InputParameters validParams<NodalMassMatrix>();

/**
 * Calculates the inertial force and mass proportional damping for a nodal mass
 */
class NodalMassMatrix : public NodalKernel
{
public:
  NodalMassMatrix(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  /// Mass associated with the node
  const Real _mass;

  /// Map between boundary nodes and nodal mass
  std::map<dof_id_type, Real> _node_id_to_mass;
};

#endif /* NodalMassMatrix_H */
