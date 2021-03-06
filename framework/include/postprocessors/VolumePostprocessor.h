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

#ifndef VOLUMEPOSTPROCESSOR_H
#define VOLUMEPOSTPROCESSOR_H

#include "ElementIntegralPostprocessor.h"

//Forward Declarations
class VolumePostprocessor;

template<>
InputParameters validParams<VolumePostprocessor>();

/**
 * This postprocessor computes the volume of a specified block.
 */
class VolumePostprocessor : public ElementIntegralPostprocessor
{
public:
  VolumePostprocessor(const InputParameters & parameters);
  VolumePostprocessor(const std::string & deprecated_name, InputParameters parameters); // DEPRECATED CONSTRUCTOR
  virtual void threadJoin(const UserObject & y);

protected:
  virtual Real computeQpIntegral();
};

#endif
