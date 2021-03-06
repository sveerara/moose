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

#include "SetupPreconditionerAction.h"
#include "Factory.h"
#include "PetscSupport.h"
#include "MoosePreconditioner.h"
#include "FEProblem.h"
#include "CreateExecutionerAction.h"

unsigned int SetupPreconditionerAction::_count = 0;

template<>
InputParameters validParams<SetupPreconditionerAction>()
{
  InputParameters params = validParams<MooseObjectAction>();
  params += Moose::PetscSupport::getPetscValidParams();
  return params;
}

SetupPreconditionerAction::SetupPreconditionerAction(InputParameters params) :
    MooseObjectAction(params)
{
}

void
SetupPreconditionerAction::act()
{
  if (_problem.get() != NULL)
  {
    // build the preconditioner
    _moose_object_pars.set<FEProblem *>("_fe_problem") = _problem.get();
    MooseSharedPointer<MoosePreconditioner> pc = MooseSharedNamespace::static_pointer_cast<MoosePreconditioner>(_factory.create(_type, _name, _moose_object_pars));
    if (!pc.get())
      mooseError("Failed to build the preconditioner.");

    _problem->getNonlinearSystem().setPreconditioner(pc);

    /**
     * Go ahead and set common precondition options here.  The child classes will still be called
     * through the action warehouse
     */
    Moose::PetscSupport::storePetscOptions(*_problem, _pars);
  }
}

// DEPRECATED CONSTRUCTOR
SetupPreconditionerAction::SetupPreconditionerAction(const std::string & deprecated_name, InputParameters params) :
    MooseObjectAction(deprecated_name, params)
{
}
