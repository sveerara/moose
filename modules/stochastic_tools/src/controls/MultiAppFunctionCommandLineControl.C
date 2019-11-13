//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "MultiAppFunctionCommandLineControl.h"
#include "Function.h"
#include "MultiApp.h"

registerMooseObject("MooseApp", MultiAppFunctionCommandLineControl);

template <>
InputParameters
validParams<MultiAppFunctionCommandLineControl>()
{
  InputParameters params = validParams<Control>();
  params.addClassDescription("Control for modifying the command line arguments of MultiApps using values from function.");

  // Set and suppress the 'execute_on' flag, it doesn't work with any other flag
  params.set<ExecFlagEnum>("execute_on") = {EXEC_PRE_MULTIAPP_SETUP};
  params.suppressParameter<ExecFlagEnum>("execute_on");

  params.addRequiredParam<MultiAppName>("multi_app", "The name of the MultiApp to control.");
  params.addRequiredParam<FunctionName>(
      "function",
      "The function object to utilize for altering the command line options of the MultiApp.");
  params.addRequiredParam<std::string>(
      "param_name", "The names of the command line parameters to set via the function data.");

  return params;
}

MultiAppFunctionCommandLineControl::MultiAppFunctionCommandLineControl(const InputParameters & parameters)
  : Control(parameters),
    _multi_app(_fe_problem.getMultiApp(getParam<MultiAppName>("multi_app"))),
    _func(getFunction("function")),
    _param_name(getParam<std::string>("param_name"))
{
}

void
MultiAppFunctionCommandLineControl::initialSetup()
{
  // Do not put anything here, this method is being called after execute because the execute_on
  // is set to PRE_MULTIAPP_SETUP for this class. It won't work any other way.
}

void
MultiAppFunctionCommandLineControl::execute()
{
  std::vector<std::string> cli_args;

  for (unsigned int i = 0; i < _multi_app->numGlobalApps(); ++i)
  {
    std::ostringstream oss;
    oss << _param_name << "=" << Moose::stringify(_func.value(0, _multi_app->position(i)));

    cli_args.push_back(oss.str());
  }

  setControllableValueByName<std::vector<std::string>>(
      "MultiApp", _multi_app->name(), "cli_args", cli_args);
}
