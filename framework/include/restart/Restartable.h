//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef RESTARTABLE_H
#define RESTARTABLE_H

// MOOSE includes
#include "MooseTypes.h"
#include "RestartableData.h"

// Forward declarations
class PostprocessorData;
class SubProblem;
class InputParameters;
class MooseObject;
class MooseApp;

/**
 * A class for creating restricted objects
 * \see BlockRestartable BoundaryRestartable
 */
class Restartable
{
public:
  Restartable(const MooseObject * moose_object, const std::string & system_name);

  Restartable(const MooseObject * moose_object, const std::string & system_name, THREAD_ID tid);

  Restartable(MooseApp & moose_app,
              const std::string & name,
              const std::string & system_name,
              THREAD_ID tid);

  /**
   * Class constructor
   * Populates the SubProblem and MooseMesh pointers
   * @param parameters The InputParameters for the object.
   * @param system_name The name of the MOOSE system.  ie "Kernel", "BCs", etc.  Should roughly
   * correspond to the section in the input file so errors are easy to understand.
   * @param subproblem An optional method for inputting the SubProblem object, this is used by
   * FEProblemBase, otherwise the SubProblem comes from the parameters
   */
  //  Restartable(const InputParameters & parameters,
  //              std::string system_name,
  //              SubProblem * subproblem = nullptr);

  /**
   * Constructor for objects that don't have "parameters"
   *
   * @param name The name of the object
   * @param system_name The name of the MOOSE system.  ie "Kernel", "BCs", etc.  Should roughly
   * correspond to the section in the input file so errors are easy to understand.
   * @param subproblem A reference to the subproblem for this object
   * @param tid Optional thread id (will default to zero)
   */
  //  Restartable(const std::string & name,
  //              std::string system_name,
  //              SubProblem & subproblem,
  //              THREAD_ID tid = 0);

  /**
   * Emtpy destructor
   */
  virtual ~Restartable() = default;

protected:
  /**
   * Declare a piece of data as "restartable".
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   */
  template <typename T>
  T & declareRestartableData(std::string data_name);

  /**
   * Declare a piece of data as "restartable" and initialize it.
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param init_value The initial value of the data
   */
  template <typename T>
  T & declareRestartableData(std::string data_name, const T & init_value);

  /**
   * Declare a piece of data as "restartable".
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param context Context pointer that will be passed to the load and store functions
   */
  template <typename T>
  T & declareRestartableDataWithContext(std::string data_name, void * context);

  /**
   * Declare a piece of data as "restartable" and initialize it.
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param init_value The initial value of the data
   * @param context Context pointer that will be passed to the load and store functions
   */
  template <typename T>
  T &
  declareRestartableDataWithContext(std::string data_name, const T & init_value, void * context);

  /**
   * NOTE: These are used internally in MOOSE.  NOT FOR PUBLIC CONSUMPTION!
   *
   * Declare a piece of data as "recoverable".
   * This means that in the event of a recovery this piece of data
   * will be restored back to its previous value.
   *
   * Note - this data will NOT be restored on _Restart_!
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   */
  template <typename T>
  T & declareRecoverableData(std::string data_name);

  /**
   * NOTE: These are used internally in MOOSE.  NOT FOR PUBLIC CONSUMPTION!
   *
   * Declare a piece of data as "restartable" and initialize it.
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * Note - this data will NOT be restored on _Restart_!
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param init_value The initial value of the data
   */
  template <typename T>
  T & declareRecoverableData(std::string data_name, const T & init_value);

  /**
   * Note: This is only used internally in MOOSE.  DO NOT use this function!
   *
   * Declare a piece of data as "restartable".
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param object_name A supplied name for the object that is declaring this data.
   */
  template <typename T>
  T & declareRestartableDataWithObjectName(std::string data_name, std::string object_name);

  /**
   * Note: This is only used internally in MOOSE.  DO NOT use this function!
   *
   * Declare a piece of data as "restartable".
   * This means that in the event of a restart this piece of data
   * will be restored back to its previous value.
   *
   * NOTE: This returns a _reference_!  Make sure you store it in a _reference_!
   *
   * @param data_name The name of the data (usually just use the same name as the member variable)
   * @param object_name A supplied name for the object that is declaring this data.
   * @param context Context pointer that will be passed to the load and store functions
   */
  template <typename T>
  T & declareRestartableDataWithObjectNameWithContext(std::string data_name,
                                                      std::string object_name,
                                                      void * context);

private:
  /// Helper function for actually registering the restartable data.
  void registerRestartableDataOnApp(std::string name, RestartableDataValue * data, THREAD_ID tid);

  /// Helper function for actually registering the restartable data.
  void registerRecoverableDataOnApp(std::string name);

  /// Reference to the application
  MooseApp & _restartable_app;

  /// The name of the object
  std::string _restartable_name;

  /// The system name this object is in
  std::string _restartable_system_name;

  /// The thread ID for this object
  THREAD_ID _restartable_tid;
};

template <typename T>
T &
Restartable::declareRestartableData(std::string data_name)
{
  return declareRestartableDataWithContext<T>(data_name, nullptr);
}

template <typename T>
T &
Restartable::declareRestartableData(std::string data_name, const T & init_value)
{
  return declareRestartableDataWithContext<T>(data_name, init_value, nullptr);
}

template <typename T>
T &
Restartable::declareRestartableDataWithContext(std::string data_name, void * context)
{
  std::string full_name = _restartable_system_name + "/" + _restartable_name + "/" + data_name;
  RestartableData<T> * data_ptr = new RestartableData<T>(full_name, context);

  registerRestartableDataOnApp(full_name, data_ptr, _restartable_tid);

  return data_ptr->get();
}

template <typename T>
T &
Restartable::declareRestartableDataWithContext(std::string data_name,
                                               const T & init_value,
                                               void * context)
{
  std::string full_name = _restartable_system_name + "/" + _restartable_name + "/" + data_name;
  RestartableData<T> * data_ptr = new RestartableData<T>(full_name, context);

  data_ptr->set() = init_value;

  registerRestartableDataOnApp(full_name, data_ptr, _restartable_tid);

  return data_ptr->get();
}

template <typename T>
T &
Restartable::declareRestartableDataWithObjectName(std::string data_name, std::string object_name)
{
  return declareRestartableDataWithObjectNameWithContext<T>(data_name, object_name, nullptr);
}

template <typename T>
T &
Restartable::declareRestartableDataWithObjectNameWithContext(std::string data_name,
                                                             std::string object_name,
                                                             void * context)
{
  std::string old_name = _restartable_name;

  _restartable_name = object_name;

  T & value = declareRestartableDataWithContext<T>(data_name, context);

  _restartable_name = old_name;

  return value;
}

template <typename T>
T &
Restartable::declareRecoverableData(std::string data_name)
{
  std::string full_name = _restartable_system_name + "/" + _restartable_name + "/" + data_name;

  registerRecoverableDataOnApp(full_name);

  return declareRestartableDataWithContext<T>(data_name, nullptr);
}

template <typename T>
T &
Restartable::declareRecoverableData(std::string data_name, const T & init_value)
{
  std::string full_name = _restartable_system_name + "/" + _restartable_name + "/" + data_name;

  registerRecoverableDataOnApp(full_name);

  return declareRestartableDataWithContext<T>(data_name, init_value, nullptr);
}

#endif // RESTARTABLE
