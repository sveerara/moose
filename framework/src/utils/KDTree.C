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

#include "KDTree.h"
#include "MooseError.h"

#include "libmesh/nanoflann.hpp"

KDTree::KDTree(std::vector<Point> & master_points)
  : _kd_tree(nullptr), _point_list_adaptor(master_points)
{
}

void
KDTree::buildTree()
{
  // Maximum number of points in each leaf of the Kd tree
  unsigned int max_leaf_size = 10;
  if (_kd_tree == nullptr)
    _kd_tree.reset(new KdTreeT(LIBMESH_DIM,
                               _point_list_adaptor,
                               nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size)));

  mooseAssert(_kd_tree != nullptr, "KDTree was not properly initalized.");

  _kd_tree->buildIndex();
}

void
KDTree::neighborSearch(Point & query_point,
                       unsigned int patch_size,
                       std::vector<unsigned int> & return_index)
{
  // The query point has to be converted from a C++ array to a C array because nanoflann library
  // expects C arrays.
  const Real query_pt[] = {query_point(0), query_point(1), query_point(2)};

  // Again size_t is used for ret_index because nanoflann expects that
  std::vector<std::size_t> ret_index(patch_size);
  std::vector<Real> ret_dist_sqr(patch_size);

  for (unsigned int i = 0; i < patch_size; ++i)
  {
    ret_dist_sqr[i] = std::numeric_limits<Real>::max();
    ret_index[i] = std::numeric_limits<std::size_t>::max();
  }

  _kd_tree->knnSearch(&query_pt[0], patch_size, &ret_index[0], &ret_dist_sqr[0]);

  if (ret_dist_sqr[0] == std::numeric_limits<Real>::max() ||
      ret_index[0] == std::numeric_limits<std::size_t>::max())
    mooseError("Unable to find closest node!");

  // ret_index to return_index conversion is basically conversion from std::size_t to unsigned int
  for (unsigned int i = 0; i < patch_size; ++i)
    return_index[i] = ret_index[i];
}
