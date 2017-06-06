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

#ifndef KDTREE_H
#define KDTREE_H

// Moose includes
#include "MooseMesh.h"

// Libmesh includes
#include "libmesh/nanoflann.hpp"

class KDTree
{
public:
  KDTree(std::vector<Point> & master_points);

  virtual ~KDTree() = default;

  void buildTree();

  void neighborSearch(Point & query_point,
                      unsigned int patch_size,
                      std::vector<unsigned int> & return_index);

  /**
   * PointListAdaptor is required to use libMesh Point coordinate type with
   * nanoflann KDTree library. The member functions within the PointListAdaptor
   * are used by nanoflann library.
   */
  template <unsigned int KDDim>
  class PointListAdaptor
  {
  private:
    const std::vector<Point> & _pts;

  public:
    PointListAdaptor(const std::vector<Point> & pts) : _pts(pts) {}

    /**
     * libMesh \p Point coordinate type
     */
    using coord_t = Real;
    /**
     * Must return the number of data points
     */
    inline size_t kdtree_get_point_count() const { return _pts.size(); }

    /**
     * Returns the distance between the vector "p1[0:size-1]"
     * and the data point with index "idx_p2" stored in the class
     */
    inline coord_t kdtree_distance(const coord_t * p1, const size_t idx_p2, size_t size) const
    {
      mooseAssert(size == LIBMESH_DIM, "Size of point should equal libmesh dimension.");
      mooseAssert(idx_p2 <= _pts.size(),
                  "The point index should be less than"
                  "total number of points used to build"
                  "the KDTree.");

      const Point & p2(_pts[idx_p2]);

      switch (size)
      {
        case 3:
        {
          const coord_t d0 = p1[0] - p2(0);
          const coord_t d1 = p1[1] - p2(1);
          const coord_t d2 = p1[2] - p2(2);

          return d0 * d0 + d1 * d1 + d2 * d2;
        }

        case 2:
        {
          const coord_t d0 = p1[0] - p2(0);
          const coord_t d1 = p1[1] - p2(1);

          return d0 * d0 + d1 * d1;
        }

        case 1:
        {
          const coord_t d0 = p1[0] - p2(0);

          return d0 * d0;
        }

        default:
          mooseError("Unknown size");
      }

      return -1.;
    }

    /**
     * Returns the dim'th component of the idx'th point in the class:
     * Since this is inlined and the "dim" argument is typically an immediate
     * value, the
     *  "if's" are actually solved at compile time.
     */
    inline coord_t kdtree_get_pt(const size_t idx, int dim) const
    {
      mooseAssert(dim < (int)LIBMESH_DIM,
                  "The required component number should be less than the LIBMESH_DIM.");
      mooseAssert(idx < _pts.size(),
                  "The index of the point should be less"
                  "than total number of points used to"
                  "construct the KDTree.");

      const Point & p(_pts[idx]);

      if (dim == 0)
        return p(0);
      if (dim == 1)
        return p(1);
      return p(2);
    }

    /**
     * Optional bounding-box computation: return false to default to a standard
     * bbox computation loop.
     * Return true if the BBOX was already computed by the class and returned in
     * "bb" so it can be
     * avoided to redo it again. Look at bb.size() to find out the expected
     * dimensionality
     * (e.g. 2 or 3 for point clouds)
     */
    template <class BBOX>
    bool kdtree_get_bbox(BBOX & /* bb */) const
    {
      return false;
    }
  };
  typedef nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<Real, PointListAdaptor<LIBMESH_DIM>>,
      PointListAdaptor<LIBMESH_DIM>,
      LIBMESH_DIM>
      KdTreeT;

protected:
  mutable UniquePtr<KdTreeT> _kd_tree;
  PointListAdaptor<LIBMESH_DIM> _point_list_adaptor;
};

#endif // KDTREE_H
