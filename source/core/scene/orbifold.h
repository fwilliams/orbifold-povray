//******************************************************************************
///
/// @file core/scene/orbifold.h
///
/// Declarations for routines to compute radiance in scenes with certain symmetry structures
///
/// @copyright
/// @parblock
///
/// Persistence of Vision Ray Tracer ('POV-Ray') version 3.7.
/// Copyright 1991-2015 Persistence of Vision Raytracer Pty. Ltd.
///
/// POV-Ray is free software: you can redistribute it and/or modify
/// it under the terms of the GNU Affero General Public License as
/// published by the Free Software Foundation, either version 3 of the
/// License, or (at your option) any later version.
///
/// POV-Ray is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU Affero General Public License for more details.
///
/// You should have received a copy of the GNU Affero General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.
///
/// ----------------------------------------------------------------------------
///
/// POV-Ray is based on the popular DKB raytracer version 2.12.
/// DKBTrace was originally written by David K. Buck.
/// DKBTrace Ver 2.0-2.12 were written by David K. Buck & Aaron A. Collins.
///
/// @endparblock
///
//******************************************************************************

#include "core/math/vector.h"
#include "core/render/ray.h"
#include "core/coretypes.h"
#include "core/scene/scenedata.h"

#ifndef SOURCE_CORE_SCENE_ORBIFOLD_H
#define SOURCE_CORE_SCENE_ORBIFOLD_H

namespace pov {

struct OrbifoldData {
  OrbifoldData() :
    scale(Vector3d(1, 1, 1)), r1(0), r2(0), r3(0), r4(0),
    num_kernel_tiles(1),
    attenuate_callback(&OrbifoldData::attenuateTrivial),
    collapse_callback(&OrbifoldData::collapseTrivial) {}

  /*
   * A scaling factor for the scense
   */
  Vector3d scale;

  /*
   * The reflectivity of each mirror. Depending on the orbifold type,
   * not all of these are used
   */
  double r1, r2, r3, r4;

  /*
   * The number of tiles in the kernel being rendered
   */
  unsigned num_kernel_tiles;

  /*
   * Computes the attenuation along a ray which may intersect one or more mirror boundaries
   * in the unfolded scene
   */
  inline double attenuate(const Ray& r, const Intersection& i) const {
    return (this->*attenuate_callback)(r, i);
  }

  /*
   * Transforms (collapses) a point in the covering space into the fundamental domain.
   * This function returns the transformed point. out_fixed_dir is mutated to point
   * in the direction corresponding to the transformation
   */
  inline Vector3d collapse(const Vector3d& point, Vector3d& out_fixed_dir) const {
    return (this->*collapse_callback)(point, out_fixed_dir);
  }

  /*
   * Transforms (collapses) a point in the covering space into the fundamental domain.
   * This function returns the transformed point.
   */
  inline Vector3d collapse(const Vector3d& point) const {
    Vector3d dummy;
    return (this->*collapse_callback)(point, dummy);
  }

  /*
   * Initialize this structure based on the type of orbifold present
   */
  void InitX333OrbifoldData();
  void InitX2222OrbifoldData();
  void InitXXOrbifoldData();
  void InitX442OrbifoldData();
  void InitX632OrbifoldData();

  /*
   * Convert a kernel radius to a number of tiles for each group
   */
  static inline unsigned X333_KernelRadToTileNum(unsigned kernel_rad) {
    unsigned ret = 1;
    for(unsigned i = 1; i <= kernel_rad; i++) {
      ret += 6 * i;
    }
    return ret * 6;
  }
  static inline unsigned X2222_KernelRadToTileNum(unsigned kernel_rad) {
    unsigned ret = 1;
    for(unsigned i = 1; i <= kernel_rad; i++) {
      ret += 8 * i;
    }
    return ret * 4;
  }
  static inline unsigned XX_KernelRadToTileNum(unsigned kernel_rad) {
    return 4 * kernel_rad + 1;
  }
  static inline unsigned X442_KernelRadToTileNum(unsigned kernel_rad) {
    return X2222_KernelRadToTileNum(kernel_rad) * 2;
  }
  static inline unsigned X632_KernelRadToTileNum(unsigned kernel_rad) {
    return X333_KernelRadToTileNum(kernel_rad) * 2;
  }

private:

  /*
   * Static information about mirror directions
   */
  struct StaticDirectionInfo {
    Vector2d o, m;
    unsigned index_array[3];
    static const unsigned NUM_MIRRORS = 3;
  };

  /*
   * Runtime information about mirror directions
   */
  struct DynamicDirectionInfo {
    Vector2d base;
    double v_sep, h_sep;
    double offset_array[2];
    static const unsigned OFFSET_ARRAY_LEN = 2;
  };

  /*
   * Mirror attenuation callback type
   */
  typedef double (OrbifoldData::*OrbifoldAttenuationCallback)(const Ray& r, const Intersection& isect) const;

  /*
   * Point collapse callback (see collapse)
   */
  typedef Vector3d (OrbifoldData::*OrbifoldCollapseCallback)(const Vector3d& pt, Vector3d& ofd) const;

  /*
   * Information about adjacent mirror directions set at runtime
   */
  DynamicDirectionInfo direction_info[2];

  /*
   * The callback function which calculates attenuation from mirror boundary intersections
   * This gets set at runtime to one of the attenuate*() functions based on the type of orbifold
   * specified.
   */
  OrbifoldAttenuationCallback attenuate_callback;

  /*
   * The callback function which collapses a virtual point into the fundamental domain.
   * This gets set at runtime to one of the collapse*() functions based on the type of orbifold specified.
   */
  OrbifoldCollapseCallback collapse_callback;

  /*
   * Called when each mirror plane consists of a repeating pattern of different mirrors.
   * Counts which mirrors are intersected by the straight path starting at S, and ending at E
   * where the mirrors are configured according to s_dir and d_dir.
   * The result is stored in mirror_count.
   */
  void countMirrorsHeterogeneous(const Vector2d& S, const Vector2d& E,
                          const Vector2d& dir, const StaticDirectionInfo& s_dir,
                          const DynamicDirectionInfo& d_dir, unsigned* mirrors) const;


  void countMirrorsHomogeneous(const Vector2d& S, const Vector2d& E,
                          const Vector2d& dir, const StaticDirectionInfo& s_dir,
                          const DynamicDirectionInfo& d_dir, unsigned* mirrors) const;

  unsigned countMirrorsHomogeneous(const Vector2d& S, const Vector2d& E,
                          const Vector2d& dir, const StaticDirectionInfo& s_dir,
                          const DynamicDirectionInfo& d_dir) const;

  double attenuateTrivial(const Ray& ray, const Intersection& isect) const {
    return 1.0;
  }
  double attenuateX333(const Ray& ray, const Intersection& isect) const;
  double attenuateX2222(const Ray& ray, const Intersection& isect) const;
  double attenuateXX(const Ray& ray, const Intersection& isect) const;
  double attenuateX632(const Ray& ray, const Intersection& isect) const;
  double attenuateX442(const Ray& ray, const Intersection& isect) const;

  Vector3d collapseTrivial(const Vector3d& pt, Vector3d& o) const { return pt; }
  Vector3d collapseX333(const Vector3d& pt, Vector3d& o) const;
  Vector3d collapseX632(const Vector3d& pt, Vector3d& o) const;
  Vector3d collapseX442(const Vector3d& pt, Vector3d& o) const;
  Vector3d collapseX2222(const Vector3d& pt, Vector3d& o) const;
  Vector3d collapseXX(const Vector3d& pt, Vector3d& o) const;
};

}


#endif /* SOURCE_CORE_SCENE_ORBIFOLD_H */
