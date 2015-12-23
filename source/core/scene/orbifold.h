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
    scale(Vector3d(1, 1, 1)), r1(0), r2(0),
    r3(0), r4(0), attenuate_callback(&OrbifoldData::attenuateTrivial) {}

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
   * Computes the attenuation along a ray which may intersect one or more mirror boundaries
   * in the unfolded scene
   */
  inline double attenuate(const Ray& r, const Intersection& i) const {
    return (this->*attenuate_callback)(r, i);
  }

  /*
   * Initialize this structure based on the type of orbifold present
   */
  void InitX333OrbifoldData();
  void InitX2222OrbifoldData();
  void InitXXOrbifoldData();
  void InitX442OrbifoldData();
  void InitX632OrbifoldData();

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

  double attenuateTrivial(const Ray& ray, const Intersection& isect) const { return 1.0; }
  double attenuateX333(const Ray& ray, const Intersection& isect) const;
  double attenuateX2222(const Ray& ray, const Intersection& isect) const;
  double attenuateXX(const Ray& ray, const Intersection& isect) const;
  double attenuateX632(const Ray& ray, const Intersection& isect) const;
  double attenuateX442(const Ray& ray, const Intersection& isect) const;
};

}


#endif /* SOURCE_CORE_SCENE_ORBIFOLD_H */
