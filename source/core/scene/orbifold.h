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

enum OrbifoldType {
  XX, X2222, X333, X442, X632, TRIVIAL
};

struct OrbifoldDirection {
  Vector2d o, m, base; //8*2 * 3 = 48
  double v_sep, h_sep; // + 16 = 64
  double offsetArray[2]; // + 16 = 80
  unsigned indexArray[3]; // + 12 = 82
  unsigned numOffsets, numMirrors; // + 8 = 90
};


struct OrbifoldData {
  OrbifoldType type;
  Vector3d scale;
  double r1, r2, r3, r4;
  OrbifoldDirection mirrorDirs[4];

  shared_ptr<SceneData> sceneData;

  typedef double (OrbifoldData::*OrbifoldAttenuationCallback)(const Ray& r, const Intersection& isect) const;

  OrbifoldAttenuationCallback attenuateFunction;

  inline double attenuate(const Ray& r, const Intersection& i) const {
    return (this->*attenuateFunction)(r, i);
  }

  void countMirrorsHetero(const Vector2d& S, const Vector2d& E,
                          const Vector2d& dir, const OrbifoldDirection& d,
                          unsigned* mirrors) const;

  double attenuateX333(const Ray& ray, const Intersection& isect) const;
  double attenuateX2222(const Ray& ray, const Intersection& isect) const;
  double attenuateXX(const Ray& ray, const Intersection& isect) const;
  double attenuateX632(const Ray& ray, const Intersection& isect) const;
  double attenuateX442(const Ray& ray, const Intersection& isect) const;

  void InitX333OrbifoldData();

  static OrbifoldData MakeX2222OrbifoldData(Vector3d scale);
  static OrbifoldData MakeXXOrbifoldData(Vector3d scale);
  static OrbifoldData MakeX442OrbifoldData(Vector3d scale);
  static OrbifoldData MakeX6322OrbifoldData(Vector3d scale);
};

struct OrbifoldInfo {
  OrbifoldType type;
  Vector3d scale;
  float r1, r2, r3, r4;
  OrbifoldDirection mirrorDirs[4];
};

void initOrbifoldX333(OrbifoldInfo& info);

}


#endif /* SOURCE_CORE_SCENE_ORBIFOLD_H */
