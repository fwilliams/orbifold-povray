#include "core/scene/orbifold.h"

namespace pov {

void OrbifoldData::InitX333OrbifoldData() {
  attenuateFunction = &OrbifoldData::attenuateX333;

  const double SQRT_3_OVER_2 = 0.8660254037844386;
  const double SQRT_3_OVER_4 = 0.4330127018922193;

  mirrorDirs[0].base = Vector2d(-0.5, -SQRT_3_OVER_4) * scale.x();
  mirrorDirs[0].o = Vector2d(0, 1);
  mirrorDirs[0].m = Vector2d(1, 0);

  mirrorDirs[0].numOffsets = 2;

  mirrorDirs[0].offsetArray[0] = 0;
  mirrorDirs[0].offsetArray[1] = 1.5 * scale.x();

  mirrorDirs[0].numMirrors = 3;
  mirrorDirs[0].indexArray[0] = 0;
  mirrorDirs[0].indexArray[1] = 1;
  mirrorDirs[0].indexArray[2] = 2;

  mirrorDirs[0].v_sep = SQRT_3_OVER_2 * scale.x();
  mirrorDirs[0].h_sep = scale.x();


  // Right tilted mirrors
  mirrorDirs[1].base = Vector2d(-0.5, -SQRT_3_OVER_4) * scale.x();
  mirrorDirs[1].o = Vector2d(-SQRT_3_OVER_2, 0.5);
  mirrorDirs[1].m = Vector2d(0.5, SQRT_3_OVER_2);

  mirrorDirs[1].numOffsets = 2;

  mirrorDirs[1].offsetArray[0] = 0;
  mirrorDirs[1].offsetArray[1] = 1.5 * scale.x();

  mirrorDirs[1].numMirrors = 3;
  mirrorDirs[1].indexArray[0] = 2;
  mirrorDirs[1].indexArray[1] = 1;
  mirrorDirs[1].indexArray[2] = 0;

  mirrorDirs[1].v_sep = SQRT_3_OVER_2 * scale.x();
  mirrorDirs[1].h_sep = scale.x();


  // Left tilted mirrors
  mirrorDirs[2].base = Vector2d(-0.5, -SQRT_3_OVER_4) * scale.x();
  mirrorDirs[2].o = Vector2d(SQRT_3_OVER_2, 0.5);
  mirrorDirs[2].m = Vector2d(-0.5, SQRT_3_OVER_2);

  mirrorDirs[2].numOffsets = 2;

  mirrorDirs[2].offsetArray[0] = 0;
  mirrorDirs[2].offsetArray[1] = 1.5 * scale.x();

  mirrorDirs[2].numMirrors = 3;
  mirrorDirs[2].indexArray[0] = 0;
  mirrorDirs[2].indexArray[1] = 1;
  mirrorDirs[2].indexArray[2] = 2;

  mirrorDirs[2].v_sep = SQRT_3_OVER_2 * scale.x();
  mirrorDirs[2].h_sep = scale.x();

//  mirrorDirs[0].base = Vector2d(-0.5, -SQRT_3_OVER_4) * scale.x();
//  mirrorDirs[0].o = Vector2d(0, 1);
//  mirrorDirs[0].m = Vector2d(1, 0);
//  mirrorDirs[0].numOffsets = 2;
//  mirrorDirs[0].offsetArray[0] = 0;
//  mirrorDirs[0].offsetArray[1] = 1.5 * scale.x();
//  mirrorDirs[0].numMirrors = 3;
//  mirrorDirs[0].indexArray[0] = 0;
//  mirrorDirs[0].indexArray[1] = 1;
//  mirrorDirs[0].indexArray[2] = 2;
//
//  mirrorDirs[1].base = Vector2d(-0.5, -SQRT_3_OVER_4) * scale.x();
//  mirrorDirs[1].o = Vector2d(-SQRT_3_OVER_2, 0.5);
//  mirrorDirs[1].m = Vector2d(0.5, SQRT_3_OVER_2);
//  mirrorDirs[1].numOffsets = 2;
//  mirrorDirs[1].offsetArray[0] = 0;
//  mirrorDirs[1].offsetArray[1] = 1.5 * scale.x();
//  mirrorDirs[1].numMirrors = 3;
//  mirrorDirs[1].indexArray[0] = 0;
//  mirrorDirs[1].indexArray[1] = 1;
//  mirrorDirs[1].indexArray[2] = 2;
//
//  mirrorDirs[2].base = Vector2d(-0.5, -SQRT_3_OVER_4) * scale.x();
//  mirrorDirs[2].o = Vector2d(SQRT_3_OVER_2, 0.5);
//  mirrorDirs[2].m = Vector2d(-0.5, SQRT_3_OVER_2);
//  mirrorDirs[2].numOffsets = 2;
//  mirrorDirs[2].offsetArray[0] = 0;
//  mirrorDirs[2].offsetArray[1] = 1.5 * scale.x();
//  mirrorDirs[2].numMirrors = 3;
//  mirrorDirs[2].indexArray[0] = 0;
//  mirrorDirs[2].indexArray[1] = 1;
//  mirrorDirs[2].indexArray[2] = 2;

}

inline unsigned fastMax(unsigned x, unsigned y) {
  return (x ^ ((x ^ y) & -(x < y)));
}

void OrbifoldData::countMirrorsHetero(const Vector2d& S, const Vector2d& E,
    const Vector2d& dir, const OrbifoldDirection& d, unsigned* mirrors) const {
  const double So = dot(S, d.o);
  const double dir_dot_o = dot(dir, d.o);

  int m0 = 0, mE = 0;

  if(dir_dot_o > 0) {
    m0 = ceil(So / d.v_sep);
    mE = floor(dot(E, d.o) / d.v_sep);
  } else if (dir_dot_o < 0) {
    m0 = floor(So / d.v_sep);
    mE = ceil(dot(E, d.o) / d.v_sep);
  } else {
    return; // Disregard this direction of mirrors
  }

  const double k0 = (m0 * d.v_sep - So) / dir_dot_o;
  const double kE = (mE * d.v_sep - So) / dir_dot_o;

  // Bail out if the path doesn't cross a mirror
  if(kE < k0) return;


  const unsigned numMirrors = abs(mE - m0);
  const double  dK = (kE - k0) / fastMax(numMirrors, 1);//max<unsigned>(numMirrors, 1);

  // TODO: What if point lies on a mirror?
  double k = k0;
  unsigned i = abs(m0) % d.numOffsets;

  for(unsigned j = 0; j <= numMirrors; j++) {
    // Where the ray intersects the mirror line
    const Vector2d P = S + k * dir; // 4 FLOPS

    // Project the intersection onto m
    const int l = (dot(P, d.m) + d.offsetArray[i]) / d.h_sep; // 5 FLOPS
    mirrors[d.indexArray[l % d.numMirrors]] += 1; // 4 INT OPS

    k += dK;
    i = (i + 1) % d.numOffsets; // 2 INT OPS
  }
}

/*
 * Orbifold attenuation functions
 */

double OrbifoldData::attenuateX333(const Ray& ray, const Intersection& isect) const {
    const Vector2d B = mirrorDirs[0].base;
    const Vector2d S = Vector2d(ray.Origin.x(), ray.Origin.z()) - B;
    const Vector2d E = Vector2d(isect.IPoint.x(), isect.IPoint.z()) - B;

    unsigned mirrorCount[3] = {0, 0, 0};

    countMirrorsHetero(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), mirrorDirs[0], mirrorCount);
    countMirrorsHetero(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), mirrorDirs[1], mirrorCount);
    countMirrorsHetero(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), mirrorDirs[2], mirrorCount);

    return pow(r1, mirrorCount[0]) *
           pow(r2, mirrorCount[1]) *
           pow(r3, mirrorCount[2]);
}

double OrbifoldData::attenuateX2222(const Ray& ray, const Intersection& isect) const {
  return 0;
}

double OrbifoldData::attenuateXX(const Ray& ray, const Intersection& isect) const {
  return 0;
}

double OrbifoldData::attenuateX632(const Ray& ray, const Intersection& isect) const {
  return 0;
}

double OrbifoldData::attenuateX442(const Ray& ray, const Intersection& isect) const {
  return 0;
}


void initOrbifoldX333(OrbifoldInfo& info) {
  const double SQRT_3_OVER_2 = 0.8660254037844386;
  const double SQRT_3_OVER_4 = 0.4330127018922193;

  // Horizontal mirrors
  info.mirrorDirs[0].base = Vector2d(-0.5, -SQRT_3_OVER_4) * info.scale.x();
  info.mirrorDirs[0].o = Vector2d(0, 1);
  info.mirrorDirs[0].m = Vector2d(1, 0);

  info.mirrorDirs[0].numOffsets = 2;

  info.mirrorDirs[0].offsetArray[0] = 0;
  info.mirrorDirs[0].offsetArray[1] = 1.5 * info.scale.x();

  info.mirrorDirs[0].numMirrors = 3;
  info.mirrorDirs[0].indexArray[0] = 0;
  info.mirrorDirs[0].indexArray[1] = 1;
  info.mirrorDirs[0].indexArray[2] = 2;

  info.mirrorDirs[0].v_sep = SQRT_3_OVER_2 * info.scale.x();
  info.mirrorDirs[0].h_sep = info.scale.x();


  // Right tilted mirrors
  info.mirrorDirs[1].base = Vector2d(-0.5, -SQRT_3_OVER_4) * info.scale.x();
  info.mirrorDirs[1].o = Vector2d(-SQRT_3_OVER_2, 0.5);
  info.mirrorDirs[1].m = Vector2d(0.5, SQRT_3_OVER_2);

  info.mirrorDirs[1].numOffsets = 2;

  info.mirrorDirs[1].offsetArray[0] = 0;
  info.mirrorDirs[1].offsetArray[1] = 1.5 * info.scale.x();

  info.mirrorDirs[1].numMirrors = 3;
  info.mirrorDirs[1].indexArray[0] = 2;
  info.mirrorDirs[1].indexArray[1] = 1;
  info.mirrorDirs[1].indexArray[2] = 0;

  info.mirrorDirs[1].v_sep = SQRT_3_OVER_2 * info.scale.x();
  info.mirrorDirs[1].h_sep = info.scale.x();


  // Left tilted mirrors
  info.mirrorDirs[2].base = Vector2d(-0.5, -SQRT_3_OVER_4) * info.scale.x();
  info.mirrorDirs[2].o = Vector2d(SQRT_3_OVER_2, 0.5);
  info.mirrorDirs[2].m = Vector2d(-0.5, SQRT_3_OVER_2);

  info.mirrorDirs[2].numOffsets = 2;

  info.mirrorDirs[2].offsetArray[0] = 0;
  info.mirrorDirs[2].offsetArray[1] = 1.5 * info.scale.x();

  info.mirrorDirs[2].numMirrors = 3;
  info.mirrorDirs[2].indexArray[0] = 0;
  info.mirrorDirs[2].indexArray[1] = 1;
  info.mirrorDirs[2].indexArray[2] = 2;

  info.mirrorDirs[2].v_sep = SQRT_3_OVER_2 * info.scale.x();
  info.mirrorDirs[2].h_sep = info.scale.x();
}

}
