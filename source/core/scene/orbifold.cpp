#include "core/scene/orbifold.h"

namespace pov {

void OrbifoldData::InitX333OrbifoldData() {
  attenuateFunction = &OrbifoldData::attenuateX333;

  const double SQRT_3_OVER_2 = 0.8660254037844386;
  const double SQRT_3_OVER_4 = 0.4330127018922193;

  // Horizontal mirrors
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
}

void OrbifoldData::InitX2222OrbifoldData() {
  attenuateFunction = &OrbifoldData::attenuateX2222;
}

void OrbifoldData::InitXXOrbifoldData() {
  attenuateFunction = &OrbifoldData::attenuateXX;
}

static inline unsigned fastMax(unsigned x, unsigned y) {
  return (x ^ ((x ^ y) & -(x < y)));
}

static inline int fastFloor(double x) {
  int i = (int)x; /* truncate */
  return i - ( i > x ); /* convert trunc to floor */
}

static inline int fastCeil(double x) {
  int i = (int)x; /* truncate and add 1 */
  return i + 1 - ( i > x ); /* convert trunc to ceil */
}

void OrbifoldData::countMirrorsHetero(const Vector2d& S, const Vector2d& E,
    const Vector2d& dir, const OrbifoldDirection& d, unsigned* mirrors) const {
  const double So = dot(S, d.o);
  const double dir_dot_o = dot(dir, d.o);

  int m0 = 0, mE = 0;

  if(dir_dot_o > 0) {
    m0 = fastCeil(So / d.v_sep);
    mE = fastFloor(dot(E, d.o) / d.v_sep);
  } else if (dir_dot_o < 0) {
    m0 = fastFloor(So / d.v_sep);
    mE = fastCeil(dot(E, d.o) / d.v_sep);
  } else {
    // Disregard rays parallel to this mirror direction since they won't intersect any mirrors
    return;
  }

  // Bail out if ray is parallel to mirrors
//  if(dir_dot_o == 0) { return; }
//
//  const int signbit = static_cast<int>(dir_dot_o < 0);
//  const int m0 = fastCeil(So / d.v_sep) - signbit;
//  const int mE = fastFloor(dot(E, d.o) / d.v_sep) + signbit;

  const double k0 = (m0 * d.v_sep - So) / dir_dot_o;
  const double kE = (mE * d.v_sep - So) / dir_dot_o;

  // Bail out if the path doesn't cross a mirror
  if(kE < k0) { return; }


  const unsigned numMirrors = abs(mE - m0);
  const double  dK = (kE - k0) / fastMax(numMirrors, 1);//max<unsigned>(numMirrors, 1);

  double k = k0;
  unsigned i = abs(m0) % d.numOffsets;

  for(unsigned j = 0; j <= numMirrors; j++) {
    // Where the ray intersects the mirror line
    const Vector2d P = S + k * dir;

    // Project the intersection onto m, and figure out which type of mirror was intersected
    const int l = floor((dot(P, d.m) + d.offsetArray[i]) / d.h_sep);
    mirrors[d.indexArray[abs(l % (int)(d.numMirrors))]] += 1;

    k += dK;
    i = (i + 1) % d.numOffsets;
  }
}

/*
 * Orbifold attenuation functions
 */


double OrbifoldData::attenuateX2222(const Ray& ray, const Intersection& isect) const {
  typedef GenericVector2d<POV_INT64> IntVector2d;

  const Vector3d rayStart = ray.Origin / scale;
  const Vector3d rayEnd = isect.IPoint / scale;

  const IntVector2d latticeStart(floor(rayStart.x() + 0.5), floor(rayStart.z() + 0.5));
  const IntVector2d latticeEnd(floor(rayEnd.x() + 0.5), floor(rayEnd.z() + 0.5));
  const IntVector2d latticeDist(abs(latticeStart.x() - latticeEnd.x()), abs(latticeStart.y() - latticeEnd.y()));

  // Bail out if the ray doesn't attenuate
  if(latticeDist.x() == 0 && latticeDist.y() == 0) {
    return 1.0;
  }

  // Whether we're starting from even and odd indices in x and y
  const IntVector2d latticeMod(latticeStart.x() % 2, latticeStart.y() % 2);

  // Mirror vector is (left, top, right, bottom)
  POV_INT64 n0 = 0, n1 = 0, n2 = 0, n3 = 0;
  POV_INT64 x1 = latticeDist.x() / 2;
  POV_INT64 x2 = latticeDist.x() - x1;
  POV_INT64 y1 = latticeDist.y() / 2;
  POV_INT64 y2 = latticeDist.y() - y1;

  // Even and going left (-) or odd and going right (+)
  if((latticeMod.x() == 0 && latticeStart.x() > latticeEnd.x()) ||
     (latticeMod.x() == 1 && latticeStart.x() < latticeEnd.x())) {
    n0 = x2;
    n2 = x1;
  } else if(latticeDist.x() != 0) { // Even and going right (+) or Odd and going left (-)
    n0 = x1;
    n2 = x2;
  }

  // Even and going down (-) or odd and going up (+)
  if((latticeMod.y() == 0 && latticeStart.y() > latticeEnd.y()) ||
     (latticeMod.y() == 1 && latticeStart.y() < latticeEnd.y())) {
    n3 = y2;
    n1 = y1;
  } else if(latticeDist.y() != 0) { // Even and going up (+) or Odd and going down (-)
    n3 = y1;
    n1 = y2;
  }

  return pow(r1, (double)n0) * pow(r2, (double)n1) *
         pow(r3, (double)n2) * pow(r4, (double)n3);
}

double OrbifoldData::attenuateXX(const Ray& ray, const Intersection& isect) const {
  typedef GenericVector2d<POV_INT64> IntVector2d;

  const Vector3d rayStart = ray.Origin / scale;
  const Vector3d rayEnd = isect.IPoint / scale;

  const POV_INT64 latticeStart = floor(rayStart.z() + 0.5);
  const POV_INT64 latticeEnd = floor(rayEnd.z() + 0.5);
  const POV_INT64 latticeDist = abs(latticeStart - latticeEnd);

  // Bail out if the ray doesn't attenuate
  if(latticeDist == 0) { return 1.0; }

  // Whether we're starting from even and odd indices in x and y
  const POV_INT64 latticeMod = latticeStart % 2;

  // Mirror vector is (left, top, right, bottom)
  POV_INT64 n0 = 0, n1 = 0, n2 = 0, n3 = 0;
  POV_INT64 x1 = latticeDist / 2;
  POV_INT64 x2 = latticeDist - x1;

  // Even and going backwards (-) or odd and going forwards (+)
  if((latticeMod == 0 && latticeStart > latticeEnd) ||
     (latticeMod == 1 && latticeStart < latticeEnd)) {
    n0 = x2;
    n1 = x1;
  } else if(latticeDist != 0) { // Even and going right (+) or Odd and going left (-)
    n0 = x1;
    n1 = x2;
  }

  return pow(r1, (double)n0) * pow(r2, (double)n1);
}

double OrbifoldData::attenuateX632(const Ray& ray, const Intersection& isect) const {
  return 0;
}

double OrbifoldData::attenuateX333(const Ray& ray, const Intersection& isect) const {
  const Vector2d B = mirrorDirs[0].base;
  const Vector2d S = Vector2d(ray.Origin.x(), ray.Origin.z()) - B;
  const Vector2d E = Vector2d(isect.IPoint.x(), isect.IPoint.z()) - B;

  unsigned mirrorCount[3] = {0, 0, 0};

  countMirrorsHetero(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), mirrorDirs[0], mirrorCount);
  countMirrorsHetero(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), mirrorDirs[1], mirrorCount);
  countMirrorsHetero(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), mirrorDirs[2], mirrorCount);

  return pow(r1, mirrorCount[0]) * pow(r2, mirrorCount[1]) * pow(r3, mirrorCount[2]);
  //  return att333(ray, isect, *this);
}

double OrbifoldData::attenuateX442(const Ray& ray, const Intersection& isect) const {
  return 0;
}

}
