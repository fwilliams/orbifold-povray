#include "core/scene/orbifold.h"

namespace pov {


/*
 * Custom math functions which perform better than the ones provided by the standard library
 * max, floor, and ceil all perform expensive branching, these ones don't branch and give
 * a significant speedup to functions in this module which are called for every ray
 */

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



/*
 * OrbifoldData initialization functions
 */

void OrbifoldData::InitX333OrbifoldData() {
  attenuate_callback = &OrbifoldData::attenuateX333;

  const double SQRT_3_OVER_2 = 0.8660254037844386;
  const double SQRT_3_OVER_4 = 0.4330127018922193;

  direction_info[0].h_sep = scale.x();
  direction_info[0].v_sep = SQRT_3_OVER_2 * scale.x();
  direction_info[0].offset_array[0] = 0;
  direction_info[0].offset_array[1] = 1.5 * scale.x();
  direction_info[0].base = Vector2d(-0.5, -SQRT_3_OVER_4) * scale.x();
}

void OrbifoldData::InitX2222OrbifoldData() {
  attenuate_callback = &OrbifoldData::attenuateX2222;
}

void OrbifoldData::InitXXOrbifoldData() {
  attenuate_callback = &OrbifoldData::attenuateXX;
}

/*
 * Counts which mirrors are intersected by the straight path starting at S, and ending at E
 * where the mirrors are configured according to s_dir and d_dir.
 * The result is stored in mirror_count.
 */
void OrbifoldData::countMirrorsHeterogeneous(const Vector2d& S, const Vector2d& E,
    const Vector2d& dir, const StaticDirectionInfo& s_dir, const DynamicDirectionInfo& d_dir,
    unsigned* mirror_count) const {

  // The projection of the start position onto the direction perpendicular to the mirror plane
  const double So = dot(S, s_dir.o);

  // The projection of the path direction onto the direction perpendicular to the mirror plane
  const double dir_dot_o = dot(dir, s_dir.o);

  // Compute the index of the first and last mirror planes intersected by the path
  int m0 = 0, mE = 0;
  if(dir_dot_o > 0) {
    m0 = fastCeil(So / d_dir.v_sep);
    mE = fastFloor(dot(E, s_dir.o) / d_dir.v_sep);
  } else if (dir_dot_o < 0) {
    m0 = fastFloor(So / d_dir.v_sep);
    mE = fastCeil(dot(E, s_dir.o) / d_dir.v_sep);
  } else {
    // Disregard rays parallel to this mirror direction since they won't intersect any mirrors
    return;
  }

  // The distance along the path of the first and last mirror intersections
  const double k0 = (m0 * d_dir.v_sep - So) / dir_dot_o;
  const double kE = (mE * d_dir.v_sep - So) / dir_dot_o;

  // Bail out if the path doesn't cross a mirror
  if(kE < k0) { return; }

  // The total number of mirrors intersected by the path
  const unsigned num_mirrors = abs(mE - m0);

  // dK is the distance along the path between mirrors
  const double  dK = (kE - k0) / fastMax(num_mirrors, 1);

  // k is the distance along the path of each mirror plane intersected
  double k = k0;
  // i keeps track of the type of mirror boundary
  unsigned i = abs(m0) % DynamicDirectionInfo::OFFSET_ARRAY_LEN;

  for(unsigned j = 0; j <= num_mirrors; j++) {
    // Where the ray intersects the mirror line
    const Vector2d P = S + k * dir;

    // Project the intersection onto the mirror plane, and figure out which type of mirror was intersected
    const int l = fastFloor((dot(P, s_dir.m) + d_dir.offset_array[i]) / d_dir.h_sep);
    mirror_count[s_dir.index_array[abs(l % static_cast<int>(StaticDirectionInfo::NUM_MIRRORS))]] += 1;

    k += dK;
    i = (i + 1) % DynamicDirectionInfo::OFFSET_ARRAY_LEN;
  }
}


unsigned OrbifoldData::countMirrorsHomogeneous(const Vector2d& S, const Vector2d& E,
                        const Vector2d& dir, const StaticDirectionInfo& s_dir,
                        const DynamicDirectionInfo& d_dir) const {
  // The projection of the start position onto the direction perpendicular to the mirror plane
  const double So = dot(S, s_dir.o);

  // The projection of the path direction onto the direction perpendicular to the mirror plane
  const double dir_dot_o = dot(dir, s_dir.o);

  // Compute the index of the first and last mirror planes intersected by the path
  int m0 = 0, mE = 0;
  if(dir_dot_o > 0) {
    m0 = fastCeil(So / d_dir.v_sep);
    mE = fastFloor(dot(E, s_dir.o) / d_dir.v_sep);
  } else if (dir_dot_o < 0) {
    // Note: Here the sign is inversed so that m0 < mE if the ray path intersects no mirrors
    m0 = -fastFloor(So / d_dir.v_sep);
    mE = -fastCeil(dot(E, s_dir.o) / d_dir.v_sep);
  } else {
    // Disregard rays parallel to this mirror direction since they won't intersect any mirrors
    return 0;
  }

  // If the ray intersects no mirrors, bail out
  if(m0 < mE) { return 0; }

  // The total number of mirrors intersected by the path
  return abs(mE - m0);

  // The distance along the path of the first and last mirror intersections
//  const double k0 = (m0 * d_dir.v_sep - So) / dir_dot_o;
//  const double kE = (mE * d_dir.v_sep - So) / dir_dot_o;
//
//  // Bail out if the path doesn't cross a mirror
//  if(kE < k0) { return 0; }


}

/*
 * Ray attenuation functions. The attenuate callback pointer (attenuate_callback) is set to one of these methods
 */
double OrbifoldData::attenuateX333(const Ray& ray, const Intersection& isect) const {
  const double SQRT_3_OVER_2 = 0.8660254037844386;
  const double SQRT_3_OVER_4 = 0.4330127018922193;

//  const Vector2d B = Vector2d(-0.5, -SQRT_3_OVER_4) * scale.x();
  const Vector2d S = Vector2d(ray.Origin.x(), ray.Origin.z()) - direction_info[0].base;
  const Vector2d E = Vector2d(isect.IPoint.x(), isect.IPoint.z()) - direction_info[0].base;

  const StaticDirectionInfo d1 = { Vector2d{0, 1},
                                 Vector2d{1, 0},
                                 { 0, 1, 2 } };

  const StaticDirectionInfo d2 = { Vector2d{-SQRT_3_OVER_2, 0.5},
                                 Vector2d{0.5, SQRT_3_OVER_2},
                                 { 2, 1, 0 } };

  const StaticDirectionInfo d3 = { Vector2d{SQRT_3_OVER_2, 0.5},
                                 Vector2d{-0.5, SQRT_3_OVER_2},
                                 { 0, 1, 2 } };

  unsigned mirrorCount[3] = {0, 0, 0};

  countMirrorsHeterogeneous(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), d1, direction_info[0], mirrorCount);
  countMirrorsHeterogeneous(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), d2, direction_info[0], mirrorCount);
  countMirrorsHeterogeneous(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), d3, direction_info[0], mirrorCount);

  return pow(r1, mirrorCount[0]) * pow(r2, mirrorCount[1]) * pow(r3, mirrorCount[2]);
}

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

double OrbifoldData::attenuateX442(const Ray& ray, const Intersection& isect) const {
  return 0;
}

}
