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
  direction_info[0].v_sep = scale.z();
  direction_info[1].v_sep = scale.x();
  direction_info[0].base = Vector2d(-0.5* scale.x(), -0.5 * scale.z());
}

void OrbifoldData::InitXXOrbifoldData() {
  attenuate_callback = &OrbifoldData::attenuateXX;
  direction_info[0].v_sep = scale.z();
  direction_info[0].base = Vector2d(-0.5* scale.x(), -0.5 * scale.z());
}

void OrbifoldData::InitX632OrbifoldData() {
  attenuate_callback = &OrbifoldData::attenuateX632;

  const double SQRT_3_OVER_2 = 0.8660254037844386;
  const double SQRT_3_OVER_4 = 0.4330127018922193;
  const double SQRT_3_OVER_6 = 0.28867513459481287;

  direction_info[0].h_sep = scale.x();
  direction_info[0].v_sep = SQRT_3_OVER_2 * scale.x();
  direction_info[0].offset_array[0] = 0;
  direction_info[0].offset_array[1] = 1.5 * scale.x();
  direction_info[0].base = Vector2d(0.5 - (1.0/3.0), SQRT_3_OVER_2 - SQRT_3_OVER_6) * scale.x();

  direction_info[1].v_sep = 1.5 * scale.x();
}

void OrbifoldData::InitX442OrbifoldData() {
  const double SQRT_2 = 1.4142135623730951;

  attenuate_callback = &OrbifoldData::attenuateX442;
  direction_info[0].v_sep = scale.x();
  direction_info[1].v_sep = scale.x() * SQRT_2;
  direction_info[0].base = Vector2d(2.0/3.0, -1.0/3.0) * scale.x();
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

void OrbifoldData::countMirrorsHomogeneous(const Vector2d& S, const Vector2d& E,
                        const Vector2d& dir, const StaticDirectionInfo& s_dir,
                        const DynamicDirectionInfo& d_dir, unsigned* mirrors) const {
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
    // Note: Here the sign is inversed so that m0 > mE if the ray path intersects no mirrors
    m0 = -fastFloor(So / d_dir.v_sep);
    mE = -fastCeil(dot(E, s_dir.o) / d_dir.v_sep);
  } else {
    // Disregard rays parallel to this mirror direction since they won't intersect any mirrors
    return;
  }

  // If the ray intersects no mirrors, bail out
  if(m0 > mE) { return; }

  // The index in the mirror array of the first mirror intersected
  const unsigned m0_type = m0 & 1;

  // The total number of mirrors intersected by the path
  const unsigned num_mirrors = abs(mE - m0) + 1;

  // Used to add one in odd cases
  const unsigned mirror_parity = num_mirrors & 1;

  mirrors[s_dir.index_array[m0_type]] += (num_mirrors/2) + mirror_parity;
  mirrors[s_dir.index_array[(m0_type + 1) % 2]] += (num_mirrors/2);
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
    // Note: Here the sign is inversed so that m0 > mE if the ray path intersects no mirrors
    m0 = -fastFloor(So / d_dir.v_sep);
    mE = -fastCeil(dot(E, s_dir.o) / d_dir.v_sep);
  } else {
    // Disregard rays parallel to this mirror direction since they won't intersect any mirrors
    return 0;
  }

  // If the ray intersects no mirrors, bail out
  if(m0 > mE) { return 0; }

  // The number of mirrors intersected
  return abs(mE - m0) + 1;
}


/*
 * Ray attenuation functions. The attenuate callback pointer (attenuate_callback) is set to one of these methods
 */

double OrbifoldData::attenuateX333(const Ray& ray, const Intersection& isect) const {
  const double SQRT_3_OVER_2 = 0.8660254037844386;
  const double SQRT_3_OVER_4 = 0.4330127018922193;

  const Vector2d S = Vector2d(ray.Origin.x(), ray.Origin.z()) - direction_info[0].base;
  const Vector2d E = Vector2d(isect.IPoint.x(), isect.IPoint.z()) - direction_info[0].base;
  const Vector2d dir = Vector2d(ray.Direction.x(), ray.Direction.z());

  const StaticDirectionInfo sdi1 = { Vector2d{0, 1},
                                     Vector2d{1, 0},
                                     { 0, 1, 2 } };

  const StaticDirectionInfo sdi2 = { Vector2d{-SQRT_3_OVER_2, 0.5},
                                     Vector2d{0.5, SQRT_3_OVER_2},
                                     { 2, 1, 0 } };

  const StaticDirectionInfo sdi3 = { Vector2d{SQRT_3_OVER_2, 0.5},
                                     Vector2d{-0.5, SQRT_3_OVER_2},
                                     { 0, 1, 2 } };

  unsigned mirrorCount[3] = {0, 0, 0};

  countMirrorsHeterogeneous(S, E, dir, sdi1, direction_info[0], mirrorCount);
  countMirrorsHeterogeneous(S, E, dir, sdi2, direction_info[0], mirrorCount);
  countMirrorsHeterogeneous(S, E, dir, sdi3, direction_info[0], mirrorCount);

  return pow(r1, mirrorCount[0]) * pow(r2, mirrorCount[1]) * pow(r3, mirrorCount[2]);
}

double OrbifoldData::attenuateX2222(const Ray& ray, const Intersection& isect) const {
  const Vector2d S = Vector2d(ray.Origin.x(), ray.Origin.z()) - direction_info[0].base;
  const Vector2d E = Vector2d(isect.IPoint.x(), isect.IPoint.z()) - direction_info[0].base;
  const Vector2d dir = Vector2d(ray.Direction.x(), ray.Direction.z());

  const StaticDirectionInfo sdi1 = { Vector2d{0, 1},
                                     Vector2d{1, 0},
                                     { 0, 1, 0 } };

  const StaticDirectionInfo sdi2 = { Vector2d{1, 0},
                                     Vector2d{0, 1},
                                     { 2, 3, 0 } };

  unsigned mirrorCount[4] = {0, 0, 0, 0};

  countMirrorsHomogeneous(S, E, dir, sdi1, direction_info[0], mirrorCount);
  countMirrorsHomogeneous(S, E, dir, sdi2, direction_info[1], mirrorCount);

  return pow(r1, mirrorCount[0]) * pow(r2, mirrorCount[1]) * pow(r3, mirrorCount[2]) * pow(r4, mirrorCount[3]);
}

double OrbifoldData::attenuateXX(const Ray& ray, const Intersection& isect) const {
  const Vector2d S = Vector2d(ray.Origin.x(), ray.Origin.z()) - direction_info[0].base;
  const Vector2d E = Vector2d(isect.IPoint.x(), isect.IPoint.z()) - direction_info[0].base;
  const Vector2d dir = Vector2d(ray.Direction.x(), ray.Direction.z());

  const StaticDirectionInfo sdi = { Vector2d{0, 1},
                                    Vector2d{1, 0},
                                    { 0, 1, 0 } };

  unsigned mirrorCount[2] = {0, 0};

  countMirrorsHomogeneous(S, E, dir, sdi, direction_info[0], mirrorCount);

  return pow(r1, mirrorCount[0]) * pow(r2, mirrorCount[1]);
}

double OrbifoldData::attenuateX632(const Ray& ray, const Intersection& isect) const {
  const double SQRT_3_OVER_2 = 0.8660254037844386;
  const double SQRT_3_OVER_4 = 0.4330127018922193;

  const Vector2d S = Vector2d(ray.Origin.x(), ray.Origin.z()) - direction_info[0].base;
  const Vector2d E = Vector2d(isect.IPoint.x(), isect.IPoint.z()) - direction_info[0].base;
  const Vector2d dir = Vector2d(ray.Direction.x(), ray.Direction.z());

  const StaticDirectionInfo sdi1 = { Vector2d{0, 1},
                                     Vector2d{1, 0},
                                     { 2, 0, 2 } };

  // Homogeneous
  const StaticDirectionInfo sdi2 = { Vector2d{-0.5, SQRT_3_OVER_2},
                                     Vector2d{SQRT_3_OVER_2, 0.5},
                                     { 1, 1, 1 } };

  const StaticDirectionInfo sdi3 = { Vector2d{-SQRT_3_OVER_2, 0.5},
                                     Vector2d{0.5, SQRT_3_OVER_2},
                                     { 2, 0, 2 } };

  // Homogeneous
  const StaticDirectionInfo sdi4 = { Vector2d{-1, 0},
                                     Vector2d{0, 1},
                                     { 1, 1, 1 } };


  const StaticDirectionInfo sdi5 = { Vector2d{SQRT_3_OVER_2, 0.5},
                                     Vector2d{0.5, -SQRT_3_OVER_2},
                                     { 2, 0, 2 } };

  // Homogeneous
  const StaticDirectionInfo sdi6 = { Vector2d{0.5, SQRT_3_OVER_2},
                                     Vector2d{SQRT_3_OVER_2, -0.5},
                                     { 1, 1, 1 } };

  unsigned mirror_count[3] = {0, 0, 0};

  mirror_count[1] = countMirrorsHomogeneous(S, E, dir, sdi2, direction_info[1]) +
                    countMirrorsHomogeneous(S, E, dir, sdi4, direction_info[1]) +
                    countMirrorsHomogeneous(S, E, dir, sdi6, direction_info[1]);

  countMirrorsHeterogeneous(S, E, dir, sdi1, direction_info[0], mirror_count);
  countMirrorsHeterogeneous(S, E, dir, sdi3, direction_info[0], mirror_count);
  countMirrorsHeterogeneous(S, E, dir, sdi5, direction_info[0], mirror_count);

  return pow(r1, mirror_count[0]) * pow(r2, mirror_count[1]) * pow(r3, mirror_count[2]);
}


double OrbifoldData::attenuateX442(const Ray& ray, const Intersection& isect) const {
  const Vector2d S = Vector2d(ray.Origin.x(), ray.Origin.z()) - direction_info[0].base;
  const Vector2d E = Vector2d(isect.IPoint.x(), isect.IPoint.z()) - direction_info[0].base;
  const Vector2d dir = Vector2d(ray.Direction.x(), ray.Direction.z());

  const double ONE_OVER_SQRT_2 = 0.7071067811865475;

  const StaticDirectionInfo sdi1 = { Vector2d{0, 1},
                                     Vector2d{1, 0},
                                     { 0, 1, 0 } };

  const StaticDirectionInfo sdi2 = { Vector2d{1, 0},
                                     Vector2d{0, 1},
                                     { 1, 0, 0 } };

  const StaticDirectionInfo sdi3 = { Vector2d{ONE_OVER_SQRT_2, ONE_OVER_SQRT_2},
                                     Vector2d{-ONE_OVER_SQRT_2, ONE_OVER_SQRT_2},
                                     { 2, 0, 0 } };

  const StaticDirectionInfo sdi4 = { Vector2d{-ONE_OVER_SQRT_2, ONE_OVER_SQRT_2},
                                     Vector2d{ONE_OVER_SQRT_2, ONE_OVER_SQRT_2},
                                     { 2, 0, 0 } };

  unsigned mirror_count[3] = {0, 0, 0};

  countMirrorsHomogeneous(S, E, dir, sdi1, direction_info[0], mirror_count);
  countMirrorsHomogeneous(S, E, dir, sdi2, direction_info[0], mirror_count);
  mirror_count[2] += countMirrorsHomogeneous(S, E, dir, sdi3, direction_info[1]);
  mirror_count[2] += countMirrorsHomogeneous(S, E, dir, sdi4, direction_info[1]);

  return pow(r1, mirror_count[0]) * pow(r2, mirror_count[1]) * pow(r3, mirror_count[2]);
}

}
