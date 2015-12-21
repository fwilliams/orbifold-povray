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
    // Disregard rays parallel to this mirror direction since they won't intersect any mirrors
    return;
  }

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

  return pow(r1, (double)n0) *
         pow(r2, (double)n1) *
         pow(r3, (double)n2) *
         pow(r4, (double)n3);
}

double OrbifoldData::attenuateXX(const Ray& ray, const Intersection& isect) const {
  typedef GenericVector2d<POV_INT64> IntVector2d;

  const Vector3d rayStart = ray.Origin / scale;
  const Vector3d rayEnd = isect.IPoint / scale;

  const POV_INT64 latticeStart = floor(rayStart.z() + 0.5);
  const POV_INT64 latticeEnd = floor(rayEnd.z() + 0.5);
  const POV_INT64 latticeDist = abs(latticeStart - latticeEnd);

  // Bail out if the ray doesn't attenuate
  if(latticeDist == 0) {
    return 1.0;
  }

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

  return pow(r1, (double)n0) *
         pow(r2, (double)n1);
}

double OrbifoldData::attenuateX632(const Ray& ray, const Intersection& isect) const {
  return 0;
}

Vector3d roundHexCoordinate(const Vector3d& hexCoord)
{
  Vector3d r = vround(hexCoord);
  Vector3d diff = vabs(hexCoord - r);

  if(diff.x() > diff.y() && diff.x() > diff.z()) {
    r.x() = -r.y() - r.z();
  } else if (diff.y() > diff.z()) {
    r.y() = -r.x() - r.z();
  } else {
    r.z() = -r.x() - r.y();
  }

  return r;
}

Vector3d triCoordForPoint(const Vector3d& point, unsigned& type, Vector2d& base) {
  const double SQRT_3 = 1.7320508075688772;
  const double HALF_SQRT_3 = 0.8660254037844386;

  const Vector3d scale = Vector3d(3.0/SQRT_3, 1.0, SQRT_3);
  const Vector3d tPoint = point / scale;
#define __HEX_COORD__(pt) Vector3d(-(pt).z() - ((pt).x()/SQRT_3), 2.0 * (pt).x() / SQRT_3, (pt).z() - ((pt).x()/SQRT_3))
  const Vector3d hexCoord = roundHexCoordinate(__HEX_COORD__(tPoint));
#undef __HEX_COORD__
  const Vector2d hexXYPos(hexCoord.y() * HALF_SQRT_3, hexCoord.z() + hexCoord.y() * 0.5);
  base = hexXYPos * Vector2d(scale.x(), scale.z());

  const Vector2d omega = Vector2d(-0.5, HALF_SQRT_3);
  const Vector2d omegaBar = Vector2d(-0.5, -HALF_SQRT_3);

  const Vector3d baseTriCoord = hexCoord.y() * Vector3d(-1, 1, -2) - hexCoord.x() * Vector3d(2, 1, 1);
  Vector3d retTriCoord = baseTriCoord;

  if(tPoint.z() > hexXYPos.y()) {
    retTriCoord.x() += 1;
  }

  if((tPoint.x() - hexXYPos.x())*omegaBar.y() - (tPoint.z() - hexXYPos.y())*omegaBar.x() > 0) {
    retTriCoord.z() += 1;
  }

  if((tPoint.x() - hexXYPos.x())*omega.y() - (tPoint.z() - hexXYPos.y())*omega.x() > 0) {
    retTriCoord.y() += 1;
  }

  const Vector3d bt = retTriCoord - baseTriCoord;
  if(bt.x() == 0 && bt.y() == 0 && bt.z() == 0) { type = 0; }
  else if(bt.x() == 0 && bt.y() == 1 && bt.z() == 0) { type = 1; }
  else if(bt.x() == 1 && bt.y() == 1 && bt.z() == 0) { type = 2; }
  else if(bt.x() == 1 && bt.y() == 1 && bt.z() == 1) { type = 3; }
  else if(bt.x() == 1 && bt.y() == 0 && bt.z() == 1) { type = 4; }
  else if(bt.x() == 0 && bt.y() == 0 && bt.z() == 1) { type = 5; }

  return retTriCoord;
}

bool intersectLines(const Vector2d& P1, const Vector2d& v1, const Vector2d& P2, const Vector2d& v2, double& alpha1, double& alpha2) {
  const double num = (v2.y()*(P1.x()-P2.x()) - v2.x()*(P1.y()-P2.y()));
  const double denom = (v2.x()*v1.y() - v2.y()*v1.x());

  if(denom == 0) {
    alpha1 = alpha2 = INFINITY;
    return num == 0;
  }

  alpha1 = num/denom;
  if(v2.x() != 0 & v2.y() != 0) {
    alpha2 = max((P1.x() - P2.x() + v1.x()*alpha1)/v2.x(), (P1.y() - P2.y() + v1.y()*alpha1)/v2.y());
  } else if(v2.x() != 0) {
    alpha2 = (P1.x() - P2.x() + v1.x()*alpha1)/v2.x();
  } else {
    alpha2 = (P1.y() - P2.y() + v1.y()*alpha1)/v2.y();
  }

  return true;

}

double att333(const Ray& ray, const Intersection& isect, const OrbifoldData& odata) {
    const double SQRT_3 = 1.7320508075688772;
    const double HALF_SQRT_3 = 0.8660254037844386;

    struct ltri {
      Vector2d corners[3];
      Vector3d coorddiffs[3];
      Vector2d posdiffs[3];
      Vector3d mAdd[3];
      unsigned next_type[3];
    };

    const unsigned next_m[3] = {0, 2, 1}; // The next mirror index based on the current index
    // Lookup table for mirror intersection information
    // Given a mirror and tile, the jump table tells us which tile is adjacent and how to update the position
    // and triangle coordinates as we march along a path
    const ltri jmptable[6] = {
        ltri{ // 0 - (0, 0, 0)
          {Vector2d(0, 0), Vector2d(-0.5, -HALF_SQRT_3), Vector2d(0.5, -HALF_SQRT_3)},
          {Vector3d(0, 0, 1), Vector3d(-1, 0, 0), Vector3d(0, 1, 0)},
          {Vector2d(0, 0), Vector2d(0, -SQRT_3), Vector2d(0, 0)},
          {Vector3d(0, 0, 1), Vector3d(1, 0, 0), Vector3d(0, 1, 0)}, // Blue, Orange, Green
          {5, 3, 1},
        },
        ltri{ // 1 - (0, 1, 0)
          {Vector2d(0, 0), Vector2d(0.5, -HALF_SQRT_3),  Vector2d(1, 0)},
          {Vector3d(0, -1, 0), Vector3d(0, 0, -1), Vector3d(1, 0, 0)},
          {Vector2d(0, 0), Vector2d(1.5, -HALF_SQRT_3), Vector2d(0, 0)},
          {Vector3d(0, 1, 0), Vector3d(1, 0, 0), Vector3d(0, 0, 1)},
          {0, 4, 2}
        },
        ltri{ // 2 - (1, 1, 0)
          {Vector2d(0, 0), Vector2d(1, 0), Vector2d(0.5, HALF_SQRT_3)},
          {Vector3d(-1, 0, 0), Vector3d(0, 1, 0), Vector3d(0, 0, 1)},
          {Vector2d(0, 0), Vector2d(1.5, HALF_SQRT_3), Vector2d(0, 0)},
          {Vector3d(0, 0, 1), Vector3d(1, 0, 0), Vector3d(0, 1, 0)},
          {1, 5, 3},
        },
        ltri{ // 3 - (1, 1, 1)
          {Vector2d(0, 0), Vector2d(0.5, HALF_SQRT_3), Vector2d(-0.5, HALF_SQRT_3)},
          {Vector3d(0, 0, -1), Vector3d(1, 0, 0), Vector3d(0, -1, 0)},
          {Vector2d(0, 0), Vector2d(0, SQRT_3), Vector2d(0, 0)},
          {Vector3d(0, 1, 0), Vector3d(1, 0, 0), Vector3d(0, 0, 1)},
          {2, 0, 4},
        },
        ltri{ // 4 - (1, 0, 1)
          {Vector2d(0, 0), Vector2d(-0.5, HALF_SQRT_3), Vector2d(-1, 0)},
          {Vector3d(0, 1, 0), Vector3d(0, 0, 1), Vector3d(-1, 0, 0)},
          {Vector2d(0, 0), Vector2d(-1.5, HALF_SQRT_3), Vector2d(0, 0)},
          {Vector3d(0, 0, 1), Vector3d(1, 0, 0), Vector3d(0, 1, 0)},
          {3, 1, 5},
        },
        ltri{ // 5 - (0, 0, 1)
          {Vector2d(0, 0), Vector2d(-1, 0), Vector2d(-0.5, -HALF_SQRT_3)},
          {Vector3d(1, 0, 0), Vector3d(0, -1, 0), Vector3d(0, 0, -1)},
          {Vector2d(0, 0), Vector2d(-1.5, -HALF_SQRT_3), Vector2d(0, 0)},
          {Vector3d(0, 1, 0), Vector3d(1, 0, 0), Vector3d(0, 0, 1)},
          {4, 2, 0},
        }};

    // The mirror intersection counters. Count the number of each type of mirror a path intersects
    unsigned m1 = 0, m2 = 0, m3 = 0;

    // Make the center at the middle of the fundamental domain, and normalize the ray to unit length mirrors
    const Vector3d ctrTx(0, 0, -SQRT_3/4);
    const Vector3d scale = Vector3d(odata.scale.x(), odata.scale.y(), odata.scale.x());
    const Vector3d rayStart = ctrTx + ray.Origin / scale;
    const Vector3d rayEnd = ctrTx + isect.IPoint / scale;

    //throw runtime_error(string("raystart = ") + vec3dToString(rayStart) + string(" rayend = ") + vec3dToString(rayEnd) + string("\n"));

    unsigned type = 0; // Index of the current triangle type in the jump table
    unsigned m = 0; // Index of the base vertex of the edge under consideration
    Vector2d basePos(0, 0); // Offset of the mirror under consideration
    const Vector3d triEnd = triCoordForPoint(rayEnd, type, basePos); // Triangle coordinates of the end of the path
    const Vector3d triStart = triCoordForPoint(rayStart, type, basePos); // Triangle coordinates of the beginning of the path
    const Vector2d d(rayEnd.x() - rayStart.x(), rayEnd.z() - rayStart.z()); // The direction of the ray
    const Vector2d P(rayStart.x(), rayStart.z()); // The start position of the ray
    Vector3d current =  triStart; // Triangle coordinate of the current triangle



    // FIXME: Heisenbug when running radiosity requires bailout to breakout of an infinite loop
    unsigned bailout = 0;
    const unsigned MAX_BAILOUT = 100000000;

    // While we haven't reached the triangle coordinate for the end of the path
    while((bailout < MAX_BAILOUT) && !(current.x() == triEnd.x() && current.y() == triEnd.y() && current.z() == triEnd.z())) {
      for(unsigned i = 0; i < 3; i++) { // Check each mirror for intersection
        // Note that the jump table is set up so m never corresponds to the mirror we just came from
        const Vector2d A = jmptable[type].corners[m] + basePos;
        const Vector2d B = jmptable[type].corners[(m+1)%3] + basePos;
        const Vector2d v = B-A;

        // Find the intersection of the mirror and the path
        double alpha, beta;
        bool isect = intersectLines(A, v, P, d, alpha, beta);

        retry_alpha:
        if(!isect) {
          continue;
        } else if(beta > 0.0) { // If the intersection lies on the path in front of the start point
          if(alpha > 0.0 && alpha < 1.0) { // If the path intersects the mirror not at one of the endpoints
            // Use the jump table to update the intersection counters and move to the next tile
            unsigned mm = m;
            m = next_m[mm];
            current += jmptable[type].coorddiffs[mm];
            basePos += jmptable[type].posdiffs[mm];
            type = jmptable[type].next_type[mm];
            m1 += jmptable[type].mAdd[mm].x();
            m2 += jmptable[type].mAdd[mm].y();
            m3 += jmptable[type].mAdd[mm].z();
            break;
          } else if(alpha == 0.0 || alpha == 1.0) { // The path intersects one of the endpoints of the mirror
            // Perturb the ray direction a little and try again
            isect = intersectLines(A, v, P, d - v.x()*0.001/beta, alpha, beta);
            goto retry_alpha;
          }

        } else if(alpha == INFINITY) { // The path and the share more than one point in common
          throw runtime_error("TODO: handle intersection of the same line!");
        }
        m = (m + 1) % 3;
      }
      bailout += 1;
    }

    return pow(odata.r1, m1) * pow(odata.r2, m2) * pow(odata.r3, m3);
  }

double OrbifoldData::attenuateX333(const Ray& ray, const Intersection& isect) const {
    const Vector2d B = mirrorDirs[0].base;
    const Vector2d S = Vector2d(ray.Origin.x(), ray.Origin.z()) - B;
    const Vector2d E = Vector2d(isect.IPoint.x(), isect.IPoint.z()) - B;
//    Vector2d dir(ray.Direction.x(), ray.Direction.z());
//    dir.normalize();

    unsigned mirrorCount[3] = {0, 0, 0};

    countMirrorsHetero(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), mirrorDirs[0], mirrorCount);
    countMirrorsHetero(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), mirrorDirs[1], mirrorCount);
    countMirrorsHetero(S, E, Vector2d(ray.Direction.x(), ray.Direction.z()), mirrorDirs[2], mirrorCount);

    return pow(r1, mirrorCount[0]) *
           pow(r2, mirrorCount[1]) *
           pow(r3, mirrorCount[2]);
//  return att333(ray, isect, *this);
}

double OrbifoldData::attenuateX442(const Ray& ray, const Intersection& isect) const {
  return 0;
}

}
