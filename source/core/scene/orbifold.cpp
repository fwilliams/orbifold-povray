#include "core/scene/orbifold.h"

namespace pov {

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
  info.mirrorDirs[2].base = Vector2d(0.5, -SQRT_3_OVER_4) * info.scale.x();
  info.mirrorDirs[2].o = Vector2d(SQRT_3_OVER_2, 0.5);
  info.mirrorDirs[2].m = Vector2d(-0.5, SQRT_3_OVER_2);

  info.mirrorDirs[2].numOffsets = 2;

  info.mirrorDirs[2].offsetArray[0] = 0;
  info.mirrorDirs[2].offsetArray[1] = 1.5 * info.scale.x();

  info.mirrorDirs[2].numMirrors = 3;
  info.mirrorDirs[2].indexArray[0] = 1;
  info.mirrorDirs[2].indexArray[1] = 0;
  info.mirrorDirs[2].indexArray[2] = 2;

  info.mirrorDirs[2].v_sep = SQRT_3_OVER_2 * info.scale.x();
  info.mirrorDirs[2].h_sep = info.scale.x();
}

}
