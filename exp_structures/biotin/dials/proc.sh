#!/bin/bash

cat > mask.phil <<+
untrusted {
  panel = 0
  polygon = 2 1006 982 1001 999 985 1023 973 1049 976 1069 991 1078 1008 1079 \
            1039 1070 1057 1049 1074 1021 1073 1004 1069 989 1053 981 1047 7 \
            1057 2 1006
}
+
dials.import biotin_xtal1/biotin_xtal1_*.mrc\
  fast_slow_beam_centre=1035.7,1029.5\
  panel.pedestal=-20\
  geometry.goniometer.axis=1,0,0
dials.generate_mask imported.expt mask.phil
dials.apply_mask imported.expt mask=pixels.mask
dials.find_spots masked.expt\
  threshold.algorithm=radial_profile blur=wide
dials.index masked.expt strong.refl\
  detector.fix=distance space_group=P212121
dials.refine indexed.expt indexed.refl\
  detector.fix=distance unit_cell.force_static=True
dials.integrate refined.expt refined.refl
