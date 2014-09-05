# Pointwise V17.2 Journal file - Fri Sep  5 19:57:48 2014

package require PWI_Glyph 2.17.2

pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}

pw::Application clearModified

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint {-13.6679210966552 -0.42578057424073 -6.39767373336465}
  $_TMP(PW_1) addPoint {-19.1100902032015 0.366385188412155 -9.28319628548404}
  $_TMP(PW_1) addPoint {-26.9121675812476 -0.311124704419976 -12.7824226902626}
  $_TMP(PW_1) addPoint {-32.5545006487571 -1.56388098389656 -15.0447905619338}
  $_TMP(PW_1) addPoint {-36.7995405094373 -5.14254553567917 -15.8199677990177}
  $_TMP(PW_1) setSlope Akima
  set _TMP(curve_1) [pw::Curve create]
  $_TMP(curve_1) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create Curve}

unset _TMP(curve_1)
set _TMP(mode_2) [pw::Application begin Create]
  set _TMP(PW_2) [pw::SegmentSpline create]
  set _DB(1) [pw::DatabaseEntity getByName "curve-1"]
  $_TMP(PW_2) addPoint [list 1 0 $_DB(1)]
  $_TMP(PW_2) addPoint {-31.3953143230786 -5.53182204125138 -13.0942859457522}
  $_TMP(PW_2) addPoint {-25.5516420345817 -5.14268044867222 -10.4318038667855}
  $_TMP(PW_2) addPoint {-20.0881174879519 -2.68633569074124 -8.67829956064628}
  $_TMP(PW_2) addPoint [list 0 0 $_DB(1)]
  $_TMP(PW_2) setSlope Akima
  set _TMP(curve_1) [pw::Curve create]
  $_TMP(curve_1) addSegment $_TMP(PW_2)
  unset _TMP(PW_2)
$_TMP(mode_2) end
unset _TMP(mode_2)
pw::Application markUndoLevel {Create Curve}

unset _TMP(curve_1)
set _TMP(mode_3) [pw::Application begin Create]
  set _TMP(PW_3) [pw::SegmentSpline create]
  $_TMP(PW_3) addPoint {-16.2560105291966 -13.1419517671157 -3.16616714081344}
  $_TMP(PW_3) addPoint {-11.5966676064922 -10.1396849906219 -1.98985077611253}
  $_TMP(PW_3) addPoint {-5.20801419534629 -8.35509318230117 0.443024152700044}
  $_TMP(PW_3) addPoint {1.84774283906344 -8.20941881831319 3.77174350033723}
  $_TMP(PW_3) addPoint {5.52644613701109 -9.27471500285603 5.90854613355536}
  $_TMP(PW_3) setSlope Akima
  set _TMP(curve_1) [pw::Curve create]
  $_TMP(curve_1) addSegment $_TMP(PW_3)
  unset _TMP(PW_3)
$_TMP(mode_3) end
unset _TMP(mode_3)
pw::Application markUndoLevel {Create Curve}

unset _TMP(curve_1)
set _TMP(mode_4) [pw::Application begin Create]
  set _TMP(PW_4) [pw::SegmentSpline create]
  set _DB(2) [pw::DatabaseEntity getByName "curve-3"]
  $_TMP(PW_4) addPoint [list 0 0 $_DB(2)]
  $_TMP(PW_4) addPoint {-9.84365115791878 -13.9259228962649 0.181227317506465}
  $_TMP(PW_4) addPoint {-3.78250888836354 -13.1815696337657 2.82298444514737}
  $_TMP(PW_4) addPoint {1.84441195906608 -11.1859503056515 4.81676225972823}
  $_TMP(PW_4) addPoint [list 1 0 $_DB(2)]
  $_TMP(PW_4) setSlope Akima
  set _TMP(curve_1) [pw::Curve create]
  $_TMP(curve_1) addSegment $_TMP(PW_4)
  unset _TMP(PW_4)
$_TMP(mode_4) end
unset _TMP(mode_4)
pw::Application markUndoLevel {Create Curve}

unset _TMP(curve_1)
set _TMP(mode_5) [pw::Application begin Create]
  set _TMP(PW_5) [pw::SegmentSpline create]
  $_TMP(PW_5) addPoint {9.9865385388185 -18.0578225374557 11.1334133856151}
  $_TMP(PW_5) addPoint {12.9898444734503 -15.6585013888335 11.7284421275916}
  $_TMP(PW_5) addPoint {15.546585825916 -15.1204351785573 12.7640100591936}
  $_TMP(PW_5) addPoint {20.0166708677155 -14.973948251142 14.8538210330803}
  $_TMP(PW_5) addPoint {23.2540342635097 -15.8708152102644 16.7199837685651}
  $_TMP(PW_5) setSlope Akima
  set _TMP(curve_1) [pw::Curve create]
  $_TMP(curve_1) addSegment $_TMP(PW_5)
  unset _TMP(PW_5)
$_TMP(mode_5) end
unset _TMP(mode_5)
pw::Application markUndoLevel {Create Curve}

unset _TMP(curve_1)
set _TMP(mode_6) [pw::Application begin Create]
  set _TMP(PW_6) [pw::SegmentSpline create]
  set _DB(3) [pw::DatabaseEntity getByName "curve-5"]
  $_TMP(PW_6) addPoint [list 1 0 $_DB(3)]
  $_TMP(PW_6) addPoint {20.6631375908869 -17.8996868106036 16.1922544821735}
  $_TMP(PW_6) addPoint {17.4295617951232 -18.8965144986249 14.993771090753}
  $_TMP(PW_6) addPoint {13.4278534956064 -18.8034621015923 13.0441007847159}
  $_TMP(PW_6) addPoint [list 0 0 $_DB(3)]
  $_TMP(PW_6) setSlope Akima
  set _TMP(curve_1) [pw::Curve create]
  $_TMP(curve_1) addSegment $_TMP(PW_6)
  unset _TMP(PW_6)
$_TMP(mode_6) end
unset _TMP(mode_6)
pw::Application markUndoLevel {Create Curve}

unset _TMP(curve_1)
set _TMP(mode_7) [pw::Application begin Create]
  set _TMP(PW_7) [pw::SegmentSpline create]
  $_TMP(PW_7) delete
  unset _TMP(PW_7)
$_TMP(mode_7) abort
unset _TMP(mode_7)
set _TMP(mode_8) [pw::Application begin Create]
  set _TMP(PW_8) [pw::Surface create]
  $_TMP(PW_8) interpolate -orient Best $_DB(2) $_DB(1)
  $_TMP(PW_8) interpolate -orient Same $_DB(2) $_DB(1)
  $_TMP(PW_8) interpolate -orient Opposite $_DB(2) $_DB(1)
$_TMP(mode_8) end
unset _TMP(mode_8)
unset _TMP(PW_8)
pw::Application markUndoLevel {Interpolate}

set _DB(4) [pw::DatabaseEntity getByName "curve-4"]
set _DB(5) [pw::DatabaseEntity getByName "curve-2"]
set _TMP(mode_9) [pw::Application begin Create]
  set _TMP(PW_9) [pw::Surface create]
  $_TMP(PW_9) interpolate -orient Opposite $_DB(4) $_DB(5)
  $_TMP(PW_9) interpolate -orient Same $_DB(4) $_DB(5)
$_TMP(mode_9) end
unset _TMP(mode_9)
unset _TMP(PW_9)
pw::Application markUndoLevel {Interpolate}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_10) [pw::Surface create]
  $_TMP(PW_10) interpolate -orient Same $_DB(2) $_DB(3)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(PW_10)
pw::Application markUndoLevel {Interpolate}

set _DB(6) [pw::DatabaseEntity getByName "curve-6"]
set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_11) [pw::Surface create]
  $_TMP(PW_11) interpolate -orient Same $_DB(6) $_DB(4)
  $_TMP(PW_11) interpolate -orient Opposite $_DB(6) $_DB(4)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(PW_11)
pw::Application markUndoLevel {Interpolate}

set _DB(7) [pw::DatabaseEntity getByName "surface-3"]
set _TMP(boundary_1) [$_DB(7) getBoundary 1]
set _TMP(PW_12) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_1)]]
unset _TMP(unused)
unset _TMP(PW_12)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_1)
set _DB(8) [pw::DatabaseEntity getByName "surface-1"]
set _TMP(boundary_2) [$_DB(8) getBoundary 2]
set _TMP(PW_13) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_2)]]
unset _TMP(unused)
unset _TMP(PW_13)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_2)
set _TMP(boundary_3) [$_DB(8) getBoundary 3]
set _TMP(PW_14) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_3)]]
unset _TMP(unused)
unset _TMP(PW_14)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_3)
set _TMP(boundary_4) [$_DB(8) getBoundary 4]
set _TMP(PW_15) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_4)]]
unset _TMP(unused)
unset _TMP(PW_15)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_4)
set _TMP(boundary_5) [$_DB(7) getBoundary 3]
set _TMP(PW_16) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_5)]]
unset _TMP(unused)
unset _TMP(PW_16)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_5)
set _TMP(boundary_6) [$_DB(7) getBoundary 2]
set _TMP(PW_17) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_6)]]
unset _TMP(unused)
unset _TMP(PW_17)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_6)
set _TMP(boundary_7) [$_DB(7) getBoundary 4]
set _TMP(PW_18) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_7)]]
unset _TMP(unused)
unset _TMP(PW_18)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_7)
set _DB(9) [pw::DatabaseEntity getByName "surface-4"]
set _TMP(boundary_8) [$_DB(9) getBoundary 1]
set _TMP(PW_19) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_8)]]
unset _TMP(unused)
unset _TMP(PW_19)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_8)
set _DB(10) [pw::DatabaseEntity getByName "surface-2"]
set _TMP(boundary_9) [$_DB(10) getBoundary 1]
set _TMP(PW_20) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_9)]]
unset _TMP(unused)
unset _TMP(PW_20)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_9)
set _TMP(boundary_10) [$_DB(10) getBoundary 3]
set _TMP(PW_21) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_TMP(boundary_10)]]
unset _TMP(unused)
unset _TMP(PW_21)
pw::Application markUndoLevel {Connectors On DB Entities}

unset _TMP(boundary_10)
