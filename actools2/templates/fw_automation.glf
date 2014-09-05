set _TMP(mode_2) [pw::Application begin Create]
  set _TMP(PW_2) [pw::SegmentSpline create]
  set _CN(1) [pw::GridEntity getByName "con-5"]
  set _CN(2) [pw::GridEntity getByName "con-7"]
  set _CN(3) [pw::GridEntity getByName "con-8"]
  $_TMP(PW_2) addPoint [$_CN(1) getPosition -arc 0]
  $_TMP(PW_2) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(1) getPosition -arc 0]] {0.05 0 0.05}]
  set _TMP(con_1) [pw::Connector create]
  $_TMP(con_1) addSegment $_TMP(PW_2)
  unset _TMP(PW_2)
  $_TMP(con_1) calculateDimension
$_TMP(mode_2) end
unset _TMP(mode_2)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_3) [pw::Application begin Create]
  set _TMP(PW_3) [pw::SegmentSpline create]
  set _CN(4) [pw::GridEntity getByName "con-6"]
  $_TMP(PW_3) addPoint [$_CN(1) getPosition -arc 1]
  $_TMP(PW_3) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(1) getPosition -arc 1]] {0 0 0.05}]
  unset _TMP(con_1)
  set _TMP(con_2) [pw::Connector create]
  $_TMP(con_2) addSegment $_TMP(PW_3)
  unset _TMP(PW_3)
  $_TMP(con_2) calculateDimension
$_TMP(mode_3) end
unset _TMP(mode_3)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_4) [pw::Application begin Create]
  set _TMP(PW_4) [pw::SegmentSpline create]
  set _CN(5) [pw::GridEntity getByName "con-15"]
  set _CN(6) [pw::GridEntity getByName "con-16"]
  $_TMP(PW_4) addPoint [$_CN(5) getPosition -arc 1]
  $_TMP(PW_4) addPoint [$_CN(6) getPosition -arc 1]
  unset _TMP(con_2)
  set _TMP(con_3) [pw::Connector create]
  $_TMP(con_3) addSegment $_TMP(PW_4)
  unset _TMP(PW_4)
  $_TMP(con_3) calculateDimension
$_TMP(mode_4) end
unset _TMP(mode_4)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_5) [pw::Application begin Create]
  set _TMP(PW_5) [pw::SegmentSpline create]
  set _CN(7) [pw::GridEntity getByName "con-17"]
  set _CN(8) [pw::GridEntity getByName "con-2"]
  set _CN(9) [pw::GridEntity getByName "con-10"]
  set _CN(10) [pw::GridEntity getByName "con-3"]
  $_TMP(PW_5) addPoint [$_CN(8) getPosition -arc 1]
  $_TMP(PW_5) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(8) getPosition -arc 1]] {120 0 0}]
  unset _TMP(con_3)
  set _TMP(con_4) [pw::Connector create]
  $_TMP(con_4) addSegment $_TMP(PW_5)
  unset _TMP(PW_5)
  $_TMP(con_4) calculateDimension
$_TMP(mode_5) end
unset _TMP(mode_5)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_6) [pw::Application begin Create]
  set _TMP(PW_6) [pw::SegmentSpline create]
  set _CN(11) [pw::GridEntity getByName "con-18"]
  $_TMP(PW_6) addPoint [$_CN(11) getPosition -arc 1]
  $_TMP(PW_6) addPoint [pwu::Vector3 add [$_CN(11) getPosition -arc 1] {0 100 0}]
  unset _TMP(con_4)
  set _TMP(con_5) [pw::Connector create]
  $_TMP(con_5) addSegment $_TMP(PW_6)
  unset _TMP(PW_6)
  $_TMP(con_5) calculateDimension
$_TMP(mode_6) end
unset _TMP(mode_6)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_7) [pw::Application begin Create]
  set _TMP(PW_7) [pw::SegmentSpline create]
  set _CN(12) [pw::GridEntity getByName "con-19"]
  $_TMP(PW_7) addPoint [$_CN(11) getPosition -arc 1]
  $_TMP(PW_7) addPoint [pwu::Vector3 add [$_CN(11) getPosition -arc 1] {0 -100 0}]
  unset _TMP(con_5)
  set _TMP(con_6) [pw::Connector create]
  $_TMP(con_6) addSegment $_TMP(PW_7)
  unset _TMP(PW_7)
  $_TMP(con_6) calculateDimension
$_TMP(mode_7) end
unset _TMP(mode_7)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_8) [pw::Application begin Create]
  set _TMP(PW_8) [pw::SegmentSpline create]
  $_TMP(PW_8) addPoint [$_CN(12) getPosition -arc 1]
  $_TMP(PW_8) addPoint [pwu::Vector3 add [$_CN(12) getPosition -arc 1] {-110 0 0}]
  unset _TMP(con_6)
  set _TMP(con_7) [pw::Connector create]
  $_TMP(con_7) addSegment $_TMP(PW_8)
  unset _TMP(PW_8)
  $_TMP(con_7) calculateDimension
$_TMP(mode_8) end
unset _TMP(mode_8)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_9) [pw::Application begin Create]
  set _TMP(PW_9) [pw::SegmentSpline create]
  set _CN(13) [pw::GridEntity getByName "con-20"]
  $_TMP(PW_9) addPoint [$_CN(13) getPosition -arc 1]
  $_TMP(PW_9) addPoint [pwu::Vector3 add [$_CN(13) getPosition -arc 1] {-110 0 0}]
  unset _TMP(con_7)
  set _TMP(con_8) [pw::Connector create]
  $_TMP(con_8) addSegment $_TMP(PW_9)
  unset _TMP(PW_9)
  $_TMP(con_8) calculateDimension
$_TMP(mode_9) end
unset _TMP(mode_9)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_10) [pw::SegmentSpline create]
  $_TMP(PW_10) delete
  unset _TMP(PW_10)
$_TMP(mode_10) abort
unset _TMP(mode_10)
unset _TMP(con_8)
pw::Display resetView -Z
set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_11) [pw::SegmentCircle create]
  set _CN(14) [pw::GridEntity getByName "con-21"]
  set _CN(15) [pw::GridEntity getByName "con-1"]
  set _CN(16) [pw::GridEntity getByName "con-22"]
  $_TMP(PW_11) addPoint [$_CN(14) getPosition -arc 1]
  $_TMP(PW_11) addPoint [$_CN(16) getPosition -arc 1]
  $_TMP(PW_11) setAngle 180 {0 0 1}
  set _TMP(con_9) [pw::Connector create]
  $_TMP(con_9) addSegment $_TMP(PW_11)
  $_TMP(con_9) calculateDimension
  unset _TMP(PW_11)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create Connector}

unset _TMP(con_9)
set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_12) [pw::SegmentCircle create]
  $_TMP(PW_12) delete
  unset _TMP(PW_12)
$_TMP(mode_10) abort
unset _TMP(mode_10)
set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_13) [pw::SegmentSpline create]
  set _CN(17) [pw::GridEntity getByName "con-9"]
  set _CN(18) [pw::GridEntity getByName "con-4"]
  $_TMP(PW_13) addPoint [$_CN(11) getPosition -arc 1]
  $_TMP(PW_13) addPoint [pwu::Vector3 add [$_CN(11) getPosition -arc 1] {0 0 4.5}]
  set _TMP(con_10) [pw::Connector create]
  $_TMP(con_10) addSegment $_TMP(PW_13)
  unset _TMP(PW_13)
  $_TMP(con_10) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_14) [pw::SegmentSpline create]
  $_TMP(PW_14) addPoint [$_CN(13) getPosition -arc 1]
  $_TMP(PW_14) addPoint [pwu::Vector3 add [$_CN(13) getPosition -arc 1] {0 0 4.5}]
  unset _TMP(con_10)
  set _TMP(con_11) [pw::Connector create]
  $_TMP(con_11) addSegment $_TMP(PW_14)
  unset _TMP(PW_14)
  $_TMP(con_11) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_15) [pw::SegmentSpline create]
  $_TMP(PW_15) addPoint [$_CN(12) getPosition -arc 1]
  $_TMP(PW_15) addPoint [pwu::Vector3 add [$_CN(12) getPosition -arc 1] {0 0 4.5}]
  unset _TMP(con_11)
  set _TMP(con_12) [pw::Connector create]
  $_TMP(con_12) addSegment $_TMP(PW_15)
  unset _TMP(PW_15)
  $_TMP(con_12) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_16) [pw::SegmentSpline create]
  set _CN(19) [pw::GridEntity getByName "con-24"]
  set _CN(20) [pw::GridEntity getByName "con-26"]
  $_TMP(PW_16) addPoint [$_CN(19) getPosition -arc 1]
  $_TMP(PW_16) addPoint [$_CN(20) getPosition -arc 1]
  unset _TMP(con_12)
  set _TMP(con_13) [pw::Connector create]
  $_TMP(con_13) addSegment $_TMP(PW_16)
  unset _TMP(PW_16)
  $_TMP(con_13) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_17) [pw::SegmentSpline create]
  set _CN(21) [pw::GridEntity getByName "con-27"]
  set _CN(22) [pw::GridEntity getByName "con-25"]
  $_TMP(PW_17) addPoint [$_CN(19) getPosition -arc 1]
  $_TMP(PW_17) addPoint [$_CN(22) getPosition -arc 1]
  unset _TMP(con_13)
  set _TMP(con_14) [pw::Connector create]
  $_TMP(con_14) addSegment $_TMP(PW_17)
  unset _TMP(PW_17)
  $_TMP(con_14) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_18) [pw::SegmentSpline create]
  set _CN(23) [pw::GridEntity getByName "con-28"]
  $_TMP(PW_18) delete
  unset _TMP(PW_18)
$_TMP(mode_10) abort
unset _TMP(mode_10)
unset _TMP(con_14)
pw::Application clearClipboard
set _CN(24) [pw::GridEntity getByName "con-23"]
pw::Application setClipboard [list $_CN(11) $_CN(13) $_CN(24) $_CN(16) $_CN(12) $_CN(14)]
pw::Application markUndoLevel {Copy}

set _TMP(mode_10) [pw::Application begin Paste]
  set _TMP(PW_19) [$_TMP(mode_10) getEntities]
  set _TMP(mode_11) [pw::Application begin Modify $_TMP(PW_19)]
    pw::Entity transform [pwu::Transform translation {0 0 90}] [$_TMP(mode_11) getEntities]
  $_TMP(mode_11) end
  unset _TMP(mode_11)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Paste}

unset _TMP(PW_19)
set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_20) [pw::SegmentSpline create]
  set _CN(25) [pw::GridEntity getByName "con-29"]
  set _CN(26) [pw::GridEntity getByName "con-31"]
  set _CN(27) [pw::GridEntity getByName "con-30"]
  $_TMP(PW_20) addPoint [$_CN(19) getPosition -arc 1]
  $_TMP(PW_20) addPoint [$_CN(25) getPosition -arc 1]
  set _TMP(con_15) [pw::Connector create]
  $_TMP(con_15) addSegment $_TMP(PW_20)
  unset _TMP(PW_20)
  $_TMP(con_15) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_21) [pw::SegmentSpline create]
  set _CN(28) [pw::GridEntity getByName "con-35"]
  set _CN(29) [pw::GridEntity getByName "con-32"]
  $_TMP(PW_21) addPoint [$_CN(20) getPosition -arc 1]
  $_TMP(PW_21) addPoint [$_CN(27) getPosition -arc 1]
  unset _TMP(con_15)
  set _TMP(con_16) [pw::Connector create]
  $_TMP(con_16) addSegment $_TMP(PW_21)
  unset _TMP(PW_21)
  $_TMP(con_16) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_22) [pw::SegmentSpline create]
  set _CN(30) [pw::GridEntity getByName "con-36"]
  set _CN(31) [pw::GridEntity getByName "con-33"]
  $_TMP(PW_22) addPoint [$_CN(22) getPosition -arc 1]
  $_TMP(PW_22) addPoint [$_CN(26) getPosition -arc 1]
  unset _TMP(con_16)
  set _TMP(con_17) [pw::Connector create]
  $_TMP(con_17) addSegment $_TMP(PW_22)
  unset _TMP(PW_22)
  $_TMP(con_17) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_23) [pw::SegmentSpline create]
  set _CN(32) [pw::GridEntity getByName "con-37"]
  $_TMP(PW_23) addPoint [$_CN(25) getPosition -arc 0]
  $_TMP(PW_23) addPoint [pwu::Vector3 add [$_CN(25) getPosition -arc 0] {-2 0 0}]
  unset _TMP(con_17)
  set _TMP(con_18) [pw::Connector create]
  $_TMP(con_18) addSegment $_TMP(PW_23)
  unset _TMP(PW_23)
  $_TMP(con_18) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_24) [pw::SegmentSpline create]
  $_TMP(PW_24) delete
  unset _TMP(PW_24)
$_TMP(mode_10) abort
unset _TMP(mode_10)
unset _TMP(con_18)
set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_25) [pw::SegmentSpline create]
  pw::Display resetView -Z
  set _CN(33) [pw::GridEntity getByName "con-38"]
  $_TMP(PW_25) addPoint [$_CN(6) getPosition -arc 1]
  $_TMP(PW_25) addPoint [$_CN(25) getPosition -arc 0]
  set _TMP(con_19) [pw::Connector create]
  $_TMP(con_19) addSegment $_TMP(PW_25)
  unset _TMP(PW_25)
  $_TMP(con_19) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_26) [pw::SegmentSpline create]
  set _CN(34) [pw::GridEntity getByName "con-39"]
  $_TMP(PW_26) addPoint [$_CN(5) getPosition -arc 1]
  $_TMP(PW_26) addPoint [$_CN(33) getPosition -arc 1]
  unset _TMP(con_19)
  set _TMP(con_20) [pw::Connector create]
  $_TMP(con_20) addSegment $_TMP(PW_26)
  unset _TMP(PW_26)
  $_TMP(con_20) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(PW_27) [pw::SegmentSpline create]
  set _CN(35) [pw::GridEntity getByName "con-40"]
  $_TMP(PW_27) delete
  unset _TMP(PW_27)
$_TMP(mode_10) abort
unset _TMP(mode_10)
unset _TMP(con_20)
pw::Display resetView -Z
set _TMP(PW_28) [pw::Collection create]
$_TMP(PW_28) set [list $_CN(9) $_CN(3) $_CN(7) $_CN(1) $_CN(33) $_CN(10) $_CN(15) $_CN(17)]
$_TMP(PW_28) do setDimension 65
$_TMP(PW_28) delete
unset _TMP(PW_28)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_29) [pw::Collection create]
$_TMP(PW_29) set [list $_CN(18) $_CN(8)]
$_TMP(PW_29) do setDimension 30
$_TMP(PW_29) delete
unset _TMP(PW_29)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_30) [pw::Collection create]
$_TMP(PW_30) set [list $_CN(2) $_CN(4)]
$_TMP(PW_30) do setDimension 31
$_TMP(PW_30) delete
unset _TMP(PW_30)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_31) [pw::Collection create]
$_TMP(PW_31) set [list $_CN(5) $_CN(6)]
$_TMP(PW_31) do setDimension 8
$_TMP(PW_31) delete
unset _TMP(PW_31)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_32) [pw::Collection create]
$_TMP(PW_32) set [list $_CN(22) $_CN(19) $_CN(20)]
$_TMP(PW_32) do setDimension 67
$_TMP(PW_32) delete
unset _TMP(PW_32)
pw::Application markUndoLevel {Dimension}

pw::Display resetView -Z
set _TMP(PW_33) [pw::Collection create]
$_TMP(PW_33) set [list $_CN(11) $_CN(16) $_CN(29) $_CN(31) $_CN(25) $_CN(14)]
$_TMP(PW_33) do setDimension 51
$_TMP(PW_33) delete
unset _TMP(PW_33)
pw::Application markUndoLevel {Dimension}

pw::Display resetView -Z
set _TMP(PW_34) [pw::Collection create]
$_TMP(PW_34) set [list $_CN(23) $_CN(27) $_CN(13) $_CN(21) $_CN(12) $_CN(26)]
$_TMP(PW_34) do setDimension 100
$_TMP(PW_34) delete
unset _TMP(PW_34)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_35) [pw::Collection create]
$_TMP(PW_35) set [list $_CN(32) $_CN(30) $_CN(28)]
$_TMP(PW_35) do setDimension 55
$_TMP(PW_35) delete
unset _TMP(PW_35)
pw::Application markUndoLevel {Dimension}

set _CN(36) [pw::GridEntity getByName "con-34"]
set _TMP(PW_36) [pw::Collection create]
$_TMP(PW_36) set [list $_CN(36) $_CN(24)]
$_TMP(PW_36) do setDimension 129
$_TMP(PW_36) delete
unset _TMP(PW_36)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_37) [pw::Collection create]
$_TMP(PW_37) set [list $_CN(34) $_CN(35)]
$_TMP(PW_37) do setDimension 55
$_TMP(PW_37) delete
unset _TMP(PW_37)
pw::Application markUndoLevel {Dimension}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(13) $_CN(12)]]
  set _TMP(PW_38) [$_CN(12) getDistribution 1]
  $_TMP(PW_38) setBeginSpacing 6.9999999999999999e-06
  unset _TMP(PW_38)
  set _TMP(PW_39) [$_CN(13) getDistribution 1]
  $_TMP(PW_39) setBeginSpacing 6.9999999999999999e-06
  unset _TMP(PW_39)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(21) $_CN(23)]]
  set _TMP(PW_40) [$_CN(21) getDistribution 1]
  $_TMP(PW_40) setBeginSpacing 1.5e-06
  unset _TMP(PW_40)
  set _TMP(PW_41) [$_CN(23) getDistribution 1]
  $_TMP(PW_41) setBeginSpacing 1.5e-06
  unset _TMP(PW_41)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(27) $_CN(26)]]
  set _TMP(PW_42) [$_CN(26) getDistribution 1]
  $_TMP(PW_42) setBeginSpacing 0.001
  unset _TMP(PW_42)
  set _TMP(PW_43) [$_CN(27) getDistribution 1]
  $_TMP(PW_43) setBeginSpacing 0.001
  unset _TMP(PW_43)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(25) $_CN(11)]]
  set _TMP(PW_44) [$_CN(11) getDistribution 1]
  $_TMP(PW_44) setBeginSpacing 0.001
  unset _TMP(PW_44)
  set _TMP(PW_45) [$_CN(25) getDistribution 1]
  $_TMP(PW_45) setBeginSpacing 0.001
  unset _TMP(PW_45)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(31) $_CN(14) $_CN(29) $_CN(16)]]
  set _TMP(PW_46) [$_CN(14) getDistribution 1]
  $_TMP(PW_46) setEndSpacing 0.01
  unset _TMP(PW_46)
  set _TMP(PW_47) [$_CN(16) getDistribution 1]
  $_TMP(PW_47) setEndSpacing 0.01
  unset _TMP(PW_47)
  set _TMP(PW_48) [$_CN(31) getDistribution 1]
  $_TMP(PW_48) setEndSpacing 0.01
  unset _TMP(PW_48)
  set _TMP(PW_49) [$_CN(29) getDistribution 1]
  $_TMP(PW_49) setEndSpacing 0.01
  unset _TMP(PW_49)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(24) $_CN(36)]]
  set _TMP(PW_50) [$_CN(24) getDistribution 1]
  $_TMP(PW_50) setBeginSpacing 0.01
  unset _TMP(PW_50)
  set _TMP(PW_51) [$_CN(24) getDistribution 1]
  $_TMP(PW_51) setEndSpacing 0.01
  unset _TMP(PW_51)
  set _TMP(PW_52) [$_CN(36) getDistribution 1]
  $_TMP(PW_52) setBeginSpacing 0.01
  unset _TMP(PW_52)
  set _TMP(PW_53) [$_CN(36) getDistribution 1]
  $_TMP(PW_53) setEndSpacing 0.01
  unset _TMP(PW_53)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(30) $_CN(35) $_CN(34) $_CN(32) $_CN(28)]]
  set _TMP(PW_54) [$_CN(28) getDistribution 1]
  $_TMP(PW_54) setBeginSpacing 0.0011999999999999999
  unset _TMP(PW_54)
  set _TMP(PW_55) [$_CN(30) getDistribution 1]
  $_TMP(PW_55) setBeginSpacing 0.0011999999999999999
  unset _TMP(PW_55)
  set _TMP(PW_56) [$_CN(32) getDistribution 1]
  $_TMP(PW_56) setBeginSpacing 0.0011999999999999999
  unset _TMP(PW_56)
  set _TMP(PW_57) [$_CN(34) getDistribution 1]
  $_TMP(PW_57) setBeginSpacing 0.0011999999999999999
  unset _TMP(PW_57)
  set _TMP(PW_58) [$_CN(35) getDistribution 1]
  $_TMP(PW_58) setBeginSpacing 0.0011999999999999999
  unset _TMP(PW_58)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(10) $_CN(9)]]
  set _TMP(PW_59) [$_CN(10) getDistribution 1]
  $_TMP(PW_59) setBeginSpacing 0.0070000000000000001
  unset _TMP(PW_59)
  set _TMP(PW_60) [$_CN(10) getDistribution 1]
  $_TMP(PW_60) setEndSpacing 0.0070000000000000001
  unset _TMP(PW_60)
  set _TMP(PW_61) [$_CN(9) getDistribution 1]
  $_TMP(PW_61) setBeginSpacing 0.0070000000000000001
  unset _TMP(PW_61)
  set _TMP(PW_62) [$_CN(9) getDistribution 1]
  $_TMP(PW_62) setEndSpacing 0.0070000000000000001
  unset _TMP(PW_62)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(15) $_CN(17)]]
  set _TMP(PW_63) [$_CN(15) getDistribution 1]
  $_TMP(PW_63) setBeginSpacing 0.0044999999999999997
  unset _TMP(PW_63)
  set _TMP(PW_64) [$_CN(15) getDistribution 1]
  $_TMP(PW_64) setEndSpacing 0.0044999999999999997
  unset _TMP(PW_64)
  set _TMP(PW_65) [$_CN(17) getDistribution 1]
  $_TMP(PW_65) setBeginSpacing 0.0044999999999999997
  unset _TMP(PW_65)
  set _TMP(PW_66) [$_CN(17) getDistribution 1]
  $_TMP(PW_66) setEndSpacing 0.0044999999999999997
  unset _TMP(PW_66)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(7) $_CN(3) $_CN(1)]]
  set _TMP(PW_67) [$_CN(1) getDistribution 1]
  $_TMP(PW_67) setBeginSpacing 0.0015
  unset _TMP(PW_67)
  set _TMP(PW_68) [$_CN(1) getDistribution 1]
  $_TMP(PW_68) setEndSpacing 0.0015
  unset _TMP(PW_68)
  set _TMP(PW_69) [$_CN(3) getDistribution 1]
  $_TMP(PW_69) setBeginSpacing 0.0015
  unset _TMP(PW_69)
  set _TMP(PW_70) [$_CN(3) getDistribution 1]
  $_TMP(PW_70) setEndSpacing 0.0015
  unset _TMP(PW_70)
  set _TMP(PW_71) [$_CN(7) getDistribution 1]
  $_TMP(PW_71) setBeginSpacing 0.0015
  unset _TMP(PW_71)
  set _TMP(PW_72) [$_CN(7) getDistribution 1]
  $_TMP(PW_72) setEndSpacing 0.0015
  unset _TMP(PW_72)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(33)]]
  set _TMP(PW_73) [$_CN(33) getDistribution 1]
  $_TMP(PW_73) setEndSpacing 0.0015
  unset _TMP(PW_73)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing(s)}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(10)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(18)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(15)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(8)
  set _TMP(dom_1) [pw::DomainStructured create]
  $_TMP(dom_1) addEdge $_TMP(edge_1)
  $_TMP(dom_1) addEdge $_TMP(edge_2)
  $_TMP(dom_1) addEdge $_TMP(edge_3)
  $_TMP(dom_1) addEdge $_TMP(edge_4)
  unset _TMP(edge_4)
  unset _TMP(edge_3)
  unset _TMP(edge_2)
  unset _TMP(edge_1)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_1)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_5) [pw::Edge create]
  $_TMP(edge_5) addConnector $_CN(9)
  set _TMP(edge_6) [pw::Edge create]
  $_TMP(edge_6) addConnector $_CN(8)
  set _TMP(edge_7) [pw::Edge create]
  $_TMP(edge_7) addConnector $_CN(17)
  set _TMP(edge_8) [pw::Edge create]
  $_TMP(edge_8) addConnector $_CN(18)
  set _TMP(dom_2) [pw::DomainStructured create]
  $_TMP(dom_2) addEdge $_TMP(edge_5)
  $_TMP(dom_2) addEdge $_TMP(edge_6)
  $_TMP(dom_2) addEdge $_TMP(edge_7)
  $_TMP(dom_2) addEdge $_TMP(edge_8)
  unset _TMP(edge_8)
  unset _TMP(edge_7)
  unset _TMP(edge_6)
  unset _TMP(edge_5)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_2)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_9) [pw::Edge create]
  $_TMP(edge_9) addConnector $_CN(15)
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(4)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(1)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(2)
  set _TMP(dom_3) [pw::DomainStructured create]
  $_TMP(dom_3) addEdge $_TMP(edge_9)
  $_TMP(dom_3) addEdge $_TMP(edge_10)
  $_TMP(dom_3) addEdge $_TMP(edge_11)
  $_TMP(dom_3) addEdge $_TMP(edge_12)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
  unset _TMP(edge_9)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_3)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(17)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(4)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(3)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(2)
  set _TMP(dom_4) [pw::DomainStructured create]
  $_TMP(dom_4) addEdge $_TMP(edge_10)
  $_TMP(dom_4) addEdge $_TMP(edge_11)
  $_TMP(dom_4) addEdge $_TMP(edge_12)
  $_TMP(dom_4) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_4)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(1)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(6)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(7)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(5)
  set _TMP(dom_5) [pw::DomainStructured create]
  $_TMP(dom_5) addEdge $_TMP(edge_10)
  $_TMP(dom_5) addEdge $_TMP(edge_11)
  $_TMP(dom_5) addEdge $_TMP(edge_12)
  $_TMP(dom_5) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_5)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(3)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(5)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(7)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(6)
  set _TMP(dom_6) [pw::DomainStructured create]
  $_TMP(dom_6) addEdge $_TMP(edge_10)
  $_TMP(dom_6) addEdge $_TMP(edge_11)
  $_TMP(dom_6) addEdge $_TMP(edge_12)
  $_TMP(dom_6) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_6)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(7)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(34)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(33)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(35)
  set _TMP(dom_7) [pw::DomainStructured create]
  $_TMP(dom_7) addEdge $_TMP(edge_10)
  $_TMP(dom_7) addEdge $_TMP(edge_11)
  $_TMP(dom_7) addEdge $_TMP(edge_12)
  $_TMP(dom_7) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_7)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(14)
  $_TMP(edge_10) addConnector $_CN(24)
  $_TMP(edge_10) addConnector $_CN(16)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(13)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(11)
  $_TMP(edge_12) addConnector $_CN(9)
  $_TMP(edge_12) addConnector $_CN(10)
  $_TMP(edge_12) addConnector $_CN(8)
  $_TMP(edge_12) removeLastConnector
  $_TMP(edge_12) addConnector $_CN(11)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(12)
  set _TMP(dom_8) [pw::DomainStructured create]
  $_TMP(dom_8) addEdge $_TMP(edge_10)
  $_TMP(dom_8) addEdge $_TMP(edge_11)
  $_TMP(dom_8) addEdge $_TMP(edge_12)
  $_TMP(dom_8) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_8)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(13)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(19)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(23)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(22)
  set _TMP(dom_9) [pw::DomainStructured create]
  $_TMP(dom_9) addEdge $_TMP(edge_10)
  $_TMP(dom_9) addEdge $_TMP(edge_11)
  $_TMP(dom_9) addEdge $_TMP(edge_12)
  $_TMP(dom_9) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_9)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(12)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(20)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(21)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(19)
  set _TMP(dom_10) [pw::DomainStructured create]
  $_TMP(dom_10) addEdge $_TMP(edge_10)
  $_TMP(dom_10) addEdge $_TMP(edge_11)
  $_TMP(dom_10) addEdge $_TMP(edge_12)
  $_TMP(dom_10) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_10)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(21)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(30)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(27)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(28)
  set _TMP(dom_11) [pw::DomainStructured create]
  $_TMP(dom_11) addEdge $_TMP(edge_10)
  $_TMP(dom_11) addEdge $_TMP(edge_11)
  $_TMP(dom_11) addEdge $_TMP(edge_12)
  $_TMP(dom_11) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_11)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(23)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(32)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(26)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(28)
  set _TMP(dom_12) [pw::DomainStructured create]
  $_TMP(dom_12) addEdge $_TMP(edge_10)
  $_TMP(dom_12) addEdge $_TMP(edge_11)
  $_TMP(dom_12) addEdge $_TMP(edge_12)
  $_TMP(dom_12) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_12)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(11)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(19)
  $_TMP(edge_11) addConnector $_CN(28)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(25)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(34)
  $_TMP(edge_13) addConnector $_CN(6)
  $_TMP(edge_13) addConnector $_CN(4)
  $_TMP(edge_13) addConnector $_CN(8)
  set _TMP(dom_13) [pw::DomainStructured create]
  $_TMP(dom_13) addEdge $_TMP(edge_10)
  $_TMP(dom_13) addEdge $_TMP(edge_11)
  $_TMP(dom_13) addEdge $_TMP(edge_12)
  $_TMP(dom_13) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_13)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(27)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(29)
  $_TMP(edge_11) addConnector $_CN(36)
  $_TMP(edge_11) addConnector $_CN(31)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(26)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(25)
  $_TMP(edge_13) addConnector $_CN(33)
  $_TMP(edge_13) addConnector $_CN(35)
  $_TMP(edge_13) addConnector $_CN(5)
  $_TMP(edge_13) addConnector $_CN(2)
  $_TMP(edge_13) removeLastConnector
  $_TMP(edge_13) removeLastConnector
  $_TMP(edge_13) removeLastConnector
  $_TMP(edge_13) addConnector $_CN(33)
  $_TMP(edge_13) addConnector $_CN(25)
  set _TMP(dom_14) [pw::DomainStructured create]
  $_TMP(dom_14) addEdge $_TMP(edge_10)
  $_TMP(dom_14) addEdge $_TMP(edge_11)
  $_TMP(dom_14) addEdge $_TMP(edge_12)
  $_TMP(dom_14) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_14)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(edge_10) [pw::Edge create]
  $_TMP(edge_10) addConnector $_CN(30)
  $_TMP(edge_10) addConnector $_CN(20)
  set _TMP(edge_11) [pw::Edge create]
  $_TMP(edge_11) addConnector $_CN(14)
  $_TMP(edge_11) addConnector $_CN(24)
  $_TMP(edge_11) addConnector $_CN(16)
  set _TMP(edge_12) [pw::Edge create]
  $_TMP(edge_12) addConnector $_CN(22)
  $_TMP(edge_12) addConnector $_CN(32)
  set _TMP(edge_13) [pw::Edge create]
  $_TMP(edge_13) addConnector $_CN(31)
  $_TMP(edge_13) addConnector $_CN(36)
  $_TMP(edge_13) addConnector $_CN(29)
  set _TMP(dom_15) [pw::DomainStructured create]
  $_TMP(dom_15) addEdge $_TMP(edge_10)
  $_TMP(dom_15) addEdge $_TMP(edge_11)
  $_TMP(dom_15) addEdge $_TMP(edge_12)
  $_TMP(dom_15) addEdge $_TMP(edge_13)
  unset _TMP(edge_13)
  unset _TMP(edge_12)
  unset _TMP(edge_11)
  unset _TMP(edge_10)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_15)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_10) [pw::Application begin Create]
$_TMP(mode_10) abort
unset _TMP(mode_10)
set _TMP(mode_10) [pw::Application begin Create]
  set _TMP(block_1) [pw::BlockStructured create]
  set _DM(1) [pw::GridEntity getByName "dom-8"]
  set _TMP(face_1) [pw::FaceStructured create]
  $_TMP(face_1) addDomain $_DM(1)
  $_TMP(block_1) addFace $_TMP(face_1)
  set _DM(2) [pw::GridEntity getByName "dom-15"]
  set _TMP(face_2) [pw::FaceStructured create]
  $_TMP(face_2) addDomain $_DM(2)
  $_TMP(block_1) addFace $_TMP(face_2)
  set _DM(3) [pw::GridEntity getByName "dom-11"]
  set _TMP(face_3) [pw::FaceStructured create]
  $_TMP(face_3) addDomain $_DM(3)
  set _DM(4) [pw::GridEntity getByName "dom-10"]
  $_TMP(face_3) addDomain $_DM(4)
  $_TMP(block_1) addFace $_TMP(face_3)
  set _DM(5) [pw::GridEntity getByName "dom-12"]
  set _TMP(face_4) [pw::FaceStructured create]
  $_TMP(face_4) addDomain $_DM(5)
  set _DM(6) [pw::GridEntity getByName "dom-9"]
  $_TMP(face_4) addDomain $_DM(6)
  $_TMP(block_1) addFace $_TMP(face_4)
  set _DM(7) [pw::GridEntity getByName "dom-14"]
  set _TMP(face_5) [pw::FaceStructured create]
  $_TMP(face_5) addDomain $_DM(7)
  $_TMP(block_1) addFace $_TMP(face_5)
  set _DM(8) [pw::GridEntity getByName "dom-13"]
  set _TMP(face_6) [pw::FaceStructured create]
  $_TMP(face_6) addDomain $_DM(8)
  set _DM(9) [pw::GridEntity getByName "dom-1"]
  $_TMP(face_6) addDomain $_DM(9)
  set _DM(10) [pw::GridEntity getByName "dom-2"]
  $_TMP(face_6) addDomain $_DM(10)
  set _DM(11) [pw::GridEntity getByName "dom-3"]
  $_TMP(face_6) addDomain $_DM(11)
  set _DM(12) [pw::GridEntity getByName "dom-4"]
  $_TMP(face_6) addDomain $_DM(12)
  set _DM(13) [pw::GridEntity getByName "dom-5"]
  $_TMP(face_6) addDomain $_DM(13)
  set _DM(14) [pw::GridEntity getByName "dom-7"]
  $_TMP(face_6) addDomain $_DM(14)
  set _DM(15) [pw::GridEntity getByName "dom-6"]
  $_TMP(face_6) addDomain $_DM(15)
  $_TMP(face_6) addDomain -linkage [list 4 1 1 1 4 1 0] $_DM(14)
  $_TMP(face_6) addDomain -linkage [list 4 1 4 4 4 1 0] $_DM(8)
  $_TMP(block_1) addFace $_TMP(face_6)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(face_6)
unset _TMP(face_5)
unset _TMP(face_4)
unset _TMP(face_3)
unset _TMP(face_2)
unset _TMP(face_1)
unset _TMP(block_1)
pw::Application markUndoLevel {Assemble Block}

# delete temporal connectors
set _CN(1) [pw::GridEntity getByName "con-11"]
set _CN(2) [pw::GridEntity getByName "con-12"]
set _CN(3) [pw::GridEntity getByName "con-13"]
set _CN(4) [pw::GridEntity getByName "con-14"]
pw::Entity delete [list $_CN(1) $_CN(2) $_CN(3) $_CN(4)]
pw::Application markUndoLevel {Delete}

# export

pw::Application setCAESolver {ANSYS FLUENT} 3
pw::Application markUndoLevel {Select Solver}

set _DM(1) [pw::GridEntity getByName "dom-7"]
set _DM(2) [pw::GridEntity getByName "dom-13"]
set _TMP(PW_1) [pw::BoundaryCondition getByName "Unspecified"]
set _DM(3) [pw::GridEntity getByName "dom-1"]
set _DM(4) [pw::GridEntity getByName "dom-2"]
set _DM(5) [pw::GridEntity getByName "dom-3"]
set _DM(6) [pw::GridEntity getByName "dom-4"]
set _DM(7) [pw::GridEntity getByName "dom-5"]
set _DM(8) [pw::GridEntity getByName "dom-6"]
set _DM(9) [pw::GridEntity getByName "dom-8"]
set _DM(10) [pw::GridEntity getByName "dom-9"]
set _DM(11) [pw::GridEntity getByName "dom-10"]
set _DM(12) [pw::GridEntity getByName "dom-11"]
set _DM(13) [pw::GridEntity getByName "dom-12"]
set _DM(14) [pw::GridEntity getByName "dom-14"]
set _DM(15) [pw::GridEntity getByName "dom-15"]
set _BL(1) [pw::GridEntity getByName "blk-1"]
set _TMP(PW_2) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_3) [pw::BoundaryCondition getByName "bc-2"]
unset _TMP(PW_2)
$_TMP(PW_3) setName "ff"
pw::Application markUndoLevel {Name BC}

$_TMP(PW_3) setPhysicalType {Pressure Far Field}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_3) apply [list [list $_BL(1) $_DM(13)] [list $_BL(1) $_DM(10)] [list $_BL(1) $_DM(11)] [list $_BL(1) $_DM(12)] [list $_BL(1) $_DM(14)] [list $_BL(1) $_DM(15)]]
pw::Application markUndoLevel {Set BC}

set _TMP(PW_4) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_5) [pw::BoundaryCondition getByName "bc-3"]
unset _TMP(PW_4)
$_TMP(PW_5) setName "sym"
pw::Application markUndoLevel {Name BC}

$_TMP(PW_5) apply [list [list $_BL(1) $_DM(9)]]
pw::Application markUndoLevel {Set BC}

set _TMP(PW_6) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_7) [pw::BoundaryCondition getByName "bc-4"]
unset _TMP(PW_6)
$_TMP(PW_7) setName "wall"
pw::Application markUndoLevel {Name BC}

$_TMP(PW_7) apply [list [list $_BL(1) $_DM(3)] [list $_BL(1) $_DM(4)] [list $_BL(1) $_DM(6)] [list $_BL(1) $_DM(5)] [list $_BL(1) $_DM(7)] [list $_BL(1) $_DM(8)]]
pw::Application markUndoLevel {Set BC}

unset _TMP(PW_1)
unset _TMP(PW_3)
unset _TMP(PW_5)
unset _TMP(PW_7)
set _TMP(mode_1) [pw::Application begin CaeExport [pw::Entity sort [list $_BL(1)]]]
  $_TMP(mode_1) initialize -type CAE {E:/ucav_tmp_sym.cas}
  if {![$_TMP(mode_1) verify]} {
    error "Data verification failed."
  }
  $_TMP(mode_1) write
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application clearClipboard
pw::Application setClipboard [list $_BL(1)]
pw::Application markUndoLevel {Copy}

set _TMP(mode_2) [pw::Application begin Paste]
  set _TMP(PW_8) [$_TMP(mode_2) getEntities]
  set _TMP(mode_3) [pw::Application begin Modify $_TMP(PW_8)]
    pw::Entity transform [pwu::Transform mirroring {0 0 1} 0] [$_TMP(mode_3) getEntities]
  $_TMP(mode_3) end
  unset _TMP(mode_3)
$_TMP(mode_2) end
unset _TMP(mode_2)
pw::Application markUndoLevel {Paste}

unset _TMP(PW_8)
set _BL(2) [pw::GridEntity getByName "blk-2"]
set _TMP(mode_4) [pw::Application begin CaeExport [pw::Entity sort [list $_BL(1) $_BL(2)]]]
  $_TMP(mode_4) initialize -type CAE {E:/ucav_tmp_sym2.cas}
  if {![$_TMP(mode_4) verify]} {
    error "Data verification failed."
  }
  $_TMP(mode_4) write
$_TMP(mode_4) end
unset _TMP(mode_4)

pw::Application exit 0