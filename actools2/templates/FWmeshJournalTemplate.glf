# Pointwise V16.03R2 Journal file - Tue Jul  8 15:38:17 2014

package require PWI_Glyph 2.3

pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}

pw::Application clearModified

set _TMP(mode_1) [pw::Application begin DatabaseImport]
$_TMP(mode_1) initialize -type Automatic {E:/codes/pyCodes/actools2/temp/fw.igs}
$_TMP(mode_1) setAttribute EntityVisibility ShowAllAndHideSupports
$_TMP(mode_1) read
$_TMP(mode_1) convert
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Import Database}


# Appended by Pointwise V16.03R2 - Tue Jul  8 15:41:32 2014

set _TMP(mode_1) [pw::Application begin Create]
set _TMP(PW_1) [pw::SegmentSpline create]
set _DB(1) [pw::DatabaseEntity getByName {quilt-7}]
set _DB(2) [pw::DatabaseEntity getByName {quilt-6}]
set _DB(3) [pw::DatabaseEntity getByName {model-7}]
set _DB(4) [pw::DatabaseEntity getByName {model-6}]
$_TMP(PW_1) delete
unset _TMP(PW_1)
set _TMP(PW_2) [pw::SegmentSurfaceSpline create]
$_TMP(PW_2) addPoint [list 1 1 1 $_DB(1)]
$_TMP(PW_2) addPoint [list 1 0 1 $_DB(1)]
set _TMP(con_1) [pw::Connector create]
$_TMP(con_1) addSegment $_TMP(PW_2)
unset _TMP(PW_2)
$_TMP(con_1) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_2) [pw::Application begin Create]
set _TMP(PW_3) [pw::SegmentSpline create]
set _CN(1) [pw::GridEntity getByName {con-1}]
$_TMP(PW_3) delete
unset _TMP(PW_3)
set _TMP(PW_4) [pw::SegmentSurfaceSpline create]
$_TMP(PW_4) addPoint [list 1 0 1 $_DB(2)]
$_TMP(PW_4) addPoint [list 1 1 1 $_DB(2)]
unset _TMP(con_1)
set _TMP(con_2) [pw::Connector create]
$_TMP(con_2) addSegment $_TMP(PW_4)
unset _TMP(PW_4)
$_TMP(con_2) calculateDimension
$_TMP(mode_2) end
unset _TMP(mode_2)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_3) [pw::Application begin Create]
set _TMP(PW_5) [pw::SegmentSpline create]
set _CN(2) [pw::GridEntity getByName {con-2}]
set _DB(5) [pw::DatabaseEntity getByName {model-5}]
set _DB(6) [pw::DatabaseEntity getByName {quilt-5}]
set _DB(7) [pw::DatabaseEntity getByName {model-8}]
set _DB(8) [pw::DatabaseEntity getByName {quilt-8}]
$_TMP(PW_5) delete
unset _TMP(PW_5)
set _TMP(PW_6) [pw::SegmentSurfaceSpline create]
$_TMP(PW_6) addPoint [list 0 1 1 $_DB(1)]
$_TMP(PW_6) addPoint [list 0 0 1 $_DB(1)]
unset _TMP(con_2)
set _TMP(con_3) [pw::Connector create]
$_TMP(con_3) addSegment $_TMP(PW_6)
unset _TMP(PW_6)
$_TMP(con_3) calculateDimension
$_TMP(mode_3) end
unset _TMP(mode_3)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_4) [pw::Application begin Create]
set _TMP(PW_7) [pw::SegmentSpline create]
set _CN(3) [pw::GridEntity getByName {con-3}]
$_TMP(PW_7) delete
unset _TMP(PW_7)
set _TMP(PW_8) [pw::SegmentSurfaceSpline create]
$_TMP(PW_8) addPoint [list 0 0 1 $_DB(2)]
$_TMP(PW_8) addPoint [list 0 1 1 $_DB(2)]
unset _TMP(con_3)
set _TMP(con_4) [pw::Connector create]
$_TMP(con_4) addSegment $_TMP(PW_8)
unset _TMP(PW_8)
$_TMP(con_4) calculateDimension
$_TMP(mode_4) end
unset _TMP(mode_4)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_5) [pw::Application begin Create]
set _TMP(PW_9) [pw::SegmentSpline create]
set _CN(4) [pw::GridEntity getByName {con-4}]
set _DB(9) [pw::DatabaseEntity getByName {model-1}]
set _DB(10) [pw::DatabaseEntity getByName {quilt-1}]
set _DB(11) [pw::DatabaseEntity getByName {quilt-4}]
set _DB(12) [pw::DatabaseEntity getByName {model-4}]
$_TMP(PW_9) delete
unset _TMP(PW_9)
set _TMP(PW_10) [pw::SegmentSurfaceSpline create]
$_TMP(PW_10) addPoint [list 0 1 1 $_DB(11)]
$_TMP(PW_10) addPoint [list 0 0 1 $_DB(11)]
unset _TMP(con_4)
set _TMP(con_5) [pw::Connector create]
$_TMP(con_5) addSegment $_TMP(PW_10)
unset _TMP(PW_10)
$_TMP(con_5) calculateDimension
$_TMP(mode_5) end
unset _TMP(mode_5)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_6) [pw::Application begin Create]
set _TMP(PW_11) [pw::SegmentSpline create]
set _CN(5) [pw::GridEntity getByName {con-5}]
$_TMP(PW_11) delete
unset _TMP(PW_11)
set _TMP(PW_12) [pw::SegmentSurfaceSpline create]
$_TMP(PW_12) addPoint [list 0 0 1 $_DB(10)]
$_TMP(PW_12) addPoint [list 0 1 1 $_DB(10)]
unset _TMP(con_5)
set _TMP(con_6) [pw::Connector create]
$_TMP(con_6) addSegment $_TMP(PW_12)
unset _TMP(PW_12)
$_TMP(con_6) calculateDimension
$_TMP(mode_6) end
unset _TMP(mode_6)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_7) [pw::Application begin Create]
set _TMP(PW_13) [pw::SegmentSpline create]
set _CN(6) [pw::GridEntity getByName {con-6}]
set _DB(13) [pw::DatabaseEntity getByName {quilt-2}]
set _DB(14) [pw::DatabaseEntity getByName {quilt-3}]
set _DB(15) [pw::DatabaseEntity getByName {model-3}]
set _DB(16) [pw::DatabaseEntity getByName {model-2}]
$_TMP(PW_13) delete
unset _TMP(PW_13)
set _TMP(PW_14) [pw::SegmentSurfaceSpline create]
$_TMP(PW_14) addPoint [list 0 0 1 $_DB(13)]
$_TMP(PW_14) addPoint [list 0 1 1 $_DB(13)]
unset _TMP(con_6)
set _TMP(con_7) [pw::Connector create]
$_TMP(con_7) addSegment $_TMP(PW_14)
unset _TMP(PW_14)
$_TMP(con_7) calculateDimension
$_TMP(mode_7) end
unset _TMP(mode_7)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_8) [pw::Application begin Create]
set _TMP(PW_15) [pw::SegmentSpline create]
set _CN(7) [pw::GridEntity getByName {con-7}]
$_TMP(PW_15) delete
unset _TMP(PW_15)
set _TMP(PW_16) [pw::SegmentSurfaceSpline create]
$_TMP(PW_16) addPoint [list 1 1 1 $_DB(11)]
$_TMP(PW_16) addPoint [list 1 0 1 $_DB(11)]
unset _TMP(con_7)
set _TMP(con_8) [pw::Connector create]
$_TMP(con_8) addSegment $_TMP(PW_16)
unset _TMP(PW_16)
$_TMP(con_8) calculateDimension
$_TMP(mode_8) end
unset _TMP(mode_8)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_9) [pw::Application begin Create]
set _TMP(PW_17) [pw::SegmentSpline create]
set _CN(8) [pw::GridEntity getByName {con-8}]
$_TMP(PW_17) delete
unset _TMP(PW_17)
set _TMP(PW_18) [pw::SegmentSurfaceSpline create]
$_TMP(PW_18) addPoint [list 1 0 1 $_DB(13)]
$_TMP(PW_18) addPoint [list 1 1 1 $_DB(13)]
unset _TMP(con_8)
set _TMP(con_9) [pw::Connector create]
$_TMP(con_9) addSegment $_TMP(PW_18)
unset _TMP(PW_18)
$_TMP(con_9) calculateDimension
$_TMP(mode_9) end
unset _TMP(mode_9)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_19) [pw::SegmentSpline create]
set _CN(9) [pw::GridEntity getByName {con-9}]
$_TMP(PW_19) delete
unset _TMP(PW_19)
set _TMP(PW_20) [pw::SegmentSurfaceSpline create]
$_TMP(PW_20) addPoint [list 1 1 1 $_DB(14)]
$_TMP(PW_20) addPoint [list 1 0 1 $_DB(14)]
unset _TMP(con_9)
set _TMP(con_10) [pw::Connector create]
$_TMP(con_10) addSegment $_TMP(PW_20)
unset _TMP(PW_20)
$_TMP(con_10) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_21) [pw::SegmentSpline create]
set _CN(10) [pw::GridEntity getByName {con-10}]
$_TMP(PW_21) delete
unset _TMP(PW_21)
set _TMP(PW_22) [pw::SegmentSurfaceSpline create]
$_TMP(PW_22) addPoint [list 0 1 1 $_DB(8)]
$_TMP(PW_22) addPoint [list 1 1 1 $_DB(8)]
unset _TMP(con_10)
set _TMP(con_11) [pw::Connector create]
$_TMP(con_11) addSegment $_TMP(PW_22)
unset _TMP(PW_22)
$_TMP(con_11) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_23) [pw::SegmentSpline create]
set _CN(11) [pw::GridEntity getByName {con-11}]
$_TMP(PW_23) delete
unset _TMP(PW_23)
set _TMP(PW_24) [pw::SegmentSurfaceSpline create]
$_TMP(PW_24) addPoint [$_CN(3) getPosition -arc 0]
$_TMP(PW_24) addPoint [$_CN(1) getPosition -arc 0]
unset _TMP(con_11)
set _TMP(con_12) [pw::Connector create]
$_TMP(con_12) addSegment $_TMP(PW_24)
unset _TMP(PW_24)
$_TMP(con_12) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_25) [pw::SegmentSpline create]
set _CN(12) [pw::GridEntity getByName {con-12}]
$_TMP(PW_25) delete
unset _TMP(PW_25)
set _TMP(PW_26) [pw::SegmentSurfaceSpline create]
$_TMP(PW_26) addPoint [list 0 0 1 $_DB(8)]
$_TMP(PW_26) addPoint [list 1 0 1 $_DB(8)]
unset _TMP(con_12)
set _TMP(con_13) [pw::Connector create]
$_TMP(con_13) addSegment $_TMP(PW_26)
unset _TMP(PW_26)
$_TMP(con_13) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_27) [pw::SegmentSpline create]
set _CN(13) [pw::GridEntity getByName {con-13}]
$_TMP(PW_27) delete
unset _TMP(PW_27)
set _TMP(PW_28) [pw::SegmentSurfaceSpline create]
$_TMP(PW_28) addPoint [$_CN(3) getPosition -arc 1]
$_TMP(PW_28) addPoint [$_CN(1) getPosition -arc 1]
unset _TMP(con_13)
set _TMP(con_14) [pw::Connector create]
$_TMP(con_14) addSegment $_TMP(PW_28)
unset _TMP(PW_28)
$_TMP(con_14) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_29) [pw::SegmentSpline create]
set _CN(14) [pw::GridEntity getByName {con-14}]
$_TMP(PW_29) delete
unset _TMP(PW_29)
set _TMP(PW_30) [pw::SegmentSurfaceSpline create]
$_TMP(PW_30) addPoint [$_CN(5) getPosition -arc 0]
$_TMP(PW_30) addPoint [list 1 1 1 $_DB(11)]
unset _TMP(con_14)
set _TMP(con_15) [pw::Connector create]
$_TMP(con_15) addSegment $_TMP(PW_30)
unset _TMP(PW_30)
$_TMP(con_15) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_31) [pw::SegmentSpline create]
set _CN(15) [pw::GridEntity getByName {con-15}]
$_TMP(PW_31) delete
unset _TMP(PW_31)
set _TMP(PW_32) [pw::SegmentSurfaceSpline create]
$_TMP(PW_32) addPoint [$_CN(7) getPosition -arc 0]
$_TMP(PW_32) addPoint [$_CN(9) getPosition -arc 0]
unset _TMP(con_15)
set _TMP(con_16) [pw::Connector create]
$_TMP(con_16) addSegment $_TMP(PW_32)
unset _TMP(PW_32)
$_TMP(con_16) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_33) [pw::SegmentSpline create]
set _CN(16) [pw::GridEntity getByName {con-16}]
$_TMP(PW_33) delete
unset _TMP(PW_33)
set _TMP(PW_34) [pw::SegmentSurfaceSpline create]
$_TMP(PW_34) addPoint [$_CN(5) getPosition -arc 1]
$_TMP(PW_34) addPoint [list 1 0 1 $_DB(11)]
unset _TMP(con_16)
set _TMP(con_17) [pw::Connector create]
$_TMP(con_17) addSegment $_TMP(PW_34)
unset _TMP(PW_34)
$_TMP(con_17) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_35) [pw::SegmentSpline create]
set _CN(17) [pw::GridEntity getByName {con-17}]
$_TMP(PW_35) delete
unset _TMP(PW_35)
set _TMP(PW_36) [pw::SegmentSurfaceSpline create]
$_TMP(PW_36) addPoint [$_CN(7) getPosition -arc 1]
$_TMP(PW_36) addPoint [$_CN(9) getPosition -arc 1]
unset _TMP(con_17)
set _TMP(con_18) [pw::Connector create]
$_TMP(con_18) addSegment $_TMP(PW_36)
unset _TMP(PW_36)
$_TMP(con_18) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_37) [pw::SegmentSpline create]
set _CN(18) [pw::GridEntity getByName {con-18}]
$_TMP(PW_37) delete
unset _TMP(PW_37)
$_TMP(mode_10) abort
unset _TMP(mode_10)
unset _TMP(con_18)
