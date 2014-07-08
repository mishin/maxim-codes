# Pointwise V16.03R2 Journal file - Wed Jun 27 15:45:01 2012

package require PWI_Glyph 2.3

pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}

pw::Application clearModified

set _TMP(mode_10) [pw::Application begin DatabaseImport]
$_TMP(mode_10) initialize -type Automatic {D:/codes/actools/pyAC/actools/temp/7y8P20U.igs}
$_TMP(mode_10) read
$_TMP(mode_10) convert
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Import Database}

pw::Application setCAESolver {ANSYS FLUENT} 3
pw::Application markUndoLevel {Select Solver}

pw::Application setCAESolver {ANSYS FLUENT} 2
pw::Application markUndoLevel {Set Dimension 2D}

set _DB(1) [pw::DatabaseEntity getByName {Spline.2-3}]
set _DB(2) [pw::DatabaseEntity getByName {Spline.1-1}]
set _TMP(PW_56) [pw::Connector createOnDatabase -merge 0 -reject _TMP(unused) [list $_DB(1) $_DB(2)]]
unset _TMP(unused)
unset _TMP(PW_56)
pw::Application markUndoLevel {Connectors On DB Entities}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_57) [pw::SegmentSpline create]
set _CN(1) [pw::GridEntity getByName {con-1}]
set _CN(2) [pw::GridEntity getByName {con-2}]
$_TMP(PW_57) addPoint [$_CN(1) getPosition -arc 1]
$_TMP(PW_57) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(1) getPosition -arc 1]] {20 0 0}]
set _TMP(con_8) [pw::Connector create]
$_TMP(con_8) addSegment $_TMP(PW_57)
unset _TMP(PW_57)
$_TMP(con_8) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_58) [pw::SegmentSpline create]
set _CN(3) [pw::GridEntity getByName {con-3}]
$_TMP(PW_58) addPoint [$_CN(3) getPosition -arc 1]
$_TMP(PW_58) addPoint [pwu::Vector3 add [$_CN(3) getPosition -arc 1] {0 15 0}]
unset _TMP(con_8)
set _TMP(con_9) [pw::Connector create]
$_TMP(con_9) addSegment $_TMP(PW_58)
unset _TMP(PW_58)
$_TMP(con_9) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_59) [pw::SegmentSpline create]
set _CN(4) [pw::GridEntity getByName {con-4}]
$_TMP(PW_59) addPoint [$_CN(4) getPosition -arc 1]
$_TMP(PW_59) addPoint [pwu::Vector3 add [$_CN(4) getPosition -arc 1] {-18 0 0}]
unset _TMP(con_9)
set _TMP(con_10) [pw::Connector create]
$_TMP(con_10) addSegment $_TMP(PW_59)
unset _TMP(PW_59)
$_TMP(con_10) calculateDimension
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_60) [pw::SegmentSpline create]
$_TMP(PW_60) delete
unset _TMP(PW_60)
$_TMP(mode_10) abort
unset _TMP(mode_10)
unset _TMP(con_10)
set _CN(5) [pw::GridEntity getByName {con-5}]
pw::Application setClipboard [list $_CN(4) $_CN(5)]
pw::Application markUndoLevel {Copy}

set _TMP(mode_10) [pw::Application begin Paste]
set _TMP(PW_61) [$_TMP(mode_10) getEntities]
set _TMP(mode_11) [pw::Application begin Modify $_TMP(PW_61)]
pw::Entity transform [pwu::Transform mirroring {0 1 0} 0] [$_TMP(mode_11) getEntities]
$_TMP(mode_11) end
unset _TMP(mode_11)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Paste}

unset _TMP(PW_61)
set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_62) [pw::SegmentSpline create]
$_TMP(PW_62) delete
unset _TMP(PW_62)
set _TMP(PW_63) [pw::SegmentCircle create]
set _CN(6) [pw::GridEntity getByName {con-7}]
$_TMP(PW_63) addPoint [$_CN(5) getPosition -arc 1]
$_TMP(PW_63) addPoint [$_CN(6) getPosition -arc 1]
$_TMP(PW_63) setAngle 180.0 {0 0 1}
set _TMP(con_11) [pw::Connector create]
$_TMP(con_11) addSegment $_TMP(PW_63)
$_TMP(con_11) calculateDimension
unset _TMP(PW_63)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create Connector}

unset _TMP(con_11)
set _CN(7) [pw::GridEntity getByName {con-6}]
set _TMP(PW_64) [pw::Collection create]
$_TMP(PW_64) set [list $_CN(4) $_CN(7)]
$_TMP(PW_64) do setDimension 51
$_TMP(PW_64) delete
unset _TMP(PW_64)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_65) [pw::Collection create]
$_TMP(PW_65) set [list $_CN(5) $_CN(6) $_CN(3)]
$_TMP(PW_65) do setDimension 65
$_TMP(PW_65) delete
unset _TMP(PW_65)
pw::Application markUndoLevel {Dimension}

set _CN(8) [pw::GridEntity getByName {con-8}]
set _TMP(PW_66) [pw::Collection create]
$_TMP(PW_66) set [list $_CN(8)]
$_TMP(PW_66) do setDimension 201
$_TMP(PW_66) delete
unset _TMP(PW_66)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_67) [pw::Collection create]
$_TMP(PW_67) set [list $_CN(1) $_CN(2)]
$_TMP(PW_67) do setDimension 101
$_TMP(PW_67) delete
unset _TMP(PW_67)
pw::Application markUndoLevel {Dimension}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(edge_9) [pw::Edge create]
$_TMP(edge_9) addConnector $_CN(4)
set _TMP(edge_10) [pw::Edge create]
$_TMP(edge_10) addConnector $_CN(5)
$_TMP(edge_10) addConnector $_CN(8)
$_TMP(edge_10) addConnector $_CN(6)
set _TMP(edge_11) [pw::Edge create]
$_TMP(edge_11) addConnector $_CN(7)
set _TMP(edge_12) [pw::Edge create]
$_TMP(edge_12) addConnector $_CN(3)
$_TMP(edge_12) addConnector $_CN(2)
$_TMP(edge_12) addConnector $_CN(1)
$_TMP(edge_12) addConnector $_CN(3)
set _TMP(dom_3) [pw::DomainStructured create]
$_TMP(dom_3) addEdge $_TMP(edge_9)
$_TMP(dom_3) addEdge $_TMP(edge_10)
$_TMP(dom_3) addEdge $_TMP(edge_11)
$_TMP(dom_3) addEdge $_TMP(edge_12)
unset _TMP(edge_11)
unset _TMP(edge_10)
unset _TMP(edge_9)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_3)
pw::Application markUndoLevel {Assemble Domain}

unset _TMP(edge_12)
set _TMP(mode_10) [pw::Application begin Create]
$_TMP(mode_10) abort
unset _TMP(mode_10)
set _TMP(mode_10) [pw::Application begin Modify [list $_CN(7) $_CN(4)]]
set _TMP(PW_68) [$_CN(4) getDistribution 1]
$_TMP(PW_68) setBeginSpacing 1.9999999999999999e-006
unset _TMP(PW_68)
set _TMP(PW_69) [$_CN(7) getDistribution 1]
$_TMP(PW_69) setBeginSpacing 1.9999999999999999e-006
unset _TMP(PW_69)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(7) $_CN(6) $_CN(5) $_CN(4) $_CN(3)]]
set _TMP(PW_70) [$_CN(3) getDistribution 1]
$_TMP(PW_70) setEndSpacing 1
unset _TMP(PW_70)
set _TMP(PW_71) [$_CN(4) getDistribution 1]
$_TMP(PW_71) setEndSpacing 1
unset _TMP(PW_71)
set _TMP(PW_72) [$_CN(5) getDistribution 1]
$_TMP(PW_72) setBeginSpacing 1
unset _TMP(PW_72)
set _TMP(PW_73) [$_CN(7) getDistribution 1]
$_TMP(PW_73) setEndSpacing 1
unset _TMP(PW_73)
set _TMP(PW_74) [$_CN(6) getDistribution 1]
$_TMP(PW_74) setBeginSpacing 1
unset _TMP(PW_74)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(6) $_CN(5) $_CN(3)]]
set _TMP(PW_75) [$_CN(3) getDistribution 1]
$_TMP(PW_75) setBeginSpacing 0.001
unset _TMP(PW_75)
set _TMP(PW_76) [$_CN(5) getDistribution 1]
$_TMP(PW_76) setEndSpacing 0.001
unset _TMP(PW_76)
set _TMP(PW_77) [$_CN(6) getDistribution 1]
$_TMP(PW_77) setEndSpacing 0.001
unset _TMP(PW_77)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(8)]]
set _TMP(PW_78) [$_CN(8) getDistribution 1]
$_TMP(PW_78) setBeginSpacing 0.01
unset _TMP(PW_78)
set _TMP(PW_79) [$_CN(8) getDistribution 1]
$_TMP(PW_79) setEndSpacing 0.01
unset _TMP(PW_79)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(2) $_CN(1)]]
set _TMP(PW_80) [$_CN(1) getDistribution 1]
$_TMP(PW_80) setBeginSpacing 0.001
unset _TMP(PW_80)
set _TMP(PW_81) [$_CN(2) getDistribution 1]
$_TMP(PW_81) setBeginSpacing 0.001
unset _TMP(PW_81)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(2) $_CN(1)]]
set _TMP(PW_82) [$_CN(1) getDistribution 1]
$_TMP(PW_82) setEndSpacing 0.001
unset _TMP(PW_82)
set _TMP(PW_83) [$_CN(2) getDistribution 1]
$_TMP(PW_83) setEndSpacing 0.001
unset _TMP(PW_83)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(PW_84) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

$_TMP(PW_84) setPhysicalType {Wall}
pw::Application markUndoLevel {Change BC Type}

set _TMP(PW_85) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

$_TMP(PW_85) setPhysicalType {Pressure Far Field}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_84) apply [list $_CN(1) $_CN(2)]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_85) apply [list $_CN(4) $_CN(5) $_CN(8) $_CN(6) $_CN(7)]
pw::Application markUndoLevel {Set BC}

unset _TMP(PW_84)
unset _TMP(PW_85)
set _DM(1) [pw::GridEntity getByName {dom-1}]
pw::Application export [list $_DM(1)] {D:/codes/actools/pyAC/actools/temp/case1.cas}
pw::Application exit 0
