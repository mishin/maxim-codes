package require PWI_Glyph 2.3
pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}
pw::Application clearModified

set _TMP(mode_1) [pw::Application begin Create]
set _TMP(PW_1) [pw::SegmentSpline create]
$_TMP(PW_1) addPoint {0 0 0}
$_TMP(PW_1) addPoint {1 0 0}
$_TMP(PW_1) setSlope Akima
set _TMP(con_1) [pw::Connector create]
$_TMP(con_1) addSegment $_TMP(PW_1)
$_TMP(con_1) calculateDimension
unset _TMP(PW_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create Connector}
unset _TMP(con_1)
set _TMP(mode_2) [pw::Application begin Create]
set _TMP(PW_2) [pw::SegmentSpline create]
set _CN(1) [pw::GridEntity getByName {con-1}]
$_TMP(PW_2) addPoint [$_CN(1) getPosition -arc 0]
$_TMP(PW_2) addPoint {0.354753 -1.24569 -0}
$_TMP(PW_2) addPoint [$_CN(1) getPosition -arc 1]
$_TMP(PW_2) setSlope Akima
set _TMP(con_2) [pw::Connector create]
$_TMP(con_2) addSegment $_TMP(PW_2)
$_TMP(con_2) calculateDimension
unset _TMP(PW_2)
$_TMP(mode_2) end
unset _TMP(mode_2)
pw::Application markUndoLevel {Create Connector}

unset _TMP(con_2)
set _TMP(mode_3) [pw::Application begin Create]
set _TMP(PW_3) [pw::SegmentSpline create]
$_TMP(PW_3) delete
unset _TMP(PW_3)
$_TMP(mode_3) abort
unset _TMP(mode_3)
set _TMP(mode_4) [pw::Application begin Create]
set _TMP(PW_4) [pw::SegmentSpline create]
set _CN(2) [pw::GridEntity getByName {con-2}]
$_TMP(PW_4) addPoint [$_CN(1) getPosition -arc 1]
$_TMP(PW_4) addPoint [pwu::Vector3 add [$_CN(1) getPosition -arc 1] {10 0 0}]
set _TMP(con_3) [pw::Connector create]
$_TMP(con_3) addSegment $_TMP(PW_4)
unset _TMP(PW_4)
$_TMP(con_3) calculateDimension
$_TMP(mode_4) end
unset _TMP(mode_4)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_5) [pw::Application begin Create]
set _TMP(PW_5) [pw::SegmentSpline create]
set _CN(3) [pw::GridEntity getByName {con-3}]
$_TMP(PW_5) addPoint [$_CN(3) getPosition -arc 1]
$_TMP(PW_5) addPoint [pwu::Vector3 add [$_CN(3) getPosition -arc 1] {0 8 0}]
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
set _CN(4) [pw::GridEntity getByName {con-4}]
$_TMP(PW_6) addPoint [$_CN(4) getPosition -arc 1]
$_TMP(PW_6) addPoint [pwu::Vector3 add [$_CN(4) getPosition -arc 1] {-10 0 0}]
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
$_TMP(PW_7) addPoint [$_CN(3) getPosition -arc 1]
$_TMP(PW_7) addPoint [pwu::Vector3 add [$_CN(3) getPosition -arc 1] {0 -8 0}]
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
set _CN(5) [pw::GridEntity getByName {con-6}]
$_TMP(PW_8) addPoint [$_CN(5) getPosition -arc 1]
$_TMP(PW_8) addPoint [pwu::Vector3 add [$_CN(5) getPosition -arc 1] {-10 0 0}]
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
$_TMP(PW_9) delete
unset _TMP(PW_9)
$_TMP(mode_9) abort
unset _TMP(mode_9)
unset _TMP(con_7)
set _TMP(mode_10) [pw::Application begin Create]
set _TMP(PW_10) [pw::SegmentSpline create]
$_TMP(PW_10) delete
unset _TMP(PW_10)
set _TMP(PW_11) [pw::SegmentCircle create]
set _CN(6) [pw::GridEntity getByName {con-5}]
set _CN(7) [pw::GridEntity getByName {con-7}]
$_TMP(PW_11) addPoint [$_CN(6) getPosition -arc 1]
$_TMP(PW_11) addPoint [$_CN(7) getPosition -arc 1]
$_TMP(PW_11) setAngle 180 {0 0 1}
set _TMP(con_8) [pw::Connector create]
$_TMP(con_8) addSegment $_TMP(PW_11)
$_TMP(con_8) calculateDimension
unset _TMP(PW_11)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Create Connector}

unset _TMP(con_8)
set _TMP(PW_12) [pw::Collection create]
$_TMP(PW_12) set [list $_CN(6) $_CN(3) $_CN(7)]
$_TMP(PW_12) do setDimension 22
$_TMP(PW_12) delete
unset _TMP(PW_12)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_13) [pw::Collection create]
$_TMP(PW_13) set [list $_CN(4) $_CN(5)]
$_TMP(PW_13) do setDimension 33
$_TMP(PW_13) delete
unset _TMP(PW_13)
pw::Application markUndoLevel {Dimension}

set _TMP(PW_14) [pw::Collection create]
$_TMP(PW_14) set [list $_CN(2) $_CN(1)]
$_TMP(PW_14) do setDimension 44
$_TMP(PW_14) delete
unset _TMP(PW_14)
pw::Application markUndoLevel {Dimension}

set _CN(8) [pw::GridEntity getByName {con-8}]
set _TMP(PW_15) [pw::Collection create]
$_TMP(PW_15) set [list $_CN(8)]
$_TMP(PW_15) do setDimension 87
$_TMP(PW_15) delete
unset _TMP(PW_15)
pw::Application markUndoLevel {Dimension}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(2) $_CN(1)]]
set _TMP(PW_16) [$_CN(1) getDistribution 1]
$_TMP(PW_16) setBeginSpacing 0.001
unset _TMP(PW_16)
set _TMP(PW_17) [$_CN(2) getDistribution 1]
$_TMP(PW_17) setBeginSpacing 0.001
unset _TMP(PW_17)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(3) $_CN(2) $_CN(1)]]
set _TMP(PW_18) [$_CN(1) getDistribution 1]
$_TMP(PW_18) setEndSpacing 0.0001
unset _TMP(PW_18)
set _TMP(PW_19) [$_CN(2) getDistribution 1]
$_TMP(PW_19) setEndSpacing 0.0001
unset _TMP(PW_19)
set _TMP(PW_20) [$_CN(3) getDistribution 1]
$_TMP(PW_20) setBeginSpacing 0.0001
unset _TMP(PW_20)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(7) $_CN(6) $_CN(3)]]
set _TMP(PW_21) [$_CN(3) getDistribution 1]
$_TMP(PW_21) setBeginSpacing 0.00050000000000000001
unset _TMP(PW_21)
set _TMP(PW_22) [$_CN(6) getDistribution 1]
$_TMP(PW_22) setEndSpacing 0.00050000000000000001
unset _TMP(PW_22)
set _TMP(PW_23) [$_CN(7) getDistribution 1]
$_TMP(PW_23) setEndSpacing 0.00050000000000000001
unset _TMP(PW_23)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(8)]]
set _TMP(PW_24) [$_CN(8) getDistribution 1]
$_TMP(PW_24) setBeginSpacing 0.01
unset _TMP(PW_24)
set _TMP(PW_25) [$_CN(8) getDistribution 1]
$_TMP(PW_25) setEndSpacing 0.01
unset _TMP(PW_25)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(5) $_CN(4)]]
set _TMP(PW_26) [$_CN(4) getDistribution 1]
$_TMP(PW_26) setBeginSpacing 0.00012
unset _TMP(PW_26)
set _TMP(PW_27) [$_CN(5) getDistribution 1]
$_TMP(PW_27) setBeginSpacing 0.00012
unset _TMP(PW_27)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(7) $_CN(6) $_CN(5) $_CN(4)]]
set _TMP(PW_28) [$_CN(4) getDistribution 1]
$_TMP(PW_28) setEndSpacing 1.05
unset _TMP(PW_28)
set _TMP(PW_29) [$_CN(6) getDistribution 1]
$_TMP(PW_29) setBeginSpacing 1.05
unset _TMP(PW_29)
set _TMP(PW_30) [$_CN(5) getDistribution 1]
$_TMP(PW_30) setEndSpacing 1.05
unset _TMP(PW_30)
set _TMP(PW_31) [$_CN(7) getDistribution 1]
$_TMP(PW_31) setBeginSpacing 1.05
unset _TMP(PW_31)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_10) [pw::Application begin Modify [list $_CN(3)]]
set _TMP(PW_32) [$_CN(3) getDistribution 1]
$_TMP(PW_32) setEndSpacing 1.05
unset _TMP(PW_32)
$_TMP(mode_10) end
unset _TMP(mode_10)
pw::Application markUndoLevel {Change Spacing}

set _TMP(mode_10) [pw::Application begin Create]
set _TMP(edge_1) [pw::Edge create]
$_TMP(edge_1) addConnector $_CN(4)
set _TMP(edge_2) [pw::Edge create]
$_TMP(edge_2) addConnector $_CN(6)
$_TMP(edge_2) addConnector $_CN(8)
$_TMP(edge_2) addConnector $_CN(7)
set _TMP(edge_3) [pw::Edge create]
$_TMP(edge_3) addConnector $_CN(5)
set _TMP(edge_4) [pw::Edge create]
$_TMP(edge_4) addConnector $_CN(3)
$_TMP(edge_4) addConnector $_CN(2)
$_TMP(edge_4) addConnector $_CN(1)
$_TMP(edge_4) addConnector $_CN(3)
set _TMP(dom_1) [pw::DomainStructured create]
$_TMP(dom_1) addEdge $_TMP(edge_1)
$_TMP(dom_1) addEdge $_TMP(edge_2)
$_TMP(dom_1) addEdge $_TMP(edge_3)
$_TMP(dom_1) addEdge $_TMP(edge_4)
unset _TMP(edge_3)
unset _TMP(edge_2)
unset _TMP(edge_1)
$_TMP(mode_10) end
unset _TMP(mode_10)
unset _TMP(dom_1)
pw::Application markUndoLevel {Assemble Domain}

unset _TMP(edge_4)
set _TMP(mode_10) [pw::Application begin Create]
$_TMP(mode_10) abort
unset _TMP(mode_10)
pw::Application setCAESolver {CGNS} 2
pw::Application markUndoLevel {Set Dimension 2D}

pw::Application setCAESolver {ANSYS FLUENT} 2
pw::Application markUndoLevel {Select Solver}

set _TMP(PW_33) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

$_TMP(PW_33) setPhysicalType {Wall}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_33) apply [list $_CN(1) $_CN(2)]
pw::Application markUndoLevel {Set BC}

set _TMP(PW_34) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

$_TMP(PW_34) setPhysicalType {Pressure Far Field}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_34) apply [list $_CN(6) $_CN(8) $_CN(7) $_CN(5) $_CN(4)]
pw::Application markUndoLevel {Set BC}

unset _TMP(PW_33)
unset _TMP(PW_34)
set _DM(1) [pw::GridEntity getByName {dom-1}]
pw::Application export [list $_DM(1)] {D:/codes/pyCodes/VCM/temp/output_ansys_file.cas}
