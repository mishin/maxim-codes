package require PWI_Glyph 2.3
pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}
pw::Application clearModified

set _TMP(mode_1) [pw::Application begin Create]
set _TMP(PW_1) [pw::SegmentSpline create]
$_TMP(PW_1) addPoint {0.000000 0.000000 0}
$_TMP(PW_1) addPoint {1.000000 0.000000 0}
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
$_TMP(PW_2) addPoint {0.444430 -0.053500 0}
$_TMP(PW_2) addPoint [$_CN(1) getPosition -arc 1]
$_TMP(PW_2) setSlope Akima
set _TMP(con_2) [pw::Connector create]
$_TMP(con_2) addSegment $_TMP(PW_2)
$_TMP(con_2) calculateDimension
unset _TMP(PW_2)
$_TMP(mode_2) end
unset _TMP(mode_2)
pw::Application markUndoLevel {Create Connector}

set _CN(1) [pw::GridEntity getByName {con-2}]
set _CN(2) [pw::GridEntity getByName {con-1}]
set _TMP(PW_3) [pw::Collection create]
$_TMP(PW_3) set [list $_CN(1) $_CN(2)]
$_TMP(PW_3) do setDimension 75
$_TMP(PW_3) delete
unset _TMP(PW_3)
pw::Application markUndoLevel {Dimension}

set _TMP(mode_3) [pw::Application begin Modify [list $_CN(2) $_CN(1)]]
set _TMP(PW_4) [$_CN(2) getDistribution 1]
$_TMP(PW_4) setBeginSpacing 0.001
unset _TMP(PW_4)
set _TMP(PW_5) [$_CN(1) getDistribution 1]
$_TMP(PW_5) setBeginSpacing 0.001
unset _TMP(PW_5)
$_TMP(mode_3) end
unset _TMP(mode_3)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_4) [pw::Application begin Modify [list $_CN(2) $_CN(1)]]
set _TMP(PW_6) [$_CN(2) getDistribution 1]
$_TMP(PW_6) setEndSpacing 0.00050000000000000001
unset _TMP(PW_6)
set _TMP(PW_7) [$_CN(1) getDistribution 1]
$_TMP(PW_7) setEndSpacing 0.00050000000000000001
unset _TMP(PW_7)
$_TMP(mode_4) end
unset _TMP(mode_4)
pw::Application markUndoLevel {Change Spacings}

set _TMP(mode_5) [pw::Application begin Create]
set _TMP(PW_8) [pw::Edge createFromConnectors [list $_CN(1) $_CN(2)]]
set _TMP(edge_2) [lindex $_TMP(PW_8) 0]
unset _TMP(PW_8)
$_TMP(edge_2) reverse
set _TMP(dom_2) [pw::DomainStructured create]
$_TMP(dom_2) addEdge $_TMP(edge_2)
$_TMP(mode_5) end
unset _TMP(mode_5)
set _TMP(mode_6) [pw::Application begin ExtrusionSolver [list $_TMP(dom_2)]]
set _DM(1) [pw::GridEntity getByName {dom-1}]
$_DM(1) setExtrusionSolverAttribute NormalInitialStepSize 3.2e-4
$_TMP(mode_6) run -keepFailingStep 85
$_TMP(mode_6) end
unset _TMP(mode_6)
unset _TMP(dom_2)
unset _TMP(edge_2)
pw::Application markUndoLevel {Extrude, Normal}

pw::Application setCAESolver {ANSYS FLUENT} 3
pw::Application markUndoLevel {Select Solver}

pw::Application setCAESolver {ANSYS FLUENT} 2
pw::Application markUndoLevel {Set Dimension 2D}

set _TMP(PW_9) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

$_TMP(PW_9) setPhysicalType {Wall}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_9) apply [list $_CN(2) $_CN(1)]
pw::Application markUndoLevel {Set BC}

set _TMP(PW_10) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

$_TMP(PW_10) setPhysicalType {Pressure Far Field}
pw::Application markUndoLevel {Change BC Type}

set _CN(3) [pw::GridEntity getByName {con-4}]
$_TMP(PW_10) apply [list $_CN(3)]
pw::Application markUndoLevel {Set BC}

unset _TMP(PW_9)
unset _TMP(PW_10)
pw::Application export [list $_DM(1)] {D:/codes/pyCodes/VCM/temp/tmp_debug_file.cas}
