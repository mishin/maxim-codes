# -*- coding: utf-8 -*-
"""
Created on Tue May 14 19:39:12 2013

@author: Maxim
"""
from numpy import zeros, array, vstack
from math import sin, cos, tan, radians
from win32com.client import Dispatch
from airfoil_cad import Airfoil3D
from os.path import abspath

class CadDesigner():
    def __init__(self,visible=True):
        self.catia = Dispatch('catia.application')
        self.catia.Visible = visible
        self.tmpReferences = list()
        
#    def measure_points(self):
#        spaw = self.catia.ActiveDocument.GetWorkbench("SPAWorkbench")
#        for ref in self.tmpReferences:
#            spaw.GetMeasurable(ref).GetPoint(coord)
#            print coord

    def add_new_part(self):
        self.partDocument1 = self.catia.Documents.Add("Part")
        self.partDocument1 = self.catia.ActiveDocument
        self.part = self.partDocument1.Part
        self.hSF = self.part.HybridShapeFactory
        self.hBodies = self.part.HybridBodies
        self.hBody = self.hBodies.Item("Geometrical Set.1")
        self.refPoint       = zeros(3)
        self.chordRoot      = 0.0
        self.secRefRoot     = None
        self.secRefTip      = None
        self.wingSurfaces   = list()
        self.wingSurface    = None
        self.elementsToHide = list()
        self.surfacesToHide = list()
        self.wingTipContour = None

    def add_new_product(self):
        documents = self.catia.Documents
        self.productDocument1 = documents.Add("Product")
        self.productDocument1 = self.catia.ActiveDocument
        self.products = self.productDocument1.Product.Products
        n = 11 # number of parts
        self.partsLocation   = zeros([n,3],dtype=float)
    
    def add_part_to_assembly(self,path,i):
        path = [abspath(path)]
        self.products.AddComponentsFromFiles(path,"All")

    def move_part(self,i,newLocation):
        product = self.products.Item(i)
        move = product.Move.MovableObject
        increment = newLocation - self.partsLocation[i-1]
        dist = zeros(12)
        dist[0] = 1
        dist[4] = 1
        dist[8] = 1
        dist[9] = increment[0]
        dist[10] = increment[1]
        dist[11] = increment[2]
        move.Apply(dist)
        self.partsLocation[i-1] = newLocation
    
    def create_canopy(self,upCurve,sideCurve,csLocation,csConic):
        # top curve
        refUpCurve = self.create_spline(upCurve,show=False)
        # side curve
        
        refSideCurve = self.create_spline(sideCurve,show=False,
                                          pointRef=[refUpCurve[0],refUpCurve[1]])
        # side curve symmetry
        refSideSym = self.symmetry(refSideCurve[2],'xz',show=False)
        # cross section conic
        crossSection = list()
        for loc,rho in zip(csLocation,csConic):
            refCs = self.create_cross_section(refUpCurve[2], refSideCurve[2], loc,rho)
            refSymCs = self.symmetry(refCs,show=False)
            cs = self.create_join([refCs,refSymCs],show=False)
            crossSection.append(cs)
        
        # surface
        surf = self.create_surface(crossSection, [refSideSym, refUpCurve[2], refSideCurve[2]], show=True, relimitation=2)
        # split
    
    def create_cross_section(self,refUpCurve, refSideCurve, offset, rho=0.5,show=False):
        originElements = self.part.OriginElements
        refPlane = originElements.PlaneYZ
        refPlane = self.part.CreateReferenceFromObject(refPlane)
        hSPlaneOffset = self.hSF.AddNewPlaneOffset(refPlane,offset,False)
        hSPlaneOffset = self.part.CreateReferenceFromObject(hSPlaneOffset)
        startPt = self.create_intersection(hSPlaneOffset,refUpCurve)
        endPt = self.create_intersection(hSPlaneOffset,refSideCurve)
        startDir = self.create_direction(0,1,0)
        endDir = self.create_direction(0,0,-1)
        conic = self.hSF.AddNewConic(hSPlaneOffset, startPt, endPt)
        conic.SetStartAndEndTangentsPlusConicParameter(startDir[0], endDir[0], rho)
        if show:
            self.hBody.AppendHybridShape(conic)
            self.part.update()
        return self.part.CreateReferenceFromObject(conic)

    def create_direction(self,xdir=1.0, ydir=0.0, zdir=0.0):
        hsDirection = self.hSF.AddNewDirectionByCoord(xdir,ydir,zdir)
        return hsDirection, self.part.CreateReferenceFromObject(hsDirection)
    
    def create_intersection(self, ref1, ref2, show=False):
        hsIntersect = self.hSF.AddNewIntersection(ref1,ref2)
        hsIntersect.PointType = 0
        if show:
            self.hBody.AppendHybridShape(hsIntersect)
            self.part.update()
        return self.part.CreateReferenceFromObject(hsIntersect)

    def create_airfoil(self,af2d,chord,incidence,rotAxis,spline=True):
        self.af3d = Airfoil3D(af2d,'xz',chord)
        self.chordRoot = chord
        upPts = self.af3d.upPts + self.refPoint
        loPts = self.af3d.loPts + self.refPoint
        if spline:
            if af2d.leRad != 0.0:
                refCirc = self.create_circle(False)
                refStartPt1,refEndPt1,refSpline1 = self.create_spline(upPts,curvature=[refCirc,None],show=False)
                refStartPt2,refEndPt2,refSpline2 = self.create_spline(loPts,curvature=[refCirc,None],
                                                orientation=[-1,1],show=False)
            else:
                refStartPt1,refEndPt1,refSpline1 = self.create_spline(upPts,show=False)
                refStartPt2,refEndPt2,refSpline2 = self.create_spline(loPts,show=False)
        else:
                refStartPt1,refEndPt1,refSpline1 = self.create_polyline(upPts,show=False)
                refStartPt2,refEndPt2,refSpline2 = self.create_polyline(loPts,show=False)
        #print self.af3d.teGap
        if self.af3d.teGap>0.0:
            refTEline = self.line_ptpt_ref(refEndPt1,refEndPt2,show=False)
            refSection = self.create_join([refSpline1,refSpline2,refTEline],False)
        else:
            refSection = self.create_join([refSpline1,refSpline2],False)
        rotAxisPt = self.refPoint + array([1,0,0])*rotAxis*chord
        rotAxisLine = self.line_pt_dir(rotAxisPt,[0,1,0])
        refStartPt1 = self.rotate(refStartPt1,rotAxisLine,incidence,False)
        refStartPt2 = self.rotate(refStartPt2,rotAxisLine,incidence,False)
        refEndPt1 = self.rotate(refEndPt1,rotAxisLine,incidence,False)
        refEndPt2 = self.rotate(refEndPt2,rotAxisLine,incidence,False)
        refSection = self.rotate(refSection,rotAxisLine,incidence,True,True)
        self.secRefRoot = self.secRefTip
        self.secRefTip = refStartPt1,refStartPt2,refEndPt1,refEndPt2,refSection

    def create_tip_section(self,af2d,chord,incidence,rotAxis,span,sweep,dihedral,c4sweep=False,spline=True):
        if c4sweep:
            self.refPoint[0] += (self.chordRoot-chord)/4.0 + span*tan(radians(sweep))
        else:
            self.refPoint[0] += span*tan(radians(sweep))
        self.refPoint[1] += span
        self.refPoint[2] += span*tan(radians(dihedral))
        self.create_airfoil(af2d,chord,incidence,rotAxis,spline)
        lineLE = self.line_ptpt_ref(self.secRefRoot[0],self.secRefTip[0],False)
        lineTE1 = self.line_ptpt_ref(self.secRefRoot[2],self.secRefTip[2],False)
        lineTE2 = self.line_ptpt_ref(self.secRefRoot[3],self.secRefTip[3],False)
        if self.af3d.teGap!=0.0:
            guide = [lineTE1,lineLE,lineTE2]
        else:
            guide = [lineTE1,lineLE]
        surface = self.create_surface([self.secRefRoot[4],self.secRefTip[4]],guide,True)
        self.wingSurfaces.append(surface)
    
    def close_wing_tip(self):
        tipContour = self.secRefTip[4]
        newFill = self.hSF.AddNewFill()
        newFill.AddBound(tipContour)
        newFill.Continuity = 0
        self.wingSurfaces.append(newFill)
        self.surfacesToHide.append(newFill)
        self.hBody.AppendHybridShape(newFill)
        self.part.update()

    def create_xtail(self,angle):
        axis = self.line_pt_dir([0,0,0],[1,0,0])
        rotate1 = self.rotate(self.wingSurface,axis,angle,True)
        for i in range(3):
            self.rotate(rotate1,axis,(i+1)*90.0,True)
    
    def hide_curves(self):
        self.hide_object(self.elementsToHide)
    def hide_surfaces(self):
        self.hide_object(self.surfacesToHide)

    def join_wing(self,hide=False):
        if self.wingSurface==None:
            if len(self.wingSurfaces)>1:
                self.wingSurface = self.create_join(self.wingSurfaces,hideSurf=hide)
            elif len(self.wingSurfaces)==1:
                self.wingSurface = self.wingSurfaces[0]
            if hide:
                self.hide_surfaces()
    
    def rotate_wing(self,angle,axisRoot,axisDir=[0,0,-1]):
        axis = self.line_pt_dir([axisRoot[0],axisRoot[1],0],axisDir)
        self.wingSurface = self.rotate(self.wingSurface,axis,angle)
    
    def rot_sym_cfd(self):
        axis = self.line_pt_dir([0,0,0],[1,0,0])
        angle = -90
        self.join_wing(True)
        self.wingSurface = self.rotate(self.wingSurface,axis,angle,show=False)
        self.symmetry(self.wingSurface,'xy',True)
        self.hide_curves()
    
    def symmetry_wing(self,surf='xz'):
        self.join_wing(True)
        self.symmetry(self.wingSurface,surf,True)
        self.hide_curves()
    
    def symmetry_wing_cfd(self):
        self.rotate_wing(90.0,array([1.,0,0]))
        self.symmetry_wing('xy')

    def close_part(self):
        self.partDocument1.close()
    def close_product(self):
        self.productDocument1.close()
    def save_part(self,path):
        self.partDocument1.SaveAs(path)
    def save_product(self,path):
        self.productDocument1.SaveAs(path)
    def save_part_igs(self,path):
        self.partDocument1.ExportData(path,'igs')
    
    def create_point(self,point):
        hspc = self.hSF.AddNewPointCoord(point[0],point[1],point[2])
        return self.part.CreateReferenceFromObject(hspc)
    
    def create_polyline(self,curve,show=True):
        polyline = self.hSF.AddNewPolyline()
        tmpRef = list()
        for i,point in enumerate(curve):
            ptRef = self.create_point(point)
            tmpRef.append(ptRef)
            polyline.InsertElement(ptRef,i+1)
        polyline.Closure = False
        if show:
            self.hBody.AppendHybridShape(polyline)
            self.part.update()
        return tmpRef[0],tmpRef[-1],polyline

    def create_spline(self,curve,tangency=[None,None],curvature=[None,None],
                      orientation=[1,1],show=True,pointRef=[None,None]):
        spline = self.hSF.AddNewSpline()
        if pointRef[0]==None:
            ptStart = self.create_point(curve[0])
        else:
            ptStart = pointRef[0]
        if pointRef[1]==None:
            ptEnd = self.create_point(curve[-1])
        else:
            ptEnd = pointRef[1]
        if not tangency[0]==None:
            spline.AddPointWithConstraintFromCurve(ptStart,tangency[0],1,orientation[0],1)
        elif curvature[0]!=None:
            spline.AddPointWithConstraintFromCurve(ptStart,curvature[0],1,orientation[0],2)
        else:
            spline.AddPointWithConstraintExplicit(ptStart, None, -1, 1, None, 0)
        
        if not curve==None:
            for point in curve[1:-1]:
                ptRef = self.create_point(point)
                #FIXME tmp reader
                self.tmpReferences.append(ptRef)
                spline.AddPointWithConstraintExplicit(ptRef, None, -1, 1, None, 0)
        
        if tangency[1]!=None:
            spline.AddPointWithConstraintFromCurve(ptEnd,tangency[1],1,orientation[1],1)
        elif curvature[1]!=None:
            spline.AddPointWithConstraintFromCurve(ptEnd,curvature[1],1,orientation[1],2)
        else:
            spline.AddPointWithConstraintExplicit(ptEnd, None, -1, 1, None, 0)
        if show:
            self.hBody.AppendHybridShape(spline)
            self.part.update()
        refSpline = self.part.CreateReferenceFromObject(spline)
        return ptStart, ptEnd, refSpline
    
    def create_circle(self,show=True):
        originElements = self.part.OriginElements
        if self.af3d.plane=='xy':
            refPlane = originElements.PlaneXY
        elif self.af3d.plane=='xz':
            refPlane = originElements.planeZX
        tmp = self.af3d.leCenter + self.refPoint
        centerPt = self.hSF.AddNewPointCoord(tmp[0],0,tmp[2])
        tmp = self.af3d.upPts[0] + self.refPoint
        lePt = self.hSF.AddNewPointCoord(tmp[0],0,tmp[2])
        ref1 = self.part.CreateReferenceFromObject(centerPt)
        ref2 = self.part.CreateReferenceFromObject(lePt)
        ref3 = self.part.CreateReferenceFromObject(refPlane)
        leCircle = self.hSF.AddNewCircleCtrPt(ref1,ref2,ref3,False)
        leCircle.SetLimitation(1)
        ref4 = self.part.CreateReferenceFromObject(leCircle)
        if show:
            self.hBody.AppendHybridShape(leCircle)
            self.part.update()
        return ref4

    def translate(self):
        pass
    
    def symmetry(self,refElemToRotate, symPlane='xz',show=True):
        originElements = self.part.OriginElements
        if symPlane=='xz':
            refPlane = originElements.PlaneZX
        elif symPlane=='xy':
            refPlane = originElements.PlaneXY
        refPlane = self.part.CreateReferenceFromObject(refPlane)
        sym = self.hSF.AddNewSymmetry(refElemToRotate,refPlane)
        sym.VolumeResult = False
        if show:
            self.hBody.AppendHybridShape(sym)
            self.part.update()
        return self.part.CreateReferenceFromObject(sym)


    def line_pt_dir(self,point=[0,0,0],direction=[0,1,0],length=20.0):
        originPt = self.create_point(point)
        refDir = self.hSF.AddNewDirectionByCoord(direction[0],direction[1],direction[2])
        line = self.hSF.AddNewLinePtDir(originPt,refDir,0,length,False)
        return self.part.CreateReferenceFromObject(line)
    
    def line_ptpt_ref(self,point1,point2,show=True):
        line = self.hSF.AddNewLinePtPt(point1,point2)
        if show:
            self.hBody.AppendHybridShape(line)
            self.part.update()
        return self.part.CreateReferenceFromObject(line)
    
    def create_surface(self,sections,guides=None,show=True, relimitation=1):
        newLoft = self.hSF.AddNewLoft()
        newLoft.SectionCoupling = 1
        newLoft.Relimitation = relimitation
        newLoft.CanonicalDetection = 2
        for sec in sections:
            newLoft.AddSectionToLoft(sec,1,None)
        if not guides==None:
            for guide in guides:
                newLoft.AddGuide(guide)
        if show:
            self.hBody.AppendHybridShape(newLoft)
            self.part.update()
        self.surfacesToHide.append(newLoft)
        return self.part.CreateReferenceFromObject(newLoft)

    def rotate(self,refElementToRotate,rotAxisLine,angle,show=True,hide=False):
        newRotate = self.hSF.AddNewEmptyRotate()
        newRotate.ElemToRotate = refElementToRotate
        newRotate.VolumeResult = False
        newRotate.RotationType = 0
        newRotate.Axis = rotAxisLine
        newRotate.AngleValue = angle
        if show:
            self.hBody.AppendHybridShape(newRotate)
            self.part.update()
        if hide:
            self.elementsToHide.append(newRotate)
        return self.part.CreateReferenceFromObject(newRotate)
    
    def create_join(self,elementsToJoin,show=True,hide=False,hideSurf=False):
        newJoin = self.hSF.AddNewJoin(elementsToJoin[0],elementsToJoin[1])
        if len(elementsToJoin)>2:
            for element in elementsToJoin[2:]:
                newJoin.AddElement(element)
        newJoin.SetConnex(1)
        newJoin.SetManifold(1)
        newJoin.SetSimplify(0)
        newJoin.SetSuppressMode(0)
        newJoin.SetDeviation(0.001)
        newJoin.SetAngularToleranceMode(0)
        newJoin.SetAngularTolerance(0.5)
        newJoin.SetFederationPropagation(0)
        if show:
            self.hBody.AppendHybridShape(newJoin)
            self.part.update()
        if hide:
            self.elementsToHide.append(newJoin)
        if hideSurf:
            self.surfacesToHide.append(newJoin)
        return self.part.CreateReferenceFromObject(newJoin)

    def hide_object(self,elementsToHide):
        partDocument = self.catia.ActiveDocument
        selection = partDocument.Selection
        for element in elementsToHide:
            selection.Add(element)
        visProperty = selection.VisProperties
        visProperty.SetShow(1)
        selection.clear()
    
    def circle_ptpt_rad(self,point1,point2,radius,show=True):
        xyPlane = self.part.OriginElements.PlaneZX
        refPlane = self.part.CreateReferenceFromObject(xyPlane)
        newCircle = self.hSF.AddNewCircle2PointsRad(point1,point2,refPlane,False,radius,1)
        if show:
            self.hBody.AppendHybridShape(newCircle)
            self.part.update()
        return self.part.CreateReferenceFromObject(newCircle)
    
    def create_revolute_x(self,refProfile,show=True):
        axis = self.line_pt_dir(direction=[1,0,0])
        newRevolute = self.hSF.AddNewRevol(refProfile,360,0,axis)
        if show:
            self.hBody.AppendHybridShape(newRevolute)
            self.part.update()
        return self.part.CreateReferenceFromObject(newRevolute)

    def create_body(self, radiusNose,angleNose,contNose,radiusTail,angleTail,
                            contTail,xoffset,coord):
        profile = list()
        tangency  = [None,None]
        curvature = [None,None]
        edgePoints = [None,None]
        if radiusNose!=0.0:
            lePt = self.create_point([0,0,0])
            curveStartPtCoord = array([1.0-cos(radians(angleNose)),0,sin(radians(angleNose))]) * radiusNose
            curveStartPt = self.create_point(curveStartPtCoord)
            noseCircle = self.circle_ptpt_rad(lePt,curveStartPt,radiusNose)
            profile.append(noseCircle)
            edgePoints[0] = curveStartPt
            if contNose=='C1':
                tangency[0] = noseCircle
            elif contNose=='C2':
                curvature[0] = noseCircle
        else:
            if coord[0]!=[0,0]:
                coord = vstack([zeros(3),coord])
        
        if radiusTail!=0.0:
            tePt = self.create_point([xoffset+radiusNose+radiusTail, 0,0])
            curveEndPt = array([cos(radians(angleTail)),0,sin(radians(angleTail))]) * radiusTail
            curveEndPt[0] += radiusNose + xoffset
            curveEndPt = self.create_point(curveEndPt)
            tailCircle = self.circle_ptpt_rad(curveEndPt,tePt,radiusTail)
            edgePoints[1] = curveEndPt
            if contTail=='C1':
                tangency[1] = tailCircle
            elif contTail=='C2':
                curvature[1] = tailCircle
        splineRef = self.create_spline(coord,tangency,curvature,pointRef=edgePoints)
        profile.append(splineRef[2])
        if radiusTail==0.0 and coord[-1,2]!=0.0:
            pointEnd = coord[-1]
            pointEnd[2] = 0.0
            pointEnd = self.create_point(pointEnd)
            lineEnd = self.line_ptpt_ref(splineRef[1],pointEnd)
            profile.append(lineEnd)
        else:
            profile.append(tailCircle)
        profileCurve = self.create_join(profile)
        self.create_revolute_x(profileCurve)
    
    def nose_conic(self,pts,isTail):
        pt1 = self.create_point(pts[0])
        pt2 = self.create_point(pts[1])
        line = self.line_ptpt_ref(pt1,pt2)
        return line

    def nose_conic_blunt(self,pts,Rnose,isTail):
        refPt = [self.create_point(pt) for pt in pts]
        if isTail:
            line = self.line_ptpt_ref(refPt[0],refPt[1],show=False)
            circle = self.circle_ptpt_rad(refPt[1],refPt[2],Rnose,show=False)
        else:
            circle = self.circle_ptpt_rad(refPt[0],refPt[1],Rnose,show=False)
            line = self.line_ptpt_ref(refPt[1],refPt[2],show=False)
        profile = self.create_join([circle,line],hide=True)
        return profile
    
    def nose_conic_flat(self,pts,isTail):
        refPt = [self.create_point(pt) for pt in pts]
        line1 = self.line_ptpt_ref(refPt[0],refPt[1],show=False)
        line2 = self.line_ptpt_ref(refPt[1],refPt[2],show=False)
        profile = self.create_join([line1,line2],hide=True)
        return profile

    def nose_ogival(self,pts,rho):
        refPt = [self.create_point(pt) for pt in pts]
        circle = self.circle_ptpt_rad(refPt[0],refPt[1],rho)
        return circle
    
    def nose_ogival_blunt(self,pts,rho,Rnose,isTail):
        refPt = [self.create_point(pt) for pt in pts]
        if isTail:
            r = [rho,Rnose]
        else:
            r = [Rnose,rho]
        circle1 = self.circle_ptpt_rad(refPt[0],refPt[1],r[0],show=False)
        circle2 = self.circle_ptpt_rad(refPt[1],refPt[2],r[1],show=False)
        profile = self.create_join([circle1,circle2],hide=True)
        return profile
        
    def nose_ogival_flat(self,pts,rho,isTail):
        refPt = [self.create_point(pt) for pt in pts]
        if isTail:
            circle = self.circle_ptpt_rad(refPt[0],refPt[1],rho,show=False)
            line = self.line_ptpt_ref(refPt[1],refPt[2],show=False)
        else:
            line = self.line_ptpt_ref(refPt[0],refPt[1],show=False)
            circle = self.circle_ptpt_rad(refPt[1],refPt[2],rho,show=False)
        profile = self.create_join([circle,line],hide=True)
        return profile
            
    
    def create_pt_on_curve_along_direction(self,refCurve,distance,vector=[1,0,0],show=False):
        #FIXME: works strange
        direction, refDirection = self.create_direction(vector[0],vector[1],vector[2])
        pt = self.hSF.AddNewPointOnCurveAlongDirection(refCurve,distance,False,direction)
        if show:
            self.hBody.AppendHybridShape(pt)
            self.part.update()
        self.part.CreateReferenceFromObject(pt)
    
    def create_pt_on_curve_percentage(self,refCurve,ratio,show=False):
        pt = self.hSF.AddNewPointOnCurveFromPercent(refCurve,ratio,False)
        if show:
            self.hBody.AppendHybridShape(pt)
            self.part.update()
        return self.part.CreateReferenceFromObject(pt)
    
    def split_curve(self,refCurve,refPt,side=1,show=False):
        split = self.hSF.AddNewHybridSplit(refCurve, refPt, side)
        if show:
            self.hBody.AppendHybridShape(split)
            self.part.update()
        return self.part.CreateReferenceFromObject(split)