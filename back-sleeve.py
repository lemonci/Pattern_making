# Back Centre Length
A = 365
# Breast Circumference - 85mm
B = 925
 
BS = 0
WL = 0
BC = A
BL = A-B/12-137
CFL = B/2+42.5
BWL = B/8+74
BY = A-80
GL = (BL+BY)/2-5
FC = BL+B/5+83
CWL = CFL-(B/8+62)
FAW = CWL-B/32
SS = (BWL+FAW)/2
 
ZeroX = SS - BS
ZeroY = BL - WL
OneX = BWL - BS
OneY = BY - BL
TwoX = BWL - BS
TwoY = BC - BY
ThreeX = FAW - BWL
ThreeY = GL - BL
FourX = CFL - SS
FourY = BL - WL
FiveX = CFL - CWL
FiveY = FC - BL
Gui.activateWorkbench("DraftWorkbench")
App.newDocument("back-sleeve")
App.setActiveDocument("back-sleeve")
App.ActiveDocument=App.getDocument("back-sleeve")
Gui.ActiveDocument=Gui.getDocument("back-sleeve")
Gui.activeDocument().activeView().viewDefaultOrientation()
Gui.activateWorkbench("DraftWorkbench")
import Draft
from math import atan, sin, cos, tan, pi
pl = FreeCAD.Placement()
pl.Rotation.Q = (0.0,0.0,0.0,1.0)

'''#Rectangle
pl.Base = FreeCAD.Vector(BS,WL,0.0)
rec = Draft.makeRectangle(length=ZeroX,height=ZeroY,placement=pl,face=False,support=None)
Draft.autogroup(rec)
#Rectangle001
pl.Base = FreeCAD.Vector(BS,BL,0.0)
rec = Draft.makeRectangle(length=OneX,height=OneY,placement=pl,face=False,support=None)
Draft.autogroup(rec)
#Rectangle002
pl.Base = FreeCAD.Vector(BS,BY,0.0)
rec = Draft.makeRectangle(length=TwoX,height=TwoY,placement=pl,face=False,support=None)
Draft.autogroup(rec)
#Rectangle003
pl.Base = FreeCAD.Vector(BWL,BL,0.0)
rec = Draft.makeRectangle(length=ThreeX,height=ThreeY,placement=pl,face=False,support=None)
#Rectangle004
pl.Base = FreeCAD.Vector(SS,WL,0.0)
rec = Draft.makeRectangle(length=FourX,height=FourY,placement=pl,face=False,support=None)
#Rectangle005
pl.Base = FreeCAD.Vector(CWL,BL,0.0)
rec = Draft.makeRectangle(length=FiveX,height=FiveY,placement=pl,face=False,support=None)'''
#point E
pointEX = BWL/2+1
#point = Draft.makePoint(pointEX, BY)
#Draft.autogroup(point)
#point BP
pointBPX = (CWL+CFL)/2-7
#point = Draft.makePoint(pointBPX, BL)
#Draft.autogroup(point)
#Front collar
collarRef = B/24 + 34
frontCollarX = CFL - collarRef
frontCollarY = FC - (collarRef + 5)
#points = [FreeCAD.Vector(frontCollarX,FC),FreeCAD.Vector(frontCollarX,frontCollarY),FreeCAD.Vector(CFL,frontCollarY)]
#bez = Draft.makeBezCurve(points,closed=False,support=None)
#Draft.autogroup(bez)
#Front shoulder line
frontShoulderLine = (CFL-CWL-collarRef)/cos(22/360*2*pi)+18
frontShoulderLineX0 = CWL - 18*cos(22/360*2*pi)
frontShoulderLineY0 = FC - frontShoulderLine*sin(22/360*2*pi)
#pl.Base = FreeCAD.Vector(frontCollarX,FC,0.0)
#points = [FreeCAD.Vector(frontCollarX,FC,0.0),FreeCAD.Vector(frontShoulderLineX0,frontShoulderLineY0,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#back collar
backCollarX1 = BS + (collarRef + 2)/3*2
backCollarX2 = BS + (collarRef + 2)
backCollarY2 = BC + (collarRef + 2)/3
points = [FreeCAD.Vector(BS,BC),FreeCAD.Vector(backCollarX1,BC),FreeCAD.Vector(backCollarX2,backCollarY2)]
bez = Draft.makeBezCurve(points,closed=False,support=None)
#back shoulder line
backShoulderLine = frontShoulderLine + (B/32 - 8)
backShoulderLineX0 = backCollarX2 + backShoulderLine*cos(18/360*2*pi)
backShoulderLineY0 = backCollarY2 - backShoulderLine*sin(18/360*2*pi)
#pl.Base = FreeCAD.Vector(backCollarX2,backCollarY2,0.0)
#points = [FreeCAD.Vector(backCollarX2,backCollarY2,0.0),FreeCAD.Vector(backShoulderLineX0,backShoulderLineY0,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#armhole
backArmHoleX0 = BWL + ((SS-BWL)/3+8)/(2**0.5)
backArmHoleY0 = BL + ((SS-BWL)/3+8)/(2**0.5)
frontArmHoleX0 = FAW - ((FAW-SS)/3+5)/(2**0.5)
frontArmHoleY0 = BL + ((FAW-SS)/3+5)/(2**0.5)
#full armhole
#points = [FreeCAD.Vector(backShoulderLineX0,backShoulderLineY0),FreeCAD.Vector(BWL+0.1,GL+7),FreeCAD.Vector(BWL+0.05,GL+5),FreeCAD.Vector(BWL,GL),FreeCAD.Vector(backArmHoleX0,backArmHoleY0),FreeCAD.Vector(SS,BL),FreeCAD.Vector(frontArmHoleX0,frontArmHoleY0),FreeCAD.Vector(FAW,GL-0.1),FreeCAD.Vector(FAW,GL),FreeCAD.Vector(FAW,GL+0.1)]
#back half armhole
#points = [FreeCAD.Vector(backShoulderLineX0,backShoulderLineY0),FreeCAD.Vector(BWL+0.1,GL+7),FreeCAD.Vector(BWL+0.05,GL+5),FreeCAD.Vector(BWL,GL),FreeCAD.Vector(backArmHoleX0,backArmHoleY0),FreeCAD.Vector(SS,BL)]
#spline = Draft.makeBSpline(points,closed=False,face=True,support=None)
#breast dart lower line
#pl.Base = FreeCAD.Vector(FAW,GL,0.0)
#points = [FreeCAD.Vector(FAW,GL,0.0),FreeCAD.Vector(pointBPX, BL,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#breast dart upper line
breastDart = ((pointBPX-FAW)**2+(BL-GL)**2)**0.5
breastDartAngle = atan((GL-BL)/(pointBPX-FAW)+(B/4-25)/3600*2*pi)
breastDartX0 = pointBPX - breastDart*sin(breastDartAngle)
breastDartY0 = BL + breastDart*cos(breastDartAngle)
#pl.Base = FreeCAD.Vector(breastDartX0,breastDartY0,0.0)
#points = [FreeCAD.Vector(breastDartX0,breastDartY0,0.0),FreeCAD.Vector(pointBPX, BL,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#back dart upperline
backDartUpperX = pointEX + 15*cos(18/360*2*pi)
backDartUpperY = backCollarY2 - (backDartUpperX-backCollarX2)*sin(18/360*2*pi)
pl.Base = FreeCAD.Vector(backDartUpperX,backDartUpperY,0.0)
points = [FreeCAD.Vector(backDartUpperX,backDartUpperY,0.0),FreeCAD.Vector(0.6*pointEX + 0.4*backDartUpperX,0.6*BY + 0.4*backDartUpperY,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#back shoulder line upper part
pl.Base = FreeCAD.Vector(backDartUpperX,backDartUpperY,0.0)
points = [FreeCAD.Vector(backDartUpperX,backDartUpperY,0.0),FreeCAD.Vector(backCollarX2,backCollarY2,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#back dart lowerline
backDartLowerX = backDartUpperX + (B/32-8)*cos(18/360*2*pi)
backDartLowerY = backDartUpperY - (B/32-8)*sin(18/360*2*pi)
pl.Base = FreeCAD.Vector(backDartLowerX,backDartLowerY,0.0)
points = [FreeCAD.Vector(backDartLowerX,backDartLowerY,0.0),FreeCAD.Vector(0.6*pointEX + 0.4*backDartLowerX, 0.6*BY + 0.4*backDartLowerY,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#back shoulder line lower part
pl.Base = FreeCAD.Vector(backDartLowerX,backDartLowerY,0.0)
points = [FreeCAD.Vector(backDartLowerX,backDartLowerY,0.0),FreeCAD.Vector(0.5*backDartLowerX + 0.5*backShoulderLineX0, 0.5*backDartLowerY + 0.5*backShoulderLineY0, 0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)

#back frames
'''pl.Base = FreeCAD.Vector(BS,WL,0.0)
points = [FreeCAD.Vector(BS,WL,0.0),FreeCAD.Vector(BS,BC,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
pl.Base = FreeCAD.Vector(BS,WL,0.0)
points = [FreeCAD.Vector(BS,WL,0.0),FreeCAD.Vector(SS,WL,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
pl.Base = FreeCAD.Vector(SS,WL,0.0)
points = [FreeCAD.Vector(SS,WL,0.0),FreeCAD.Vector(SS,BL,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
'''

#arm hole front
controlPointRadius = ((pointBPX-FAW)**2+(BL-breastDartY0)**2)**0.5
controlPointAngle = atan((breastDartY0-BL)/(pointBPX-FAW))
armholecontrolPointX = pointBPX - controlPointRadius*cos(controlPointAngle+(B/4-25)/3600*2*pi)
armholecontrolPointY = BL + controlPointRadius*sin(controlPointAngle+(B/4-25)/3600*2*pi)
#points = [FreeCAD.Vector(breastDartX0,breastDartY0,0.0),FreeCAD.Vector(armholecontrolPointX,armholecontrolPointY,0.0),FreeCAD.Vector(frontShoulderLineX0,frontShoulderLineY0,0.0)]
#bez = Draft.makeBezCurve(points,closed=False,support=None)
#Draft.autogroup(bez)
#point frontShoulderLineR
frontShoulderLineRadius = ((pointBPX-frontShoulderLineX0)**2+(frontShoulderLineY0-BL)**2)**0.5
frontShoulderAngle = atan((frontShoulderLineY0-BL)/(pointBPX-frontShoulderLineX0))
frontShoulderLineXR = pointBPX - frontShoulderLineRadius*cos(frontShoulderAngle-(B/4-25)/3600*2*pi)
frontShoulderLineYR = BL + frontShoulderLineRadius*sin(frontShoulderAngle-(B/4-25)/3600*2*pi)
#point = Draft.makePoint(frontShoulderLineXR, frontShoulderLineYR)
#Draft.autogroup(point)

#Sleeve top
sleeveTop6 = (backShoulderLineY0 + frontShoulderLineYR)/2
sleeveTopX = SS
sleeveTopY = 5/6*(sleeveTop6-BL)+BL
#point = Draft.makePoint(sleeveTopX, sleeveTopY)
#Draft.autogroup(point)

#(PENDING!!!)get front half length of sleeve hole
frontHalfSleeveHole = 201.6

#draw the sleeve triangle
triangleFrontHalfSleeveX = sleeveTopX + (frontHalfSleeveHole**2 - (sleeveTopY-BL)**2)**0.5
triangleFrontHalfSleeveY = BL
#points = [FreeCAD.Vector(triangleFrontHalfSleeveX,triangleFrontHalfSleeveY,0.0),FreeCAD.Vector(sleeveTopX, sleeveTopY,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)

triangleBackHalfSleeveX = sleeveTopX - ((frontHalfSleeveHole + 11)**2 - (sleeveTopY-BL)**2)**0.5
triangleBackHalfSleeveY = BL
#points = [FreeCAD.Vector(triangleBackHalfSleeveX,triangleBackHalfSleeveY,0.0),FreeCAD.Vector(sleeveTopX, sleeveTopY,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)

twooverthreeFrontSleeveX = triangleFrontHalfSleeveX-(FAW-SS)*2/3
#An estimation for twooverthreeFrontSleeveY
twooverthreeFrontSleeveY = BL + 13.3*ThreeY/58.5
#points = [FreeCAD.Vector(twooverthreeFrontSleeveX,BL,0.0),FreeCAD.Vector(twooverthreeFrontSleeveX, twooverthreeFrontSleeveY,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)

oneoverthreeBackSleeveX = (SS-BWL)*2/3+triangleBackHalfSleeveX
#An estimation for oneoverthreeBackSleeveY
oneoverthreeBackSleeveY = BL + 18*ThreeY/58.5
#points = [FreeCAD.Vector(oneoverthreeBackSleeveX,BL,0.0),FreeCAD.Vector(oneoverthreeBackSleeveX, oneoverthreeBackSleeveY,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)

oneoverfourFrontSleeveX = sleeveTopX + (triangleFrontHalfSleeveX-sleeveTopX)*0.25
oneoverfourFrontSleeveY = sleeveTopY - (sleeveTopY-triangleFrontHalfSleeveY)*0.25
#point = Draft.makePoint(oneoverfourFrontSleeveX, oneoverfourFrontSleeveY)
#Draft.autogroup(point)
cutlength = (((triangleFrontHalfSleeveX-sleeveTopX)*0.25)**2 + ((sleeveTopY-triangleFrontHalfSleeveY)*0.25)**2)**0.5

frontSleeveAngle = atan((sleeveTopY-triangleFrontHalfSleeveY)/(triangleFrontHalfSleeveX-sleeveTopX))
frontSleevePeakX = oneoverfourFrontSleeveX + 18.5*sin(frontSleeveAngle)
frontSleevePeakY = oneoverfourFrontSleeveY + 18.5*cos(frontSleeveAngle)
#points = [FreeCAD.Vector(oneoverfourFrontSleeveX,oneoverfourFrontSleeveY,0.0),FreeCAD.Vector(frontSleevePeakX, frontSleevePeakY,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)

backSleeveAngle = atan((sleeveTopY-triangleBackHalfSleeveY)/(sleeveTopX-triangleBackHalfSleeveX))
cutBackSleeveX = sleeveTopX - cutlength*cos(backSleeveAngle)
cutBackSleeveY = sleeveTopY - cutlength*sin(backSleeveAngle)
backSleevePeakX = cutBackSleeveX - 19.5*sin(backSleeveAngle)
backSleevePeakY = cutBackSleeveY + 19.5*cos(backSleeveAngle)
#points = [FreeCAD.Vector(cutBackSleeveX,cutBackSleeveY,0.0),FreeCAD.Vector(backSleevePeakX, backSleevePeakY,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)

frontSleeveIntersectXAdj = triangleFrontHalfSleeveX - (GL - BL)/tan(frontSleeveAngle) - 10*cos(frontSleeveAngle)
frontSleeveIntersectYAdj = GL + 10*sin(frontSleeveAngle)
#point = Draft.makePoint(frontSleeveIntersectXAdj, frontSleeveIntersectYAdj)
#Draft.autogroup(point)

backSleeveIntersectXAdj = triangleBackHalfSleeveX + (GL - BL)/tan(backSleeveAngle) - 10*sin(backSleeveAngle)
backSleeveIntersectYAdj = GL - 10*cos(frontSleeveAngle)
#point = Draft.makePoint(backSleeveIntersectXAdj, backSleeveIntersectYAdj)
#Draft.autogroup(point)
#Draw the sleeve curve
#points = [FreeCAD.Vector(triangleBackHalfSleeveX,triangleBackHalfSleeveY),FreeCAD.Vector(oneoverthreeBackSleeveX,oneoverthreeBackSleeveY),FreeCAD.Vector(backSleeveIntersectXAdj,backSleeveIntersectYAdj),FreeCAD.Vector(backSleevePeakX,backSleevePeakY),FreeCAD.Vector(sleeveTopX,sleeveTopY),FreeCAD.Vector(frontSleevePeakX,frontSleevePeakY),FreeCAD.Vector(frontSleeveIntersectXAdj,frontSleeveIntersectYAdj),FreeCAD.Vector(twooverthreeFrontSleeveX,twooverthreeFrontSleeveY),FreeCAD.Vector(triangleFrontHalfSleeveX,triangleFrontHalfSleeveY)]
#spline = Draft.makeBSpline(points,closed=False,face=True,support=None)
'''
#Draw the back half sleeve curv
points = [FreeCAD.Vector(triangleBackHalfSleeveX,triangleBackHalfSleeveY),FreeCAD.Vector(oneoverthreeBackSleeveX,oneoverthreeBackSleeveY),FreeCAD.Vector(backSleeveIntersectXAdj,backSleeveIntersectYAdj),FreeCAD.Vector(backSleevePeakX,backSleevePeakY),FreeCAD.Vector(sleeveTopX,sleeveTopY)]
spline = Draft.makeBSpline(points,closed=False,face=True,support=None)'''
#FrameforSleeve
sleeveLenX = sleeveTopX - triangleBackHalfSleeveX
sleeveLenY = 130
sleeveDownY = triangleBackHalfSleeveY-sleeveLenY

#Draw the shifted back half sleeve curve
shiftX = backShoulderLineX0 + 10 * cos(18/360*2*pi) - sleeveTopX
shiftY = backShoulderLineY0 - 10 * sin(18/360*2*pi) - sleeveTopY
sleeveTopXS = sleeveTopX + shiftX
sleeveTopYS = sleeveTopY + shiftY
triangleBackHalfSleeveXS = triangleBackHalfSleeveX + shiftX
triangleBackHalfSleeveYS = triangleBackHalfSleeveY + shiftY
oneoverthreeBackSleeveXS = oneoverthreeBackSleeveX + shiftX
oneoverthreeBackSleeveYS = oneoverthreeBackSleeveY + shiftY
oneoverfourFrontSleeveXS = oneoverfourFrontSleeveX + shiftX
oneoverfourFrontSleeveYS = oneoverfourFrontSleeveY + shiftY
backSleevePeakXS = backSleevePeakX + shiftX
backSleevePeakYS = backSleevePeakY + shiftY
backSleeveIntersectXAdjS = backSleeveIntersectXAdj + shiftX
backSleeveIntersectYAdjS = backSleeveIntersectYAdj + shiftY
sleeveDownLeftXS = triangleBackHalfSleeveX + shiftX
sleeveDownRightXS = sleeveDownLeftXS + sleeveLenX
sleeveDownYS = sleeveDownY + shiftY

#Draw the shifted back half sleeve curve
#points = [FreeCAD.Vector(triangleBackHalfSleeveXS,triangleBackHalfSleeveYS),FreeCAD.Vector(oneoverthreeBackSleeveXS,oneoverthreeBackSleeveYS),FreeCAD.Vector(backSleeveIntersectXAdjS,backSleeveIntersectYAdjS),FreeCAD.Vector(backSleevePeakXS,backSleevePeakYS),FreeCAD.Vector(sleeveTopXS,sleeveTopYS)]
#spline = Draft.makeBSpline(points,closed=False,face=True,support=None)

#Rotate the sleeve curve with central point sleeveTop
triangleBackHalfSleeveRadius = ((triangleBackHalfSleeveXS - sleeveTopXS)**2 + (triangleBackHalfSleeveYS - sleeveTopYS)**2)**0.5
triangleBackHalfSleeveAngle = atan((sleeveTopXS - triangleBackHalfSleeveXS)/(sleeveTopYS - triangleBackHalfSleeveYS))
oneoverthreeBackSleeveRadius = ((oneoverthreeBackSleeveXS - sleeveTopXS)**2 + (oneoverthreeBackSleeveYS - sleeveTopYS)**2)**0.5
oneoverthreeBackSleeveAngle = atan((sleeveTopXS - oneoverthreeBackSleeveXS)/(sleeveTopYS - oneoverthreeBackSleeveYS))
backSleevePeakRadius = ((backSleevePeakXS - sleeveTopXS)**2 + (backSleevePeakYS - sleeveTopYS)**2)**0.5
backSleevePeakAngle = atan((sleeveTopXS - backSleevePeakXS)/(sleeveTopYS - backSleevePeakYS))
backSleeveIntersectRadius = ((backSleeveIntersectXAdjS - sleeveTopXS)**2 + (backSleeveIntersectYAdjS - sleeveTopYS)**2)**0.5
backSleeveIntersectAngle = atan((sleeveTopXS - backSleeveIntersectXAdjS)/(sleeveTopYS - backSleeveIntersectYAdjS))
sleeveDownLeftRadius = ((sleeveTopXS - sleeveDownLeftXS)**2 + (sleeveTopYS - sleeveDownYS)**2)**0.5
sleeveDownLeftAngle = atan((sleeveTopXS - sleeveDownLeftXS)/(sleeveTopYS - sleeveDownYS))
sleeveDownRightRadius = ((sleeveTopXS - sleeveDownRightXS)**2 + (sleeveTopYS - sleeveDownYS)**2)**0.5
sleeveDownRightAngle = atan((sleeveTopXS - sleeveDownRightXS)/(sleeveTopYS - sleeveDownYS))

#pi is 180 degree, we do some rotations here
backSleeveRotateAngle = pi/6
triangleBackHalfSleeveXR = sleeveTopXS - triangleBackHalfSleeveRadius * sin(triangleBackHalfSleeveAngle - backSleeveRotateAngle)
triangleBackHalfSleeveYR = sleeveTopYS - triangleBackHalfSleeveRadius * cos(triangleBackHalfSleeveAngle - backSleeveRotateAngle)
oneoverthreeBackSleeveXR = sleeveTopXS - oneoverthreeBackSleeveRadius * sin(oneoverthreeBackSleeveAngle - backSleeveRotateAngle)
oneoverthreeBackSleeveYR = sleeveTopYS - oneoverthreeBackSleeveRadius * cos(oneoverthreeBackSleeveAngle - backSleeveRotateAngle)
backSleevePeakXR = sleeveTopXS - backSleevePeakRadius * sin(backSleevePeakAngle - backSleeveRotateAngle)
backSleevePeakYR = sleeveTopYS - backSleevePeakRadius * cos(backSleevePeakAngle - backSleeveRotateAngle)
backSleeveIntersectXAdjR = sleeveTopXS - backSleeveIntersectRadius * sin(backSleeveIntersectAngle - backSleeveRotateAngle)
backSleeveIntersectYAdjR = sleeveTopYS - backSleeveIntersectRadius * cos(backSleeveIntersectAngle - backSleeveRotateAngle)
sleeveDownLeftXR = sleeveTopXS - sleeveDownLeftRadius * sin(sleeveDownLeftAngle - backSleeveRotateAngle)
sleeveDownLeftYR = sleeveTopYS - sleeveDownLeftRadius * cos(sleeveDownLeftAngle - backSleeveRotateAngle)
sleeveDownRightXR = sleeveTopXS - sleeveDownRightRadius * sin(sleeveDownRightAngle - backSleeveRotateAngle)
sleeveDownRightYR = sleeveTopYS - sleeveDownRightRadius * cos(sleeveDownRightAngle - backSleeveRotateAngle)


#draw rotated backsleeve line
#points = [FreeCAD.Vector(triangleBackHalfSleeveXR,triangleBackHalfSleeveYR),FreeCAD.Vector(oneoverthreeBackSleeveXR,oneoverthreeBackSleeveYR),FreeCAD.Vector(backSleeveIntersectXAdjR,backSleeveIntersectYAdjR),FreeCAD.Vector(backSleevePeakXR,backSleevePeakYR),FreeCAD.Vector(sleeveTopXS,sleeveTopYS)]
#spline = Draft.makeBSpline(points,closed=False,face=True,support=None)
#draw rotated backsleeve line (part)
points = [FreeCAD.Vector(triangleBackHalfSleeveXR,triangleBackHalfSleeveYR),FreeCAD.Vector(oneoverthreeBackSleeveXR,oneoverthreeBackSleeveYR),FreeCAD.Vector(backSleeveIntersectXAdjR,backSleeveIntersectYAdjR)]
spline = Draft.makeBSpline(points,closed=False,face=True,support=None)


#draw back central cut line
backcentralcut0X = 1/3 * BS + 2/3 * backCollarX2
backcentralcut0Y = 1/2 * BC + 1/2 * backCollarY2
backcentralcut1X = 3/5 * pointEX + 2/5 * backDartLowerX
backcentralcut1Y = 3/5 * BY + 2/5 * backDartLowerY
backcentralcut2X = 1/2 * oneoverthreeBackSleeveXR + 1/2 * backSleeveIntersectXAdjR
backcentralcut2Y = 1/2 * oneoverthreeBackSleeveYR + 1/2 * backSleeveIntersectYAdjR
points = [FreeCAD.Vector(backcentralcut0X,backcentralcut0Y), FreeCAD.Vector(backcentralcut1X,backcentralcut1Y), FreeCAD.Vector(backcentralcut2X,backcentralcut2Y), FreeCAD.Vector(backArmHoleX0,backArmHoleY0)]
spline = Draft.makeBSpline(points,closed=False,face=True,support=None)

#draw rotated sleeve lower frame
sleeveDownLeftXR = sleeveDownLeftXR * 0.9 + sleeveDownRightXR*0.1
sleeveDownLeftYR = sleeveDownLeftYR * 0.9 + sleeveDownRightYR*0.1
points = [FreeCAD.Vector(sleeveDownLeftXR,sleeveDownLeftYR,0.0),FreeCAD.Vector(sleeveDownRightXR, sleeveDownRightYR,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
points = [FreeCAD.Vector(sleeveDownLeftXR,sleeveDownLeftYR,0.0),FreeCAD.Vector(triangleBackHalfSleeveXR,triangleBackHalfSleeveYR,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
points = [FreeCAD.Vector(0.8*sleeveTopXS + 0.2*sleeveDownRightXR,0.8*sleeveTopYS + 0.2*sleeveDownRightYR,0.0),FreeCAD.Vector(sleeveDownRightXR, sleeveDownRightYR,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)

#connect the shoulder line and the sleeve frame
points = [FreeCAD.Vector(0.5*backDartLowerX + 0.5*backShoulderLineX0, 0.5*backDartLowerY + 0.5*backShoulderLineY0,0.0), FreeCAD.Vector(sleeveTopXS,sleeveTopYS), FreeCAD.Vector(0.8*sleeveTopXS + 0.2*sleeveDownRightXR,0.8*sleeveTopYS + 0.2*sleeveDownRightYR)]
bez = Draft.makeBezCurve(points,closed=False,support=None)

Gui.SendMsgToActiveView("ViewFit")