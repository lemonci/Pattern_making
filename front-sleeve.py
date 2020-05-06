# Back Centre Length
A = 360
# Breast Circumference - 85mm
B = 1000
 
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
App.newDocument("front-sleeve")
App.setActiveDocument("front-sleeve")
App.ActiveDocument=App.getDocument("front-sleeve")
Gui.ActiveDocument=Gui.getDocument("front-sleeve")
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
points = [FreeCAD.Vector(frontCollarX,FC),FreeCAD.Vector(frontCollarX,frontCollarY),FreeCAD.Vector(CFL,frontCollarY)]
bez = Draft.makeBezCurve(points,closed=False,support=None)
Draft.autogroup(bez)
#Front shoulder line
frontShoulderLine = (CFL-CWL-collarRef)/cos(22/360*2*pi)+18
frontShoulderLineX0 = CWL - 18*cos(22/360*2*pi)
frontShoulderLineY0 = FC - frontShoulderLine*sin(22/360*2*pi)
pl.Base = FreeCAD.Vector(frontCollarX,FC,0.0)
points = [FreeCAD.Vector(frontCollarX,FC,0.0),FreeCAD.Vector(0.25*frontCollarX+0.75*frontShoulderLineX0,0.25*FC + 0.75*frontShoulderLineY0,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#back collar
backCollarX1 = BS + (collarRef + 2)/3*2
backCollarX2 = BS + (collarRef + 2)
backCollarY2 = BC + (collarRef + 2)/3
#points = [FreeCAD.Vector(BS,BC),FreeCAD.Vector(backCollarX1,BC),FreeCAD.Vector(backCollarX2,backCollarY2)]
#bez = Draft.makeBezCurve(points,closed=False,support=None)
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

#front half armhole
#points = [FreeCAD.Vector(SS,BL),FreeCAD.Vector(frontArmHoleX0,frontArmHoleY0),FreeCAD.Vector(FAW,GL-0.1),FreeCAD.Vector(FAW,GL),FreeCAD.Vector(FAW,GL+0.1)]
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
'''
#back dart upperline
backDartUpperX = pointEX + 15*cos(18/360*2*pi)
backDartUpperY = backCollarY2 - (backDartUpperX-backCollarX2)*sin(18/360*2*pi)
#pl.Base = FreeCAD.Vector(backDartUpperX,backDartUpperY,0.0)
#points = [FreeCAD.Vector(backDartUpperX,backDartUpperY,0.0),FreeCAD.Vector(pointEX,BY,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#back shoulder line upper part
#pl.Base = FreeCAD.Vector(backDartUpperX,backDartUpperY,0.0)
#points = [FreeCAD.Vector(backDartUpperX,backDartUpperY,0.0),FreeCAD.Vector(backCollarX2,backCollarY2,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#back dart lowerline
backDartLowerX = backDartUpperX + (B/32-8)*cos(18/360*2*pi)
backDartLowerY = backDartUpperY - (B/32-8)*sin(18/360*2*pi)
#pl.Base = FreeCAD.Vector(backDartLowerX,backDartLowerY,0.0)
#points = [FreeCAD.Vector(backDartLowerX,backDartLowerY,0.0),FreeCAD.Vector(pointEX,BY,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
#back shoulder line lower part
#pl.Base = FreeCAD.Vector(backDartLowerX,backDartLowerY,0.0)
#points = [FreeCAD.Vector(backDartLowerX,backDartLowerY,0.0),FreeCAD.Vector(backShoulderLineX0,backShoulderLineY0,0.0)]
#line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)


#front frames
pl.Base = FreeCAD.Vector(CFL,WL,0.0)
points = [FreeCAD.Vector(CFL,WL,0.0),FreeCAD.Vector(CFL, frontCollarY,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
pl.Base = FreeCAD.Vector(CFL,WL,0.0)
points = [FreeCAD.Vector(CFL,WL,0.0),FreeCAD.Vector(SS,WL,0.0)]
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
frontHalfSleeveHole = 500/2-100
#201.6
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

#Draw the front half sleeve curve
#points = [FreeCAD.Vector(sleeveTopX,sleeveTopY),FreeCAD.Vector(frontSleevePeakX,frontSleevePeakY),FreeCAD.Vector(frontSleeveIntersectXAdj,frontSleeveIntersectYAdj),FreeCAD.Vector(twooverthreeFrontSleeveX,twooverthreeFrontSleeveY),FreeCAD.Vector(triangleFrontHalfSleeveX,triangleFrontHalfSleeveY)]
#spline = Draft.makeBSpline(points,closed=False,face=True,support=None)

#FrameforSleeve
sleeveLenX = triangleFrontHalfSleeveX - sleeveTopX
sleeveLenY = 130
sleeveDownY = triangleFrontHalfSleeveY - sleeveLenY

#draw shifted sleeve lower frame
'''points = [FreeCAD.Vector(sleeveTopX,sleeveDownY,0.0),FreeCAD.Vector(sleeveTopX + sleeveLenX, sleeveDownY,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
points = [FreeCAD.Vector(sleeveTopX,sleeveDownY,0.0),FreeCAD.Vector(sleeveTopX,sleeveTopY,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
points = [FreeCAD.Vector(triangleFrontHalfSleeveX, triangleFrontHalfSleeveY,0.0),FreeCAD.Vector(sleeveTopX + sleeveLenX, sleeveDownY,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)'''


#Draw the shifted front half sleeve curve
shiftY = frontShoulderLineY0 + 10 * sin(45/360*2*pi) - sleeveTopY
shiftX = frontShoulderLineX0 - 10 * cos(45/360*2*pi) - sleeveTopX
sleeveTopXS = sleeveTopX + shiftX
sleeveTopYS = sleeveTopY + shiftY
frontSleevePeakXS = frontSleevePeakX + shiftX
frontSleevePeakYS = frontSleevePeakY + shiftY
frontSleeveIntersectXAdjS = frontSleeveIntersectXAdj + shiftX 
frontSleeveIntersectYAdjS = frontSleeveIntersectYAdj + shiftY
twooverthreeFrontSleeveXS = twooverthreeFrontSleeveX + shiftX
twooverthreeFrontSleeveYS = twooverthreeFrontSleeveY + shiftY
triangleFrontHalfSleeveXS = triangleFrontHalfSleeveX + shiftX
triangleFrontHalfSleeveYS = triangleFrontHalfSleeveY + shiftY
sleeveDownLeftXS = sleeveTopX + shiftX
sleeveDownRightXS = sleeveDownLeftXS + sleeveLenX
sleeveDownYS = sleeveDownY + shiftY

'''
#draw shified front sleeve curve
points = [FreeCAD.Vector(sleeveTopXS,sleeveTopYS),FreeCAD.Vector(frontSleevePeakXS,frontSleevePeakYS),FreeCAD.Vector(frontSleeveIntersectXAdjS,frontSleeveIntersectYAdjS),FreeCAD.Vector(twooverthreeFrontSleeveXS,twooverthreeFrontSleeveYS),FreeCAD.Vector(triangleFrontHalfSleeveXS,triangleFrontHalfSleeveYS)]
splinespline = Draft.makeBSpline(points,closed=False,face=True,support=None)

#draw shifted sleeve lower frame
points = [FreeCAD.Vector(sleeveDownLeftXS,sleeveDownYS,0.0),FreeCAD.Vector(sleeveDownRightXS, sleeveDownYS,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
points = [FreeCAD.Vector(sleeveDownLeftXS,sleeveDownYS,0.0),FreeCAD.Vector(sleeveTopXS,sleeveTopYS,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
points = [FreeCAD.Vector(triangleFrontHalfSleeveXS, triangleFrontHalfSleeveYS,0.0),FreeCAD.Vector(sleeveDownRightXS, sleeveDownYS,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
'''

#Rotate the sleeve curve with central point sleeveTop
frontSleevePeakRadius = ((frontSleevePeakXS - sleeveTopXS)**2 + (frontSleevePeakYS - sleeveTopYS)**2)**0.5
frontSleevePeakAngle = atan((sleeveTopXS - frontSleevePeakXS)/(sleeveTopYS - frontSleevePeakYS))
frontSleeveIntersectRadius = ((frontSleeveIntersectXAdjS - sleeveTopXS)**2 + (frontSleeveIntersectYAdjS - sleeveTopYS)**2)**0.5
frontSleeveIntersectAngle = atan((sleeveTopXS - frontSleeveIntersectXAdjS)/(sleeveTopYS - frontSleeveIntersectYAdjS))
twooverthreeFrontSleeveRadius = ((twooverthreeFrontSleeveXS - sleeveTopXS)**2 + (twooverthreeFrontSleeveYS - sleeveTopYS)**2)**0.5
twooverthreeFrontSleeveAngle = atan((sleeveTopXS - twooverthreeFrontSleeveXS)/(sleeveTopYS - twooverthreeFrontSleeveYS))
triangleFrontHalfSleeveRadius = ((triangleFrontHalfSleeveXS - sleeveTopXS)**2 + (triangleFrontHalfSleeveYS - sleeveTopYS)**2)**0.5
triangleFrontHalfSleeveAngle = atan((sleeveTopXS - triangleFrontHalfSleeveXS)/(sleeveTopYS - triangleFrontHalfSleeveYS))
sleeveDownLeftRadius = ((sleeveTopXS - sleeveDownLeftXS)**2 + (sleeveTopYS - sleeveDownYS)**2)**0.5
sleeveDownLeftAngle = atan((sleeveTopXS - sleeveDownLeftXS)/(sleeveTopYS - sleeveDownYS))
sleeveDownRightRadius = ((sleeveTopXS - sleeveDownRightXS)**2 + (sleeveTopYS - sleeveDownYS)**2)**0.5
sleeveDownRightAngle = atan((sleeveTopXS - sleeveDownRightXS)/(sleeveTopYS - sleeveDownYS))

#pi is 180 degree, we do some rotations here
backSleeveRotateAngle = -pi/6
frontSleevePeakXR = sleeveTopXS - frontSleevePeakRadius * sin(frontSleevePeakAngle - backSleeveRotateAngle)
frontSleevePeakYR = sleeveTopYS - frontSleevePeakRadius * cos(frontSleevePeakAngle - backSleeveRotateAngle)
frontSleeveIntersectXAdjR = sleeveTopXS - frontSleeveIntersectRadius * sin(frontSleeveIntersectAngle - backSleeveRotateAngle)
frontSleeveIntersectYAdjR = sleeveTopYS - frontSleeveIntersectRadius * cos(frontSleeveIntersectAngle - backSleeveRotateAngle)
twooverthreeFrontSleeveXR = sleeveTopXS - twooverthreeFrontSleeveRadius * sin(twooverthreeFrontSleeveAngle - backSleeveRotateAngle)
twooverthreeFrontSleeveYR = sleeveTopYS - twooverthreeFrontSleeveRadius * cos(twooverthreeFrontSleeveAngle - backSleeveRotateAngle)
triangleFrontHalfSleeveXR = sleeveTopXS - triangleFrontHalfSleeveRadius * sin(triangleFrontHalfSleeveAngle - backSleeveRotateAngle)
triangleFrontHalfSleeveYR = sleeveTopYS - triangleFrontHalfSleeveRadius * cos(triangleFrontHalfSleeveAngle - backSleeveRotateAngle)
sleeveDownLeftXR = sleeveTopXS - sleeveDownLeftRadius * sin(sleeveDownLeftAngle - backSleeveRotateAngle)
sleeveDownLeftYR = sleeveTopYS - sleeveDownLeftRadius * cos(sleeveDownLeftAngle - backSleeveRotateAngle)
sleeveDownRightXR = sleeveTopXS - sleeveDownRightRadius * sin(sleeveDownRightAngle - backSleeveRotateAngle)
sleeveDownRightYR = sleeveTopYS - sleeveDownRightRadius * cos(sleeveDownRightAngle - backSleeveRotateAngle)

points = [FreeCAD.Vector(sleeveTopXS,sleeveTopYS),FreeCAD.Vector(frontSleevePeakXR,frontSleevePeakYR),FreeCAD.Vector(frontSleeveIntersectXAdjR,frontSleeveIntersectYAdjR),FreeCAD.Vector(twooverthreeFrontSleeveXR,twooverthreeFrontSleeveYR),FreeCAD.Vector(triangleFrontHalfSleeveXR,triangleFrontHalfSleeveYR)]
#spline = Draft.makeBSpline(points,closed=False,face=True,support=None)

#draw rotated sleeve lower frame

sleeveDownRightXR = 0.9*sleeveDownRightXR + 0.1*sleeveDownLeftXR
sleeveDownRightYR = 0.9*sleeveDownRightYR + 0.1*sleeveDownLeftYR
points = [FreeCAD.Vector(sleeveDownLeftXR,sleeveDownLeftYR,0.0),FreeCAD.Vector(sleeveDownRightXR, sleeveDownRightYR,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
points = [FreeCAD.Vector(sleeveDownLeftXR,sleeveDownLeftYR,0.0),FreeCAD.Vector(0.2*sleeveDownLeftXR + 0.8*sleeveTopXS, 0.2*sleeveDownLeftYR + 0.8*sleeveTopYS,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)
points = [FreeCAD.Vector(triangleFrontHalfSleeveXR, triangleFrontHalfSleeveYR,0.0),FreeCAD.Vector(sleeveDownRightXR, sleeveDownRightYR,0.0)]
line = Draft.makeWire(points,placement=pl,closed=False,face=True,support=None)

#draw the front sleeve-trunk cut line
points = [FreeCAD.Vector(0.9*frontCollarX + 0.1*CFL, 0.45*frontCollarY+0.55*FC,0.0),
FreeCAD.Vector(breastDartX0,breastDartY0,0.0),
FreeCAD.Vector(twooverthreeFrontSleeveXR,twooverthreeFrontSleeveYR,0.0),
FreeCAD.Vector(triangleFrontHalfSleeveXR,triangleFrontHalfSleeveYR,0.0)
]

#
spline = Draft.makeBSpline(points,closed=False,face=True,support=None)

#draw shoulder connection
points = [FreeCAD.Vector(0.2*sleeveDownLeftXR + 0.8*sleeveTopXS, 0.2*sleeveDownLeftYR + 0.8*sleeveTopYS,0.0),
FreeCAD.Vector(frontShoulderLineX0,frontShoulderLineY0,0.0),
FreeCAD.Vector(0.75*frontShoulderLineX0 + 0.25*frontCollarX, 0.75*frontShoulderLineY0 + 0.25*FC,0.0)
]
bez = Draft.makeBezCurve(points,closed=False,support=None)

Gui.SendMsgToActiveView("ViewFit")
