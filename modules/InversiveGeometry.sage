print "InversiveGeometry(IG) module loading ..."

# O is center of inversion
# OB is diameter of circle to invert
# q1 is quadrance of OA
# q2 is quadrance of AB
def InvertCircle2Line(pO,pB,q1,q2):
    # compute q = k^2 = q1 - q2
    q = q1 - q2
    # compute diameter d = 2r where r is radius of circle to be inverted
    q3 = RT_Quadrance(pO,pB)
    # compute D image of B
    qOD = q^2/q3
    pD = Extend(pO,pB,qOD)
    return pD

# InverseRadius returns radius of circle (O,r) after inversion in circle (P,k) where q is quadrance from P to O
def InverseRadius(r,q,k):
    return r*k^2/abs(q - r^2)

# RadicalAxisFoot returns radical axis foot for two circles centered pG1,pG2 with squared radius q1,q2
def RadicalAxisFoot(pG1,q1,pG2,q2):
    d12 = RT_Distance(pG1,pG2)
    q12 = d12^2
    xm = q12 - q2 + q1
    x = xm/(2*d12)
    pX = LineBarycentric(pG1,pG2,1-(x/d12),(x/d12))
    print "Dist(G1,X) = ",RT_Distance(pG1,pX)," vs ",x
    return pX

# LimitingPoints return two limiting points for two disjoint circles
# when making inversion with center one of these two points and any radius, than images of two circles (G1,r1) and (G2,r2)
# are concentric circles

def LimitingPoints(pG1,r1,pG2,r2):
    q1 = r1^2 ; q2 = r2^2
    # Set the middle point of segment G1-G2 an get altitude line
    # the line is the locus of points with same distance to both centers G1,G2
    p1 = Middle(pG1,pG2)
    l = RT_Altitude_Line(pG1,pG2,p1)
    # Compute point on altitude line far away from centers
    p2 = RT_Any_Point(l,p1)
    p3 = Extend(p1,p2,q1 + q2)
    # Intersect circle from p3 with two other circles
    r3 = RT_Distance(p3,pG1) ; q3 = r3^2
    [p4,p5] = CirclesIntersect(p3,q3,pG1,q1)
    [p6,p7] = CirclesIntersect(p3,q3,pG2,q2)
    # Get radical center point p8 of three circles, intersecting two radical axis p4-p5 and p6-p7
    p8 = Intersect(p4,p5,p6,p7)
    # Compute squared radius of orthogonal circle R^2 = d^2 - r^2
    q8 = RT_Quadrance(p8,pG1) - q1 ; r8 = sqrt(q8)
    # Compute chord line, intersection of line of centers with circle ...are the limiting points
    [p9,p10] = Chord(p8,q8,pG1,pG2)
    return [p9,p10]

# InvertPoint2Point returns point image of one point by circle inversion
# p1 is center of inversion O
# q1 is quadrance of inversion OM' OM  = q1
# p2 is the point to invert M
def InvertPoint2Point(p1,q1,p2):
    q = (q1/RT_Distance(p1,p2))^2
    p3 = Extend(p1,p2,q)
    return p3

def RT_InvertPoint2Point(p1,q1,p2):
    q = (q1/RT_Distance(p1,p2))^2
    p3 = RT_Extend(p1,p2,q)
    return p3

# InvertCircle2Circle returns circle image of one circle by circle inversion
# p1 is center of inversion O
# q1 is quadrance of inversion OM' OM  = q1
# p2 is the center of circle and r2 its radius
def InvertCircle2Circle(p1,q1,p2,r2):
    # Compute a chord of circle : intersection of line of centers with the circle
    [p3,p4] = Chord(p2,r2^2,p1,p2)  # p3,p4 are two chords points of diameter circle (p2,r2) thru p1
    # Compute images of chord limits points and their Middle
    p5 = InvertPoint2Point(p1,q1,p3) # p5 is image of p3
    p6 = InvertPoint2Point(p1,q1,p4) # p6 is image of p4
    p7 = Middle(p5,p6)               # p7 is middle of p5p6 or the center of image circle
    # Power of a point theorem, to get distance on tangent
    q7 = RT_Distance(p1,p5)*RT_Distance(p1,p6)
    # Get radius with formula r/r' = q/qT
    r7 = (r2*q7)/q1
    return [p3,p4,p5,p6,p7,r7]

def RT_InvertCircle2Circle(p1,q1,p2,r2):
    # Compute a chord of circle : intersection of line of centers with the circle
    [p3,p4] = Chord(p2,r2^2,p1,p2)
    # Compute images of chord limits points and their Middle
    p5 = RT_InvertPoint2Point(p1,q1,p3)
    p6 = RT_InvertPoint2Point(p1,q1,p4)
    p7 = Middle(p5,p6)
    # Power of a point theorem, to get distance on tangent
    q7 = RT_Distance(p1,p5)*RT_Distance(p1,p6)
    # Get radius with formula r/r' = q/qT
    r7 = (r2*q7)/q1
    return [p3,p4,p5,p6,p7,r7]

# RT_InvertedCircleRadius returns radius of inverted circle
# q1 : squared radius of circle inversion OM OM' = q1
# q is quadrance between center of inversion circle and center of circle radius R1
# to be inverted
def RT_InvertedCircleRadius(q1,q,R1):
    R = R1*q1/(R1^2 - q)
    return R

# RT_LimitingPointQuadrances returns quadrances between image circles radius R and r to center
# or circle inversion
# j = (1 - sn)/(1 + sn)  Example  if n = pi/6 => sn = 1/2 => j = (1-1/2)/(1+1/2) = 1/3
# q1 : squared radius of circle inversion OM OM' = q1
# q is quadrance between center of inversion circle and center of circle radius R1
#     to be inverted
def RT_LimitingPointQuadrances(j,q1,q,R1):
    q3 = q*q1^2/(R1^2 - q)^2
    q4 = q*q1^2/(R1^2*j^2 - q)^2
    return [q3,q4]


# RT_SpreadLimitingPointQuadrances returns quadrances between image circles radius R and r to center
# or circle inversion
# j = (1 - sn)/(1 + sn)  Example  if n = pi/6 => sn = 1/2 => j = (1-1/2)/(1+1/2) = 1/3
# q1 : squared radius of circle inversion OM OM' = q1
# q is quadrance between L center of inversion circle and O center of circle radius R1
#     to be inverted
# s is spread of lines OL,AB
#
# q = x^2 + y^2 ; s = y^2/q => y^2 = q*s and x^2 = q*(1-s)
def RT_SpreadLimitingPointQuadrances(j,q1,q,R1,s):
    q3 = q*q1^2/((R1^2 + q)^2 - 4*R1^2*q*(1 - s))
    q4 = q*q1^2/((R1^2*j^2 + q)^2 - 4*R1^2*j^2*q*(1 - s))
    return [q3,q4]


# RT_HexletQuadrance returns quadrance between centers of image circles
# origin concentric circles is outer circle radius R1 and inner circle radius r1 = R1(1 - sin pi/6)/(1 + sin pi/6)
# q1 : squared radius of circle inversion OM OM' = q1
# q is quadrance between center of inversion circle and center of circle radius R1
# to be inverted
def RT_HexletQuadrance(q1,q,R1):
    q2 = 64*q*R1^4*q1^2/((R1^2 - 9*q)^2*(R1^2 - q)^2)
    return q2

# RT_QuadranceToInvertedCircleCenter returns quadrance from p2 to center of inverted circle
# p1 : center of circle inversion
# q1 : squared radius of circle inversion OM OM' = q1
# p2 : center of circle
# r2 : radius of circle
# q12 is quadrance between p1 and p2
def RT_QuadranceToInvertedCircleCenter(q1,q12,r2):
    q =  (r2^2 - q12 + q1)^2*q12/(r2^2 - q12)^2
    return q


# InvertRatioRadii returns the ratio rR/rA between radii of image circle and circle (pA,rA) thru cricle inversion limiting point A0
# and quadrance q.
# Circle (R,rR) is image of circle (A,rA) and Circle (S,rS) is image of circle (B,rB)
# Point I is intersection of radical axis of both circles with line of centers
def InvertRatioRadii(pA,pB,rA,rB,q):
    dAB = RT_Distance(pA,pB)
    dAI = (dAB^2 - rB^2 + rA^2 )/(2*dAB)  # distance from A to middle point I of limiting points A0 and B0
    qIA0 = dAI^2 - rA^2
    m = sqrt(qIA0)  # the square root factor
    xA0 = dAI - m
    dA0R1 = q/(-xA0 - rA)
    dA0R2 = q/(-xA0 + rA)
    rR = (dA0R2 - dA0R1)/2
    ratio_RA = rR/rA
    dA0S1 = q/(dAB - rB - xA0)
    dA0S2 = q/(dAB + rB - xA0)
    rS = (dA0S1 - dA0S2)/2
    ratio_SB = rS/rB
    return [ratio_RA,ratio_SB]

# OrthogonalCircle returns a circle orthogonal at p2 from circle center p1 thru p2
# parameter p3 is any point different of p1,p2...can be used to produce various range of circles
# Being orthogonal, the returned circle is invariant under inversion with respect to p1 circle
def OrthogonalCircle(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2)
    [p4,p9] = Chord(p1,q12,p1,p3) # intersect with radial to get p4 another point on circle...p1p2p4 isoscele
    l24 = RT_Line(p2,p4)
    p5 = RT_Line_Reflection(l24,p3) # reflection of point p3 thru chord p2p4
    p6 = Intersect(p2,p5,p1,p4) # image of p3 under circle inversion circle center p1
    p7 = TriangleCircumcenter(p2,p3,p6)
    R = TriangleCircumradius(p2,p3,p6)
    return [p7,R]

print "...EuclideanGeometry module loaded"