print "EuclideanGeometry (EG) module loading ..."

# Factorisation of symbolic cartesian coordinates
def CartesianFactor(p1):
    p2 = p1
    if  p1[0] != 0:
        p2[0] = p1[0].factor()
    if  p2[0] != 0:
        p2[1] = p1[1].factor()
    return p2

# Numeric Evaluation of cartesian coordinates
def CartesianEvalReal(p1):
    p2 = p1
    p2[0] = p1[0].n() ; p2[1] = p1[1].n()
    return p2

# Cartesian2Barycentric returns barycentric coordinates for a point p given three points p1,p2,p3
# coordinates are not homogenous (sum 1)...use next function to get them
def Cartesian2Barycentric(p1,p2,p3,p):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3 ; [x,y] = p
    w1 = (y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)
    w2 = (y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)
    w3 = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
    return [w1,w2,w3]

# Cartesian2HomogeneousBarycentric returns homogenous (sum is 1) barycentric coordinates for a point p given three points p1,p2,p3
def Cartesian2HomogeneousBarycentric(p1,p2,p3,p):
    [w1,w2,w3] = Cartesian2Barycentric(p1,p2,p3,p); w = w1 + w2 + w3
    return [w1/w,w2/w,w3/w]

# Barycentric2Cartesian returns carteseian coordinates from barycentric
def Barycentric2Cartesian(p1,w1,p2,w2,p3,w3):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3
    x = (w1*x1 + w2*x2 + w3*x3)/(w1 + w2 + w3) ; y = (w1*y1 + w2*y2 + w3*y3)/(w1 + w2 + w3)
    return [x,y]

# Midpoint for segment
# p1 first point of segment and p2 second point
def Middle(p1,p2):
    [x1,y1] = p1; [x2,y2] = p2
    return vector( [ (x1+x2)/2, (y1+y2)/2] )

# Returns point symetric to another point (midpoint)
# p1 mid point of segment and p2 second point
def Symetric(p1,p2):
    [x1,y1] = p1; [x2,y2] = p2
    return [ 2*x1-x2, 2*y1-y2 ]

# Segment Isotomic point returns isotomic point of p3 on segment p1p2
# the isotomic point is symetric with respect to middle of p1p2
# Function is used to get "isotomic conjugate" of interior point of triangle
def SegmentIsotomic(p1,p2,p3):
    p4 = Middle(p1,p2)
    p5 = Symetric(p4,p3)
    return p5

# Returns point at quadrance q from point point p1 on line p1-p2
def Extend(p1,p2,q):
    [x1,y1] = p1; [x2,y2] = p2
    q12 = RT_Quadrance(p1,p2)
    # n(digits=nd) is mandatory ore we have an endless sqrt() or real() call
    k = real(sqrt(q/q12).n(digits=nd))   # must take real part, because sometimes sqrt return approx of sqrt(0) as I.exp-N
    x3 = x1 + k*(x2-x1) ; y3 = y1 + k*(y2-y1)
    return [ x3, y3 ]

# Same as extend but symbolic
def RT_Extend(p1,p2,q):
    [x1,y1] = p1; [x2,y2] = p2
    q12 = RT_Quadrance(p1,p2)
    k = sqrt(q/q12)
    #if k!= 0:
    #    k = k.factor()
    x3 = x1 + k*(x2-x1) ; y3 = y1 + k*(y2-y1)
    return [ x3, y3 ]

# Same as extend but rational if distance between two points is rational
def ExtendPoint(p1,p2,d):
    [x1,y1] = p1; [x2,y2] = p2
    d12 = RT_Distance(p1,p2)  # compute distance in order to have vector p1p2 norm 1 when divided by d12
    k = d/d12
    x3 = x1 + k*(x2-x1) ; y3 = y1 + k*(y2-y1)
    return [ x3, y3 ]

# OppositeExtend returns point at quadrance q from point point p1 on line p1-p2 butin opposite direction
def OppositeExtend(p1,p2,q):
    [x1,y1] = p1; [x2,y2] = p2
    q12 = RT_Quadrance(p1,p2) ; k = sqrt(q/q12)
    x3 = x1 - k*(x2-x1) ; y3 = y1 - k*(y2-y1)
    return [ x3, y3 ]

# Intersect returns point at intsection of two lines : one thru p1 and p2 and one thru p3 and p4
def Intersect(p1,p2,p3,p4):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3; [x4,y4] = p4
    u = ( (x3 - x1)*(y4-y3) - (y3 - y1)*(x4-x3) ) / ( (x2-x1)*(y4-y3) - (y2-y1)*(x4-x3) )
    x = x1 + u*(x2-x1) ; y = y1 + u*(y2-y1)
    return [ x, y ]

# Perpendicular bissector of segment p1-p2, returns middle point of p1-p2 and any point on altitude
def PerpendicularBissector(p1,p2):
    p3 = Middle(p1,p2)
    l34 = RT_Altitude_Line(p1,p2,p3)
    p4 = RT_Any_Point(l34,p3)
    return [p3,p4]

# CirclesIntersect returns intersections of two circles
def CirclesIntersect(p1,q1,p2,q2):
    [x1,y1] = p1; [x2,y2] = p2
    e = x2 - x1
    f = y2 - y1
    q12 = e^2 + f^2
    d12 = sqrt(q12) #.n(digits=nd)
    k = (q12 + q1 - q2)/(2.0*d12)
    ql = abs(q1 - k^2)
    l = real(sqrt(ql)) # must take real part, because sometimes sqrt return approx of sqrt(0) as I.exp-N
    #
    xm = x1 + (e*k + f*l) / d12
    ym = y1 + (f*k - e*l) / d12
    xp = x1 + (e*k - f*l) / d12
    yp = y1 + (f*k + e*l) / d12
    #
    return [ [xm,ym], [xp,yp] ]

# RT_CirlcesIntersect is same than CirclesIntersect for symbolic computations, not using abs,real func, and real numbers
def RT_CirclesIntersect(p1,q1,p2,q2):
    [x1,y1] = p1; [x2,y2] = p2
    e = x2 - x1
    f = y2 - y1
    q12 = e^2 + f^2
    d12 = sqrt(q12)
    k = (q12 + q1 - q2)/(2*d12)
    ql = q1 - k^2
    l = sqrt(ql) # must take real part, because sometimes sqrt return approx of sqrt(0) as I.exp-N
    #
    xm = x1 + (e*k + f*l) / d12
    ym = y1 + (f*k - e*l) / d12
    xp = x1 + (e*k - f*l) / d12
    yp = y1 + (f*k + e*l) / d12
    #
    return [ [xm,ym], [xp,yp] ]


# Chord return chords points on line p2-p3 intersecting circle p1 quadrance q1
def Chord(p1,q1,p2,p3):
    p4 = RT_Foot(p2,p3,p1)   # p4 is foot of p1 on line p2-p3
    q2 = q1 - RT_Quadrance(p1,p4)
    # Handle when foot p4 is one of parameter p2 or p3 by choosing max altitude point
    if RT_Quadrance(p3,p4) < RT_Quadrance(p2,p4):
        p5 = Extend(p4,p2,q2)
    else:
        p5 = Extend(p4,p3,q2)
    p6 = Symetric(p4,p5)
    return [p5,p6]

# Bissect angle p2-p1 and p2-p3, returning point on bissector on bissector line
# at quadrance q from p2
def Bissect(p1,p2,p3,q):
    q13 = RT_Quadrance(p1,p3)
    p4 = Arc(p2,p1,p2,1) ; p5 = Arc(p2,p3,p2,1)
    [ p6,p7 ] = CirclesIntersect(p4,q13,p5,q13)
    # Choose nearest point
    if RT_Quadrance(p2,p7) < RT_Quadrance(p2,p6):
        p8 = p6
    else:
        p8 = p7
    p9 = Extend(p2,p8,q)
    return p9

# Same as Bisect but symbolic an choose one point
def RT_Bissect(p1,p2,p3,q):
    q13 = RT_Quadrance(p1,p3)
    p4 = RT_Arc(p2,p1,p2,1) ; p5 = RT_Arc(p2,p3,p2,1)
    [ p6,p7 ] = RT_CirclesIntersect(p4,q13,p5,q13)
    p8 = p6
    p9 = RT_Extend(p2,p8,q)
    return p9

# CircleTangents returns two tangent points T on circle with tangent lines intesecting at outside point P
# Use of the "power of a point theorem" : PT^2 = PA PA' where AA' is any secant segment to circle
# p1 is center of triangle and q1 its squared radius
# p2 is outside of tirangle point
def CircleTangents(p1,q1,p2):
    # If p2 is outside cicle then we get points p3 and p4...secant goes thru center
    p3 = Extend(p1,p2,q1) ; p4 = Symetric(p1,p3)
    # Power of point theorem for last arg
    [p5,p6] = CirclesIntersect(p1,q1,p2,RT_Distance(p2,p3)*RT_Distance(p2,p4))
    return [p5,p6]

# CircleTangencyPoints is the same than CircleTangents with geometric constuction
def CircleTangencyPoints(p1,q1,p2):
    p3 = Middle(p1,p2)  # Middle point, defining center of orthogonal circle with radius half distance outside point from center circle
    # Intersect both circles to get tangency points
    [p4,p5] = CirclesIntersect(p1,q1,p3,RT_Quadrance(p1,p2)/4)
    return [p4,p5]

# CirclesHomotheticPoints returns the external and internal homothetic points of two circles
def CirclesHomotheticPoints(p1,q1,p2,q2):   # p1 = pA q1 = qA p2 = pB q2 = qB
    # Get circles intersections with line of centers
    p3 = RT_Arc(p1,p2,p1,q1)   # p3 = pA0
    p4 = RT_Arc(p1,p2,p2,q2)   # p4 = pB0
    # Compute last point triangle on circle
    #k = 1    # use equilateral triangle
    k = 2/5  # use Pythagorean triangle 3-4-5
    q35 = k*q1 ; q46 = k*q2
    [p9,p5] = CirclesIntersect(p1,q1,p3,q35) # p5 = pA1
    [p9,p6] = CirclesIntersect(p2,q2,p4,q46) # p6 = pB1
    # Get diameters points
    p7 = RT_Arc(p5,p1,p1,q1) # p7 = pA2
    p8 = RT_Arc(p6,p2,p2,q2) # p8 = pB2
    # Get homothetical centers as intersection of lines
    p9 = Intersect(p5,p6,p7,p8) # p9 = pE_AB
    p10 = Intersect(p5,p8,p7,p6)
    #
    return [p9,p10,p5,p6,p7,p8]

# CirclesAntihomologousPoints returns the antihomologous points of two circles
def CirclesAntihomologousPoints(p1,q1,p2,q2):
    # Get external homothetic point p3 and two rays p5-p6, p7-p8
    [p3,p4,p5,p6,p7,p8] = CirclesHomotheticPoints(p1,q1,p2,q2)
    # Get other rays intersection points, rays from external homothetic point
    [p9,p0] = Chord(p1,q1,p3,p5)
    [p10,p0] = Chord(p1,q1,p3,p7)
    [p0,p11] = Chord(p2,q2,p3,p6)
    [p0,p12] = Chord(p2,q2,p3,p8)
    #
    return [p9,p10,p11,p12]

# CirclesAntihomologousCenter returns center of circumcircle of antihomologous points quadrilateral
def CirclesAntihomologousCenter(p1,q1,p2,q2):
    # Get quadrilateral of antihomologous points
    [p3,p4,p5,p6] = CirclesAntihomologousPoints(p1,q1,p2,q2)
    # Intersecting two perpendicular bissector, choosing two edges not parallel
    [p35,p53] = PerpendicularBissector(p3,p5)
    [p46,p64] = PerpendicularBissector(p4,p6)
    p7 = Intersect(p35,p53,p46,p64)
    return p7

# CirclesRadicalAxis of two circles returns the radical axis of two circles
def CirclesRadicalAxis(p1,q1,p2,q2):
    # Get quadrilateral of antihomologous points
    [p3,p4,p5,p6] = CirclesAntihomologousPoints(p1,q1,p2,q2)
    # Intersect radical axis of couple of circles : circle p1 (or p2) and circumcircle of antihomologous points quadrilateral
    p7 = Intersect(p3,p4,p5,p6)
    p8 = RT_Foot(p1,p2,p7)
    return [p8,p7]

# CirclesRadicalCenter returns radical circle , centered at radical center intersection of the three radical axis
def CirclesRadicalCircle(p1,q1,p2,q2,p3,q3):
    # Get Radical center of three circles as intersection of two radical axis for two couples of circles
    [p12,p21] = CirclesRadicalAxis(p1,q1,p2,q2)
    [p23,p32] = CirclesRadicalAxis(p2,q2,p3,q3)
    p4 = Intersect(p12,p21,p23,p32)
    [p5,p6] = CircleTangencyPoints(p1,q1,p4)
    q45 = RT_Quadrance(p4,p5)
    return [p4,q45]

# Returns two points tangent to two circles center p1 squared radius q1 and center p2 squared radius q2
def CirclesCommonTangent(p1,q1,p2,q2):
    r1 = sqrt(q1) ; r2 = sqrt(q2)
    p12 = Middle(p1,p2) ; r12 = (r1+r2)/2 ; q12 = r12^2 # center center middle point
    if r1 > r2:
        r3 = r1 - r2; q3 = r3^2
        # Get tangency point as intersection of circles one radius r1-r2 and the other r12 = (r1+r2)/2
        [p7,p8] = CirclesIntersect(p1,q3,p12,q12)
        # Intersect radius vector from center of biggest radius circle thru tangency point to get first tangent point
        [p5,p9] = Chord(p1,q1,p1,p8)
        # Use parallel (or rectangular shape) to get tangent to lowest radius circle
        [p7,p6] = CirclesIntersect(p5,RT_Quadrance(p2,p8),p2,RT_Quadrance(p8,p5))
    else:
        r3 = r2 - r1; q3 = r3^2
        # Get tangency point as intersection of circles one radius r2-r1 and the other r12 = (r1+r2)/2
        [p8,p7] = CirclesIntersect(p2,q3,p12,q12)
        # Intersect radius vector from center of biggest radius circle thru tangency point to get first tangent point
        [p5,p9] = Chord(p2,q2,p2,p8)
        # Use parallel (or rectangular shape) to get tangent to lowest radius circle
        [p6,p7] = CirclesIntersect(p5,RT_Quadrance(p1,p8),p1,RT_Quadrance(p8,p5))
    return [p5,p6]

# ScalarProduct returns product between vectors p1p2 and p1p3
def ScalarProduct(p1,p2,p3):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    return (x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)

# TriangleCosineAngle returns cosine angle p1p2p3  (apex at p2)
def TriangleCosineAngle(p1,p2,p3):
    cs = (RT_Quadrance(p2,p1) + RT_Quadrance(p2,p3) - RT_Quadrance(p1,p3))/(2*RT_Distance(p2,p1)*RT_Distance(p2,p3))
    return cs

# TriangleEdgesQuadrances returns quadrances for edges
def TriangleEdgesQuadrances(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2); q23 = RT_Quadrance(p2,p3); q31 = RT_Quadrance(p3,p1)
    return [q12,q23,q31]

# TriangleEdgesLengths returns lengths f dges
def TriangleEdgesLengths(p1,p2,p3):
    [q12,q23,q31] = TriangleEdgesQuadrances(p1,p2,p3)
    d12 = sqrt(q12) ; d23 = sqrt(q23) ; d31 = sqrt(q31)
    return [d12,d23,d31]

# TriangleCosinesAngles returns cosines angles p1p2p3, p2p3p1 and p3,p1p2
def TriangleCosinesAngles(p1,p2,p3):
    [q12,q23,q31] = TriangleEdgesQuadrances(p1,p2,p3)
    [d12,d23,d31] = TriangleEdgesLengths(p1,p2,p3)
    cs123 = (q12 + q23 - q31)/(2*d12*d23)
    cs231 = (q23 + q31 - q12)/(2*d23*d31)
    cs312 = (q31 + q12 - q23)/(2*d31*d12)
    return [cs123,cs231,cs312]

# Arc returns point intersection at quadrance q from center of circle p3 and on line between p1 and p2
def Arc(p1,p2,p3,q):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    # (x-x3)^2 + (y-y3)^2 = q
    # x-x1 = u*(x2-x1) ; y-y1 = u*(y2-y1)
    # (u*(x2-x1) + x1-x3)^2 + (u*(y2-y1) + y1-y3)^2 = q
    # ( (x2-x1)^2 + (y2-y1)^2 )*u^2 - 2*((x2-x1)*(x1-x3)+((y2-y1)*(y1-y3))*u + (x3-x1)^2+(y3-y1)^2-q = 0
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3)
    a = q12 ; b = - ScalarProduct(p1,p2,p3); c = q13 - q
    d = sqrt(b^2 - a*c).n(digits=nd)
    um = (-b-d)/a; xm = x1 + um*(x2-x1) ; ym = y1 + um*(y2-y1)
    up = (-b+d)/a; xp = x1 + up*(x2-x1) ; yp = y1 + up*(y2-y1)
    if RT_Quadrance(p2,[xm,ym]) < RT_Quadrance(p2,[xp,yp]):
        [x,y] = [xm,ym]
    else:
        [x,y] = [xp,yp]
    return [x,y]

# Same as Arc but symbolic and choose up rejecting um
def RT_Arc(p1,p2,p3,q):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    # (x-x3)^2 + (y-y3)^2 = q
    # x-x1 = u*(x2-x1) ; y-y1 = u*(y2-y1)
    # (u*(x2-x1) + x1-x3)^2 + (u*(y2-y1) + y1-y3)^2 = q
    # ( (x2-x1)^2 + (y2-y1)^2 )*u^2 - 2*((x2-x1)*(x1-x3)+((y2-y1)*(y1-y3))*u + (x3-x1)^2+(y3-y1)^2-q = 0
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3)
    a = q12 ; b = - ScalarProduct(p1,p2,p3); c = q13 - q
    d = sqrt(b^2 - a*c)
    up = (-b+d)/a; xp = x1 + up*(x2-x1) ; yp = y1 + up*(y2-y1)
    [x,y] = [xp,yp]
    return [x,y]


# Arcs returns the TWO intersection points
def Arcs(p1,p2,p3,q):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    # (x-x3)^2 + (y-y3)^2 = q
    # x-x1 = u*(x2-x1) ; y-y1 = u*(y2-y1)
    # (u*(x2-x1) + x1-x3)^2 + (u*(y2-y1) + y1-y3)^2 = q
    # ( (x2-x1)^2 + (y2-y1)^2 )*u^2 - 2*((x2-x1)*(x1-x3)+((y2-y1)*(y1-y3))*u + (x3-x1)^2+(y3-y1)^2-q = 0
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3)
    a = q12 ; b = - ScalarProduct(p1,p2,p3); c = q13 - q
    d = sqrt(b^2 - a*c).n(digits=nd)
    um = (-b-d)/a; xm = x1 + um*(x2-x1) ; ym = y1 + um*(y2-y1)
    up = (-b+d)/a; xp = x1 + up*(x2-x1) ; yp = y1 + up*(y2-y1)
    return [[xm,ym],[xp,yp]]

# EulerLine2Vertices returns triangle vertices given three lines: Euler line (l3) and two edge lines (l1,l2)
# line can be obtained with l12 = RT_Line(p1,p2) where p1 and p2 are two points on the line
def EulerLine2Vertices(p1,p2,p3,p4,p5,p6):
    l1 = RT_Line(p1,p2) ; l2 =  RT_Line(p3,p4)   # side lines
    l3 = RT_Line(p5,p6)  # Euler line
    # Intersecting side lines to get one vertice A
    p7 = Intersect(p1,p2,p3,p4)
    # Reflecting two points of Euler line and intersecting the two reflected lines
    p51 = RT_Line_Reflection(l1,p5) ; p52 = RT_Line_Reflection(l2,p5)
    p61 = RT_Line_Reflection(l1,p6) ; p62 = RT_Line_Reflection(l2,p6)
    p8 = Intersect(p51,p61,p52,p62)  # Euler reflection point
    # Intersecting bisector with Euler line to ge circumcenter
    [p9,p10] = PerpendicularBissector(p8,p7)
    p11 = Intersect(p9,p10,p5,p6)
    # Intersecting circumcircle with side lines
    q11 = RT_Quadrance(p11,p7)
    [p12,p13] = Chord(p11,q11,p1,p2) # B vertice
    [p14,p15] = Chord(p11,q11,p3,p4) # C vertice
    # Return A,B and C vertices of triangle
    return [p7,p12,p14]

# CevianTriangle returns the cevian triangle of a point with respect to a triangle
# of edges with  medians
def CevianTriangle(p1,p2,p3,p4):
    p14 = Intersect(p1,p4,p2,p3)
    p24 = Intersect(p2,p4,p1,p3)
    p34 = Intersect(p3,p4,p1,p2)
    return [p14,p24,p34]

# MedianTriangle returns the median triangle of a triangle : triangle whose vertices are intersections
# of edges with  medians
def MedianTriangle(p1,p2,p3):
    p12 = Middle(p1,p2)
    p13 = Middle(p1,p3)
    p23 = Middle(p2,p3)
    return [p12,p13,p23]

# MedianQuadrance returns the quadrance of the median from p1 to middle of edge p2-p3
def MedianQuadrance(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2) ; q13 = RT_Quadrance(p1,p3) ; q23 = RT_Quadrance(p2,p3)
    return (2*q12 + 2*q13 - q23)/4

# Orthic triangle returns the orthic triangle of a triangle : triangle whose vertices are intersections
# of vertices altitudes crossing at orthocenter
def OrthicTriangle(p1,p2,p3):
    p4 = TriangleOrthocenter(p1,p2,p3)
    p5 = Intersect(p1,p4,p2,p3)
    p6 = Intersect(p2,p4,p1,p3)
    p7 = Intersect(p3,p4,p1,p2)
    return [p5,p6,p7]

# Contact triangle returns the contact triangle of a triangle : triangle whose vertices are intersections
# of incircle with edges or intersections of three mutually tangent circles centered at vertices
def ContactTriangle(p1,p2,p3):
    d12 = RT_Distance(p1,p2) ; d13 = RT_Distance(p1,p3) ; d23 = RT_Distance(p2,p3)
    s = (d12 + d13 + d23)/2
    p6 = ExtendPoint(p1,p2,s-d23)
    p4 = ExtendPoint(p2,p3,s-d13)
    p5 = ExtendPoint(p3,p1,s-d12)
    return [p4,p5,p6]

# TriangleAngleBisector returns the ratio BD/DC = AB/AC as given by angle bisector theorem where AD is segment
# bisecting angle A in triangle ABC A = p1, B = p2, C = p3
# BD/DC = r and BD + DC = a => DC(1 + rq) = a => DC^2 = a^2/(1 + r)^2 and BD^2 = rq DC^2
def TriangleAngleBisectorRatio(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2) ; q13 = RT_Quadrance(p1,p3)
    rq = q13/q12
    return rq

def TriangleAngleBisector(p1,p2,p3):
    rq = TriangleAngleBisectorRatio(p1,p2,p3)
    q23 = RT_Quadrance(p2,p3)
    q24 = q23/(1 + sqrt(rq))^2 ; q34 = rq*q24
    [p4,p5] = CirclesIntersect(p2,q24,p3,q34)
    return p4

# TriangleIncircle returns incenter and contact triangle
def TriangleIncircle(p1,p2,p3):
    [p5,p6,p7] = ContactTriangle(p1,p2,p3)
    p45 = RT_HalfTurn(p2,p3,p5)
    p46 = RT_HalfTurn(p3,p1,p6)
    p47 = RT_HalfTurn(p1,p2,p7)
    p4 = Intersect(p5,p45,p6,p46)
    return [p4,p5,p6,p7]

# TriangleSymmedianLine returns a point on the symmedian line for vertice p3
# It's the locus of midpoints of the antiparallels to the side p1-p2 of triangle
def TriangleSymmedianLine(p1,p2,p3):
    p4 = TriangleCircumcenter(p1,p2,p3)
    p5 = RT_HalfTurn(p4,p1,p1)
    p6 = RT_HalfTurn(p4,p2,p2)
    p7 = Intersect(p1,p5,p2,p6)
    return p7

# TriangleSymmedianPoint returns symmedian point as intersection of two symmedian lines
def TriangleSymmedianPoint(p1,p2,p3):
    p12 = TriangleSymmedianLine(p1,p2,p3)
    p23 = TriangleSymmedianLine(p2,p3,p1)
    p4 = Intersect(p3,p12,p1,p23)
    return p4

# TriangleGergonnePoint returns Gergonne point, which is symmedian of contact triangle
def TriangleGergonnePoint(p1,p2,p3):
    [p4,p5,p6] = ContactTriangle(p1,p2,p3)
    p7 = TriangleSymmedianPoint(p4,p5,p6)
    return p7

# TriangleNinePointCircleCenter returns center of Nine Point circle
# Nine point circles goes thru the three vertices of median triangles (median points of edges) that's why
# we compute the circumcircle of these points
def TriangleNinePointCircleCenter(p1,p2,p3):
    [p12,p13,p23] = MedianTriangle(p1,p2,p3)
    p4 = TriangleCircumcenter(p12,p13,p23)
    return [p4,p12]

# TriangleSpiekerCircleCenter returns Spieker center of triangle
def TriangleSpiekerCircleCenter(p1,p2,p3):
    [p12,p13,p23] = MedianTriangle(p1,p2,p3)
    [p4,p5,p6,p7] = TriangleIncircle(p23,p13,p12)
    return [p4,p5]

# TriangleApolloniusCircle returns center of Apollonius circle and one point of circle
def TriangleApolloniusCircle(p1,p2,p3):
    p4 = TriangleCircumcenter(p1,p2,p3)
    p5 = TriangleSymmedianPoint(p1,p2,p3)
    [p6,p0] = TriangleNinePointCircleCenter(p1,p2,p3)
    [p7,p0] = TriangleSpiekerCircleCenter(p1,p2,p3)
    p8 = Intersect(p4,p5,p6,p7)
    p0 = Middle(p2,p3)
    p9 = ExtendPoint(p7,p0,-RT_Distance(p7,p0)*RT_Distance(p7,p8)/RT_Distance(p7,p6))
    return [p8,p9]

# ExtouchTriangle returns contact points between excircles and edges of triangle
# by dropping center of circle thru excenters on edges
def ExtouchTriangle(p1,p2,p3):
    p4 = TriangleBevanCircleCenter(p1,p2,p3)
    p5 = RT_Foot(p2,p3,p4)
    p6 = RT_Foot(p1,p3,p4)
    p7 = RT_Foot(p1,p2,p4)
    return [p5,p6,p7]

# TriangleExcenters returns excenters of triangle
def TriangleExcenters(p1,p2,p3):
    p4 = TriangleBevanCircleCenter(p1,p2,p3)
    [p5,p6,p7] = ExtouchTriangle(p1,p2,p3)
    p8 = TriangleIncenter(p1,p2,p3)
    p9 = Intersect(p1,p8,p4,p5)
    p10 = Intersect(p2,p8,p4,p6)
    p11 = Intersect(p3,p8,p4,p7)
    return [p9,p10,p11]

# TriangleExcirclesRadicalCircleSquaredRadius returns squared radius of radical circle of excircles
def TriangleExcirclesRadicalCircleSquaredRadius(a,b,c):
    return (1/4)*(a^2*b + a*b^2 + a^2*c + a*b*c + b^2*c + a*c^2 + b*c^2)/(a + b + c)

# TriangleExcirclesTangents returns tangency point between edge lines and excircles
def TriangleExcirclesTangents(p1,p2,p3):
    d12 = RT_Distance(p1,p2)
    d13 = RT_Distance(p1,p3)
    d23 = RT_Distance(p2,p3)
    s = TriangleSemiperimeter(p1,p2,p3)
    p12 = ExtendPoint(p1,p2,s) ; p13 = ExtendPoint(p1,p3,s)
    p21 = ExtendPoint(p2,p1,s) ; p23 = ExtendPoint(p2,p3,s)
    p31 = ExtendPoint(p3,p1,s) ; p32 = ExtendPoint(p3,p2,s)
    return [p12,p13,p21,p23,p31,p32]

# TriangleStevanovicCircleCenter returns Stevanovic circle center
# It's circle orthogonal to nine main triangle circles
def TriangleStevanovicCircleCenter(p1,p2,p3):
    # Compute trilinear polar of Gergonne point
    p4 = TriangleGergonnePoint(p1,p2,p3)
    [p41,p42,p43] = CevianTriangle(p1,p2,p3,p4)
    [p44,p45,p46] = TrianglePerspectrix(p1,p2,p3,p41,p42,p43)
    # Compute trilinear polar of Nagel point
    [p51,p52,p53] = ExtouchTriangle(p1,p2,p3)
    [p54,p55,p56] = TrianglePerspectrix(p1,p2,p3,p51,p52,p53)
    # Intersect the the trilinear polars
    p6 = Intersect(p44,p45,p54,p55)
    return p6

# TriangleStevanovicCircleRadius returns Stevanovic circle radius
def TriangleStevanovicCircleRadius(p1,p2,p3):
    d12 = RT_Distance(p1,p2)
    d13 = RT_Distance(p1,p3)
    d23 = RT_Distance(p2,p3)
    q0 = d23*d13*d12*(d23^2+d13^2+d12^2)-d23^4*(d13+d12-d23)-d13^4*(d12+d23-d13)-d12^4*(d23+d13-d12)
    q = d23*d13*d12*q0/(4*(d23-d13)^2*(d13-d12)^2*(d12-d23)^2)
    r = sqrt(q)
    return r

# Numeric equalizer for a triangle ... a CAR (Compass And Ruler) algorithm
# assuming d12 <= d23 <= d13
def Triangle_Equalizer(p1,p2,p3):
    d12 = RT_Distance(p1,p2)
    d13 = RT_Distance(p1,p3)
    d23 = RT_Distance(p2,p3)
    s = (d12 + d13 + d23)/2  # semiperimeter
    d14 = (s - sqrt(s^2 - 2*d12*d13))/2 ; q14 = d14^2
    d15 = d12*d13/(2*d14) ; q15 = d15^2
    p4 = Extend(p1,p2,q14)
    p5 = Extend(p1,p3,q15)
    return [p4,p5]

# Numeric equalizer for a triangle : find sub-triangle half area
# assuming d12 <= d23 <= d13
# assuming p1,p2,p3 numerical values with square root allowed
def Triangle_Equalizer_Fully_Simplified(p1,p2,p3):
    d12 = RT_Distance(p1,p2) # edge length
    d13 = RT_Distance(p1,p3) # edge length
    d23 = RT_Distance(p2,p3) # edge length
    s = (d12 + d13 + d23)/2  # semiperimeter
    disc = ((4*s^2 - 8*d12*d13).full_simplify()).factor()
    t = (((2*s - sqrt(disc))/4).full_simplify()).factor()  # edge length of sub-triangle
    q14 = (t^2).full_simplify() # its quadrance
    h = (((d12 + d13 + d23)*(d12 + d13 - d23)*(d12 - d13 + d23)*(-d12 + d13 + d23)).full_simplify()).factor()
    area = sqrt(h)/4 # area main triangle
    q15 = ((((-d23^2 + d12^2 + d13^2)/(4*t))^2 +  (area/t)^2).full_simplify()).factor() # quadrance other edge sub-triangle
    p4 = RT_Extend(p1,p2,q14)
    p5 = RT_Extend(p1,p3,q15)
    return [p4,p5]

# Symbolic equalizer for a triangle : find sub-triangle half area
# assuming d12 <= d23 <= d13
def RT_Triangle_Equalizer(p1,p2,p3):
    d12 = RT_Distance(p1,p2) # edge length
    d13 = RT_Distance(p1,p3) # edge length
    d23 = RT_Distance(p2,p3) # edge length
    s = (d12 + d13 + d23)/2  # semiperimeter
    disc = 4*s^2 - 8*d12*d13
    t = (2*s - sqrt(disc))/4   # edge length of sub-triangle
    q14 = t^2 # its quadrance
    h = (d12 + d13 + d23)*(d12 + d13 - d23)*(d12 - d13 + d23)*(-d12 + d13 + d23)
    area = sqrt(h)/4 # area main triangle
    q15 = ((-d23^2 + d12^2 + d13^2)/(4*t))^2 +  (area/t)^2 # quadrance other edge sub-triangle
    p4 = RT_Extend(p1,p2,q14)
    p5 = RT_Extend(p1,p3,q15)
    return [p4,p5]

# Dissect triangle into seven (next function : eight) acute isosceles triangle using Hogatt and Densen way
# returns p4 incenter, cutting points of pentagon on edges, and vertex pp3
def RT_Triangle_Dissect_To_Isosceles(p1,p2,p3):
    p4 = TriangleIncenter(p1,p2,p3)
    q = RT_Quadrance(p4,p3)
    p5 = RT_Arc(p3,p1,p4,q)
    p6 = RT_Arc(p2,p1,p4,q)
    p7 = RT_Arc(p3,p2,p4,q)
    p8 = RT_Arc(p1,p2,p4,q)
    return [p4,p5,p6,p7,p8,p3]

def Triangle_Dissect_To_Isosceles(p1,p2,p3):
    [p4,p5,p6,p7,p8,p9] = RT_Triangle_Dissect_To_Isosceles(p1,p2,p3)
    # Check internal angles
    q34 = RT_Quadrance(p3,p4); q56 = RT_Quadrance(p5,p6); q78 = RT_Quadrance(p7,p8)
    if (q56 - 2*q34 > 0) or (q78 - 2*q34 > 0):
        q13 = RT_Quadrance(p1,p3); q23 = RT_Quadrance(p2,p3)
        if q13 < q23:
            p10 = RT_Arc(p1,p2,p1,q13)
            [p4,p5,p6,p7,p8,p9] = RT_Triangle_Dissect_To_Isosceles(p2,p3,p10)
            p9 = p10
        else:
            p10 = RT_Arc(p2,p1,p2,q23)
            [p4,p5,p6,p7,p8,p9] = RT_Triangle_Dissect_To_Isosceles(p3,p1,p10)
            p9 = p10
    return [p4,p5,p6,p7,p8,p9]

# TrianglesPerspectrix returns the perspectrix line with respect to two triangles in perspective with a line
# To check if two triangles are in perspective with a line : respective intersecion points of line edges are lying
# on the same line
def TrianglePerspectrix(p11,p12,p13,p21,p22,p23):
    p1 = Intersect(p12,p13,p22,p23)
    p2 = Intersect(p11,p13,p21,p23)
    p3 = Intersect(p11,p12,p21,p22)
    return [p1,p2,p3]

# Steiner Circumellipse
def TriangleSteinerCircumellipseZ(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2) ; q13 = RT_Quadrance(p1,p3) ; q23 = RT_Quadrance(p2,p3)
    SquaredZ = q12^2 + q13^2 + q23^2 - q12*q13 - q12*q23 - q13*q23
    Z = sqrt(SquaredZ)
    return Z

def TriangleSteinerCircumellipseSemiAxislengths(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2) ; q13 = RT_Quadrance(p1,p3) ; q23 = RT_Quadrance(p2,p3)
    Z =  TriangleSteinerCircumellipseZ(p1,p2,p3)
    MajorLength = (1/3)*sqrt(q12 + q13 + q23 + 2*Z)
    MinorLength = (1/3)*sqrt(q12 + q13 + q23 - 2*Z)
    return [MinorLength,MajorLength]

def TriangleSteinerCircumellipseFocalLength(p1,p2,p3):
    Z =  TriangleSteinerCircumellipseZ(p1,p2,p3)
    c = (2/3)*sqrt(Z)
    return c

# CyclicQuad returns a cyclic quadrilateral centered at a point p0 and rotated by angle given by tangent of half this angle
# given edge lengths  (from Garza's paper FG201803.pdf)
# When t = 0, horizontal base of cyclic quad is the edge between 2nd and 3rd vertices, in order to have my usual drawing of
# ABC triangle when quad is degenerated (d=0)
def CyclicQuad(a,b,c,d,p0,t):
    [x0,y0] = p0 ; sn = 2*t/(1 + t^2); cs = (1 - t^2)/(1 + t^2)
    qa = a^2; qb = b^2; qc = c^2; qd = d^2
    s = (a + b + c + d)/2
    AreaMax2 = (s - a)*(s - b)*(s - c)*(s - d); AreaMax = sqrt(AreaMax2)
    m = 1/2*((qb + qd) - (qa + qc))
    p = 2*b*(qb - m - qc); q = 4*AreaMax*b; r = (m - qb)^2 + qb*qa + 4*AreaMax2 - qc*(qb + qa)
    u = r*p/(q^2 + p^2); v = r/q - r*p^2/(q*(q^2 + p^2))
    w = (m*(b - u) + 2*AreaMax*v)/(qb + qa - 2*b*u); z = (2*AreaMax*(b - u) - m*v)/(qb + qa - 2*b*u)
    p2 = [x2,y2] = [0,0]; p1 = [x1,y1] = [cs*u - sn*v,sn*u + cs*v]
    p4 = [x4,y4] = [cs*w - sn*z,sn*w + cs*z]; p3 = [x3,y3] = [cs*b,sn*b]
    [x9,y9] = PolygonCircumcenterOfMass([p1,p2,p3,p4,p1])
    p1 = [x1 - x9 + x0,y1 - y9 + y0]; p2 = [x2 - x9 + x0,y2 - y9 + y0]
    p3 = [x3 - x9 + x0,y3 - y9 + y0]; p4 = [x4 - x9 + x0,y4 - y9 + y0]
    return [p1,p2,p3,p4]

# Origin cyclic quad returns vertices of cyclic quad ABCD with B at origin and BC parallel x-axis
def OriginCyclicQuad(a,b,c,d):
    return CyclicQuad(a,b,c,d,[0,0],0)

# OnSphere returns condition (if zero) for five points  on a sphere(Beyer 1987)
# setting p1 = (x,y,z) for equation in x,y,z
def OnSphere(p1,p2,p3,p4,p5):
    [x1,y1,y1] = p1; [x2,y2,y2] = p2; [x3,y3,y3] = p3; [x4,y4,y4] = p4; [x5,y5,y5] = p5
    m = matrix(SR,5,5)
    m[0] = [ x1^2 + y1^2 + z1^2, x1, y1, z1, 1 ]
    m[1] = [ x2^2 + y2^2 + z2^2, x2, y2, z2, 1 ]
    m[2] = [ x3^2 + y3^2 + z3^2, x3, y3, z3, 1 ]
    m[3] = [ x4^2 + y4^2 + z4^2, x4, y4, z4, 1 ]
    m[4] = [ x5^2 + y5^2 + z5^2, x5, y5, z5, 1 ]
    delta = m.determinant()
    return delta

# RichmondPentagon returns a pentagon inscribed in circle center p1, one vertice is p2
# done using Richmond construction
def RichmondPentagon(p1,q1,p2):
    # Build two chords points
    p3 = Extend(p1,p2,q1)
    p4 = Symetric(p1,p3)
    # Bissec p3-p1-p4 angle and get another chord point
    p5 = Bissect(p4,p1,p3,q1)
    # Construct middle point of p1-p5
    p6 = Middle(p1,p5)
    # Bissect p1-p6-p3 and intersect with p1-p3
    p7 = Bissect(p1,p6,p3,q1)
    p8 = Intersect(p6,p7,p1,p3)
    # Construct chord thru p8
    l = RT_Altitude_Line(p1,p3,p8)
    p9 = RT_Any_Point(l,p8)
    [p10,p11] = Chord(p1,q1,p8,p9)
    # Intersect two circles to get last two vertices
    q3 = RT_Quadrance(p3,p10)
    [p12,p13] = CirclesIntersect(p1,q1,p10,q3)
    [p14,p15] = CirclesIntersect(p1,q1,p11,q3)
    return [p3,p10,p13,p14,p11]

# Pentagon returns a pentagon inscribed in circle center p1, one vertice is p2
def Pentagon(p1,q1,p2):
    # Constants Richmond's construction
    q2 = 1.381966011 * q1
    q3 = 3.618033989 * q1
    # Get point in direction p2 at quadrance q1 from p1
    p3 = Extend(p1,p2,q1)
    # Construct two vertices at quadrance q2 = edge length Vi Vi+1
    [p4,p5] = CirclesIntersect(p1,q1,p3,q2)
    # Construct two vertices at quadrance q3 = length Vi Vi+2
    [p6,p7] = CirclesIntersect(p1,q1,p3,q3)
    return [p3,p5,p7,p6,p4]

# Heptagon returns a heptagon inscribed in circle center p1, one vertice is p2
def Heptagon(p1,q1,p2):
    # Constants Archimede's construction
    q2 = 0.753020396279769 * q1
    q3 = 2.44504186790693  * q1
    q4 = 3.80193773580560  * q1
    # Get point in direction p2 at quadrance q1 from p1
    p3 = Extend(p1,p2,q1)
    # Construct two vertices at quadrance q2 = edge length Vi Vi+1
    [p4,p5] = CirclesIntersect(p1,q1,p3,q2)
    # Construct two vertices at quadrance q3 = length Vi Vi+2
    [p6,p7] = CirclesIntersect(p1,q1,p3,q3)
    # Construct two vertices at quadrance q3s = length Vi Vi+3
    [p8,p9] = CirclesIntersect(p1,q1,p3,q4)
    return [p3,p5,p7,p9,p8,p6,p4]


# Rectangle_QuadrancesCut returns quadrance to cut rectangle axb such as 1 <= a/b <= 2
def Rectangle_QuadrancesCut(a,b):
    # quadrance side of square
    q1 = a*b
    # quadrance for point on the border
    q2 = a*(b-a)
    # quadrance for point inside
    q3 = b*(b-a)
    return [q1,q2,q3]

# Rectangle_Split returns splitting points for rectangle axb such as 1 <= a/b <= 2
def Rectangle_Split(p1,p2,p3,p4):
    # Compute side lengths  w < h
    h = RT_Distance(p1,p2)
    w = RT_Distance(p2,p3)
    # Check ratio of side length
    r = QQ(h / w)
    assert((1 <=  r) and (r <= 2))
    # Compute quadrances for splitting and merging
    [q1,q2,q3] = Rectangle_QuadrancesCut(w,h)
    # Compute point on border and interior point
    p5 = Extend(p3,p4,q2)
    p6 = Extend(p2,p5,q3)
    return [p5,p6]

# Rectangle_Merge returns points for square with parts of rectangle axb such as 1 <= a/b <= 2
def Rectangle_Merge(p1,p2,p3,p4,p5,p6):
    # Compute points for square
    p7 = OppositeExtend(p2,p6,RT_Quadrance(p5,p6))
    p8 = Middle(p1,p7)  # orthocenter of square
    p9 = Symetric(p8,p6)
    p10 = Extend(p1,p2,RT_Quadrance(p3,p5))
    return [p7,p9,p10]

# BambooShoot returns a heptagon with shape bamboo shoot, given one vertice p1 and cutting segment
def BamboShoot(p1,q1,p2):
    # Constants for construction
    q2 = q1 / 1.81521228562209
    q3 = 2.34729814281104 * q2
    # Get point in direction p2 at quadrance q1 from p1
    p3 = Extend(p1,p2,q1)
    [p4,p5] = CirclesIntersect(p1,q3,p3,q2)
    [p6,p7] = CirclesIntersect(p3,q3,p1,q2)
    [p8,p9] = CirclesIntersect(p5,q2,p6,q2)
    return [p1,p6,p8,p5,p3,p4,p7]

# GenericBambooShoot returns a heptagon with shape bamboo shoot, given one vertice p1 and cutting segment
# parameters ratio of quadrances k1,k2,k3,k4 for shape definition
def GenericBamboShoot(k1,k2,k3,k4,p1,q1,p2):
    # Compute quadrances
    q2 = q1/k1 # qBD = qFG = qAE/k1
    q3 = k2*q2 # qBE = qAF
    q4 = k3*q2 # qAB = qEF
    q5 = k4*q2 # qDE = qAG
    # Get point in direction p2 at quadrance q1 from p1
    p3 = Extend(p1,p2,q1)
    [p9,p4] = CirclesIntersect(p1,q4,p3,q3)
    [p9,p5] = CirclesIntersect(p4,q2,p3,q5)
    [p9,p6] = CirclesIntersect(p3,q4,p1,q3)
    [p7,p9] = CirclesIntersect(p1,q5,p6,q2)
    [p9,p8] = CirclesIntersect(p4,q2,p5,q2)
    return [p1,p4,p8,p5,p3,p6,p7]

# SoddyianLastVertex returns last vertex p4 of Soddyian triangle p1-p2-p4 given one edge p1-p2 and one altitude point p3
def SoddyianLastVertex(p1,p2,p3):
    p5 = Middle(p1,p2) ; q15 = RT_Quadrance(p1,p5)
    l = RT_Altitude_Line(p1,p2,p3); p6 = RT_Any_Point(l,p3)
    p7 = Arc(p3,p6,p5,q15)
    q13 = RT_Quadrance(p1,p3) ; q23 = RT_Quadrance(p2,p3)
    [p8,p9] = RT_Altitudes(p1,p2,p1,q13)
    [p10,p11] = RT_Altitudes(p1,p2,p2,q23)
    p12 = Intersect(p1,p2,p8,p7)
    p13 = Intersect(p1,p2,p10,p7)
    p14 = TriangleCircumcenter(p7,p12,p13)
    q147 = RT_Quadrance(p14,p7)
    p15 = Arc(p7,p3,p14,q147)
    q315 = RT_Quadrance(p3,p15)
    [p16,p17] = Chord(p3,q315,p1,p2)
    q117 = RT_Quadrance(p1,p17); q216 = RT_Quadrance(p2,p16)
    [p18,p4] = CirclesIntersect(p1,q117,p2,q216)
    return p4

# SoddyianParameter returns parameter t given area and semiperimeter s
# t is solution of k = area/s^2 = t^2 (t+1)^2/(t^2 + t + 1)^3
def SoddyianParameter(area,s):
    k = area/s^2
    g1 = (1/6*sqrt(1/3)*sqrt((27*k - 4)/k)/k + 1/2/k)^(1/3)
    g2 = g1 + 1/3/(k*g1) - 1
    # solve t + 1/t = g2 or t^2 - g2*t + 1 = 0
    g3 = g2^2 - 4
    t1 = (g2 - sqrt(g3))/2
    t2 = (g2 + sqrt(g3))/2
    return [t1,t2]

# SoddyianParameteredVertices returns vertices given by parameter t and other vertices p1,p2
def SoddyianParameteredVertices(p1,p2,t):
    q12 = RT_Quadrance(p1,p2)
    q13 = q12 * ( (2*t^2 + 2*t + 1)/((t^2 + 1)*(t + 1)^2) )^2
    q23 = q12 * ( (t^2 + 2*t + 2)*t^2/((t^2 + 1)*(t + 1)^2) )^2
    [p3,p4] = CirclesIntersect(p1,q13,p2,q23)
    return [p3,p4]

# RT_SoddyianParameteredVertices is same as SoddyianParameteredVertices for symbolic computations
def RT_SoddyianParameteredVertices(p1,p2,t):
    q12 = RT_Quadrance(p1,p2)
    q13 = q12 * ( (2*t^2 + 2*t + 1)/((t^2 + 1)*(t + 1)^2) )^2
    q23 = q12 * ( (t^2 + 2*t + 2)*t^2/((t^2 + 1)*(t + 1)^2) )^2
    [p3,p4] = RT_CirclesIntersect(p1,q13,p2,q23)
    return [p3,p4]

# SoddyianTriangle returns a Soddyian triangle with length c of p1-p2 edge and given parameter t
def SoddyianTriangle(c,t):
    p1 = [0,0] ; p2 = [c,0]
    x = c*(1 + t + 2*t^2)/ ( (1 + t)*(1 + t^2)^2 )
    y = c*2*t^2*(1 + t + t^2) / ( (1 + t)^2*(1 + t^2)^2 )
    p3 = [x,y]
    return [p1,p2,p3]

# SoddyianQuadrances returns quadrance of edges
def SoddyianQuadrances(c,t):
    [p1,p2,p3] = SoddyianTriangle(c,t)
    q12 = RT_Quadrance(p1,p2)
    q13 = RT_Quadrance(p1,p3)
    q23 = RT_Quadrance(p2,p3)
    return [q12,q13,q23]

# SoddyianDiameter returns diameter of circumcircle
def SoddyianDiameter(c,t):
    D = 1/2*c*(2*t^2 + 2*t + 1)*(t^2 + 2*t + 2)/((t^2 + t + 1)*(t + 1)^2)
    return D

# SoddyianIncenter returns incenter
def SoddyianIncenter(c,t):
    x =  c/(t^2 + 1)
    y =  c*t^2/((t^2 + t + 1)*(t^2 + 1))
    p1 = [x,y]
    return p1

# SoddyianCentroid returns centroid
def SoddyianCentroid(c,t):
    [p1,p2,p3] = SoddyianTriangle(c,t)
    p4 =  TriangleCentroid(p1,p2,p3)
    return p4

# SoddyianNagelPoint returns Nagel point
# for Soddyian triangle Nagel point has barycentric coordinates (1/t + 1)^2 : (t + 1)^2 : 1
def SoddyianNagelPoint(c,t):
    [p1,p2,p3] = SoddyianTriangle(c,t)
    p4 =  TriangleBarycentric(p1,p2,p3 ,(1/t + 1)^2,(t + 1)^2,1)
    return p4

# SoddyianSpiekerCenter returns Spieker center
# for Soddyian triangle Spieker center has barycentric coordinates (t + 1/t + 1)^2 + (1/t + 1)^2 : (t + 1/t + 1)^2 + (t + 1)^2 : (t + 1/t + 1)^2 + 1
def SoddyianSpiekerCenter(c,t):
    [p1,p2,p3] = SoddyianTriangle(c,t)
    p4 =  TriangleBarycentric(p1,p2,p3 ,(t + 1/t + 1)^2 + (1/t + 1)^2,(t + 1/t + 1)^2 + (t + 1)^2,(t + 1/t + 1)^2 + 1)
    return p4

# SoddyanNinePointCenter returns center of nine point circle
def SoddyianNinePointCenter(c,t):
    x =  1/4*(t^5 + t^4 + 2*t^3 + 6*t^2 + 3*t + 3)*c/((t^2 + 1)^2*(t + 1))
    y =  1/8*(2*t^3 + 3*t^2 + 4*t + 1)*(t^3 + 4*t^2 + 3*t + 2)*c*t/((t^2 + t + 1)*(t^2 + 1)^2*(t + 1)^2)
    p1 = [x,y]
    return p1


# RevereSoddyianTriangle returns reflected Soddyian triangle
def ReverseSoddyianTriangle(c,t):
    p1 = [0,0] ; p2 = [c,0]
    x = c*(1 + t + 2*t^2)/ ( (1 + t)*(1 + t^2)^2 )
    y = -c*2*t^2*(1 + t + t^2) / ( (1 + t)^2*(1 + t^2)^2 )
    p3 = [x,y]
    return [p1,p2,p3]

# RevereSoddyianQuadrilateral returns Soddyian quadrilateral made of two Soddyian trinagles
def SoddyianQuadrilateral(c,t,u):
    [p1,p3,p4] = SoddyianTriangle(c,t)
    [p1,p3,p2] = ReverseSoddyianTriangle(c,u)
    return [p1,p2,p3,p4]

# SoddyianQuadrilateral_Tangential_Losange_Parameter returns parameter t of two congruent Soddyian
# triangles joined by one common edge. The result is a losange made of four Pythagorean 3-4-5 triangles
# It's one case where the quadrilateral is tangential (has incircle)
def SoddyianQuadrilateral_Tangential_Losange_Parameter():
    t = 1
    return t

# SoddyianQuadrilateral_Tangential_Losange returns losange tangential made of Soddyian triangles
def SoddyianQuadrilateral_Tangential_Losange(c):
    t = SoddyianQuadrilateral_Tangential_Losange_Parameter()
    [p1,p2,p3,p4] = SoddyianQuadrilateral(c,t,t)
    return [p1,p2,p3,p4]

# SoddyianQuadrilateral_Tangential_Kyte_Parameter returns parameter t of two congruent Soddyian
# triangles joined by one common edge. The result is a kyte.
# It's one case where the quadrilateral is tangential (has incircle)
def SoddyianQuadrilateral_Tangential_Kyte_Parameter():
    w = (1/9)*(36*sqrt(113) + 388)^(1/3)
    z = sqrt((81*w^2 + 9*w + 16)/w)
    t = 1/18*(z + sqrt(-81*w - 16/w + 1458/z + 18) - 9)
    return t

# SoddyianQuadrilateral_Tangential_Losange returns kyte made of Soddyian triangles
def SoddyianQuadrilateral_Tangential_Kyte(c):
    t = SoddyianQuadrilateral_Tangential_Kyte_Parameter()
    [p1,p2,p3,p4] = SoddyianQuadrilateral(c,t,t)
    return [p1,p2,p3,p4]

# Functions for Soddyian triangle given parameter z = t + 1/t
def SoddyianTriangle_Semiperimeter(c,z):
    s = c*(z + 1)^2/(z*(z + 2))
    return s

def SoddyianTriangle_Area(c,z):
    a = c^2*(z + 1)/(z^2*(z + 2))
    return a

def SoddyianTriangle_Height(c,z):
    h = c*(z + 1)/(z^2*(z + 2))
    return h

def SoddyianTriangle_Ratio_Area_SquaredSemiperimeter(c,z):
    k = (z + 2)/(z + 1)^3
    return k

def SoddyianTriangle_Inradius(c,z):
    r = c/(z*(z + 1))
    return r

def SoddyianTriangle_Circumradius(c,z):
    R = c*(2*z^2 + 6*z + 5)/(4*(z + 1)*(z + 2))
    return R

def SoddyianTriangle_Diameter(c,z):
    D = 2 * SoddyianTriangle_Circumradius(c,z)
    return D

def SoddyianTriangle_Cosine(z):
    cs = - (2*z + 3)/(2*z^2 + 6*z + 5)
    return cs

def SoddyianTriangle_Sine(z):
    sn = 2*(z + 1)*(z + 2)/(2*z^2 + 6*z + 5)
    return sn

def SoddyianTriangle_Parameter(z):
    t = (z + sqrt(z^2 - 4))/2
    return t

def SoddyianTriangle_One_Turn(z1):
    t = SoddyianTriangle_Parameter(z1)
    u = (-t-1)/t
    z2 = u + 1/u
    return z2

def SoddyianTriangle_Two_Turn(z1):
    t = SoddyianTriangle_Parameter(z1)
    v=  1/(-t-1)
    z3 = v + 1/v
    return z3

def SoddyianTriangle_Distance_CircumcircleChord_Vertice(c,z):
    d = c*(2*z^2 + 3*z + 2)/(4*z^2*(z + 1))
    return d

def SoddyianTriangle_Distance_CircumcircleChord_Incenter(c,z):
    d = c*(2*z^2 + 7*z + 8)/(4*z*(z + 1)*(z + 2))
    return d

def SoddyianTriangle_Quadrance_Circumcenter_Incenter(c,z):
    q = c^2*(2*z^3 + 6*z^2 - 3*z - 16)*(2*z^2 + 6*z + 5)/(16*z*(z + 1)^2*(z + 2)^2)
    return q

# Soddyian triangle cosines and sines parameter t displaying the nice properties of angles swap

def SoddyianTriangle_Cosine_Longest(t):
    z = t + 1/t
    cs = SoddyianTriangle_Cosine(z)  # cs =  - (2*t^2 + 3*t + 2)*t/((2*t^2 + 2*t + 1)*(t^2 + 2*t + 2))
    return cs

def SoddyianTriangle_Cosine_Shortest_Positive(t):
    t1 = - (t + 1)
    cs = SoddyianTriangle_Cosine_Longest(t1)  # cs = (2*t^2 + t + 1)*(t + 1)/((2*t^2 + 2*t + 1)*(t^2 + 1))
    return cs

def SoddyianTriangle_Cosine_Shortest_Negative(t):
    t1 = - (t + 1)/t
    cs = SoddyianTriangle_Cosine_Longest(t1) # cs = (t^2 + t + 2)*(t + 1)*t/((t^2 + 2*t + 2)*(t^2 + 1))
    return cs

def SoddyianTriangle_Sine_Longest(t):
    z = t + 1/t
    sn = SoddyianTriangle_Sine(z) # sn =  2*(t^2 + t + 1)*(t + 1)^2/((2*t^2 + 2*t + 1)*(t^2 + 2*t + 2))
    return sn

def SoddyianTriangle_Sine_Shortest_Positive(t):
    t1 = - (t + 1)
    sn = SoddyianTriangle_Sine_Longest(t1)  # sn = 2*(t^2 + t + 1)*t^2/((2*t^2 + 2*t + 1)*(t^2 + 1))
    return sn

def SoddyianTriangle_Sine_Shortest_Negative(t):
    t1 = - (t + 1)/t
    sn = SoddyianTriangle_Sine_Longest(t1)  # sn = 2*(t^2 + t + 1)/((t^2 + 2*t + 2)*(t^2 + 1))
    return sn

# SoddyianTriangle_QuadArea_Factors(t,z) returns factors in ratio of quadarea and square of quadrance longest edge
# triangle ABC longest edge AB with a = BC, b = AC, c = AB
# from Heron formula : 16 QuadArea = (a+b+c)(a+b-c)(a+c-b)(b+c-a)
# or                   16 QuadArea / c^4 = ( (a+b+c)/c )((a+b-c)/c )((a+c-b)/c )((b+c-a)/c )
# but for near_Soddyian triangle, factors are given in terms of parameters t and z = t + 1/t
#                      16 QuadArea / c^4 = ( 2(z+1)^2/z(z+2) )( 2/z(z+2))(2t/z )(2/tz)
# We returns factor (a+b+c)/c = 2(z+1)^2/z(z+2) ; (a+b-c)/c = 2/z(z+2) ; (a+c-b)/c = 2t/z ; (b+c-a)/c = 2/tz
# Note : Area / c^2 =  (z+1)/z^2(z+2)
def SoddyianTriangle_QuadArea_Factors(t,z):
    f1 = 2*(z+1)^2/(z*(z+2))  # (a+b+c)/c
    f2 = 2/(z*(z+2))  # (a+b-c)/c
    f3 = 2*t/z  # (a+c-b)/c
    f4 = 2/(t*z)  # (b+c-a)/c
    return [f1,f2,f3,f4]

# TangentialQuadrilateral returns a tangential quadrilateral p3p4p5p6 given one direction p1p2 for AB
# a tangential quadrilateral is defined by 6 numbers (a,b,c,d,e,f) and edge lengths are sum of four of these numbers
#
# AB = a + e + f + b
# BC = b + f + e + c
# CD = c + e + f + d
# AD = d + f + e + a
# AC = (a + e) + (e + c)
def TangentialQuadrilateral(p1,p2,a,b,c,d,e,f):
    l34 = a + e + f + b ; q34 = l34^2
    l45 = b + f + e + c ; q45 = l45^2
    l56 = c + e + f + d ; q56 = l56^2
    l36 = d + f + e + a ; q36 = l36^2
    l35 = (a + e) + (e + c) ; q35 = l35^2
    p3 = p1
    p4 = Extend(p1,p2,q34)
    [p9,p5] = CirclesIntersect(p3,q35,p4,q45)
    [p9,p6] = CirclesIntersect(p3,q36,p5,q56)
    return [p3,p4,p5,p6]

# Routh_Triangle returns cevians cut and the cevian triangle of p1-p2-p3 triangle
# with the cutting ratios k1 (for p2-p3 segment), k2 (for p1-p3 segment), k3 (for p1-p2 segment)
def Routh_Triangle(p1,k1,p2,k2,p3,k3):
    p4 = LineBarycentric(p2,p3,k1,1-k1)
    p5 = LineBarycentric(p3,p1,k2,1-k2)
    p6 = LineBarycentric(p1,p2,k3,1-k3)
    p7 = Intersect(p1,p4,p3,p6)
    p8 = Intersect(p2,p5,p1,p4)
    p9 = Intersect(p2,p5,p3,p6)
    return [p4,p5,p6,p7,p8,p9]

def Routh_Quadrances(p1,k1,p2,k2,p3,k3):
    k,u,v,w = var('k,u,v,w')
    psi(u,v) = u*v - u + 1
    theta(k,u,v,w) = k^2*u - k*u - k*v + k*w + v
    alpha(u,v,w) = u*v + u*w + v*w - u - v - w + 1
    gamma(u,v,w) = 2*u*v*w - alpha(u,v,w)
    #
    q23 = RT_Quadrance(p2,p3)
    q13 = RT_Quadrance(p1,p3)
    q12 = RT_Quadrance(p1,p2)
    #
    q14 = theta(k1,q23,q13,q12)
    q25 = theta(k2,q13,q12,q23)
    q36 = theta(k3,q12,q23,q13)
    #
    q78 = q14*(gamma(k1,k2,k3)/(psi(k1,k2)*psi(k3,k1)))^2
    q89 = q25*(gamma(k1,k2,k3)/(psi(k1,k2)*psi(k2,k3)))^2
    q79 = q36*(gamma(k1,k2,k3)/(psi(k3,k1)*psi(k2,k3)))^2
    #
    q18 = q14*(k2/psi(k1,k2))^2
    q29 = q25*(k3/psi(k2,k3))^2
    q37 = q36*(k1/psi(k3,k1))^2
    #
    q17 = q14*((k3-1)/psi(k3,k1))^2
    q28 = q25*((k1-1)/psi(k1,k2))^2
    q39 = q36*((k2-1)/psi(k2,k3))^2
    #
    return [q23,q13,q12,q14,q25,q36,q78,q89,q79,q18,q29,q37,q17,q28,q39]

# TCubic_Interpolation returns two control points p4,p5 given p1,p2,p3 points
# p1,p2 are two consecutive points on the curve
# p3 is point at intersection of tangents
# p4,p5 are Bezier control points and Bezier path is p1,p4,p5,p2

# Using ratios of lengths for take care of scaling
#d = (-1/2*(b^2*c + b*c^2 - sqrt(4*a^2 - 3*b^2 + 6*b*c - 3*c^2)*b*c)/(a^2 - b^2 + b*c - c^2))
#d = (-1/2*(b^2*c + b*c^2 - 2*a*b*c*sqrt(1 - 3/4*(b/a - c/a)^2))/(a^2 - b^2 + b*c - c^2))
#d = (-1/2*(b^2*c + b*c^2) + a*b*c*sqrt(1 - 3/4*(b/a - c/a)^2))/(a^2 - b^2 + b*c - c^2)
#d/a = (-1/2*((b/a)^2*(c/a) + (b/a)*(c/a)^2) + (b/a)*(c/a)*sqrt(1 - 3/4*(b/a - c/a)^2))/(1 - (b/a)^2 + (b/a)*(c/a) - (c/a)^2)
#
# L = a*r4 is length of curve

def TCubic_Parameters(q12,q13,q23):
    r1 = sqrt(q23/q12)  # r1 = b/a
    r2 = sqrt(q13/q12)  # r2 = c/a
    disc = 1 - 3/4*(r1 - r2)^2
    if disc < 0.0:
        # No interpolation, we return 0.0
        r3 = 0.0
        r4  = 0.0  # r4 = L/a
    else:
        r3 = (-1/2*(r1^2*r2 + r1*r2^2) + r1*r2*sqrt(disc))/(1 - r1^2 + r1*r2 - r2^2) # d/a
        #L =  1/2*(a^2*b^2*c - b^4*c + 2*a^2*b*c^2 + 2*b^3*c^2 + b^2*c^3 - 2*b*c^4 + a^4*d - 3*a^2*b^2*d + 2*b^4*d - 2*a^2*b*c*d - 2*b^3*c*d - 2*a^2*c^2*d - 3*b^2*c^2*d + 2*b*c^3*d + c^4*d)/((a^2 - b^2 + b*c - c^2)*b*(c - d))
        r4 =  1/2*(r1^2*r2 - r1^4*r2 + 2*r1*r2^2 + 2*r1^3*r2^2 + r1^2*r2^3 - 2*r1*r2^4 + r3 - 3*r1^2*r3 + 2*r1^4*r3 - 2*r1*r2*r3 - 2*r1^3*r2*r3 - 2*r2^2*r3 - 3*r1^2*r2^2*r3 + 2*r1*r2^3*r3 + r2^4*r3)/((1 - r1^2 + r1*r2 - r2^2)*r1*(r2 - r3))
    return [r1,r2,r3,r4]

def TCubic_Interpolation_Barycentric(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2).n(digits=nd) ; q13 = RT_Quadrance(p1,p3).n(digits=nd) ; q23 = RT_Quadrance(p2,p3).n(digits=nd)
    [r1,r2,r3,r4] = TCubic_Parameters(q12,q13,q23)
    if r4 == 0.0:
        [k41,k42,k43] = [1.0,0.0,0.0]  # No T-cubic interpolation, returns barycentric or p1,p2
        [k51,k52,k53] = [0.0,1.0,0.0]  # No T-cubic interpolation, returns barycentric or p1,p2
    else:
        # p4 has trilinear coordinates [(c-d)/a,d/b,0] with respect to (p3,p1,p2) or p4 = (c-d)/c p3 + d/c p1
        # p5 has trilinear coordinates [(b-d)/a,0,d/c] or p5 = (b-d)/b p3 + d/b p2
        [k41,k42,k43] = [r3/r2,0,1-r3/r2]
        [k51,k52,k53] = [0,r3/r1,1-r3/r1]
    return [[k41,k42,k43],[k51,k52,k53]]

def TCubic_Interpolation(p1,p2,p3):
    [[k41,k42,k43],[k51,k52,k53]] = TCubic_Interpolation_Barycentric(p1,p2,p3)
    p4 = TriangleBarycentric(p1,p2,p3,k41,k42,k43)
    p5 = TriangleBarycentric(p1,p2,p3,k51,k52,k53)
    return [p4,p5]


# RT_Quadrance3 returns quadrance for quadrance in 3D (R^3)
def RT_Quadrance3(p1,p2):
    [x1,y1,z1] = p1 ; [x2,y2,z2] = p2
    q = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
    return q

# CayleyMengerDeterminant3 returns determinant for simplex of 3 points
def CayleyMengerDeterminant3(q12,q13,q23):
    m = matrix(SR,4,4)
    m[0] = [   0, q12, q13, 1 ]
    m[1] = [ q12,   0, q23, 1 ]
    m[2] = [ q13, q23,   0, 1 ]
    m[3] = [   1,   1,   1, 0 ]
    delta = m.determinant()
    return delta

# SimplexSquaredArea returns squared area using CayleyMengerDeterminant3
def RT_SimplexSquaredArea(q12,q13,q23):
    delta = CayleyMengerDeterminant3(q12,q13,q23)
    area2 = (-1/16)*delta
    return area2

# SimplexSquaredArea returns squared area using CayleyMengerDeterminant3
def SimplexSquaredArea(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2) ; q13 = RT_Quadrance(p1,p3) ; q23 = RT_Quadrance(p2,p3)
    area2 = RT_SimplexSquaredArea(q12,q13,q23)
    return area2

# CayleyMengerDeterminant4 returns determinant for simplex of 4 points
def CayleyMengerDeterminant4(q12,q13,q14,q23,q24,q34):
    m = matrix(SR,5,5)
    m[0] = [   0, q12, q13, q14, 1 ]
    m[1] = [ q12,   0, q23, q24, 1 ]
    m[2] = [ q13, q23,   0, q34, 1 ]
    m[3] = [ q14, q24, q34,   0, 1 ]
    m[4] = [   1,   1,   1,   1, 0 ]
    delta = m.determinant()
    return delta

# SimplexSquaredVolume returns signed volume using CayleyMengerDeterminant4
def SimplexSquaredVolume(p1,p2,p3,p4):
    q12 = RT_Quadrance3(p1,p2) ; q13 = RT_Quadrance3(p1,p3) ; q14 = RT_Quadrance3(p1,p4)
    q23 = RT_Quadrance3(p2,p3) ; q24 = RT_Quadrance3(p2,p4) ; q34 = RT_Quadrance3(p3,p4)
    delta = CayleyMengerDeterminant4(q12,q13,q14,q23,q24,q34)
    volume2 = 1/288*delta
    return volume2

# QuadrilateralAffineTransformation returns coefficients of matrix M for transformation of quadrilateral P1-P2-P3-P4
# to quadrilateral A(0,0) B(0,1) C(1,0) D(d_x,d_y)
#       (-k_x*s*t + c*k_x       c*k_x*t + k_x*s      -(k_x*s*t - c*k_x)*u_x + (c*k_x*t + k_x*s)*u_y    )
# M  =  (-k_y*s                 c*k_y                -k_y*s*u_x + c*k_y*u_y                            )
#       (   0,                    0,                              1                                    )
#

def QuadAffineParameters(p1,p2,p3,p4):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3 ; [x4,y4] = p4
    u_x = -x1 ; u_y = -y1  # translation vector of P1 to origine A(0,0)
    if (u_y + y2) == 0:
        # No rotation
        s = 0 ; c = 1
    else:
        cotg =  (u_x + x2)/(u_y + y2) # cotg(theta) where theta is angle rotation at origin A
        z = - cotg - sqrt(1 + cotg^2) # compute z = tg(theta/2)
        s = 2*z/(1 + z^2) ; c = (1 - z^2)/(1 + z^2) # sine and cosine parametrization, for sin(theta) and cos(theta)
    t =  -( -c*(x1 - x4) - s*(y1 - y4) ) /  ( s*(x1 - x4) - c*(y1 - y4) )  # shear parallel with x-axis
    k_x = 1 / ( -(s*(u_x + x2) - c*(u_y + y2))*t + c*(u_x + x2) + s*(u_y + y2) )  # scale x coordinates
    k_y = 1 / ( -s*(u_x + x4) + c*(u_y + y4) )  # scale y coordinates
    d_x =  -((s*(u_x + x3) - c*(u_y + y3))*t - c*(u_x + x3) - s*(u_y + y3))*k_x   # x-coordinate of D point
    d_y = -(s*(u_x + x3) - c*(u_y + y3))*k_y  # y-coordinate of D point
    return [d_x,d_y,u_x,u_y,k_x,k_y,c,s,t]

def QuadrilateralAffineTransformation(p1,p2,p3,p4):
    [d_x,d_y,u_x,u_y,k_x,k_y,c,s,t] = QuadAffineParameters(p1,p2,p3,p4)
    Real_Matrix = MatrixSpace(RR,3,3)
    M = Real_Matrix()
    M[0] = (  -k_x*s*t + c*k_x,   c*k_x*t + k_x*s,  -(k_x*s*t - c*k_x)*u_x + (c*k_x*t + k_x*s)*u_y   )
    M[1] = (      -k_y*s,                c*k_y,                -k_y*s*u_x + c*k_y*u_y                )
    M[2] = (         0,                    0,                              1                         )
    return M

def QuadrilateralAffineMap(M,p1):
    [x1,y1] = p1
    Real_Vector = VectorSpace(RR,3)
    p1v = Real_Vector([x1, y1, 1])
    p2v = M*p1v
    [x2,y2] = [p2v[0],p2v[1]] ; p2 = [x2,y2]
    return p2

def QuadrilateralAffineUnmap(M,p1):
    [x1,y1] = p1
    Real_Vector = VectorSpace(RR,3)
    p1v = Real_Vector([x1, y1, 1])
    p2v = M.inverse()*p1v
    [x2,y2] = [p2v[0],p2v[1]] ; p2 = [x2,y2]
    return p2

# same functions but symbolic
def RT_QuadrilateralAffineTransformation(p1,p2,p3,p4):
    [d_x,d_y,u_x,u_y,k_x,k_y,c,s,t] = QuadAffineParameters(p1,p2,p3,p4)
    Symbolic_Matrix = MatrixSpace(SR,3,3)
    M = Symbolic_Matrix()
    M[0] = (  -k_x*s*t + c*k_x,   c*k_x*t + k_x*s,  -(k_x*s*t - c*k_x)*u_x + (c*k_x*t + k_x*s)*u_y   )
    M[1] = (      -k_y*s,                c*k_y,                -k_y*s*u_x + c*k_y*u_y                )
    M[2] = (         0,                    0,                              1                         )
    return M

def RT_QuadrilateralAffineMap(M,p1):
    [x1,y1] = p1
    Symbolic_Vector = VectorSpace(SR,3)
    p1v = Symbolic_Vector([x1, y1, 1])
    p2v = M*p1v
    [x2,y2] = [p2v[0],p2v[1]] ; p2 = [x2,y2]
    return p2

def RT_QuadrilateralAffineUnmap(M,p1):
    [x1,y1] = p1
    Symbolic_Vector = VectorSpace(SR,3)
    p1v = Symbolic_Vector([x1, y1, 1])
    p2v = M.inverse()*p1v
    [x2,y2] = [p2v[0],p2v[1]] ; p2 = [x2,y2]
    return p2

# returning the last point D
# D1 =  ((x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4)/(x2*y1 - x4*y1 - x1*y2 + x4*y2 + x1*y4 - x2*y4), (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3)/(x2*y1 - x4*y1 - x1*y2 + x4*y2 + x1*y4 - x2*y4))
# D1 =  (( (x3*y1 + x4*y3 + x1*y4) - (x1*y3 + x3*y4 + x4*y1)) /( (x4*y2 + x1*y4 + x2*y1) - (x2*y4 + x4*y1 + x1*y2 ), ( (x1*y3 + x2*y1 + x3*y2) - (y1*x3 + y2*x1 + y3*x2))/( ( x4*y2 + x1*y4 + x2*y1) - (y4*x2 + y1*x4 + y2*x1 ))

# TriangleShoeLaceArea returns area of a triangle following the shoelace's formula for triangle P1P2P3 counterclocwise
# It's a "signed" area : if we take vertices clockwise then we can negate the area
def TriangleShoeLaceArea(p1,p2,p3):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3
    return 1/2*( (x1*y2 + x2*y3 + x3*y1) - (x2*y1 + x3*y2 + x1*y3) )

# QuadrilateralIsConvex checks if quad if convex by computing signed areas, and if same sign quad is convex
def QuadrilateralIsConvex(p1,p2,p3,p4):
    Area_P1P2P3 = TriangleShoeLaceArea(p1,p2,p3) ; s1 = sign(Area_P1P2P3)
    Area_P2P3P4 = TriangleShoeLaceArea(p2,p3,p4) ; s2 = sign(Area_P2P3P4)
    Area_P3P4P1 = TriangleShoeLaceArea(p3,p4,p1) ; s3 = sign(Area_P3P4P1)
    Area_P4P1P2 = TriangleShoeLaceArea(p4,p1,p2) ; s4 = sign(Area_P4P1P2)
    return ((s2 == s1)and(s3 == s1)and(s4 == s1))

# QuadrilateralLastPoint returns last point D(d_x,d_y) for quadrilateral affine transformation
def QuadrilateralLastPoint(p1,p2,p3,p4):
    Area_P3P4P1 = TriangleShoeLaceArea(p3,p4,p1)
    Area_P4P1P2 = TriangleShoeLaceArea(p4,p1,p2)
    Area_P1P2P3 = TriangleShoeLaceArea(p1,p2,p3)
    d_x = Area_P3P4P1 / Area_P4P1P2 ;  d_y = Area_P1P2P3 / Area_P4P1P2
    return [d_x,d_y]

# QuadrilateralShoeLaceArea returns area of a quadrilateral following the shoelace's formula for quad P1P2P3P4 counterclocwise
# It's a "signed" area : if we take vertices clockwise then we can negate the area
def QuadrilateralShoeLaceArea(p1,p2,p3,p4):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3 ; [x4,y4] = p4
    return 1/2*( (x1*y2 + x2*y3 + x3*y4 + x4*y1) - (x2*y1 + x3*y2 + x4*y3 + x1*y4) )

# PentagonShoeLaceArea returns area of a pentagon following the shoelace's formula for quad P1P2P3P4P5 counterclocwise
# It's a "signed" area : if we take vertices clockwise then we can negate the area
def PentagonShoeLaceArea(p1,p2,p3,p4,p5):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3 ; [x4,y4] = p4 ; [x5,y5] = p5
    return 1/2*( (x1*y2 + x2*y3 + x3*y4 + x4*y5 + x5*y1) - (x2*y1 + x3*y2 + x4*y3 + x5*y4 + x1*y5) )

# QuadilateralDiagonalsIntersect returns intersection of diagonals
def QuadilateralDiagonalsIntersect(p1,p2,p3,p4):
    [x1,y1] = p1; [x2,y2] = p2 ; [x3,y3] = p3 ; [x4,y4] = p4
    Area_124 = TriangleShoeLaceArea(p1,p2,p4)
    Area_234 = TriangleShoeLaceArea(p2,p3,p4)
    Area_123 = TriangleShoeLaceArea(p1,p2,p3)
    Area_341 = TriangleShoeLaceArea(p3,p4,p1)
    return [ (x4*Area_123 + x2*Area_341)/( Area_123 + Area_341) , (y3*Area_124 + y1*Area_234)/(Area_124  + Area_234)]


# LeonAnneLine returns two points on line locus of points
#    Area(LAB) + Area(LCD) = k*( Area(LAC) + Area(LBD) )
# When k = 1, one of the line is the "Newton Line", and we have from Leon Anne's theorem, two points on that line : middle of diagonals
# Line has equation :  (d_y - 1)*x  + (1 - d_x)*y == (- d_x + d_y*k)/(k + 1)
# and is member of family of parallel lines slope (1 - d_y)/(1 - d_x)

def LeonAnneLine(k,d_x,d_y):
    u = (- d_x + d_y*k)/(k + 1)
    if d_y == 1:
        x1 = 0 ; y1 = u/(1 - d_x) ; x2 = 1 ; y2 = y1
    else:
        if d_x == 1:
            x1 = u/(d_y - 1) ; y1 = 0 ; x2 = x1 ; y2 = 1
        else:
            x1 = 0 ; y1 = u/(1 - d_x) ; x2 = u/(d_y - 1) ; y2 = 0
    return [ [x1,y1], [x2,y2] ]

# LeonAnneExtrema returns extremal values for k, when line cross the vertices A,B,C or D
def LeonAnneExtrema(p1,p2,p3,p4):
    [d_x,d_y] = QuadrilateralLastPoint(p1,p2,p3,p4)
    k1 = d_x/d_y           # plug coordinates of A(0,0) into the line equation
    k2 = d_x + d_y - 1     # plug coordinates of B(1,0) into the line equation
    k3 = 1/k2              # plug coordinates of C(0,1) into the line equation
    k4 = 1/k1              # plug coordinates of D(d_x,d_y) into the line equation
    return [k1,k2,k3,k4]

# LeonAnneMinimum returns minimal value for k  (and maximal is 1/k)
def LeonAnneMinimum(p1,p2,p3,p4):
    [k1,k2,k3,k4] = LeonAnneExtrema(p1,p2,p3,p4)
    k = min(k1,k2,k3,k4)
    return k

# RT_LeonAnneThruPoint returns value for k , for which the k-generalized goes thru a point
def LeonAnneThruPoint(p1,p2,p3,p4,p5):
    M = RT_QuadrilateralAffineTransformation(p1,p2,p3,p4)  # affine transformation matrix
    [x,y] = RT_QuadrilateralAffineMap(M,p5)
    [d_x,d_y] = QuadrilateralLastPoint(p1,p2,p3,p4)
    u = ( (d_y - 1)*x  + (1 - d_x)*y)
    k = (d_x + u) / (d_y - u)
    return k

# QuadrilateralLeonAnne(k,p1,p2,p3,p4) returns one point (if parallelogram) or two points on the diagonals
# when linking the two points, we get a segment on the "generalized" Newton line
def QuadrilateralLeonAnne(k,p1,p2,p3,p4):
    [d_x,d_y] = QuadrilateralLastPoint(p1,p2,p3,p4)
    if d_x == 1 and d_y == 1:
        # Parallelogram : return intersection of diagonals, many lines
        p7 = Intersect(p1,p3,p2,p4) ; p8 = p7
    else:
        M = QuadrilateralAffineTransformation(p1,p2,p3,p4)  # affine transformation matrix
        [q1,q2] = LeonAnneLine(k,d_x,d_y)   # line in the image plane
        p5 = QuadrilateralAffineUnmap(M,q1) # unmap point to get point in original plane
        p6 = QuadrilateralAffineUnmap(M,q2) # unmap point to get point in original plane
        p7 = Intersect(p5,p6,p1,p3)
        p8 = Intersect(p5,p6,p2,p4)
    return [p7,p8]

def RT_QuadrilateralLeonAnne(k,p1,p2,p3,p4):
    M = RT_QuadrilateralAffineTransformation(p1,p2,p3,p4)  # affine transformation matrix
    [d_x,d_y] = QuadrilateralLastPoint(p1,p2,p3,p4)
    [q1,q2] = LeonAnneLine(k,d_x,d_y)   # line in the image plane
    p5 = RT_QuadrilateralAffineUnmap(M,q1) # unmap point to get point in original plane
    p6 = RT_QuadrilateralAffineUnmap(M,q2) # unmap point to get point in original plane
    p7 = Intersect(p5,p6,p1,p3)
    p8 = Intersect(p5,p6,p2,p4)
    return [p7,p8]

# TriangleLeonAnne(k,p1,p2,p3) returns two points of segment parallel with p1p2
# when linking the two points, we get a segment on the "generalized" Newton line
def TriangleLeonAnne(k,p1,p2,p3):
    [d_x,d_y] = [0,1]  # D = C in the affine transformation image
    M = QuadrilateralAffineTransformation(p1,p2,p3,p3)  # affine transformation matrix
    [q1,q2] = LeonAnneLine(k,d_x,d_y)   # line in the image plane
    p5 = QuadrilateralAffineUnmap(M,q1) # unmap point to get point in original plane
    p6 = QuadrilateralAffineUnmap(M,q2) # unmap point to get point in original plane
    p7 = Intersect(p5,p6,p1,p3)
    p8 = Intersect(p5,p6,p2,p3)
    return [p7,p8]

# PolygonArea returns area of a 2D polygon (p1,p2,...,pn,p1)
def PolygonArea(plist):
    area = 0
    for i in xrange(len(plist) - 1):
        xi = plist[i][0] ; yi = plist[i][1]
        ui = plist[i+1][0] ; vi = plist[i+1][1]
        area += xi*vi - ui*yi
    area = area/2
    return area

# PolygonCenterOfMass returns centroid (center of mass) CM of a closed 2D polygon  (p1,p2,...,pn,p1)
# Example :  Polygon = [(0,0), (2,0), (0,1), (0,0)] ; print "Center of Mass : ", PolygonCenterOfMass(Polygon)
# Polygon is a plist : list of couple of two cartesian coordinates
# We can use more complex list. For example : Polygon = [(tan(x),x) for x in srange(-2*float(pi),2*float(pi),0.01)]
def PolygonCenterOfMass(plist):
    area = PolygonArea(plist)
    x = 0 ; y = 0
    for i in xrange(len(plist) - 1):
        xi = plist[i][0] ; yi = plist[i][1]
        ui = plist[i+1][0] ; vi = plist[i+1][1]
        cross = (xi*vi - ui*yi)
        x += (xi + ui)*cross
        y += (yi + vi)*cross
    cm = [x/(6*area),y/(6*area)]
    return cm

# PolygonCircumcenterOfMass returns circumcenter of mass CCM of a 2D polygon  (p1,p2,...,pn,p1)
def PolygonCircumcenterOfMass(plist):
    area = PolygonArea(plist)
    x = 0 ; y = 0
    for i in xrange(len(plist) - 1):
        xi = plist[i][0] ; yi = plist[i][1]
        ui = plist[i+1][0] ; vi = plist[i+1][1]
        x += - yi*vi^2 + yi^2*vi + xi^2*vi - ui^2*yi
        y += - ui*yi^2 + xi*vi^2 + xi*ui^2 - xi^2*ui
    ccm = [x/(4*area),y/(4*area)]
    return ccm

# PolygonVisibility returns -1 if two points p1,p2 are visible (no segment [Pi,Pi+1] intersecting [p1,p2] segment)
# otherwise returns first i
# we skip Pi or Pi+1 = p1 to avoid special case prod1 = prod2 = 0
def PolygonVisibility(plist,p1,p2):
    [x1,y1] = p1 ; [x2,y2] = p2
    cross = -1
    for i in xrange(len(plist) - 1):
        xi = plist[i][0] ; yi = plist[i][1]
        xj = plist[i+1][0] ; yj = plist[i+1][1]
        if ((xi != x1)or(yi != y1)) and ((xj != x1)or(yj != y1)) and ((xi != x2)or(yi != y2)) and ((xj != x2)or(yj != y2)):
            prod1 = sign((x2-x1)*(yi-y1)-(y2-y1)*(xi-x1))*sign((x2-x1)*(yj-y1)-(y2-y1)*(xj-x1))
            prod2 = sign((xj-xi)*(y1-yi)-(yj-yi)*(x1-xi))*sign((xj-xi)*(y2-yi)-(yj-yi)*(x2-xi))
            if  (prod1 == 0)or(prod2 == 0)or((prod1 == -1) and (prod2 == -1)):
                cross = i
                break
    return cross

# PolygonLeastSquaresPoint returns point minimizing sum of quadrances to edge lines of a polygon [p1,p2,...pn,p1]
# for a triangle it is the symmedian point
def PolygonLeastSquaresPoint(plist):
    ll = mm = lm = ln = mn = 0
    for i in xrange(len(plist) - 1):
        xi = plist[i][0] ; yi = plist[i][1]
        ui = plist[i+1][0] ; vi = plist[i+1][1]
        li = vi - yi ; mi = xi - ui ; qi = sqrt(li^2 + mi^2) ; li = li / qi ; mi = mi / qi ; ni = li*xi + mi*yi
        ll += li^2 ; mm += mi^2 ; lm += li*mi ; ln += li*ni ; mn += mi*ni
    x0 = (mm*ln - lm*mn)/(ll*mm - lm^2) ; y0 = (ll*mn - lm*ln)/(ll*mm - lm^2)
    return [x0,y0]

# Tripolar coordinates of P point, x = PA y = PB z = PC for reference triangle ABC  a = BC  b = CA  c = AB
# From the following Euler equation in quadrances qx = x^2  qy = y^2  qz = z^2  qa = a^2  qb = b^2 qc = c^2 :
#    (qa + qb - qc)*(qx*qy + qc*qz) + (- qa + qb + qc)*(qy*qz + qa*qx) + (qa - qb + qc)*(qz*qx + qb*qy) - (qa*qx^2 + qb*qy^2 + qc*qz^2) - qa*qb*qc = 0
# if we know two quadrances about three quadrances(qx,qy,qz) then we can get the last one

# TriEuler_Tripolar_A returns two possible quadrances qx1 or qx2 given qy,qz
def TriEuler_Tripolar_A(qa,qb,qc,qy,qz):
    disc = (qa^2 + qb^2 + qc^2 - 2*qa*qb - 2*qa*qc - 2*qb*qc)*(qa^2 - 2*qa*qy + qy^2 - 2*qa*qz - 2*qy*qz + qz^2)
    k = qa^2 - qa*qb - qa*qc - qa*qy - qb*qy + qc*qy - qa*qz + qb*qz - qc*qz
    qx1 = -1/2*(k - sqrt(disc))/qa ; qx2 = -1/2*(k + sqrt(disc))/qa
    return [qx1,qx2]

# TriEuler_Tripolar_B returns two possible quadrances qy1 or qy2 given qx,qz
def TriEuler_Tripolar_B(qa,qb,qc,qx,qz):
    disc = (qa^2 + qb^2 + qc^2 - 2*qa*qb - 2*qa*qc - 2*qb*qc)*(qb^2 - 2*qb*qx + qx^2 - 2*qb*qz - 2*qx*qz + qz^2)
    k =  qb^2 - qb*qc - qb*qa - qb*qz - qc*qz + qa*qz - qb*qx + qc*qx - qa*qx
    qy1 = -1/2*(k - sqrt(disc))/qb ; qy2 = -1/2*(k + sqrt(disc))/qb
    return [qy1,qy2]

# TriEuler_Tripolar_C returns two possible quadrances qz1 or qz2 given qx,qy
def TriEuler_Tripolar_C(qa,qb,qc,qx,qy):
    disc = (qa^2 + qb^2 + qc^2 - 2*qa*qb - 2*qa*qc - 2*qb*qc)*(qc^2 - 2*qc*qx + qx^2 - 2*qc*qy - 2*qx*qy + qy^2)
    k =  qc^2 - qc*qa - qc*qb  - qc*qx - qa*qx + qb*qx  - qc*qy + qa*qy - qb*qy
    qz1 = -1/2*(k - sqrt(disc))/qc ; qz2 = -1/2*(k + sqrt(disc))/qc
    return [qz1,qz2]


# Functions for searching n-gon (3 <= n <= 11) circumradius, using a Newton's method for finding root
# of the equation, given the cosine of the sum of internal angles
# As starting value for the Newton's method, we use formula for a regular polygon with 11 sides
# circumradius = 1/2*sidelength*csc(pi/n) = 1/2*sidelength*1/sin(pi/11)
# and we set sidelength as the arithmetic mean of li
# but then we compute number of turns thru the circle...and TBD maybe scale using number of turns
def HendecagonCircumradius(l):
    [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11] = l
    x = var('x',domain='real')
    f(x) = 1 + cos(asin(l1/(2*x)) + asin(l2/(2*x)) + asin(l3/(2*x)) + asin(l4/(2*x)) + asin(l5/(2*x)) + asin(l6/(2*x)) + asin(l7/(2*x)) + asin(l8/(2*x)) + asin(l9/(2*x)) + asin(l10/(2*x)) + asin(l11/(2*x)))
    f1(x) = f.derivative(x)
    # solving f(x) = 0 by Newton's method
    #x1 = (max([l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11])).n()
    x1 = 1/2*(sum(l)/11)*(1/sin(pi.n()/11))
    nturns = sum(l)/(2*pi.n()*x1) #; print "nturns = ",nturns
    prec = 9.9999999999999998e-131
    while True:
        x1 = x1 - f(x1).n()/f1(x1).n() #;print "x1 = ",x1," f(x1) = ",f(x1)
        if abs(f(x1)) < prec:
            break;
    return x1

def NGonCircumradius(p):
    l = [None]*11
    for i in xrange(11):
        if i < len(p)-1:
            l[i] = RT_Distance(p[i],p[i+1])
        else:
            l[i] = 0.0
    x = HendecagonCircumradius(l)
    return x

# TriangleFractal4Points returns points in triangle p1-p2-p3 (=ABC) dissection into four triangles
# given as three points p4 (on p1p2 edge) p5 (on p2p3 edge) and p6 on (p1p3 edge)
# two bezier paths exist from p1 to p2 or p3
def TriangleFractal4Quadrances(q1,q2,q3):
    k2 = q3/q2 ; x = 1/(1 + k2); y = x^2 ; z = (1-x)^2
    q4 = y*q3 ; q5 = z*q3
    q6 = z*q1 ; q7 = y*q1
    q8 = z*q2 ; q9 = y*q2
    return [q4,q5,q6,q7,q8,q9]

def TriangleFractal4Points(p1,p2,p3):
    q1 = RT_Quadrance(p2,p3) ; q2 = RT_Quadrance(p1,p3) ; q3 = RT_Quadrance(p1,p2)
    [q14,q24,q25,q35,q16,q36] = TriangleFractal4Quadrances(q1,q2,q3)
    [p4,p9] = CirclesIntersect(p1,q14,p2,q24)
    [p5,p9] = CirclesIntersect(p2,q25,p3,q35)
    [p6,p9] = CirclesIntersect(p1,q16,p3,q36)
    return [p4,p5,p6]

# ---- Projective geometry : plane space P2 of homogeneous coordinates x:y:z

# x:y:z can be interpreted as cartesian coordinates in P2  (x/z,y/z)
#                       or as trilinear coordinates of points in plane wrt (= with reference triangle) ABC


# Cross-ratio [123] for homogeneous coordinates : if zero, three points colinear
def cxrl(p1,p2,p3):
    [x1,y1,z1] = p1; [x2,y2,z2] = p2; [x3,y3,z3] = p3
    M3 = matrix(SR,3,3,[x1,y1,z1,x2,y2,z2,x3,y3,z3])
    return M3.determinant()

# Cross-ratio [123456] for homogeneous coordinates : if zero, six points lying on conic
def cxrc(p1,p2,p3,p4,p5,p6):
    [x1,y1,z1] = p1; [x2,y2,z2] = p2; [x3,y3,z3] = p3
    [x4,y4,z4] = p4; [x5,y5,z5] = p5; [x6,y6,z6] = p6
    M6 = matrix(SR,6,6,[x1^2,x1*y1,x1*z1,y1^2,y1*z1,z1^2,
                    x2^2,x2*y2,x2*z2,y2^2,y2*z2,z2^2,
                    x3^2,x3*y3,x3*z3,y3^2,y3*z3,z3^2,
                    x4^2,x4*y4,x4*z4,y4^2,y4*z4,z4^2,
                    x5^2,x5*y5,x5*z5,y5^2,y5*z5,z5^2,
                    x6^2,x6*y6,x6*z6,y6^2,y6*z6,z6^2])
    return M6.determinant()

def RT_cxrc(p1,p2,p3,p4,p5,p6,p7,p8):
    k123 = cxrl(p1,p2,p3); k145 = cxrl(p1,p4,p5); k246 = cxrl(p2,p4,p6); k356 = cxrl(p3,p5,p6)
    k124 = cxrl(p1,p2,p4); k135 = cxrl(p1,p3,p5); k236 = cxrl(p2,p3,p6); k456 = cxrl(p4,p5,p6)
    return (k123*k145*k246*k356 - k124*k135*k236*k456)

# Cross-ratio (1234)_5678 = [567813]*[567824]/([567814]*[567823])
def cxlm(p1,p2,p3,p4,p5,p6,p7,p8):
    return cxrc(p5,p6,p7,p8,p1,p3)*cxrc(p5,p6,p7,p8,p2,p4)/(cxrc(p5,p6,p7,p8,p1,p4)*cxrc(p5,p6,p7,p8,p2,p3))

def RT_cxlm(p1,p2,p3,p4,p5,p6,p7,p8):
    return RT_cxrc(p5,p6,p7,p8,p1,p3)*RT_cxrc(p5,p6,p7,p8,p2,p4)/(RT_cxrc(p5,p6,p7,p8,p1,p4)*RT_cxrc(p5,p6,p7,p8,p2,p3))


# Cross-ratio [1;2345678] for homogeneous coordinates : if zero, eight points lying on cubic singular at first point
#  singular point p1 is either a double point (2 tangents) or a cusp (inflection point)

# TBD : cubs function is not working, mistake somewhere ... RT_cubs is valid function and quicker
#def cubs(p1,p2,p3,p4,p5,p6,p7,p8):
#    [x1,y1,z1] = p1; [x2,y2,z2] = p2; [x3,y3,z3] = p3
#    [x4,y4,z4] = p4; [x5,y5,z5] = p5; [x6,y6,z6] = p6
#    [x7,y7,z7] = p7; [x8,y8,z8] = p8
#    M10 = matrix(SR,10,10,[x2^3,   x2^2*y2, x2^2*z2, x2*y2^2, x2*y2*z2, x2*z2^2, y2^3,   y2^2*z2, y2*z2^2, z2^3,
#                           x3^3,   x3^2*y3, x3^2*z3, x3*y3^2, x3*y3*z3, x3*z3^2, y3^3,   y3^2*z3, y3*z3^2, z3^3,
#                           x4^3,   x4^2*y4, x4^2*z4, x4*y4^2, x4*y4*z4, x4*z4^2, y4^3,   y4^2*z4, y4*z4^2, z4^3,
#                           x5^3,   x5^2*y5, x5^2*z5, x5*y5^2, x5*y5*z5, x5*z5^2, y5^3,   y5^2*z5, y5*z5^2, z5^3,
#                           x6^3,   x6^2*y6, x6^2*z6, x6*y6^2, x6*y6*z6, x6*z6^2, y6^3,   y6^2*z6, y6*z6^2, z6^3,
#                           x7^3,   x7^2*y7, x7^2*z7, x7*y7^2, x7*y7*z7, x7*z7^2, y7^3,   y7^2*z7, y7*z7^2, z7^3,
#                           x8^3,   x8^2*y8, x8^2*z8, x8*y8^2, x8*y8*z8, x8*z8^2, y8^3,   y8^2*z8, y8*z8^2, z8^3,
#                           3*x1^2, 2*x1*y1, 2*x1*z1, y1^2,    y1*z1,    z1^2,    0,      0,       0,       0,
#                           0,      x1^2,    0,       2*x1*y1, x1*z1,    0,       3*y1^2, 2*y1*z1, z1^2,    0,
#                           0,      0,       x1^2,    0,       x1*y1,    2*x1*z1, 0,      y2^2,    2*y1*z1, 3*z1^2])
#    return M10.determinant()

def RT_cubs_params(p1,p2,p3,p4,p5,p6,p7,p8):
    k123 = cxrl(p1,p2,p3); k126 = cxrl(p1,p2,p6); k128 = cxrl(p1,p2,p8); k173 = cxrl(p1,p7,p3) ; k176 = cxrl(p1,p7,p6)
    k178 = cxrl(p1,p7,p8); k423 = cxrl(p4,p2,p3); k426 = cxrl(p4,p2,p6); k428 = cxrl(p4,p2,p8)
    k473 = cxrl(p4,p7,p3); k476 = cxrl(p4,p7,p6); k478 = cxrl(p4,p7,p8); k523 = cxrl(p5,p2,p3) ; k526 = cxrl(p5,p2,p6)
    k528 = cxrl(p5,p2,p8); k573 = cxrl(p5,p7,p3); k576 = cxrl(p5,p7,p6); k578 = cxrl(p5,p7,p8)
    k647 = cxrl(p6,p4,p7); k657 = cxrl(p6,p5,p7); k847 = cxrl(p8,p4,p7); k857 = cxrl(p8,p5,p7)
    return [k123,k126,k128,k173,k176,k178,k423,k426,k428,k473,k476,k478,k523,k526,k528,k573,k576,k578,k647,k657,k847,k857]

def RT_cubs(p7,p1,p2,p3,p4,p5,p6,p8):
    [k123,k126,k128,k173,k176,k178,k423,k426,k428,k473,k476,k478,k523,k526,k528,k573,k576,k578,k647,k657,k847,k857] = RT_cubs_params(p1,p2,p3,p4,p5,p6,p7,p8)
    return 3*(k647*k857*k478*k128*k173*k423*k573*k526*k176 -
              k647*k857*k473*k428*k178*k123*k573*k526*k176 +
              k647*k857*k473*k428*k178*k576*k126*k173*k523 +
              k657*k847*k573*k528*k178*k123*k473*k426*k176 -
              k657*k847*k578*k128*k173*k523*k473*k426*k176 -
              k657*k847*k573*k528*k178*k476*k126*k173*k423)


# RT_RectangularTriangleInradi returns 0 if qr is squared radius of one of the four incircle/excircles of A-rectangular triangle ABC
# with edge quadrances qb,qc, qa = qb + qc. See Wildbeger video: https://www.youtube.com/watch?v=iaKLyz7vxb0&t=1819s
def RT_RectangularTriangleInradi(qb,qc,qr):
    return 16*qr^4 - 32*(qb + qc)*qr^3 + 8*(2*qb^2 + 3*qb*qc + 2*qc^2)*qr^2 - 8*qb*qc*(qb + qc)*qr + qb^2*qc^2

# --- Cyclic polygon drawing

# Defining function equation for circumcenter inside n-gon  (translating f(r))
#   s is the array of sides lengths
#   r is one value near circumradius
def InsideCircumcenter(s,r):
    eq = 0
    for l in s:
        eq += acos(1 - l^2/(2*r^2))
    eq -= 2*pi.n()
    return eq

# Defining function equation for circumcenter outside n-gon  (translating f2(r))
#   s is the array of sides lengths
#   r is one value near circumradius
def OutsideCircumcenter(s,r,n,indmax):
    eq = 0
    for i in range(n):
        l = s[i]
        if i != indmax:
            k = 1
        else:
            k = -1
        eq += k*acos(1 - l^2/(2*r^2))
    return eq

# Defining my argmax function (in order not to import it from numpy) : return index of maximum value
def MyArgmax(l):
    n = len(l)
    indmax = 0; argmax = l[0]
    for i in xrange(2,n):
        if l[i] > argmax:
            argmax = l[i]
            indmax = i
    return indmax

# CyclicPolygonVertices returns vertices for a drawing of a cyclic polygon given side lengths
def CyclicPolygonVertices(SidesLengths):
    n = len(SidesLengths)
    # Compute radius of circumcircle.
    # Quasi-Newton method start with rleft as max side/2. Function is not defined for smaller r
    indmax = MyArgmax(SidesLengths); amax = SidesLengths[indmax]
    rleft = amax/2.0
    r = rleft
    # Quasi-Newton iteration
    # First find rright such that f(rleft) has different sign from f(rright).
    # If not it means the centre of the circumcircle is outside the polygon and we have to use f2.
    flag = 0; it = 2; rright = it * rleft
    vleft = InsideCircumcenter(SidesLengths,rleft)
    vright = InsideCircumcenter(SidesLengths,rright)
    signsproduct = vleft * vright
    while it < 10 and abs(vleft) > 0.005 and signsproduct > 0:
        it += 1
        rright = it * rright
        vright = InsideCircumcenter(SidesLengths,rright)
        signsproduct = vleft * vright
    if ( signsproduct <= 0 ):
        flag = 1
    # If f did not have root look for root of f2:
    if flag == 0:
        it = 1
        rright = 2 * rleft
        vleft = OutsideCircumcenter(SidesLengths,rleft,n,indmax)
        vright = OutsideCircumcenter(SidesLengths,rright,n,indmax)
        signsproduct = vleft * vright
        while it < 30 and abs(vleft) > 0.001 and signsproduct > 0:
            it += 1
            rright = it * rright
            vright = OutsideCircumcenter(SidesLengths,rright,n,indmax)
            signsproduct = vleft * vright
    # Now do quasi newton:
    j = 1
    if flag == 1:
        v = InsideCircumcenter(SidesLengths,r)
        while j < 20 and abs(v) > 0.001 :
            rm = (rleft + rright) / 2.0
            vm = InsideCircumcenter(SidesLengths,rm)
            signsproduct = vleft * vm
            if signsproduct <= 0:
                rright = rm; vright = vm
            else:
                signsproduct = vm * vright
                if signsproduct <= 0:
                    rleft = rm; vleft = vm
                else:
                    assert false # failure f to find circumcenter
            r = rm
            j += 1
        #print "Circumcenter inside n-gon, and circumradius = ",r
    else:
        v = OutsideCircumcenter(SidesLengths,r,n,indmax)
        while j < 50 and abs(v) > 0.00000001:   # Improve accurency by increasing number of iterations (j max) and decreasing equation limit value
            rm = (rleft + rright) / 2.0
            vm = OutsideCircumcenter(SidesLengths,rm,n,indmax)
            signsproduct = vleft * vm
            if signsproduct <= 0:
                rright = rm; vright = vm
            else:
                signsproduct = vm * vright
                if signsproduct <= 0:
                    rleft = rm; vleft = vm
                else:
                    assert false # failure f2 to find circumcenter
            j += 1
        signsproduct = vleft * vright
        if  signsproduct <= 0:
            r = (rleft + rright) / 2.0
            flag = 2
        else:
            r = rm
        #print "Circumcenter outside n-gon, and circumradius = ",r
    #
    # Compute squared circumradius
    qr = r^2
    # If we used f2 and the longest side is the first then place the center below the axis
    # Otherwise place it above the axis
    x0 = SidesLengths[0]/2; y0 = sqrt(qr - x0^2)
    if flag == 2 and indmax == 0:
        y0 = - y0
    O = [x0,y0] ; pO = vector(O);
    # Compute vertices
    Vertices = n*[None]
    Vertices[0] = [0,0]
    Vertices[n - 1] = [2*x0,0]
    for k in xrange(1,n - 1):
        if flag == 2 and k == indmax:
            [Pk,P9] = CirclesIntersect(Vertices[k - 1],SidesLengths[k]^2,pO,qr)
        else:
            [P9,Pk] = CirclesIntersect(Vertices[k - 1],SidesLengths[k]^2,pO,qr)
        Vertices[k] = Pk
    return Vertices

print "...EuclideanGeometry module loaded"