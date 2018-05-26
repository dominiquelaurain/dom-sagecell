print "RationalTrigonometry (RT) module loading ..."

# Given two points p1,p1 returns quadrance Q(P1,P2)
def RT_Quadrance(p1,p2):
    [x1,y1] = p1; [x2,y2] = p2
    q = (x1-x2)^2+(y1-y2)^2
    return q

# Rational Trigonometry : distance of two points p1,p2
def RT_Distance(p1,p2):
    q = RT_Quadrance(p1,p2)
    return sqrt(q)

# Rational Trigonometry : line (a:b:c) going thru two points p1,p2
def RT_Line(p1,p2):
    [x1,y1] = p1; [x2,y2] = p2
    return [ y1-y2, x2-x1, x1*y2-x2*y1 ]

#  RT_Any_Point returns one point different from p1 on line l
def RT_Any_Point(l,p1):
    [x1,y1] = p1
    [a,b,c] = l  # line ax + by + c
    if x1 != 0:
        if b != 0:
            [x,y] = [0,-c/b]
        else:
            if a != 0:
                [x,y] = [-c/a,y1+1]
            else:
                [x,y] = [x1+1,-c/b]
    else:
        if a != 0:
            [x,y] = [-c/a,0]
        else:
            if b != 0:
                [x,y] = [x1+1,-c/b]
            else:
                [x,y] = [-c/a,y1+1]
    return [x,y]

#  RT_HalfTurn returns point on altitude line
#   p1p2 is segment of line p3 is point on segment, function returns p4 such as p3p4 orthogonal to p1p1
#   (y1-y2)x + (x1-x2)y + (x2-x1)x3 + (y2-y1)y3 = 0 is equation of line segment p1p1
#   (x1-x2)(x-x3) + (y2-y1)(y-y3) = 0 is equation of altitude thru p3
#  Setting x4 = x3 + (y1-y2)/2 and y4 = y3 + (x1-x2)/2 we have altitude point
def RT_HalfTurn(p1,p2,p3):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3
    x4 = x3 + (y1-y2)/2 ; y4 = y3 - (x1-x2)/2
    return vector( [x4,y4] )

# RT_QuarterTurn returns the point making a rectangular isosecles triangle with a segment (a turn 45 degrees)
# same idea than RT_HalfTurn where p3 = middle of p1-p2
def RT_QuarterTurn(p1,p2):
    [x1,y1] = p1 ; [x2,y2] = p2
    x4 = (x1 + x2 + y1 - y2)/2 ; y4 = (y1 + y2 - x1 + x2)/2
    return vector( [x4,y4] )

# Homogoneous barycentric coordinates of Quarter turn point N for triangle ABC
# edges lengths a = BC (base) ,b,c and triangle area S
def QuarterTurnBarycentrics(a,b,c,S):
    d = a/2 ; h = S/d  # half base length and height
    l = sqrt(c^2 - h^2) # distance between foot of A from B on base
    return [d/h, (h+l-2*d)/(2*h), (h-l)/(2*h)]

# RT_Line_Reflection returns reflection of point p1 on line l
def RT_Line_Reflection(l,p1):
    [x1,y1] = p1
    [a,b,c] = l  # line ax + by + c
    q = a^2 + b^2
    x = ((b^2 - a^2)*x1 - 2*a*b*y1 - 2*a*c)/q
    y = (-2*a*b*x1 + (a^2 - b^2)*y1 - 2*b*c)/q
    return [x,y]

# RT_Line_Symmetric returns reflection of point p3 on line l given by two points p1,p2
def RT_Line_Symmetric(p1,p2,p3):
    l = RT_Line(p1,p2)
    p4 = RT_Line_Reflection(l,p3)
    return p4

# Rational Trigonometry : rectangular projection point of p3 on line thru points p1,p2
def RT_Foot(p1,p2,p3):
    [a,b,c] =  RT_Line(p1,p2) ; [x,y] = p3
    d2 = a^2+b^2
    return [ (b^2*x -a*b*y -a*c)/d2, (-a*b*x+a^2*y-b*c)/d2 ]

# Rational Trigonometry : altitude from point p3 to line thru p1,p2
def RT_Altitude_Line(p1,p2,p3):
    [a,b,c] =  RT_Line(p1,p2) ; [x,y] = p3
    return [ -b, a, b*x-a*y ]

# Rational Trigonometry : rectangular altitudes points of point p3 on line thru points p1,p2
def RT_Altitudes(p1,p2,p3,q):
    [a,b,c] = RT_Altitude_Line(p1,p2,p3)   # line with p3 point rectangular with line thru p1-p2
    [x4,y4] = RT_Foot(p1,p2,p3) # rectangular projection point of p3 on p1-p2 line
    # search for (x,y) such as ax+by+c = 0  and (x-x4)^2 + (y-y4)^2 = q
    if (b != 0):
        # by = (-ax-c) and (bx-bx4)^2 + (ax+c+by4)^2 - qb^2 = 0
        # (a^2+b^2)*x^2 +2*(-b^2*x4 +ac+aby4)*x + b^2*x4^2 + (c+by4)^2 - q*b^2 = 0
        a1 = a^2 + b^2 ; b1 = -b^2*x4 +a*(c+b*y4); c1 = (b*x4)^2 + (c+b*y4)^2 -q*b^2 ; d1 = b1^2 -a1*c1
        xm = (-b1 -sqrt(d1))/a1 ; ym = (-a*xm -c)/b
        xp = (-b1 +sqrt(d1))/a1 ; yp = (-a*xp -c)/b
    else:
        # ax = (-by-c) and (by+c+ax4)^2 + (ay-ay4)^2 - qa^2 = 0
        # (a^2+b^2)*y^2 +2*(-a^2*y4 +bc+abx4)*y + a^2*y4^2 + (c+ax4)^2 - q*a^2 = 0
        a1 = a^2 + b^2 ; b1 = -a^2*y4 +b*(c+a*x4); c1 = (a*y4)^2 + (c+a*x4)^2 -q*a^2 ; d1 = b1^2 -a1*c1
        ym = (-b1 -sqrt(d1))/a1 ; xm = (-b*ym -c)/a
        yp = (-b1 +sqrt(d1))/a1 ; xp = (-b*yp -c)/a
    return [[xm,ym],[xp,yp]]

# Rational Trigonometry : relative polar spread s of point p2 with origin point p1
def RT_Polar_Spread(p1,p2):
    [x1,y1] = p1; [x2,y2] = p2
    return ((y1-y2)/(x1-x2))^2

# Rational Trigonometry : quadrancy of altitude point p3, given base line thru points p1,p2
def RT_Altitude_Quadrance(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3); q23 = RT_Quadrance(p2,p3)
    return (q13*q23 - ((q12-q13-q23)/2)^2)/q12

# Squared Area
# Reference : N.Wildberger  Symetric formula for edges quadrances Q1,Q2,Q3. Gives 0 iff three aligned points
# symetric formula for quadrances parameters
# 16 Area^2 = (Q1+Q2+Q3)^2 - 2*(Q1^2+Q2^2+Q2^2)
def RT_TriangleSignedSquaredArea(q1,q2,q3):
    D = (q1+q2+q3)^2 - 2*(q1^2+q2^2+q3^2)
    SquaredArea = D / 16
    return SquaredArea

def RT_TriangleSquaredArea(q1,q2,q3):
    D = (q1+q2+q3)^2 - 2*(q1^2+q2^2+q3^2)
    if D < 0:
        D = -D
    SquaredArea = D / 16
    return SquaredArea

def RT_TriangleArea(d12,d13,d23):
    q12 = d12^2 ; q13 = d13^2 ; q23 = d23^2
    S2 = RT_TriangleSquaredArea(q12,q13,q23)
    S = sqrt(S2)
    return S

# Squared Area of triangle p1-p2-p3
def TriangleQuadArea(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3); q23 = RT_Quadrance(p2,p3)
    QuadArea123 = RT_TriangleSquaredArea(q12,q13,q23)
    return QuadArea123

# Area of triangle p1-p2-p3
def TriangleArea(p1,p2,p3):
    QuadArea123 = TriangleQuadArea(p1,p2,p3)
    Area123 = sqrt(QuadArea123)
    return Area123

# Inradius of triangle p1-p2-p3
def TriangleInradius(p1,p2,p3):
    d12 = RT_Distance(p1,p2); d13 = RT_Distance(p1,p3); d23 = RT_Distance(p2,p3)
    p = d12 + d13 + d23 # Perimeter
    q = (d13 + d23 - d12)*(d23 + d12 - d13)*(d12 + d13 - d23)/p
    r = sqrt(q)
    return r

# Inradius of triangle p1-p2-p3 given area S and semi-perimeter from formula S = r s
# and s semiperimeter, half os sum of edge lengths
def RT_TriangleInradius(S,s):
    r = S / s
    return r

# Semiperimeter of triangle p1-p2-p3
def RT_TriangleSemiperimeter(d12,d13,d23):
    s = (d12 + d13 + d23)/2
    return s

def TriangleSemiperimeter(p1,p2,p3):
    d12 = RT_Distance(p1,p2); d13 = RT_Distance(p1,p3); d23 = RT_Distance(p2,p3)
    s = RT_TriangleSemiperimeter(d12,d13,d23)
    return s

# Law of cosines triangle p1-p2-p3 : cosine of angle at p3
def RT_TriangleCosineLaw(d12,d13,d23):
    cs = (d12^2 + d13^2 - d23^2)/(2*d12*d13)
    return cs

# Law of sines for triangle p1-p2-p3 : sine of angle at p3
def RT_TriangleSineLaw(d12,d13,d23):
    S = RT_TriangleArea(d12,d13,d23)   # S = rs  where r is inradius and s semi-perimeter
    R = (d12*d13*d23)/(4*S)            # circumradius R as product of edge lengths divided by 4rs = 4S
    D = 2*R                            # diameter D is two times R
    si = d23 / D                       # finally law of sines
    return si

# RT_TriangleHalfHarmonicMeans returns lengths of bisector segments, as harmonic means of two edge lengths
# it can be used for a 120 degrees triangle where one bisector length is harmonic of two edge lengths
def RT_TriTriangleHalfHarmonicMeans(a,b,c):
    return [1/b+1/c,1/c+1/a,1/a+1/b]

# Area of quadrilateral p1-p2-p3-p4
def QuadrilateralArea(p1,p2,p3,p4):
    return TriangleArea(p1,p2,p3)+TriangleArea(p1,p3,p4)

# Half perimeter of cyclic quad
def RT_CyclicQuadSemiperimeter(d12,d23,d34,d41):
    s = (d12+d23+d34+d41)/2
    return s

# Cross ratio (= anharmonic ratio) for four colinear points (p1,p2;p3,p4)
# we compute in complex numbers to get the "signed" cross-ratio
# When -1 result we have an harmonic quadrangle (p1,p2,p3,p4)
def RT_CrossRatio(p1,p2,p3,p4):
    z0 = CC(I)
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3; [x4,y4] = p4
    z1 = x1 + y1*z0; z2 = x2 + y2*z0; z3 = x3 + y3*z0; z4 = x4 + y4*z0
    return (((z3 - z1)*(z4 - z2))/((z3 - z2)*(z4 - z1))).real()

# Cyclic quad ABCD
# a = AB ; b = BC ; c = CD ; d = DA ; p = AC ; q = BD
# Ptolemy : p*q = (a*c + b*d)
# Ratios  : p/q = (a*d + b*c)/(a*b + c*d)
#
# multiplying and dividing the two equations giving p^2 and q^2

# Crossed products of cycle quad
def RT_CyclicQuadCrossProducts(d12,d23,d34,d41):
    z1 = d12*d23 + d34*d41 ; z2 = d12*d34 + d23*d41 ;z3 = d12*d41 + d23*d34
    return [z1,z2,z3]

# Squared diagonal lengths of a cyclic quadrilateral given by edges distances
def RT_CyclicQuadSquaredDiagonals(d12,d23,d34,d41):
    [z1,z2,z3] = RT_CyclicQuadCrossProducts(d12,d23,d34,d41)
    p2 = z1*z2/z3 ; q2 = z2*z3/z1 ; r2 = z3*z1/z2
    return [p2,q2,r2]

# Diagonal lengths of a cyclic quadrilateral given by edges distances
def RT_CyclicQuadDiagonals(d12,d23,d34,d41):
    [p2,q2,r2] = RT_CyclicQuadSquaredDiagonals(d12,d23,d34,d41)
    [z1,z2,z3] = RT_CyclicQuadCrossProducts(d12,d23,d34,d41)
    p = sqrt(p2) ; q = sqrt(q2) ; r = sqrt(r2)
    return [p,q,r]

# Squared circumradius a cyclic quadrilateral given by edges distances
def RT_CyclicQuadSquaredCircumradius(d12,d23,d34,d41):
    s = RT_CyclicQuadSemiperimeter(d12,d23,d34,d41)
    [z1,z2,z3] = RT_CyclicQuadCrossProducts(d12,d23,d34,d41)
    R2 = 1/16*(z1*z2*z3)/((s-d12)*(s-d23)*(s-d34)*(s-d41))
    return R2

# Squared area of a cyclic quadrilateral given by edges distances
def RT_CyclicQuadSquaredArea(d12,d23,d34,d41):
    [p2,q2,r2] = RT_CyclicQuadSquaredDiagonals(d12,d23,d34,d41)
    s = RT_CyclicQuadSemiperimeter(d12,d23,d34,d41)
    R2 = RT_CyclicQuadSquaredCircumradius(d12,d23,d34,d41)
    S2 = (p2*q2*r2)/(16*R2)
    return S2

# Area of a cyclic quadrilateral given by edges distances
def RT_CyclicQuadArea(d12,d23,d34,d41):
    S2 = RT_CyclicQuadSquaredArea(d12,d23,d34,d41)
    S = sqrt(S2)

# Formula to get edge lengths of bicentric quad (cyclic and tangential) from inradius r, semiperimeter s and circumradius R
# they are the four real roots of quartic  equation
#   y^4 - 2*s*y^3 + (s^2 + 2*r^2 + 2*r*sqrt(4*R^2 + r^2))*y^2 - 2*r*s*(sqrt(4*R^2 + r^2) + r)*y + r^2*s^2 = 0
def RT_BicentricQuadEdgeLengths(r,s,R):
    y = var('y')
    f(y) = y^4 - 2*s*y^3 + (s^2 + 2*r^2 + 2*r*sqrt(4*R^2 + r^2))*y^2 - 2*r*s*(sqrt(4*R^2 + r^2) + r)*y + r^2*s^2
    res = solve([f == 0],y)
    a = res[0].rhs()
    b = res[1].rhs()
    c = res[2].rhs()
    d = res[3].rhs()
    return [a,b,c,d]

# function from quartic equation for edge lengths bicentric quad
def RT_SquaredBicentricQuadInradius(a,b,c,d):
    s = (a + b + c + d)/2
    r2 = a*b*c*d/s^2
    return r2

# Semiperimeter of a triangle given by points p1,p2,p3
def TriangleSemiperimeter(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3); q23 = RT_Quadrance(p2,p3)
    s = (sqrt(q12) + sqrt(q13) + sqrt(q23))/2
    return s

# Semiperimeter of a triangle given by points p1,p2,p3,p4
def QuadrilateralSemiperimeter(p1,p2,p3,p4):
    q12 = RT_Quadrance(p1,p2); q23 = RT_Quadrance(p2,p3); q34 = RT_Quadrance(p3,p4); q14 = RT_Quadrance(p1,p4)
    s = (sqrt(q12) + sqrt(q23) + sqrt(q34) + sqrt(q14))/2
    return s

# Incenter of triangle given by points p1,p2,p3
def TriangleIncenter(p1,p2,p3):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    d12 = RT_Distance(p1,p2); d13 = RT_Distance(p1,p3); d23 = RT_Distance(p2,p3)
    d = d12 + d13 + d23
    x = (d23*x1 + d13*x2 + d12*x3)/d; y = (d23*y1 + d13*y2 + d12*y3)/d
    return [ x, y ]

# Inradius of triangle given by points p1,p2,p3
def TriangleInradius(p1,p2,p3):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    d12 = RT_Distance(p1,p2); d13 = RT_Distance(p1,p3); d23 = RT_Distance(p2,p3)
    s = (d12 + d13 + d23)/2
    r = sqrt((s-d12)*(s-d13)*(s-d23)/s)
    return r

# Barycentric function given two points and two weights
def LineBarycentric(p1,p2,w1,w2):
    w = w1 + w2
    return [ ((w1*p1[0])+(w2*p2[0]))/w,((w1*p1[1])+(w2*p2[1]))/w]

def LineCentroid(p1,p2):
    return LineBarycentric(p1,p2,1.0,1.0)

# Barycentric function given three points and three weights
def TriangleBarycentric(p1,p2,p3,w1,w2,w3):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    w = w1 + w2 + w3
    return vector( [ (w1*x1 + w2*x2 + w3*x3)/w, (w1*y1 + w2*y2 + w3*y3)/w ] )

def TriangleCentroid(p1,p2,p3):
    return TriangleBarycentric(p1,p2,p3,1.0,1.0,1.0)

# Barycentric function given four points and four weights
def QuadrilateralBarycentric(p1,p2,p3,p4,w1,w2,w3,w4):
    w = w1 + w2 + w3 + w4
    return [ ((w1*p1[0])+(w2*p2[0])+(w3*p3[0])+(w4*p4[0]))/w,((w1*p1[1])+(w2*p2[1])+(w3*p3[1])+(w4*p4[1]))/w]

# Centroid of quadrilateral is intersection of segments joining centroid of the four sub-triangles
# it is center of mass of quadrilateral
def QuadrilateralCentroid(p1,p2,p3,p4):
    p51 = TriangleCentroid(p2,p3,p4)
    p52 = TriangleCentroid(p1,p3,p4)
    p53 = TriangleCentroid(p1,p2,p4)
    p54 = TriangleCentroid(p1,p2,p3)
    p5 = Intersect(p51,p53,p52,p54)
    return p5

# Barycentric function given five points and five weights
def PentagonBarycentric(p1,p2,p3,p4,p5,w1,w2,w3,w4,w5):
    w = w1 + w2 + w3 + w4 + w5
    return [ ((w1*p1[0])+(w2*p2[0])+(w3*p3[0])+(w4*p4[0])+(w5*p5[0]))/w,((w1*p1[1])+(w2*p2[1])+(w3*p3[1])+(w4*p4[1])+(w5*p5[1]))/w]

def PentagonCentroid(p1,p2,p3,p4,p5):
    return PentagonBarycentric(p1,p2,p3,p4,p5,1.0,1.0,1.0,1.0,1.0)

# Orthocenter of triangle given by points p1,p2,p3
# intersection of altitudes
def TriangleOrthocenter(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3); q23 = RT_Quadrance(p2,p3)
    w1 = (-q12 + q13 + q23)*(q12 - q13 + q23)
    w2 = (-q12 + q13 + q23)*(q12 + q13 - q23)
    w3 = ( q12 - q13 + q23)*(q12 + q13 - q23)
    return TriangleBarycentric(p1,p2,p3,w1,w2,w3)

# Circumcenter of triangle given by points p1,p2,p3
def TriangleCircumcenter(p1,p2,p3):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3); q23 = RT_Quadrance(p2,p3)
    return TriangleBarycentric(p1,p2,p3,q23*(-q23+q12+q13),q13*(-q13+q12+q23),q12*(-q12+q13+q23))

# NinePoint center of triangle given by points p1,p2,p3
def TriangleNinePointCenter(p1,p2,p3):
    p4 = TriangleOrthocenter(p1,p2,p3)
    p5 = TriangleCircumcenter(p1,p2,p3)
    p6 = Middle(p4,p5)
    return p6

# Circumradius of triangle given by points p1,p2,p3
# from the propery that  product of inradius by circumradius is product of edge lengths divided by two
# times sum of edge lengths
def TriangleCircumradius(p1,p2,p3):
    #[x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    d12 = RT_Distance(p1,p2); d13 = RT_Distance(p1,p3); d23 = RT_Distance(p2,p3)
    s = (d12 + d13 + d23)/2
    ri = TriangleInradius(p1,p2,p3)
    rc = (d12*d13*d23)/(4*s*ri)
    return rc

# TriangleNinePointsCircleRadius returns radius of nine-point circle of triangle given by points p1,p2,p3
# simply half of circumradius
def TriangleNinePointsCircleRadius(p1,p2,p3):
    rc = TriangleCircumradius(p1,p2,p3)
    return rc/2

# Spieker center of triangle given by points p1,p2,p3
# Spieker is incenter of median triangle and intersection of cleavers (cleaver is segment from vertice splitting perimeter in two)
# It is center of mass of three wires with uniform mass
def TriangleSpiekerCenter(p1,p2,p3):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    d12 = RT_Distance(p1,p2); d13 = RT_Distance(p1,p3); d23 = RT_Distance(p2,p3)
    return TriangleBarycentric(p1,p2,p3 ,d12 + d13 , d12 + d23 , d13 + d23)

# Bevan circle goses thru  the three excenters of triangle

# Bevan circle center of triangle given by points p1,p2,p3
# It's name is "Bevan point"
def TriangleBevanCircleCenter(p1,p2,p3):
    p4 = TriangleCircumcenter(p1,p2,p3)
    p5 = TriangleIncenter(p1,p2,p3)
    p6 = Symetric(p4,p5)
    return p6

# Bevan circle radius of triangle given by points p1,p2,p3
# It is two times circumradius
def TriangleBevanCircleRadius(p1,p2,p3):
    R = TriangleCircumradius(p1,p2,p3)
    return 2*R

# Cleaver is segment halving perimeter, it is parallel of bissector p3 for triangle given by points p1,p2,p3
# it is also line thru middle of segment and Spieker center
# swapping arguments give three cleavers
def TriangleCleaver(p1,p2,p3):
    p4 = TriangleSpiekerCenter(p1,p2,p3)
    p5 = Middle(p1,p2)
    p6 = Intersect(p4,p5,p1,p3)
    return [p5,p6]

# Nagel point of triangle given by points p1,p2,p3
# Perimeter excepted, coordinates are factors in quadarea formula 16 area^2 = (a+b+c)(a+b-c)(a-b+c)(-a+b+c)
def TriangleNagelPoint(p1,p2,p3):
    [x1,y1] = p1; [x2,y2] = p2; [x3,y3] = p3
    d12 = RT_Distance(p1,p2); d13 = RT_Distance(p1,p3); d23 = RT_Distance(p2,p3)
    return TriangleBarycentric(p1,p2,p3 ,d12 + d13 - d23, d12 + d23 - d13, d13 + d23 - d12)

# Perimeter of quadrilateral given by four distances
def RT_QuadrilateralPerimeter(d12,d23,d34,d41):
    perimeter = d12 + d23 + d34 + d41
    return perimeter

# Perimeter of quadrilateral given by four points p1,p2,p3,p4 (in that order along perimeter)
def QuadrilateralPerimeter(p1,p2,p3,p4):
    d12 = RT_Distance(p1,p2); d23 = RT_Distance(p2,p3); d34 = RT_Distance(p3,p4); d41 = RT_Distance(p4,p1)
    perimeter = RT_QuadrilateralPerimeter(d12,d23,d34,d41)
    return perimeter

# Sum of diagonal lengths of a quadrilateral
def SumDiagonalsLengths(p1,p2,p3,p4):
    d13 = RT_Distance(p1,p3); d24 = RT_Distance(p2,p4)
    return d13 + d24

# Squared area of quadrilateral given by six quadrances
def RT_QuadrilateralSquaredArea(q12,q23,q34,q41,q13,q24):
    quadarea = 1/16*(4*q13*q24 - (q12 + q34 - q23 - q41)^2)
    return quadarea

# Squared area of parallelogram given by three quadrances (two consecutive edge lengths and one diagonal length)
def RT_ParallelogramSquaredArea(q12,q23,q13):
    quadarea = 1/4*(q13^2 - (q12 - q23)^2)
    return quadarea


# Squared area of quadrilateral given by four points p1,p2,p3,p4 (in that order along perimeter)
def QuadrilateralSquaredArea(p1,p2,p3,p4):
    q12 = RT_Quadrance(p1,p2); q23 = RT_Quadrance(p2,p3); q34 = RT_Quadrance(p3,p4); q41 = RT_Quadrance(p4,p1)
    q13 = RT_Quadrance(p1,p3); q24 = RT_Quadrance(p2,p4)
    quadarea = RT_QuadrilateralSquaredArea(q12,q23,q34,q41,q13,q24)
    return quadarea

def QuadrilateralArea(p1,p2,p3,p4):
    quadarea = QuadrilateralSquaredArea(p1,p2,p3,p4)
    area = sqrt(quadarea)
    return area

# QuadrilateralCyclicChordsRatios returns ratio of chords lengths, from Intersecting chords theorem
def QuadrilateralCyclicChordsRatios(p1,p2,p3,p4):
    p5 = Intersect(p1,p3,p2,p4)
    d15 = RT_Distance(p1,p5)
    d25 = RT_Distance(p2,p5)
    d35 = RT_Distance(p3,p5)
    d45 = RT_Distance(p4,p1)
    d14 = RT_Distance(p1,p4)
    d23 = RT_Distance(p2,p3)
    return [ d15/d25 , d45/d35 , d14/d23 ]

# Check if quadrilateral given by four points p1,p2,p3,p4 is cyclic (has a circumcircle)
# Ptolemy's theorem : pq = ac + bd
def QuadrilateralIsCyclic(p1,p2,p3,p4):
    d12 = RT_Distance(p1,p2); d23 = RT_Distance(p2,p3); d34 = RT_Distance(p3,p4); d41 = RT_Distance(p4,p1)
    d13 = RT_Distance(p1,p3); d24 = RT_Distance(p2,p4)
    iscyclic = d13*d24 == d12*d34 + d23*d41
    return iscyclic

# Anticenter of cyclic quadrilateral given by four points p1,p2,p3,p4
def QuadrilateralAnticenter(p1,p2,p3,p4):
    # Compute maltitudes (line from one midpoint to edge orthogonal to opposite edge)
    p5 = Middle(p4,p3); p6 = RT_Foot(p1,p2,p5)
    p7 = Middle(p1,p4); p8 = RT_Foot(p2,p3,p7)
    # Intersect the two maltitudes to get the anticenter
    p9 = Intersect(p5,p6,p7,p8)
    return p9

# Squared area of cyclic quadrilateral given by four points p1,p2,p3,p4
# using Brahmagupta formula

def RT_CyclicQuadrilateralQuadArea(d12,d23,d34,d41):
    s = (d12 + d23 + d34 + d41)/2
    quadarea = (s-d12)*(s-d23)*(s-d34)*(s-d41)
    return quadarea

def CyclicQuadrilateralQuadArea(p1,p2,p3,p4):
    d12 = RT_Distance(p1,p2); d23 = RT_Distance(p2,p3); d34 = RT_Distance(p3,p4); d41 = RT_Distance(p4,p1)
    quadarea = RT_CyclicQuadrilateralQuadArea(d12,d23,d34,d41)
    return quadarea

# Area of cyclic quadrilateral given by four points p1,p2,p3,p4
# using Brahmagupta formula
def CyclicQuadrilateralArea(p1,p2,p3,p4):
    quadarea = CyclicQuadrilateralQuadArea(p1,p2,p3,p4)
    area = sqrt(quadarea)
    return area

# Circumradius of cyclic quadrilateral given by four points p1,p2,p3,p4
def RT_SquaredQuadrilateralCircumradius(d12,d23,d34,d41):
    quadarea = RT_CyclicQuadrilateralQuadArea(d12,d23,d34,d41)
    q = (d12*d34 + d23*d41)*(d12*d23 + d34*d41)*(d23*d34 + d12*d41)/(16*quadarea)
    return q

def QuadrilateralCircumradius(p1,p2,p3,p4):
    d12 = RT_Distance(p1,p2); d23 = RT_Distance(p2,p3); d34 = RT_Distance(p3,p4); d41 = RT_Distance(p4,p1)
    q = RT_SquaredQuadrilateralCircumradius(d12,d23,d34,d41)
    r = sqrt(q)
    return r

# Squared cicumradius of cyclic quadrilateral given by four points p1,p2,p3,p4
# from NJ Wildberger math foundations 146 video
# p1 is chosen point origin on diameter

def RT_QuadrilateralQuadCircumradius(p1,p2,p3,p4):
    q23 = RT_Quadrance(p2,p3) ; q24 = RT_Quadrance(p2,p4) ; q34 = RT_Quadrance(p3,p4)
    quadarea =  (q23 + q24 + q34)^2 - 2*(q23^2 + q24^2 + q34^2)
    R = q23*q24*q34/quadarea
    return R

# Check if quadrilateral given by four points p1,p2,p3,p4 is cicumscriptible(=tangential)
# diagonals split quadrilateral in four triangles, we check sum of inverse inradii
def QuadrilateralIsCircumscriptible(p1,p2,p3,p4):
    p5 = Intsersect(p1,p3,p2,p4)
    r1 = TriangleInradius(p1,p2,p5) ; r3 = TriangleInradius(p3,p4,p5)
    r2 = TriangleInradius(p2,p3,p5) ; r4 = TriangleInradius(p4,p1,p5)
    iscircumscriptible = (1/r1 + 1/r3) == (1/r2 + 1/r4)
    return iscircumscriptible

# Inradius of circumscriptible quadrilateral given by four points p1,p2,p3,p4
def QuadrilateralInradius(p1,p2,p3,p4):
    d12 = RT_Distance(p1,p2); d23 = RT_Distance(p2,p3); d34 = RT_Distance(p3,p4); d41 = RT_Distance(p4,p1)
    p = d12 + d23 + d34 + d41 # Perimeter
    q = (d12*d23*d34 + d23*d34*d41 + d34*d41*d12 + d41*d12*d23)/p
    r = sqrt(q)
    return r

# Incenter of circumscriptible quadrilateral given by four points p1,p2,p3,p4
# all four bissectors intersect in single point, incenter
def QuadrilateralIncenter(p1,p2,p3,p4):
    p5 =  Bissect(p4,p1,p2,1)  # bissector angle at p1
    p6 =  Bissect(p1,p2,p3,1)  # bissector angle at p2
    p7 = Intersect(p1,p5,p2,p6)
    return p7

# Inradius of tangential quadrilateral given by four points p1,p2,p3,p4
def QuadrilateralTangentialInradius(p1,p2,p3,p4):
    s = QuadrilateralSemiperimeter(p1,p2,p3,p4)
    K = QuadrilateralArea(p1,p2,p3,p4)
    r = K / s
    return r

# Check if circumscriptible quadrilateral is cyclic  (see Hajja paper FG200814)
def QuadrilateralCircumscriptibleIscyclic(p1,p2,p3,p4):
    d12 = RT_Distance(p1,p2); d23 = RT_Distance(p2,p3); d34 = RT_Distance(p3,p4); d41 = RT_Distance(p4,p1)
    d13 = RT_Distance(p1,p3); d24 = RT_Distance(p2,p4)
    iscyclic = (d13 / d24) == (d12 + d34)/(d23 + d41)
    return iscyclic

# Returns squared circumradius for a cyclic quadrilateral given only by edge lengths - Brahmagupta's theorem extension to quad
# a1,a2,a3,a4 are edge lengths of a cyclic quadrilateral p1-p2-p3-p4 see Moritsugu's paper
def RT_QuadrilateralSquaredCircumradius(a1,a2,a3,a4):
    u = (a1*a2 + a3*a4)*(a1*a3 + a2*a4)*(a1*a4 + a2*a3)
    v = (-a1 + a2 + a3 + a4)*(a1 - a2 + a3 + a4)*(a1 + a2 - a3 + a4)*(a1 + a2 + a3 - a4)
    r2 = u/v
    return r2

def RT_TriangleSquaredCircumradius(a1,a2,a3):
    u = (a1*a2)*(a1*a3)*(a2*a3)
    v = (-a1 + a2 + a3)*(a1 - a2 + a3)*(a1 + a2 - a3)*(a1 + a2 + a3)
    r2 = u/v
    return r2

# RT_TangentialQuadrilateralSquaredArea returns squared area given edge and diagonal lengths
# when the quadrilateral is effectily "tangential" (eg  from Pitot theorem : a + c = b + d)
# area is maximal when ac = bd ie when quadrilateral is cyclic too
def RT_TangentialQuadrilateralSquaredArea(a,b,c,d,p,q):
    K2 = 1/4*(p^2*q^2 - (a*c - b*d)^2)
    return K2

def RT_TangentialQuadrilateralArea(a,b,c,d,p,q):
    K2 = RT_TangentialQuadrilateralSquaredArea(a,b,c,d,p,q)
    K = sqrt(K2)
    return K

def RT_TangentialQuadrilateralSemiPerimeter(a,b,c,d):
    s = 1/2*(a + b + c + d)
    return s

def RT_TangentialQuadrilateralInradius(a,b,c,d,p,q):
    K = RT_TangentialQuadrilateralArea(a,b,c,d,p,q)
    s = RT_TangentialQuadrilateralSemiPerimeter(a,b,c,d)
    r = K / s
    return r

def RT_IsTangentialQuadrilateralInradius(a,b,c,d):
    return (a + c) == (b + d)

# RT_IsAcute
# Reference : N.Wildberger  applications of rational trigonometry
# Given three points P1,P2,P3 and quadrances q12 = Q(P1,P2), q23 = Q(P2,P3), q13 = Q(P1,P3)
# then angle P1-P2-P3 is "strictly" acute in P2 iff q12 + q23 > q13  of (q12 + q23) / q13 > 1.0
# "strictly acute" meaning angle not 90 degrees or P1=P2, forbidding  q12 + q23 = q13
def RT_IsAcute(q12,q23,q13):
    d = ((q12.n(digits=nd) + q23.n(digits=nd)) / q13.n(digits=nd) - 1.0)
    ac = ( d > 0 )
    return ac

# 2nd parameter p2 is the point in angle
def IsAcute(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3); q23 = RT_Quadrance(p2,p3)
    ac = RT_IsAcute(q12,q23,q13)
    return ac

def RT_IsRectangular(q12,q23,q13):
    d = ((q12.n(digits=nd) + q23.n(digits=nd)) / q13.n(digits=nd) - 1.0)
    rect = ( d == 0 )
    return rect

def IsRectangular(p1,p2,p3):
    q12 = RT_Quadrance(p1,p2); q13 = RT_Quadrance(p1,p3); q23 = RT_Quadrance(p2,p3)
    rect = RT_IsRectangular(q12,q23,q13)
    return rect

# From 402976.main.pdf Dragutin Svrtan June 5 2009 On Circumradius Equations of Cyclic Polygons
# we set some functions for squared curvature (ro = 1/r^2, where r is circumradius) of cyclic polygons
def RT_CyclicTriangleSquaredCurvarture(q12,q23,q31):
    # Set the symmetric functions s1 = x + y + z ; s2 = xy + xz + yZ ; s3 = xyz
    e = SymmetricFunctions(QQ).e()
    x,y,z = var('x,y,z')
    s1(x,y,z) = e([1]).expand(3,alphabet='x,y,z')
    s2(x,y,z) = e([2]).expand(3,alphabet='x,y,z')
    s3(x,y,z) = e([3]).expand(3,alphabet='x,y,z')
    # Computing values for quadrances (square of edge lengths)
    k1 = s1(q12,q23,q31)
    k2 = s2(q12,q23,q31)
    k3 = s3(q12,q23,q31)
    # Heron's forumula for area and formula for diameter is continuation for result qa = q12, qb = q23 , qc = q31
    # ro = (4*(qa*qb + qa*qc + qb*qc) - (qa + qb + qc)^2)/(qa*qb*qc) = (-qa^2 -qb^2 -qc^2 +2*qa*qb + 2*qa*qc + 2*qb*qc)/(qa*qb*qc)
    # ro = ((qa + qb + qc)^2 - 2*(qa^2 + qb^2 + qc^2))/(qa*qb*qc) = 16*area^2/(qa*qb*qc) = 4* (4*area^2)/(qa*qb*qc)
    # but ro = 1/r^2 => r = 1/2( a*b*c/(2*area)) => diameter d =  a*b*c/(2*area)
    ro = (4*k2 - k1^2)/k3
    return ro

def RT_CyclicTriangleSquaredCicumradius(q12,q23,q31):
    ro = RT_CyclicTriangleSquaredCurvarture(q12,q23,q31)
    qr = 1/ro
    return qr

def RT_CyclicQuadSquaredCurvarture(q12,q23,q34,q41):
    # Set the symmetric functions
    e = SymmetricFunctions(QQ).e()
    x,y,z,t = var('x,y,z,t')
    s1(x,y,z,t) = e([1]).expand(4,alphabet='x,y,z,t')
    s2(x,y,z,t) = e([2]).expand(4,alphabet='x,y,z,t')
    s3(x,y,z,t) = e([3]).expand(4,alphabet='x,y,z,t')
    s4(x,y,z,t) = e([4]).expand(4,alphabet='x,y,z,t')
    # Computing values for quadrances (square of edge lengths)
    k1 = s1(q12,q23,q34,q41)
    k2 = s2(q12,q23,q34,q41)
    k3 = s3(q12,q23,q34,q41)
    k4 = s4(q12,q23,q34,q41)
    # ro is root of k1^2*k4*ro^2 - k1^4 - 2*k1^2*k3*ro - k3^2*ro^2 + 8*k1^2*k2 + 8*k2*k3*ro - 16*k1*k4*ro - 16*k2^2 + 64*k4 = 0
    # two roots : one is for convex quad and the other for non convex quad
    z1 = (k1^2 - 4*k2)*k3 + 8*k1*k4
    z2 = k1^3 - 4*k1*k2 + 8*k3
    z3 = sqrt(k4)  # k4 = qa*qb*qc*qd => z3 = a*b*c*d
    z4 = k1^2*k4 - k3^2
    # Returns two roots
    ro1 = (z1 - z2*z3)/z4
    ro2 = (z1 + z2*z3)/z4
    return [ro1,ro2]

def RT_CyclicQuadSquaredCicumradius(q12,q23,q34,q41):
    [ro1,ro2] = RT_CyclicQuadSquaredCurvarture(q12,q23,q34,q41)
    qr1 = 1/ro1 ; qr2 = 1/ro2
    return [qr1,qr2]


def RT_RationalGcd8(a):
    [a0,a1,a2,a3,a4,a5,a6,a7] = a
    a0 = QQ(a0); a1 = QQ(a1); a2 = QQ(a2); a3 = QQ(a3); a4 = QQ(a4); a5 = QQ(a5); a6 = QQ(a6); a7 = QQ(a7)
    a01 = gcd(a0,a1); a23 = gcd(a2,a3); a45 = gcd(a4,a5); a67 = gcd(a6,a7); a0123 = gcd(a01,a23); a4567 = gcd(a45,a67); a01234567 = gcd(a0123,a4567)
    return [a0/a01234567, a1/a01234567, a2/a01234567, a3/a01234567, a4/a01234567, a5/a01234567, a6/a01234567, a7/a01234567]

# RT_CyclicPentagonSquaredCurvartureEq returns squared curvature
#
# when ro = 0 is root(s) of eq5, then shift coeff result to remove these roots
# for example :
# coeff =  [0, -30984162508800, 550819182394736640, -3122062715938751456256, 4758112267355215557128109, 4926687172919212979146474144, 584189390020640128881889365248, 4349768537262367169044380188672]
# shift to
# coeff =  [-30984162508800, 550819182394736640, -3122062715938751456256, 4758112267355215557128109, 4926687172919212979146474144, 584189390020640128881889365248, 4349768537262367169044380188672, 0]
#
# try to have rational fixed precision for input quadrances, for example : 1/4 = 0.25 is ok but 1/7 = 0.14... is not

def RT_CyclicPentagonSquaredCurvartureEq(q12,q23,q34,q45,q51):
    # Set the symmetric functions
    e = SymmetricFunctions(QQ).e()
    x,y,z,t,w = var('x,y,z,t,w')
    s1(x,y,z,t,w) = e([1]).expand(5,alphabet='x,y,z,t,w')
    s2(x,y,z,t,w) = e([2]).expand(5,alphabet='x,y,z,t,w')
    s3(x,y,z,t,w) = e([3]).expand(5,alphabet='x,y,z,t,w')
    s4(x,y,z,t,w) = e([4]).expand(5,alphabet='x,y,z,t,w')
    s5(x,y,z,t,w) = e([5]).expand(5,alphabet='x,y,z,t,w')
    # Computing values for quadrances (square of edge lengths)
    k1 = s1(q12,q23,q34,q45,q51)
    k2 = s2(q12,q23,q34,q45,q51)
    k3 = s3(q12,q23,q34,q45,q51)
    k4 = s4(q12,q23,q34,q45,q51)
    k5 = s5(q12,q23,q34,q45,q51)
    # Use compact form
    ro = var('ro')
    l0(ro) =  -k5*ro^5 + 2*k4*ro^4 - 6*k3*ro^3 + 20*k2*ro^2 - 70*k1*ro + 252
    l1(ro) =  k4*ro^4 - 4*k3*ro^3 + 15*k2*ro^2 - 56*k1*ro + 210
    l2(ro) =  -k3*ro^3 + 6*k2*ro^2 - 28*k1*ro + 120
    l3(ro) =  k2*ro^2 - 8*k1*ro + 45
    l4(ro) =  -k1*ro + 10
    # compact Pell form equation for ro : alpha5^2 - beta5^2*delta5 = 0
    alpha5(ro) = l4^4 + (-3*l3 + 2*l2 + l1 -3)*l4^2 + (-2*l3 - 4*l1 + 2)*l4 + 2*l3^2 + (-2*l2 -2*l1 + 4)*l3 + l2^2 + 2*l2 -2*l1 + (l3 + 3)*l0 + 2
    beta5(ro) = -l4^3 + 2*l4^2 + (2*l3 - l2)*l4 - 2*l3 + 2*l1 - l0 -2
    delta5(ro) = l0 + 2*(l1 + l2 + l3 + l4 + 1)
    eq5(ro) = ((alpha5^2 - beta5^2*delta5)/ro^8).factor()
    # compute coefficients of equation degree 7 in ro and return them  (..could call RT_RationalGcd8 to convert them)
    coef = eq5.coefficients(sparse=False)
    return coef

def RT_s5_1(x, y, z, t, w):
    return t + w + x + y + z

def RT_s5_2(x, y, z, t, w):
    return t*w + t*x + w*x + t*y + w*y + x*y + t*z + w*z + x*z + y*z

def RT_s5_3(x, y, z, t, w):
    return t*w*x + t*w*y + t*x*y + w*x*y + t*w*z + t*x*z + w*x*z + t*y*z + w*y*z + x*y*z

def RT_s5_4(x, y, z, t, w):
    return t*w*x*y + t*w*x*z + t*w*y*z + t*x*y*z + w*x*y*z

def RT_s5_5(x, y, z, t, w):
    return t*w*x*y*z


# CyclicPentagonSquaredCurvartureEq has to been used for numerical computation
# comparing with RT_CyclicPentagonSquaredCurvartureEq : we don't divide by ro^8
# becquse for find_root it is then difficult to find root
def CyclicPentagonSquaredCurvartureEq(q12,q23,q34,q45,q51):
    # Computing values for quadrances (square of edge lengths)
    k1 =  RT_s5_1(q12,q23,q34,q45,q51).n()
    k2 =  RT_s5_2(q12,q23,q34,q45,q51).n()
    k3 =  RT_s5_3(q12,q23,q34,q45,q51).n()
    k4 =  RT_s5_4(q12,q23,q34,q45,q51).n()
    k5 =  RT_s5_5(q12,q23,q34,q45,q51).n()
    # Use compact form
    ro = var('ro')
    l0(ro) =  -k5*ro^5 + 2*k4*ro^4 - 6*k3*ro^3 + 20*k2*ro^2 - 70*k1*ro + 252
    l1(ro) =  k4*ro^4 - 4*k3*ro^3 + 15*k2*ro^2 - 56*k1*ro + 210
    l2(ro) =  -k3*ro^3 + 6*k2*ro^2 - 28*k1*ro + 120
    l3(ro) =  k2*ro^2 - 8*k1*ro + 45
    l4(ro) =  -k1*ro + 10
    # compact Pell form equation for ro : alpha5^2 - beta5^2*delta5 = 0
    alpha5(ro) = l4^4 + (-3*l3 + 2*l2 + l1 -3)*l4^2 + (-2*l3 - 4*l1 + 2)*l4 + 2*l3^2 + (-2*l2 -2*l1 + 4)*l3 + l2^2 + 2*l2 -2*l1 + (l3 + 3)*l0 + 2
    beta5(ro) = -l4^3 + 2*l4^2 + (2*l3 - l2)*l4 - 2*l3 + 2*l1 - l0 -2
    delta5(ro) = l0 + 2*(l1 + l2 + l3 + l4 + 1)
    eq5(ro) = alpha5^2 - beta5^2*delta5
    # finding root of equation, limiting to two times the circumradius of cyclid quad
    d12 = sqrt(q12) ; d23 =  sqrt(q23); d34 = sqrt(q34) ; d45 = sqrt(q45) ; d51 =  sqrt(q51)
    qO = RT_CyclicQuadSquaredCircumradius(d12,d23,d34 + d45,d51)
    ro1 = find_root(eq5,0.0,1/qO)
    return ro1


# Discriminant for a cubic a3 z^3 + a2z z^2 + a1 z + a0
def CubicDisc(a0,a1,a2,a3):
    return  a1^2*a2^2 - 4*a0*a2^3 - 4*a1^3*a3 + 18*a0*a1*a2*a3- 27*a0^2*a3^2

# Equation for u = 16 S^2 with S^2 squared area of the cyclic pentagon given by quadrance of edges
def RT_CyclicPentagonSquaredAreaEq(q12,q23,q34,q45,q51):
    # Set the symmetric functions
    e = SymmetricFunctions(QQ).e()
    x,y,z,t,w = var('x,y,z,t,w')
    s1(x,y,z,t,w) = e([1]).expand(5,alphabet='x,y,z,t,w')
    s2(x,y,z,t,w) = e([2]).expand(5,alphabet='x,y,z,t,w')
    s3(x,y,z,t,w) = e([3]).expand(5,alphabet='x,y,z,t,w')
    s4(x,y,z,t,w) = e([4]).expand(5,alphabet='x,y,z,t,w')
    s5(x,y,z,t,w) = e([5]).expand(5,alphabet='x,y,z,t,w')
    # Computing values for quadrances (square of edge lengths)
    k1 = s1(q12,q23,q34,q45,q51)
    k2 = s2(q12,q23,q34,q45,q51)
    k3 = s3(q12,q23,q34,q45,q51)
    k4 = s4(q12,q23,q34,q45,q51)
    k5 = s5(q12,q23,q34,q45,q51)
    # u = 16 S^2 with S^2 squared area of the cyclic pentagon
    u = var('u')
    t2 = u - 4*k2 + k1^2
    t3 = 8*k3 + k1*t2
    t4 = -64*k4 + t2^2
    t5 = 128*k5
    # u is root of (z^3 + 2t3 z^2 - ut4 z + 2u^2t5)/u^2
    eq(u) = (CubicDisc(2*u^2*t5,-u*t4,2*t3,1)/u^2).factor()
    # compute coefficients of equation degree 7 in u and return them
    coef = eq.coefficients(sparse=False)
    return coef

# Areas of Alhazen lunes (from Wolfram's formulas)
# a,b,c are edge lengths of a rectangular triangle with a^2 = b^2 + c^2
# last parameter t = arctan(b/c)
def RT_AlhazenLunesAreas(b,c,t):
    zr = 2*b*c ; zi = b^2*pi - 2*(b^2 + c^2)*t
    return [1/8*(zr + zi),1/8*(zr - zi)]

# RT_LastQuadVertex returns cartesian coordinates for image of C point in inverse of affine map
#    d_z = (x3*(y1 - y2) - x2*(y1 - y3) + x1*(y2 - y3))        : area of triangle p1p2p3
#    d_x = (x4*(y1 - y2) - x2*(y1 - y4) + x1*(y2 - y4))/d_z    : ratio of area triangle p1p2p4 / area of triangle p1p2p3
#    d_y = (x4*(y2 - y3) - x3*(y2 - y4) + x2*(y3 - y4))/d_z    : ratio of area triangle p2p3p4 / area of triangle p1p2p3
def RT_LastQuadVertex(p1,p2,p3,p4):
    d_z = TriangleShoeLaceArea(p1,p2,p3)
    d_x = TriangleShoeLaceArea(p1,p2,p4) / d_z
    d_y = TriangleShoeLaceArea(p2,p3,p4) / d_z
    return [d_x,d_y]

def RT_Map2QuadMatrixElements(d_x,d_y,x1,y1,x2,y2,x4,y4):
    k = -(d_y*x1 - (d_y - 1)*x2 - x4)/d_x ; l = x1 - x2
    m = -(d_y*y1 - (d_y - 1)*y2 - y4)/d_x ; n = y1 - y2
    o = x2 ; p = y2
    return [k,l,m,n,o,p]

def RT_Map2QuadMatrix(p1,p2,p3,p4):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3 ; [x4,y4] = p4
    [d_x,d_y] = RT_LastQuadVertex(p1,p2,p3,p4)
    return RT_Map2QuadMatrixElements(d_x,d_y,x1,y1,x2,y2,x4,y4)

# RT_Map2Quad returns point from the affine map of (0,1) (0,0) (1,0) (d_x,d_y) to p1p2p3p4 = ABDC quadrilateral
def RT_Map2Quad(z,p1):
    [k,l,m,n,o,p] = z ; [x,y] = p1
    return [ k*x + l*y + o , m*x + n*y + p]


print "...RationalTrigonometry module loaded"