print "QuadriplanarCoordinates (QC) module loading ..."

# Quadriplanar coordinates are coordinates with respect to a reference regular tetrahedron ABCD

# Module includes 3D functions

# ---- Quadriplanar coordinates


# ---- 3D functions

# Difference of two 3D vectors
def SpaceDifference(u,v):
    [u1,u2,u3] = u ; [v1,v2,v3] = v
    w1 = u1 - v1
    w2 = u2 - v2
    w3 = u3 - v3
    w = [w1,w2,w3]
    return w

# Cross Product of two 3D vectors
def SpaceCrossProduct(u,v):
    [u1,u2,u3] = u ; [v1,v2,v3] = v
    w1 = u2*v3 - u3*v2
    w2 = u3*v1 - u1*v3
    w3 = u1*v2 - u2*v1
    w = [w1,w2,w3]
    return w

# Dot Product of two 3D vectors
def SpaceDotProduct(u,v):
    [u1,u2,u3] = u ; [v1,v2,v3] = v
    sp = u1*v1 + u2*v2 + u3*v3
    return sp

# Norm of 3D vector
def SpaceNorm(u):
    n = sqrt(SpaceDotProduct(u,u))
    return n

# Unit normal of two 3D vectors (vector orthogonal to plane of vector and with norm 1)
def SpaceUnitNormal(u,v):
    w = SpaceCrossProduct(u,v) ; [w1,w2,w3] = w
    n = SpaceNorm(w)
    wn = [w1/n,w2/n,w3/n]
    return wn

# Distance from point x to plane defined by points u1,u2,u3
# distance is signed
def SpacePlaneDistance(x,u1,u2,u3):
    u21 = SpaceDifference(u2,u1)
    u31 = SpaceDifference(u3,u1)
    n123 = SpaceUnitNormal(u21,u31)
    x1 = SpaceDifference(x,u1)
    d =  SpaceDotProduct(x1,n123)
    return d

def SpacePlaneQuadrance(x,u1,u2,u3):
    d = SpacePlaneDistance(x,u1,u2,u3)
    q = d^2
    return q

#  Squared distance between two 3D points
def SpaceQuadrance(u,v):
    w = SpaceDifference(u,v)
    q = SpaceDotProduct(w,w)
    return q

#  Distance between two 3D points
def SpaceDistance(u,v):
    q = SpaceQuadrance(u,v)
    d = sqrt(q)
    return d

# Squared Volume of tetahedron
# Using Tartaglia's formula and the Cayley-Menger determinant
def SpaceSquaredVolume(u1,u2,u3,u4):
    q12 = SpaceQuadrance(u1,u2) ; q13 = SpaceQuadrance(u1,u3) ; q14 = SpaceQuadrance(u1,u4)
    q23 = SpaceQuadrance(u2,u3) ; q24 = SpaceQuadrance(u2,u4) ; q34 =  SpaceQuadrance(u3,u4)
    m = Matrix(5,5,[[0,q12,q13,q14,1],[q12,0,q23,q24,1],[q13,q23,0,q34,1],[q14,q24,q34,0,1],[1,1,1,1,0]])
    v2 = m.determinant()/288
    return v2

def SpaceVolume(u1,u2,u3,u4):
    v2 = SpaceSquaredVolume(u1,u2,u3,u4)
    v = sqrt(v2)
    return v

print "...Quadriplanar coordinates module loaded"









