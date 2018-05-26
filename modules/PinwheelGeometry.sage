print "Pinwheel Geometry(PWG) module loading ..."

# PWG_Tiling returns tiling such as BJ line is split into BV-VT-TW-WJ segments with lengths of same length
# point S splits segment TW in two
def PWG_Tiling(pA,pB,pC,Red,Green,Blue):
    pD = Middle(pA,pB)
    pE = Middle(pB,pC)
    pF = Middle(pA,pC)
    pG = Middle(pB,pE)
    pH = Middle(pC,pE)
    pI = Middle(pC,pF)
    pJ = Middle(pA,pF)
    pK = Middle(pE,pF)
    pM = Intersect(pD,pG,pB,pJ)
    pU = Intersect(pC,pM,pA,pE)
    pL = Middle(pB,pU)
    pN = Middle(pD,pM)
    pP = Intersect(pA,pE,pB,pJ)
    pO = Middle(pP,pU)
    pQ = Middle(pF,pO)
    pR = Middle(pK,pO)
    pT = Middle(pB,pJ)
    pW = Middle(pT,pJ)
    pS = Middle(pT,pW)
    pV = Middle(pB,pT)
    #
    vertices = []
    vertices.append(pA)
    vertices.append(pB)
    vertices.append(pC)
    vertices.append(pD)
    vertices.append(pE)
    vertices.append(pF)
    vertices.append(pG)
    vertices.append(pH)
    vertices.append(pI)
    vertices.append(pJ)
    vertices.append(pK)
    vertices.append(pL)
    vertices.append(pM)
    vertices.append(pN)
    vertices.append(pO)
    vertices.append(pP)
    vertices.append(pQ)
    vertices.append(pR)
    vertices.append(pS)
    vertices.append(pT)
    vertices.append(pU)
    vertices.append(pV)
    vertices.append(pW)
    edges1 = []
    edges1.append([pA,pB])
    edges1.append([pA,pJ])
    edges1.append([pD,pN])
    edges1.append([pN,pV])
    edges1.append([pA,pW])
    edges1.append([pB,pN])
    edges1.append([pN,pT])
    edges1.append([pD,pW])
    edges1.append([pD,pS])
    edges1.append([pB,pJ])
    edges2 = []
    edges2.append([pB,pJ])
    edges2.append([pB,pC])
    edges2.append([pC,pJ])
    edges2.append([pE,pF])
    edges2.append([pC,pK])
    edges2.append([pH,pK])
    edges2.append([pI,pK])
    edges2.append([pE,pO])
    edges2.append([pO,pW])
    edges2.append([pG,pL])
    edges2.append([pL,pV])
    edges2.append([pF,pO])
    edges2.append([pJ,pO])
    edges2.append([pK,pO])
    edges2.append([pJ,pQ])
    edges2.append([pQ,pR])
    edges2.append([pE,pR])
    edges2.append([pG,pO])
    edges2.append([pL,pO])
    edges2.append([pO,pV])
    edges2.append([pO,pS])
    edges2.append([pO,pT])
    edges2.append([pB,pL])
    edges2.append([pM,pP])
    edges2.append([pB,pM])
    colors1 = []
    pCEK = TriangleCentroid(pC,pE,pK); colors1.append([pCEK,Blue,"st"]) # st Blue
    pFOK = TriangleCentroid(pF,pO,pK); colors1.append([pFOK,Blue,"xy"]) # xy Blue
    pFJQ = TriangleCentroid(pF,pJ,pQ); colors1.append([pFJQ,Blue,"c"]) # c  Blue
    pCFK = TriangleCentroid(pC,pF,pK); colors1.append([pCFK,Red,"ab"]) # ab Red
    pJQOW = QuadrilateralCentroid(pJ,pQ,pO,pW); colors1.append([pJQOW,Green,"de"]) # de Green
    pOSW = TriangleCentroid(pO,pS,pW); colors1.append([pOSW,Red,"u"]) # u Red
    pOTS = TriangleCentroid(pO,pT,pS); colors1.append([pOTS,Blue,"v"]) # v  Blue
    pBGLV = QuadrilateralCentroid(pB,pG,pL,pV); colors1.append([pBGLV,Blue,"no"]) # no Blue
    pLVTO = QuadrilateralCentroid(pL,pV,pT,pO); colors1.append([pLVTO,Green,"wm"]) # wm Green
    pGLO = TriangleCentroid(pG,pL,pO); colors1.append([pGLO,Red,"p"]) # p Red
    pEGO = TriangleCentroid(pE,pG,pO); colors1.append([pEGO,Blue,"q"]) # q Blue
    pEKO = TriangleCentroid(pE,pK,pO); colors1.append([pEKO,Green,"rz"]) # rz Green
    colors2 = []
    pAJW = TriangleCentroid(pA,pJ,pW); colors2.append([pAJW,Blue,"f"]) # f Blue
    pADSW = QuadrilateralCentroid(pA,pD,pS,pW) ; colors2.append([pADSW,Red,"gh"]) # gh Red
    pDNTS = QuadrilateralCentroid(pD,pN,pT,pS); colors2.append([pDNTS,Blue,"i"]) # i  Blue
    pBDNTV = PentagonOCentroid(pB,pD,pN,pT,pV); colors2.append([pBDNTV,Red,"jkl"]) # jkl Red
    return [vertices,edges1,edges2,colors1,colors2]

# PWG_Convex_Tiling returns convex tiles and with only Middle operations justifying coloring vertices in three colors
# simply by coloring A,B,C and then defining color of Middle(p1,p2) as last color different from those of p1 and p2
def PWG_Convex_Tiling(pA,pB,pC,Red,Green,Blue):
    # A Red   B Green   C Blue
    pD = Middle(pA,pB)                # Blue
    pE = Middle(pB,pC)                # Red
    pF = Middle(pA,pC)                # Green
    pG = Middle(pB,pE)                # Blue
    pH = Middle(pC,pE)                # Green
    pI = Middle(pC,pF)                # Red
    pJ = Middle(pA,pF)                # Blue
    pK = Middle(pE,pF)                # Blue
    pO = Intersect(pD,pI,pA,pE)       #  White because adjacent to Red,Blue,Green vertices
    pS = Middle(pD,pO)                # Red
    pQ = Middle(pF,pO)                # Red
    pR = Middle(pK,pO)                # Green
    pP = Middle(pJ,pS)                # Green
    pM = Middle(pB,pP)                # Blue
    pT = Middle(pM,pS)                # Green
    pN = Middle(pD,pM)                # Red
    pU = Middle(pE,pO)                # Blue
    pL = Middle(pB,pU)                # Red
    pV = pM
    pW = pP
    #
    vertices = []
    vertices.append(pA)
    vertices.append(pB)
    vertices.append(pC)
    vertices.append(pD)
    vertices.append(pE)
    vertices.append(pF)
    vertices.append(pG)
    vertices.append(pH)
    vertices.append(pI)
    vertices.append(pJ)
    vertices.append(pK)
    vertices.append(pL)
    vertices.append(pM)
    vertices.append(pN)
    vertices.append(pO)
    vertices.append(pP)
    vertices.append(pQ)
    vertices.append(pR)
    vertices.append(pS)
    vertices.append(pT)
    vertices.append(pU)
    vertices.append(pV)
    vertices.append(pW)
    edges1 = []
    edges1.append([pA,pB])
    edges1.append([pA,pJ])
    edges1.append([pD,pN])
    edges1.append([pN,pV])
    edges1.append([pA,pW])
    edges1.append([pB,pN])
    edges1.append([pN,pT])
    edges1.append([pD,pW])
    edges1.append([pD,pS])
    edges1.append([pB,pJ])
    edges2 = []
    edges2.append([pB,pJ])
    edges2.append([pB,pC])
    edges2.append([pC,pJ])
    edges2.append([pE,pF])
    edges2.append([pC,pK])
    edges2.append([pH,pK])
    edges2.append([pI,pK])
    edges2.append([pE,pO])
    edges2.append([pO,pW])
    edges2.append([pG,pL])
    edges2.append([pL,pV])
    edges2.append([pF,pO])
    edges2.append([pJ,pO])
    edges2.append([pK,pO])
    edges2.append([pJ,pQ])
    edges2.append([pQ,pR])
    edges2.append([pE,pR])
    edges2.append([pG,pO])
    edges2.append([pL,pO])
    edges2.append([pO,pV])
    edges2.append([pO,pS])
    edges2.append([pO,pT])
    edges2.append([pB,pL])
    edges2.append([pM,pP])
    edges2.append([pB,pM])
    colors1 = []
    pCEK = TriangleCentroid(pC,pE,pK); colors1.append([pCEK,Blue,"st"]) # st Blue
    pFOK = TriangleCentroid(pF,pO,pK); colors1.append([pFOK,Blue,"xy"]) # xy Blue
    pFJQ = TriangleCentroid(pF,pJ,pQ); colors1.append([pFJQ,Blue,"c"]) # c  Blue
    pCFK = TriangleCentroid(pC,pF,pK); colors1.append([pCFK,Red,"ab"]) # ab Red
    pJQOW = QuadrilateralCentroid(pJ,pQ,pO,pW); colors1.append([pJQOW,Green,"de"]) # de Green
    pOSW = TriangleCentroid(pO,pS,pW); colors1.append([pOSW,Red,"u"]) # u Red
    pOTS = TriangleCentroid(pO,pT,pS); colors1.append([pOTS,Blue,"v"]) # v  Blue
    pBGLV = QuadrilateralCentroid(pB,pG,pL,pV); colors1.append([pBGLV,Blue,"no"]) # no Blue
    pLVTO = QuadrilateralCentroid(pL,pV,pT,pO); colors1.append([pLVTO,Green,"wm"]) # wm Green
    pGLO = TriangleCentroid(pG,pL,pO); colors1.append([pGLO,Red,"p"]) # p Red
    pEGO = TriangleCentroid(pE,pG,pO); colors1.append([pEGO,Blue,"q"]) # q Blue
    pEKO = TriangleCentroid(pE,pK,pO); colors1.append([pEKO,Green,"rz"]) # rz Green
    colors2 = []
    pAJW = TriangleCentroid(pA,pJ,pW); colors2.append([pAJW,Blue,"f"]) # f Blue
    pADSW = QuadrilateralCentroid(pA,pD,pS,pW) ; colors2.append([pADSW,Red,"gh"]) # gh Red
    pDNTS = QuadrilateralCentroid(pD,pN,pT,pS); colors2.append([pDNTS,Blue,"i"]) # i  Blue
    pBDNTV = PentagonCentroid(pB,pD,pN,pT,pV); colors2.append([pBDNTV,Red,"jkl"]) # jkl Red
    return [vertices,edges1,edges2,colors1,colors2]

# Pinwheel tiling in order to have a clear cut for Haberdasher
def PWG_Haberdasher_Tiling(pA,pB,pC,Red,Green,Blue):
    # A Red   B Green   C Blue
    pD = Middle(pA,pB)                # Blue
    pF = Middle(pA,pC)                # Green
    [pG,pZ] = Arcs(pB,pC,pF,TriangleArea(pA,pB,pC))  # Blue
    [pZ,pH] = Arcs(pB,pC,pF,RT_Quadrance(pD,pG)) # Green
    pI = Middle(pC,pF)                # Red
    pJ = Middle(pA,pF)                # Blue
    pM = Intersect(pB,pJ,pD,pG)       # Blue
    pK = Middle(pF,pH)                # Blue
    pO = Middle(pF,pG)                # White because adjacent to Red,Blue,Green vertices
    pE = Intersect(pJ,pO,pB,pC)       # Red
    pP = Intersect(pB,pJ,pA,pE)       # Green
    pS = Middle(pM,pP)                # Red
    pQ = Middle(pF,pO)                # Red
    pR = Middle(pK,pO)                # Green
    pT = Middle(pM,pS)                # Green
    pN = Middle(pD,pM)                # Red
    pU = Middle(pE,pO)                # Blue
    pL = Middle(pM,pG)                # Red
    pV = pM
    pW = pP
    #
    vertices = []
    vertices.append(pA)
    vertices.append(pB)
    vertices.append(pC)
    vertices.append(pD)
    vertices.append(pE)
    vertices.append(pF)
    vertices.append(pG)
    vertices.append(pH)
    vertices.append(pI)
    vertices.append(pJ)
    vertices.append(pK)
    vertices.append(pL)
    vertices.append(pM)
    vertices.append(pN)
    vertices.append(pO)
    vertices.append(pP)
    vertices.append(pQ)
    vertices.append(pR)
    vertices.append(pS)
    vertices.append(pT)
    vertices.append(pU)
    vertices.append(pV)
    vertices.append(pW)
    edges1 = []
    edges1.append([pA,pB])
    edges1.append([pA,pJ])
    edges1.append([pD,pN])
    edges1.append([pN,pV])
    edges1.append([pA,pW])
    edges1.append([pB,pN])
    edges1.append([pN,pT])
    edges1.append([pD,pW])
    edges1.append([pD,pS])
    edges1.append([pB,pJ])
    edges2 = []
    edges2.append([pB,pJ])
    edges2.append([pB,pC])
    edges2.append([pC,pJ])
    edges2.append([pE,pK])
    edges2.append([pF,pK])
    edges2.append([pC,pK])
    edges2.append([pH,pK])
    edges2.append([pI,pK])
    edges2.append([pE,pO])
    edges2.append([pO,pW])
    edges2.append([pG,pL])
    edges2.append([pL,pV])
    edges2.append([pF,pO])
    edges2.append([pJ,pO])
    edges2.append([pK,pO])
    edges2.append([pJ,pQ])
    edges2.append([pQ,pR])
    edges2.append([pE,pR])
    edges2.append([pG,pO])
    edges2.append([pL,pO])
    edges2.append([pO,pV])
    edges2.append([pO,pS])
    edges2.append([pO,pT])
    edges2.append([pB,pL])
    edges2.append([pM,pP])
    edges2.append([pB,pM])
    colors1 = []
    colors2 = []
    # pCEK = TriangleOCentroid(pC,pE,pK); colors1.append([pCEK,Blue,"st"]) # st Blue
    pEHK = TriangleCentroid(pE,pH,pK); colors1.append([pEHK,Red,"s"]) # s Red
    pCHK = TriangleCentroid(pC,pH,pK); colors1.append([pCHK,Blue,"t"]) # t Blue
    #pFOK = TriangleCentroid(pF,pO,pK); colors1.append([pFOK,Blue,"xy"]) # xy Blue
    pFOK = TriangleCentroid(pF,pO,pK); colors1.append([pFOK,Red,"xy"]) # xy Red
    #pFJQ = TriangleCentroid(pF,pJ,pQ); colors1.append([pFJQ,Blue,"c"]) # c  Blue
    pFJQ = TriangleCentroid(pF,pJ,pQ); colors1.append([pFJQ,Red,"c"]) # c  Red
    #pCFK = TriangleCentroid(pC,pF,pK); colors1.append([pCFK,Red,"ab"]) # ab Red
    pCIK = TriangleCentroid(pC,pI,pK); colors1.append([pCIK,Green,"a"]) # a Green
    pFIK = TriangleCentroid(pF,pI,pK); colors1.append([pFIK,Blue,"b"]) #  b Blue
    pJQOW = QuadrilateralCentroid(pJ,pQ,pO,pW); colors1.append([pJQOW,Green,"de"]) # de Green
    #pJOQ = TriangleOCentroid(pJ,pO,pQ); colors1.append([pJOQ,Blue,"d"]) # d Blue
    #pJOW = TriangleCentroid(pJ,pO,pW); colors1.append([pJOW,Green,"e"]) # e Green
    pOSW = TriangleCentroid(pO,pS,pW); colors1.append([pOSW,Red,"u"]) # u Red
    pOTS = TriangleCentroid(pO,pT,pS); colors1.append([pOTS,Blue,"v"]) # v  Blue
    #pBGLV = QuadrilateralCentroid(pB,pG,pL,pV); colors1.append([pBGLV,Blue,"no"]) # no Blue
    pBLV = TriangleCentroid(pB,pL,pV); colors1.append([pBLV,Blue,"n"]) # n Blue
    pBGL = TriangleCentroid(pB,pG,pL); colors1.append([pBGL,Green,"o"]) # o Green
    pLVTO = QuadrilateralCentroid(pL,pV,pT,pO); colors1.append([pLVTO,Green,"wm"]) # wm Green
    pGLO = TriangleCentroid(pG,pL,pO); colors1.append([pGLO,Red,"p"]) # p Red
    pEGO = TriangleCentroid(pE,pG,pO); colors1.append([pEGO,Blue,"q"]) # q Blue
    pEKO = TriangleCentroid(pE,pK,pO); colors1.append([pEKO,Green,"rz"]) # rz Green
    colors2 = []
    pAJW = TriangleCentroid(pA,pJ,pW); colors2.append([pAJW,Blue,"f"]) # f Blue
    pADSW = QuadrilateralCentroid(pA,pD,pS,pW) ; colors2.append([pADSW,Red,"gh"]) # gh Red
    pDNTS = QuadrilateralCentroid(pD,pN,pT,pS); colors2.append([pDNTS,Blue,"i"]) # i  Blue
    #pBDNTV = PentagonCentroid(pB,pD,pN,pT,pV); colors2.append([pBDNTV,Red,"jkl"]) # jkl Red
    pBDN = TriangleCentroid(pB,pD,pN); colors2.append([pBDN,Green,"j"]) # j Green
    pBNT = TriangleCentroid(pB,pN,pT); colors2.append([pBNT,Red,"kl"]) # kl Red
    return [vertices,edges1,edges2,colors1,colors2]

# PWG_Corners_Square_Tiling returns corners for Pinwheel tiling of square ABCD
def PWG_Corners_Square_Tiling(pA,pB,pC,pD):
    pE = Middle(pC,pD)
    pF = Middle(pA,pD)
    pG = Intersect(pA,pE,pB,pF)
    pH = Middle(pB,pG)
    pI = Middle(pB,pC)
    pJ = Intersect(pA,pE,pD,pI)
    pK = Symetric(pE,pJ)
    pL = Symetric(pK,pJ)
    return [pE,pF,pG,pH,pI,pJ,pK,pL]

# PWG_Square_Tiling returns Pinwheel tiling of square ABCD
# g1,g2,g3 are colors, for example g1 = c2 (red) ; g2 = c3 (green) ; g3 = c7 (blue)
def PWG_Square_Tiling(pA,pB,pC,pD,g1,g2,g3):
    # Set lists to empty
    vertices = edges = colors = []
    # Get corners of Pinwheel tiles
    [pE,pF,pG,pH,pI,pJ,pK,pL] = PWG_Corners_Square_Tiling(pA,pB,pC,pD)
    # Get tiling for triangle 1
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Convex_Tiling(pG,pA,pB,g1,g2,g3)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1 + Tile_edges2; colors = colors + Tile_colors1 + Tile_colors2
    # Get tiling for triangle 2
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Convex_Tiling(pH,pB,pC,g3,g2,g1)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1 + Tile_edges2; colors = colors + Tile_colors1 + Tile_colors2
    # Get tiling for triangle 3
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Convex_Tiling(pH,pG,pC,g1,g2,g3)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1 + Tile_edges2; colors = colors + Tile_colors1 + Tile_colors2
    # Get tiling for triangle 4
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Convex_Tiling(pJ,pD,pA,g3,g2,g1)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1 + Tile_edges2; colors = colors + Tile_colors1 + Tile_colors2
    # Get tiling for triangle 5 part 1
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Convex_Tiling(pK,pC,pG,g3,g2,g1)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges2; colors = colors + Tile_colors1
    # Get tiling for triangle 5 part 2
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Convex_Tiling(pJ,pD,pL,g1,g2,g3)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1; colors = colors + Tile_colors2
    #
    return [vertices,edges,colors]

# PWG_Rectangle_Tiling returns Pinwheel tiling of rectangle ABCD
def PWG_Rectangle_Tiling(pA,pB,pC,pD):
    # Set lists to empty
    vertices = edges = colors = []
    # Get middle points in order to split rectangle in two parts
    pE = Middle(pB,pC) ; pF = Middle(pA,pD)
    # Get tiling for square ABEF
    [Tile_vertices,Tile_edges,Tile_colors] = PWG_Square_Tiling(pA,pB,pE,pF,c2,c3,c7)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges; colors = colors + Tile_colors
    # Get tiling for square DCEF
    [Tile_vertices,Tile_edges,Tile_colors] = PWG_Square_Tiling(pD,pC,pE,pF,c7,c3,c2)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges; colors = colors + Tile_colors
    #
    return [vertices,edges,colors]

# PWG_Double_Rectangle_Tiling returns Pinwheel tiling of a square ABCD made of two rectangles
def PWG_Double_Rectangle_Tiling(pA,pB,pC,pD):
    # Set lists to empty
    vertices = edges = colors = []
    # Get middle points in order to split rectangle in two parts
    pE = Middle(pB,pC) ; pF = Middle(pA,pD)
    # Get tiling for rectangle ABEF
    [Tile_vertices,Tile_edges,Tile_colors] = PWG_Rectangle_Tiling(pA,pF,pE,pB)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges; colors = colors + Tile_colors
    # Get tiling for rectangle ECDF
    [Tile_vertices,Tile_edges,Tile_colors] = PWG_Rectangle_Tiling(pE,pC,pD,pF)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges; colors = colors + Tile_colors
    #
    return [vertices,edges,colors]

# PWG_Square_Haberdasher_Tiling returns Pinwheel tiling of square UBGF
def PWG_Square_Haberdasher_Tiling(pU,pB,pG,pF,Red,Green,Blue):
    # Compute quadrance edge of square
    q = RT_Quadrance(pU,pB)
    # Compute Haberdasher square
    pH = Middle(pF,pU)
    pD = Middle(pB,pG)
    pO = Middle(pF,pG)
    pQ = Middle(pF,pO)
    [pA,pZ] = CirclesIntersect(pD,q/4,pF,q)
    pJ = Middle(pA,pF)
    pA1 = Symetric(pD,pA)
    pC1 = Symetric(pF,pA)
    #[pG1,pZ1] = Arcs(pA1,pC1,pF,RT_Quadrance(pB,pG))
    #[pZ1,pH1] = Arcs(pA1,pC1,pF,RT_Quadrance(pD,pG))
    #pE1 = Middle(pG1,pH1)
    pE1 = Intersect(pJ,pO,pA1,pC1)
    pW = Intersect(pA1,pJ,pA,pE1)
    pV = Intersect(pJ,pW,pD,pG)
    pN = Middle(pD,pV)
    pL = Middle(pG,pV)
    pS = Middle(pV,pW)
    pT = Middle(pS,pV)
    pK = Middle(pF,pH)
    pX = Middle(pH,pU)
    pP = Middle(pB,pU)
    pR = Middle(pX,pP)
    pC = Middle(pU,pP)
    pE = Extend(pD,pB,RT_Quadrance(pD,pV))
    pI = Middle(pB,pE)
    pM = Middle(pD,pE)
    #
    vertices = []
    vertices.append(pA)
    vertices.append(pB)
    vertices.append(pC)
    vertices.append(pD)
    vertices.append(pE)
    vertices.append(pF)
    vertices.append(pG)
    vertices.append(pH)
    vertices.append(pI)
    vertices.append(pJ)
    vertices.append(pK)
    vertices.append(pL)
    vertices.append(pM)
    vertices.append(pN)
    vertices.append(pO)
    vertices.append(pP)
    vertices.append(pQ)
    vertices.append(pR)
    vertices.append(pS)
    vertices.append(pT)
    vertices.append(pU)
    vertices.append(pV)
    vertices.append(pW)
    vertices.append(pX)
    #
    edges = []
    edges.append([pB,pG])
    edges.append([pG,pF])
    edges.append([pF,pU])
    edges.append([pU,pB])
    edges.append([pA,pB])
    edges.append([pA,pD])
    edges.append([pA,pF])
    edges.append([pA,pH])
    edges.append([pA,pW])
    edges.append([pJ,pV])
    edges.append([pA,pI])
    edges.append([pA,pE])
    edges.append([pA,pM])
    edges.append([pA,pK])
    edges.append([pA,pX])
    edges.append([pA,pR])
    edges.append([pA,pP])
    edges.append([pD,pS])
    edges.append([pD,pW])
    edges.append([pO,pL])
    edges.append([pO,pS])
    edges.append([pO,pT])
    edges.append([pO,pV])
    edges.append([pO,pW])
    edges.append([pJ,pK])
    edges.append([pJ,pO])
    edges.append([pJ,pQ])
    edges.append([pR,pC])
    edges.append([pR,pP])
    edges.append([pR,pX])
    edges.append([pN,pT])
    #
    colors = []
    pAHX = TriangleCentroid(pA,pH,pX); colors.append([pAHX,Red,"s"])               # s   Red
    pAHK = TriangleCentroid(pA,pH,pK); colors.append([pAHK,Blue,"t"])              # t   Blue
    pUXP = TriangleCentroid(pU,pX,pP); colors.append([pUXP,Red,"xy"])              # xy  Red
    pFJQ = TriangleCentroid(pF,pJ,pQ); colors.append([pFJQ,Red,"c"])               # c   Red
    pAJK = TriangleCentroid(pA,pJ,pK); colors.append([pAJK,Green,"a"])             # a   Green
    pFJK = TriangleCentroid(pF,pJ,pK); colors.append([pFJK,Blue,"b"])              # b   Blue
    pJQOW = QuadrilateralCentroid(pJ,pQ,pO,pW); colors.append([pJQOW,Green,"de"])  # de  Green
    pOSW = TriangleCentroid(pO,pS,pW); colors.append([pOSW,Red,"u"])               # u   Red
    pOST = TriangleCentroid(pO,pS,pT); colors.append([pOST,Blue,"v"])              # v   Blue
    pAEI = TriangleCentroid(pA,pE,pI); colors.append([pAEI,Blue,"n"])              # n   Blue
    pABI = TriangleCentroid(pA,pB,pI); colors.append([pABI,Green,"o"])             # o   Green
    pLVTO = QuadrilateralCentroid(pL,pV,pT,pO); colors.append([pLVTO,Green,"wm"])     # wm  Green
    pGLO = TriangleCentroid(pG,pL,pO); colors.append([pGLO,Red,"p"])               # p   Red
    pABP = TriangleCentroid(pA,pB,pP); colors.append([pABP,Blue,"q"])              # q   Blue
    pAXP = TriangleCentroid(pA,pX,pP); colors.append([pAXP,Green,"rz"])            # rz  Green
    pAJW = TriangleCentroid(pA,pJ,pW); colors.append([pAJW,Blue,"f"])              # f   Blue
    pADSW = QuadrilateralCentroid(pA,pD,pS,pW); colors.append([pADSW,Red,"gh"])    # gh  Red
    pDNTS = QuadrilateralCentroid(pD,pN,pT,pS); colors.append([pDNTS,Blue,"i"])    # i   Blue
    pADM = TriangleCentroid(pA,pD,pM); colors.append([pADM,Green,"j"])             # j  Green
    pAEM = TriangleCentroid(pA,pE,pM); colors.append([pAEM,Red,"k"])               # k  Red
    pNTV = TriangleCentroid(pN,pT,pV); colors.append([pNTV,Red,"l"])               # l  Red
    #
    return [vertices,edges,colors]

# PWG_5PW_Square_Haberdasher_Tiling returns Pinwheel tiling of square ABCD as five Pinwheel
# g1,g2,g3 are colors, for example g1 = c2 (red) ; g2 = c3 (green) ; g3 = c7 (blue)
def PWG_5PW_Square_Haberdasher_Tiling(pA,pB,pC,pD,g1,g2,g3):
    # Set lists to empty
    vertices = edges = colors = []
    # Get corners of Pinwheel tiles
    [pE,pF,pG,pH,pI,pJ,pK,pL] = PWG_Corners_Square_Tiling(pA,pB,pC,pD)
    # Get tiling for triangle 1
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Haberdasher_Tiling(pG,pA,pB,g1,g2,g3)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1 + Tile_edges2; colors = colors + Tile_colors1 + Tile_colors2
    # Get tiling for triangle 2
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Haberdasher_Tiling(pH,pB,pC,g3,g2,g1)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1 + Tile_edges2; colors = colors + Tile_colors1 + Tile_colors2
    # Get tiling for triangle 3
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Haberdasher_Tiling(pH,pG,pC,g1,g2,g3)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1 + Tile_edges2; colors = colors + Tile_colors1 + Tile_colors2
    # Get tiling for triangle 4
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Haberdasher_Tiling(pJ,pD,pA,g3,g2,g1)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1 + Tile_edges2; colors = colors + Tile_colors1 + Tile_colors2
    # Get tiling for triangle 5 part 1
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Haberdasher_Tiling(pK,pC,pG,g3,g2,g1)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges2; colors = colors + Tile_colors1
    # Get tiling for triangle 5 part 2
    [Tile_vertices,Tile_edges1,Tile_edges2,Tile_colors1,Tile_colors2] = PWG_Haberdasher_Tiling(pJ,pD,pL,g1,g2,g3)
    vertices = vertices + Tile_vertices; edges = edges + Tile_edges1; colors = colors + Tile_colors2
    #
    return [vertices,edges,colors]

print "...Pinwheel Geometry module loaded"