module DeUni.Dim3.ReTri3D where

import Hammer.Math.Vector

import DeUni.Types



-- | Based on the papers: "Parallel dynamic and kinetic regular triangulation in three dimensions" (1993) and 
-- "A data-parallel algorithm for three-dimensional Delaunay triangulation and its implementation" (2005)

getCircumSphere::WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> (Double, Vec3)
getCircumSphere a b c d = (radius, center) 
  where
    (dist, center) = getFaceDistCenter a b c d
    radius         = (normsqr $ point a &- center) - weight a
  
getFaceDistCenter::WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> (Double, Vec3)  
getFaceDistCenter a b c d = (signDist, center)
  where
    center          = (-0.5) *& ((mux *& q1) &+ ( muy *& q2) &+ (muz *& q3))
    signDist        = if signRDet r then -dist else dist
    dist            = muz*0.5 + (q3 &. (point a))
    m               = getM a b c d
    (q,r)           = qrGramSchmidt m --qrDecomp m
    (mux, muy, muz) = solveMu ((-1) *& (getAlpha a b c d)) r
    (Mat3 q1 q2 q3) = q
    
getM :: WPoint Point3D -> WPoint  Point3D -> WPoint Point3D -> WPoint Point3D -> Mat3
getM a b c d = Mat3 (point b &- point a) (point c &- point a) (point d &- point a)

getAlpha :: WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> Vec3
getAlpha a b c d = Vec3 (fun b a) (fun c a) (fun d a)
  where fun x y = (normsqr.point) x - (normsqr.point) y - weight x + weight y

solveMu::Vec3 -> Mat3 -> (Double,Double,Double)
solveMu a r@(Mat3 r1 r2 r3) = (mux, muy, muz)
  where
    (ax, ay, az)    = unvec3 a
    (r11, r21, r31) = unvec3 r1
    (r12, r22, r32) = unvec3 r2
    (r13, r23, r33) = unvec3 r3
    mux = ax / r11
    muy = (ay - mux*r12) / r22
    muz = (az - mux*r13 - muy* r23) / r33        

signRDet::Mat3 -> Bool
signRDet (Mat3 r1 r2 r3) = r11 * r22 * r33 >= 0
  where
    (r11, r21, r31) = unvec3 r1
    (r12, r22, r32) = unvec3 r2
    (r13, r23, r33) = unvec3 r3

unvec3 (Vec3 a b c) = (a,b,c)


-- | A Householder reflection (or Householder transformation) is a transformation that
-- takes a vector and reflects it about some plane or hyperplane. We can use this operation
-- to calculate the QR factorization of an m-by-n matrix A with m ≥ n.
-- Q can be used to reflect a vector in such a way that all coordinates but one disappear.
qrDecomp::Mat3 -> (Mat3,Mat3)
qrDecomp m = (q, r)
  where
    x1 = _1 m
    a1 = let k = norm x1 in if _1 x1 > 0 then k else -k
    u1 = x1 &+ a1 *& (Vec3 1 0 0)
    q1 = householder $ mkNormal u1
    
    m' = (trimHead $ m .*. q1)::Mat2
    
    x2 = _1 m'
    a2 = let k = norm x2 in if _1 x2 > 0 then k else -k
    u2 = x2 &+ a2 *& (Vec2 1 0)
    t  = householder $ mkNormal u2 :: Mat2
    q2 = extendHeadWith 1 t
    
    q   = transpose q2 .*. transpose q1
    r   = m .*. q
    

-- | Gram-Schmidt is a very simple process to find an orthonormal basis for a matrix.
-- It takes the columns of the matrix and subtracts the part of each vector that is
-- parallel to the previous columns. After performing this process, we are left with
-- the Q matrix of our decomposition. 
orthoBasis :: Mat3 -> Mat3
orthoBasis (Mat3 a1 a2 a3) = Mat3 e1 e2 e3 
  where
    projAonB what dir = dir &* ((what &. dir) / (dir &. dir))
    e1 = normalize a1
    e2 = normalize $ a2 &- projAonB a2 e1
    e3 = normalize $ a3 &- projAonB a3 e1 &- projAonB a3 e2

qrGramSchmidt :: Mat3 -> (Mat3, Mat3)
qrGramSchmidt m = (q,r)
  where q = orthoBasis m
        r = m .*. (transpose q)



