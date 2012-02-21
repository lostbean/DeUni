module DeUni.Dim2.ReTri2D where

import Hammer.Math.Vector

import DeUni.Types
import DeUni.Dim2.Base2D


-- | Based on the papers: "Parallel dynamic and kinetic regular triangulation in three dimensions" (1993) and 
-- "A data-parallel algorithm for three-dimensional Delaunay triangulation and its implementation" (2005)

getCircumCircle::WPoint Point2D -> WPoint Point2D -> WPoint Point2D -> (Double, Vec2)
getCircumCircle a b c = (radius, center) 
  where
    (_, center) = getFaceDistCenter a b c
    radius      = (normsqr $ point a &- center) - weigth a
  
getFaceDistCenter::WPoint Point2D -> WPoint Point2D -> WPoint Point2D -> (Double, Vec2)
getFaceDistCenter a b c = (signDist, center)
  where
    center       = (-0.5) *& ((mux *& q1) &+ ( muy *& q2))
    signDist     = if signRDet r then -dist else dist
    dist         = muy*0.5 + (q2 &. (point a))
    m            = getM a b c
    (q,r)        = qrDecomp m
    (mux, muy)   = solveMu ((-1) *& (getAlpha a b c)) r
    (Mat2 q1 q2) = q
    
getM :: WPoint Point2D -> WPoint  Point2D -> WPoint Point2D -> Mat2
getM a b c = Mat2 (point b &- point a) (point c &- point a)

getAlpha :: WPoint Point2D -> WPoint Point2D -> WPoint Point2D -> Vec2
getAlpha a b c = Vec2 (fun b a) (fun c a)
  where fun x y = (normsqr.point) x - (normsqr.point) y - weigth x + weigth y

solveMu::Vec2 -> Mat2 -> (Double,Double)
solveMu a r@(Mat2 r1 r2) = (mux, muy)
  where
    (ax, ay)   = unvec2 a
    (r11, r21) = unvec2 r1
    (r12, r22) = unvec2 r2
    mux = ax / r11
    muy = (ay - mux*r12) / r22

signRDet::Mat2 -> Bool
signRDet (Mat2 r1 r2) = r11 * r22 >= 0
  where
    (r11, r21) = unvec2 r1
    (r12, r22) = unvec2 r2

unvec2 (Vec2 a b) = (a,b)


-- | A Householder reflection (or Householder transformation) is a transformation that
-- takes a vector and reflects it about some plane or hyperplane. We can use this operation
-- to calculate the QR factorization of an m-by-n matrix A with m ≥ n.
-- Q can be used to reflect a vector in such a way that all coordinates but one disappear.
-- "Stability of Householder QR Factorization for Weighted Least Squares Problems"
qrDecomp::Mat2 -> (Mat2,Mat2)
qrDecomp m = (q, r)
  where
    x  = _1 m
    a  = let k = norm x in if _1 x > 0 then k else -k
    u  = x &+ a *& (Vec2 1 0)
    q1 = householder $ mkNormal u
    q  = transpose q1
    r  = m .*. q
    

-- | Gram-Schmidt is a very simple process to find an orthonormal basis for a matrix.
-- It takes the columns of the matrix and subtracts the part of each vector that is
-- parallel to the previous columns. After performing this process, we are left with
-- the Q matrix of our decomposition. 
orthoBasis :: Mat2 -> Mat2
orthoBasis (Mat2 a1 a2) = Mat2 e1 e2 
  where
    projAonB what dir = dir &* ((what &. dir) / (dir &. dir))
    e1 = normalize a1
    e2 = normalize $ a2 &- projAonB a2 e1

qrGramSchmidt :: Mat2 -> (Mat2, Mat2)
qrGramSchmidt m = (q,r)
  where q = orthoBasis m
        r = m .*. (transpose q)



