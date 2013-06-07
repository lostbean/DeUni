module DeUni.Dim3.ReTri3D where

import Hammer.Math.Algebra

import DeUni.Types


-- | Based on the papers: "Parallel dynamic and kinetic regular triangulation in three dimensions" (1993) and 
-- "A data-parallel algorithm for three-dimensional Delaunay triangulation and its implementation" (2005)

getCircumSphere :: WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> (Double, Vec3)
getCircumSphere a b c d = (radius, center) 
  where
    (_, center) = getFaceDistCenter a b c d
    radius      = (normsqr $ point a &- center) - weight a
  
getFaceDistCenter :: WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> (Double, Vec3)  
getFaceDistCenter a b c d = (signDist, center)
  where
    center          = (-0.5) *& ((mux *& q1) &+ (muy *& q2) &+ (muz *& q3))
    signDist        = if signRDet r then -dist else dist
    dist            = muz*0.5 + (q3 &. (point a))
    (mux, muy, muz) = solveMu ((-1) *& (getAlpha a b c d)) r
    (Mat3 q1 q2 q3) = q

    m = getM a b c d
    -- hand made QR for row vector matrix
    q = orthoRowsGram  m
    r = m .*. transpose q
    
getM :: WPoint Point3D -> WPoint  Point3D -> WPoint Point3D -> WPoint Point3D -> Mat3
getM a b c d = Mat3 (point b &- point a) (point c &- point a) (point d &- point a)

getAlpha :: WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> Vec3
getAlpha a b c d = Vec3 (fun b a) (fun c a) (fun d a)
  where fun x y = (normsqr.point) x - (normsqr.point) y - weight x + weight y

solveMu :: Vec3 -> Mat3 -> (Double,Double,Double)
solveMu a (Mat3 r1 r2 r3) = (mux, muy, muz)
  where
    (ax, ay, az)    = unvec3 a
    (r11, _,   _)   = unvec3 r1
    (r12, r22, _)   = unvec3 r2
    (r13, r23, r33) = unvec3 r3
    mux = ax / r11
    muy = (ay - mux*r12) / r22
    muz = (az - mux*r13 - muy* r23) / r33        

signRDet :: Mat3 -> Bool
signRDet = (>= 0) . vecFoldr (*) . diagVec

unvec3 :: Vec3 -> (Double, Double, Double)
unvec3 (Vec3 a b c) = (a, b, c)




