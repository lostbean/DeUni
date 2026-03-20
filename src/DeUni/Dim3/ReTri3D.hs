{-# LANGUAGE FlexibleContexts #-}

module DeUni.Dim3.ReTri3D where

import Linear.Decomp
import Linear.Mat
import Linear.Vect

import DeUni.Types

{- | Based on the papers: "Parallel dynamic and kinetic regular triangulation in three dimensions" (1993) and
"A data-parallel algorithm for three-dimensional Delaunay triangulation and its implementation" (2005)
-}
getCircumSphere :: WPoint Vec3 -> WPoint Vec3 -> WPoint Vec3 -> WPoint Vec3 -> (Double, Vec3D)
getCircumSphere a b c d = (radius, center)
  where
    (_, center) = getFaceDistCenter a b c d
    radius = (normsqr $ point a &- center) - weight a

getFaceDistCenter :: WPoint Vec3 -> WPoint Vec3 -> WPoint Vec3 -> WPoint Vec3 -> (Double, Vec3D)
getFaceDistCenter a b c d
    | isDegenerate = (0, point a)
    | otherwise = (signDist, center)
  where
    center = (-0.5) *& ((mux *& q1) &+ (muy *& q2) &+ (muz *& q3))
    signDist = if signRDet r then -dist else dist
    dist = muz * 0.5 + (q3 &. (point a))
    (mux, muy, muz) = solveMu ((-1) *& getAlpha a b c d) r
    (Mat3 q1 q2 q3) = q

    m = getM a b c d
    -- hand made QR for row vector matrix
    q = orthoRowsGram m
    r = m .*. transpose q
    -- Check if the tetrahedron is degenerate (near-zero volume).
    -- Volume is proportional to |det(M)| = |scalar triple product of edges|.
    e1 = point b &- point a
    e2 = point c &- point a
    e3 = point d &- point a
    tripleProduct = abs (e1 &. (e2 &^ e3))
    maxEdgeSq = maximum [normsqr e1, normsqr e2, normsqr e3]
    -- Degenerate when volume is tiny relative to edge lengths cubed
    isDegenerate = maxEdgeSq < 1e-20 || tripleProduct < 1e-6 * maxEdgeSq * sqrt maxEdgeSq

getM :: WPoint Vec3 -> WPoint Vec3 -> WPoint Vec3 -> WPoint Vec3 -> Mat3D
getM a b c d = Mat3 (point b &- point a) (point c &- point a) (point d &- point a)

getAlpha :: WPoint Vec3 -> WPoint Vec3 -> WPoint Vec3 -> WPoint Vec3 -> Vec3D
getAlpha a b c d = Vec3 (fun b a) (fun c a) (fun d a)
  where
    fun x y = (normsqr . point) x - (normsqr . point) y - weight x + weight y

solveMu :: Vec3D -> Mat3D -> (Double, Double, Double)
solveMu a (Mat3 r1 r2 r3) = (mux, muy, muz)
  where
    (ax, ay, az) = unVec3 a
    (r11, _, _) = unVec3 r1
    (r12, r22, _) = unVec3 r2
    (r13, r23, r33) = unVec3 r3
    mux = ax / r11
    muy = (ay - mux * r12) / r22
    muz = (az - mux * r13 - muy * r23) / r33

signRDet :: Mat3D -> Bool
signRDet = (>= 0) . product . diagVec
