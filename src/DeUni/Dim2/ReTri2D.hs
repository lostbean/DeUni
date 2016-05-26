{-# LANGUAGE FlexibleContexts #-}
module DeUni.Dim2.ReTri2D where

import Linear.Vect
import Linear.Mat
import Linear.Decomp

import DeUni.Types

-- | Based on the papers: "Parallel dynamic and kinetic regular triangulation in three dimensions" (1993) and
-- "A data-parallel algorithm for three-dimensional Delaunay triangulation and its implementation" (2005)

getCircumCircle :: WPoint Vec2 -> WPoint Vec2 -> WPoint Vec2 -> (Double, Vec2D)
getCircumCircle a b c = (radius, center)
  where
    (_, center) = getFaceDistCenter a b c
    radius      = (normsqr $ point a &- center) - weight a

getFaceDistCenter :: WPoint Vec2 -> WPoint Vec2 -> WPoint Vec2 -> (Double, Vec2D)
getFaceDistCenter a b c = let
    center       = (-0.5) *& ((mux *& q1) &+ ( muy *& q2))
    dist         = muy * 0.5 + (q2 &. point a)
    (mux, muy)   = solveMu ((-1) *& getAlpha a b c) r
    (Mat2 q1 q2) = q

    m = getM a b c
    -- hand made QR for row vector matrix
    q = orthoRowsGram  m
    r = m .*. transpose q

    nd   = let (Vec2 x y) = point b &- point a in Vec2 (-y) x
    dir  = (nd &. (point c &- point a)) * (nd &. (center &- point a))
    absdist  = abs dist
    signDist = if dir > 0 then absdist else -absdist

    -- For some reason the sign determined by matrix don't work
    --signDist     = if signRDet r then -dist else dist

    in (signDist, center)

getM :: WPoint Vec2 -> WPoint  Vec2 -> WPoint Vec2 -> Mat2D
getM a b c = Mat2 (point b &- point a) (point c &- point a)

getAlpha :: WPoint Vec2 -> WPoint Vec2 -> WPoint Vec2 -> Vec2D
getAlpha a b c = Vec2 (fun b a) (fun c a)
  where fun x y = (normsqr.point) x - (normsqr.point) y - weight x + weight y

solveMu :: Vec2D -> Mat2D -> (Double, Double)
solveMu a (Mat2 r1 r2) = (mux, muy)
  where
    (ax, ay)   = unVec2 a
    (r11, _)   = unVec2 r1
    (r12, r22) = unVec2 r2
    mux = ax / r11
    muy = (ay - mux*r12) / r22

signRDet::Mat2D -> Bool
signRDet (Mat2 r1 r2) = r11 * r22 >= 0
  where
    (r11, _) = unVec2 r1
    (_, r22) = unVec2 r2
