{-# LANGUAGE FlexibleContexts #-}
module DeUni.GeometricTools where

import Prelude
import Data.List (foldl')

import Linear.Class

import DeUni.Types

-- | Projection A on B = B * (A°B)/(B°B)
projAonB::(LinearMap a v, DotProd a v, Fractional a) => v a -> v a -> v a
projAonB a b = b &* ((a &. b) / (b &. b))

-- | Normal component of A to B
normalofAtoB::(LinearMap a v, DotProd a v, Norm a v) => v a -> v a -> v a
normalofAtoB a b = normalize $ a &- (projAonB a b)

-- | retrieve the radius of a weigthed point
radius :: (PointND a) => WPoint a -> Double
radius = sqrt . weight

powerDist :: (PointND a) => WPoint a -> WPoint a -> Double
powerDist a b = (normsqr $ point a &- point b) - weight a - weight b

whichSideOfPlane :: (PointND a, Norm Double a) => Plane a -> a Double -> Position
whichSideOfPlane plane p
  | onPlane           = OnPlane
  | projection > dist = B1
  | otherwise         = B2
  where
    projection = p &. (normalize . planeNormal) plane
    dist       = planeDist plane
    onPlane    = isMainlyZero (projection - dist)

-- | Project a vector-point on the plane that goes throw the oringe.
--   It discard the distance on Plane data. It assumes that the plane pass throw the oringe
getProjOnPlane::(PointND a) => Plane a -> a Double -> a Double
getProjOnPlane plane p = projOnPlane
  where
    nd = planeNormal plane
    projOnPlane = p &- projAonB p nd

pointSetPartition::(PointND a) => (a Double -> Position) -> SetPoint a -> [PointPointer] -> PointPartition
pointSetPartition func sP ps = convert $ splitInBox ([],[],[]) ps
  where
    convert (p1,p2,pA) = PointPartition p1 p2 pA
    splitInBox (p1,p2,pA) []     = (p1,p2,pA)
    splitInBox (p1,p2,pA) (x:xs) = case func (sP!.x) of
      B1      -> splitInBox (x:p1,p2,pA) xs
      B2      -> splitInBox (p1,x:p2,pA) xs
      OnPlane -> splitInBox (p1,p2,x:pA) xs
      _       -> splitInBox (p1,p2,pA)   xs

whichBoxIsIt::(PointND a) => BoxPair a -> a Double -> Position
whichBoxIsIt pairBox p
  | inbox1           = B1
  | inbox2           = B2
  | inbox1 && inbox2 = OnPlane
  | otherwise        = None
    where
      inbox1 = isInBox (halfBox1 pairBox) p
      inbox2 = isInBox (halfBox2 pairBox) p


edgePos::(PointND a) => BoxPair a -> SetPoint a -> PointPointer -> PointPointer -> Position
edgePos pairBox sP a b = case (findPos $ sP!.a, findPos $ sP!.b) of
  (B1,B1)      -> B1
  (B2,B2)      -> B2
  (B1,OnPlane) -> B1
  (B2,OnPlane) -> B2
  (OnPlane,B1) -> B1
  (OnPlane,B2) -> B2
  (None,_)     -> None
  (_,None)     -> None
  _            -> CrossPlane
  where findPos  = whichBoxIsIt pairBox

facePos::(PointND a) => BoxPair a -> SetPoint a -> PointPointer -> PointPointer -> PointPointer -> Position
facePos pairBox sP a b c = case (findPos $ sP!.a, findPos $ sP!.b, findPos $ sP!.c) of
    (B1,B1,B1)                -> B1
    (B2,B2,B2)                -> B2
    (B1,OnPlane,OnPlane)      -> B1
    (B2,OnPlane,OnPlane)      -> B2
    (OnPlane,B1,OnPlane)      -> B1
    (OnPlane,B2,OnPlane)      -> B2
    (OnPlane,OnPlane,B1)      -> B1
    (OnPlane,OnPlane,B2)      -> B2
    (B1,B1,OnPlane)           -> B1
    (B2,B2,OnPlane)           -> B2
    (OnPlane,B1,B1)           -> B1
    (OnPlane,B2,B2)           -> B2
    (B1,OnPlane,B1)           -> B1
    (B2,OnPlane,B2)           -> B2
    (OnPlane,OnPlane,OnPlane) -> None
    (None,_,_)                -> None
    (_,None,_)                -> None
    (_,_,None)                -> None
    _                         -> CrossPlane
    where findPos  = whichBoxIsIt pairBox


-- | Performance can be improve by removing the duplicate call to "func" in "dropZero"
-- and the first "(func x, x)"
findMinimunButZero :: (PointPointer -> Double) -> [PointPointer] -> Maybe (Double, PointPointer)
findMinimunButZero func ps = let
  ds     = map dist ps
  dist i = Just (func i, i)

  foldMaybe Nothing old = old
  foldMaybe new@(Just (d, _)) old = case old of
    Just (olddist, _)
      | d == 0      -> old
      | d > olddist -> old
      | d < olddist -> new
      | otherwise   -> error $ "Multiple points on circle or sphere! " ++ show new
    Nothing -> new

  in foldl' foldMaybe Nothing ds

-- | Performance can be improve by removing the duplicate call to "func" in "dropZero"
-- and the first "(func x, x)"
findMinimunButZero' :: (PointPointer -> Maybe Double) -> [PointPointer] -> Maybe (Double, PointPointer)
findMinimunButZero' func ps = let
  ds = map dist ps

  dist i = case func i of
    Just x -> Just (x, i)
    _      -> Nothing

  foldMaybe Nothing old = old
  foldMaybe new@(Just (d, _)) old = case old of
    Just (olddist, _)
      | d == 0      -> old
      | d > olddist -> old
      | d < olddist -> new
      | otherwise   -> error $ "Multiple points on circle or sphere! " ++ show new
    Nothing -> new

  in foldl' foldMaybe Nothing ds
