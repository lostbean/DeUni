{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}

module DeUni.GeometricTools where

import Prelude hiding (null, lookup)
import Control.Applicative ((<$>))
import Control.Monad.State.Lazy
import Data.List (foldl')
import Data.Maybe

import Hammer.Math.Vector

import DeUni.Types

truncation::Double
truncation = 1e-10

-- | Projection A on B = B * (A°B)/(B°B)
projAonB::(Vector a, DotProd a) => a -> a -> a
projAonB a b = b &* ((a &. b) / (b &. b))

-- | Normal component of A to B
normalofAtoB::(Vector a, DotProd a) => a -> a -> a
normalofAtoB a b = normalize $ a &- (projAonB a b)

-- | retrieve the radius of a weigthed point
radius :: (PointND a) => WPoint a -> Double
radius = sqrt . weight

powerDist :: (PointND a) => WPoint a -> WPoint a -> Double
powerDist a b = (normsqr $ point a &- point b) - weight a - weight b

whichSideOfPlane::(PointND a) => Plane a -> a -> Position
whichSideOfPlane plane p 
  | truncation > delta = OnPlane
  | projection > dist  = B1
  | otherwise          = B2
  where
    projection = p &. (normalize.planeNormal) plane
    dist       = planeDist plane
    delta      = abs (projection - dist)
    
-- | Project a vector-point on the plane that goes throw the oringe.
--   It discard the distance on Plane data. It assumes that the plane pass throw the oringe
getProjOnPlane::(PointND a) => Plane a -> a -> a
getProjOnPlane plane p = projOnPlane
  where
    nd = planeNormal plane
    projOnPlane = p &- projAonB p nd


pointSetPartition::(PointND a) => (a -> Position) -> SetPoint a -> [PointPointer] -> PointPartition
pointSetPartition func sP ps = convert $ splitInBox ([],[],[]) ps
  where
    convert (p1,p2,pA) = PointPartition p1 p2 pA
    splitInBox (p1,p2,pA) []     = (p1,p2,pA)
    splitInBox (p1,p2,pA) (x:xs) = case func (sP!.x) of
      B1      -> splitInBox (x:p1,p2,pA) xs
      B2      -> splitInBox (p1,x:p2,pA) xs
      OnPlane -> splitInBox (p1,p2,x:pA) xs
      _       -> splitInBox (p1,p2,pA)   xs

whichBoxIsIt::(PointND a) => BoxPair a -> a -> Position
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
findMinimunButZero::(PointND a)=>(PointPointer -> Double) -> SetPoint a -> [PointPointer] -> Maybe (Double, PointPointer)
findMinimunButZero func sP ps = let
  ds     = map dist ps
  dist i = Just (func i, i)
  
  foldMaybe Nothing old = old
  foldMaybe new@(Just (dist, i)) old = case old of
    Just (olddist, oldi)
      | dist == 0      -> old
      | dist > olddist -> old
      | dist < olddist -> new
      | otherwise      -> error $ "Multiple points on circle or sphere! " ++ show new
    Nothing -> new

  in foldl' foldMaybe Nothing ds
     
-- | Performance can be improve by removing the duplicate call to "func" in "dropZero"
-- and the first "(func x, x)"
findMinimunButZero'::(PointND a)=>(PointPointer -> Maybe Double) -> SetPoint a -> [PointPointer] -> Maybe (Double, PointPointer)
findMinimunButZero' func sP ps = let
  ds = map dist ps
  
  dist i = case func i of
    Just x -> Just (x, i)
    _      -> Nothing
  
  foldMaybe Nothing old = old
  foldMaybe new@(Just (dist, i)) old = case old of
    Just (olddist, oldi)
      | dist == 0      -> old
      | dist > olddist -> old
      | dist < olddist -> new
      | otherwise      -> error $ "Multiple points on circle or sphere! " ++ show new
    Nothing -> new

  in foldl' foldMaybe Nothing ds