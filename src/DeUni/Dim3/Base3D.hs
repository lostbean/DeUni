
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE OverlappingInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}


module DeUni.Dim3.Base3D where

import Control.Applicative ((<$>))
import Control.Monad.State.Lazy
import Data.Array.Diff hiding (elems)
import Data.List (map, foldl', filter, head, (\\), minimumBy, maximumBy)
import qualified Data.List as L

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import Math.Vector


instance PointND Point3D where
  data Box Point3D     =  Box3D
    { xMax3D::Double
    , xMin3D::Double
    , yMax3D::Double
    , yMin3D::Double
    , zMax3D::Double
    , zMin3D::Double
    } deriving (Show, Eq)
  
  data Plane Point3D   = Plane3D
    { plane3DNormal::Vec3
    , plane3DDist  ::Double
    } deriving (Show, Eq)
 
  data S0 Point3D      = Edge3D
    { edge3DL :: PointPointer
    , edge3DR :: PointPointer
    } deriving (Show)
  
  data S1 Point3D      = Face3D
   { face3DPoints   :: (PointPointer, PointPointer, PointPointer)
   , face3DND       :: Point3D
   } deriving (Show)
  
  data S2 Point3D      = Tetrahedron
   { circumSphereCenter :: Point3D
   , circumRadius       :: Double  
   , tetraPoints        :: (PointPointer, PointPointer, PointPointer, PointPointer)
   } deriving (Show, Eq)

  compS0 a b = compEdge (edge3DL a) (edge3DR a) (edge3DL b) (edge3DR b)    
  
  compS1 a b = compFace (face3DPoints a) (face3DPoints b)
  
  isInBox box (Vec3 x y z) = let 
    between min max x
      --will get points on the edge of the box and store if P1 those are on the commun face
      | min < max = (x >= min) && (max >= x)
      | min > max = (x <= min) && (max <= x)
      | otherwise = error ("Zero size box: " ++ show (box))
    in between (xMin3D box) (xMax3D box) x
    && between (yMin3D box) (yMax3D box) y
    && between (zMin3D box) (zMax3D box) z


  planeNormal = plane3DNormal
  planeDist   = plane3DDist
  makePlane n dist = Plane3D n dist

  calcPlane sp face
    | nSize == 0 = Nothing
    | d >= 0     = Just $ makePlane normN d
    | d < 0      = Just $ makePlane (inv normN) (-d)
    where
      (a,b,c) = face3DPoints face
      n       = (sp!.b &- sp!.a) &^ (sp!.c &- sp!.a)
      nSize   = len n
      normN   = normalize n
      d       = normN &. (sp!.a)
      inv n   = (-1) *& n

  touchPlane refdir divPlane a b
    | nSize == 0 = Nothing
    | d >= 0     = Just $ makePlane normND  d
    | otherwise  = Just $ makePlane (neg normND) (-d)
    where
      refOnDivPlane = getProjOnPlane divPlane refdir
      fixRotDir     = normalize $ (planeNormal divPlane) &^ refOnDivPlane
      nd            = (b &- a) &^ fixRotDir
      nSize         = norm nd
      normND        = normalize nd
      d             = normND &. a
  
  cutBox box subB
    | null subB  = smartBox box box
    | otherwise = func box subB
    where
      func sub [] = smartBox box sub
      func sub (p:ps) = case p of
        B1 -> func (halfBox1.snd $ smartBox sub sub) ps
        B2 -> func (halfBox2.snd $ smartBox sub sub) ps
        _  -> func sub ps
      smartBox box subbox@Box3D{..}
        | (deltaX >= (max deltaY deltaZ)) = (Plane3D (Vec3 1 0 0) halfX, cutX)
        | (deltaY >= (max deltaX deltaZ)) = (Plane3D (Vec3 0 1 0) halfY, cutY)
        | otherwise                       = (Plane3D (Vec3 0 0 1) halfZ, cutZ)
        where
          cutX = BoxPair box { xMax3D = halfX } box { xMin3D = halfX }
          cutY = BoxPair box { yMax3D = halfY } box { yMin3D = halfY }
          cutZ = BoxPair box { zMax3D = halfZ } box { zMin3D = halfZ }
          deltaX = abs (xMax3D - xMin3D)
          deltaY = abs (yMax3D - yMin3D)
          deltaZ = abs (zMax3D - zMin3D)
          halfX = (xMax3D + xMin3D)/2
          halfY = (yMax3D + yMin3D)/2
          halfZ = (zMax3D + zMin3D)/2


getThrirdPoint :: (PointND Point3D) => SetPoint Point3D -> PointPointer -> PointPointer -> [PointPointer] -> Maybe (PointPointer, Point3D)
getThrirdPoint sP pA pB ps = scan ps
  where
    cleanList x = ps \\ [pA, pB, x]
    scan [] = Nothing
    scan (x:xs)
      | x == pA || x == pB = scan xs
      | otherwise =
        let face = Face3D { face3DPoints = (pA, pB, x), face3DND = undefined }
        in case calcPlane sP face of
          Just plane
            | (L.null.pointsOnB1) pp &&
              (L.null.pointsOnB2) pp -> Nothing
            | (L.null.pointsOnB1) pp -> return (x, nd)
            | (L.null.pointsOnB2) pp -> return (x, ind)
            | otherwise              -> scan xs
            where pp = pointSetPartition (whichSideOfPlane plane) sP (cleanList x)
                  nd = planeNormal plane
                  ind = nd &* (-1)
          Nothing    -> scan xs
          