
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


module DeUni.Dim2.Base2D where

import Control.Applicative ((<$>))
import Control.Monad.State.Lazy
import Data.Array.Diff hiding (elems)
import Data.List (map, foldl', filter, head, (\\), minimumBy, maximumBy)
import qualified Data.List as L

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import Math.Vector


instance PointND Point2D where
  data Box Point2D     = Box2D
    { xMax2D::Double
    , xMin2D::Double
    , yMax2D::Double
    , yMin2D::Double
    } deriving (Show, Eq)
  
  data Plane Point2D   = Plane2D
    { plane2DNormal::Vec2
    , plane2DDist  ::Double
    } deriving (Show, Eq)
 
  data S0 Point2D      = P2D
    { point2D :: PointPointer
    } deriving (Show, Eq) 
  
  data S1 Point2D      = Edge2D
    { edge2DL  :: PointPointer
    , edge2DR  :: PointPointer
    } deriving (Show)  
               
  data S2 Point2D      = Face2D
    { face2DPoints :: (PointPointer, PointPointer, PointPointer)
    , circleRadius :: Double
    , circleCenter :: Point2D
    } deriving (Show)  
  
  compS0 a b = undefined

  compS1 a b = compEdge (edge2DL a) (edge2DR a) (edge2DL b) (edge2DR b)    
    
  isInBox box (Vec2 x y) = let 
    between min max x
      --will get points on the edge of the box and store if P1 those are on the commun face
      | min < max = (x >= min) && (max >= x)
      | min > max = (x <= min) && (max <= x)
      | otherwise = error ("Zero size box: " ++ show (box))
    in between (xMin2D box) (xMax2D box) x
    && between (yMin2D box) (yMax2D box) y


  planeNormal = plane2DNormal
  planeDist   = plane2DDist
  makePlane n dist = Plane2D n dist

  calcPlane sp edge = let
    a = edge2DL edge
    b = edge2DR edge
    in plane2D (sp!.a) (sp!.b)

  touchPlane refdir divPlane = plane2D
  
  cutBox box subB
    | null subB  = smartBox box box
    | otherwise = func box subB
    where
      func sub [] = smartBox box sub
      func sub (p:ps) = case p of
        B1 -> func (halfBox1.snd $ smartBox sub sub) ps
        B2 -> func (halfBox2.snd $ smartBox sub sub) ps
        _  -> func sub ps
      smartBox box subbox@Box2D{..}
        | deltaX >= deltaY = (Plane2D (Vec2 1 0) halfX, cutX)
        | otherwise        = (Plane2D (Vec2 0 1) halfY, cutY)
        where
          cutX = BoxPair box { xMax2D = halfX } box { xMin2D = halfX }
          cutY = BoxPair box { yMax2D = halfY } box { yMin2D = halfY }
          deltaX = abs (xMax2D - xMin2D)
          deltaY = abs (yMax2D - yMin2D)
          halfX = (xMax2D + xMin2D)/2
          halfY = (yMax2D + yMin2D)/2

plane2D::Vec2 -> Vec2 -> Maybe (Plane Point2D)
plane2D a b
  | nSize == 0 = Nothing
  | d >= 0     = Just $ makePlane normN d
  | d < 0      = Just $ makePlane (neg normN) (-d)
  where
    (Vec2 x y) = b &- a
    n     = Vec2 (-y) x
    nSize = len n
    normN = normalize n
    d     = normN &. a


