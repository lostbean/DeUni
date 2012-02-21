-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| First Face |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module DeUni.FirstSeed (getFirstEdge) where

import Prelude hiding (null, lookup)
import Data.List (map, foldl', filter, head, (\\), minimumBy, maximumBy)
import qualified Data.List as L
import Control.Applicative ((<$>))
import Data.Maybe
import Data.Array.Diff hiding (elems)

import Hammer.Math.Vector

import DeUni.GeometricTools
import DeUni.Types

-- | Finds the first valid edge (a,b) across the division plane alpha laying
-- on the convex hull. The is a base to construct either the fist valid regular 
-- triangulation or the first face of the convex hull. Functons are generics and
-- the algorithm works for 2D and 3D, but no attempt was made for higher dimensions.
getFirstEdge::(PointND a)=> Plane a -> SetPoint a -> [PointPointer] -> [PointPointer] -> Maybe (PointPointer, PointPointer)
getFirstEdge divPlane sP ps1 ps2 = do
  (p1, d1) <- getMaxDistPoint divPlane sP ps1
  (p2, d2) <- getMaxDistPoint divPlane sP ps2
  if d1 > d2
    then do
      let refdir = normalize $ getProjOnPlane divPlane (sP!.p1)
      (p, _) <- getMaxDistPointOnDir refdir sP ps2
      climber refdir divPlane sP p1 p ps1 ps2
    else do
      let refdir = normalize $ getProjOnPlane divPlane (sP!.p2)
      (p, _) <- getMaxDistPointOnDir refdir sP ps1
      climber refdir divPlane sP p p2 ps1 ps2

-- | Finds the point which lays on the most distant (from origin) plane perpendicular
-- to the division plane.
getMaxDistPoint::(PointND a)=> Plane a -> SetPoint a -> [PointPointer] -> Maybe (PointPointer, Double)
getMaxDistPoint _ _ []         = Nothing
getMaxDistPoint divPlane sP (x:xs) = Just $ foldl' func d1 xs
  where
    d1     = (x, dist x)
    dist x = norm $ getProjOnPlane divPlane (sP!.x)
    func old x
      | d > (snd old) = (x, d)
      | otherwise     = old
      where d = dist x
              
-- | Finds the most distant point along the direction refdir.
getMaxDistPointOnDir::(PointND a)=> a -> SetPoint a -> [PointPointer] -> Maybe (PointPointer, Double)
getMaxDistPointOnDir _      _  []     = Nothing
getMaxDistPointOnDir refdir sP (x:xs) = Just $ foldl' func d1 xs
  where
    d1       = (x, dist x)
    dist x   = norm $ projAonB (sP!.x) refdir
    func old x
        | d > (snd old) = (x, d)
        | otherwise     = old
        where d = dist x              

-- | Given two initial points, it climbs outward the set of points using the touchPlane function (the plane conteins the two points and the direction perpenicular to refdir and conteined on the divison plane)
climber::(PointND a)=> a -> Plane a -> SetPoint a -> PointPointer -> PointPointer -> [PointPointer] -> [PointPointer] -> Maybe (PointPointer, PointPointer)
climber refdir divPlane sP p1 p2 ps1 ps2 = goTop p1 p2
  where
    getPP _     [] = Nothing
    getPP plane ps = return $ pointSetPartition (whichSideOfPlane plane) sP ps 
    
    goTop p1 p2 = do
      facePlane <- touchPlane refdir divPlane (sP!.p1) (sP!.p2)
      ppPS1     <- getPP facePlane ps1
      ppPS2     <- getPP facePlane ps2
      let
        ps1B1 = (pointsOnB1 ppPS1) \\ [p1, p2]
        ps1B2 = (pointsOnB2 ppPS1) \\ [p1, p2]
        ps2B1 = (pointsOnB1 ppPS2) \\ [p1, p2]
        ps2B2 = (pointsOnB2 ppPS2) \\ [p1, p2]

        okExit = return (p1, p2)

        move = goTop (nextP p1 ps1) (nextP p2 ps2)

        nextP p ps
          | mvdir > 0 = selByDist p (filter (\x -> dist x > 0 && x /= p) ps)
          | mvdir < 0 = selByDist p (filter (\x -> dist x < 0 && x /= p) ps)

        -- Find move direction
        mvdir = (planeNormal facePlane) &. (getProjOnPlane divPlane refdir)

        -- get farest plane (point contiened in) 
        selByDist x [] = x
        selByDist x ps = maximumBy (\a b -> dist a `compare` dist b) ps

        dist x =  (sP!.x &- sP!.p1) &. (planeNormal facePlane)

      case (ps1B1, ps1B2, ps2B1, ps2B2) of
        ([], _ , [], _ ) -> okExit
        (_ , [], _ , []) -> okExit
        _                -> move
          

