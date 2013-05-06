module DeUni.FirstSeed (getFirstEdge) where

import Prelude hiding (null, lookup)
import Data.List      (foldl', minimumBy, maximumBy)
import Data.Vector    ((!))

import Hammer.Math.Algebra

import DeUni.GeometricTools
import DeUni.Types

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| First Face |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-- | Finds the first valid edge (a,b) across the division plane alpha laying
-- on the convex hull. The is a base to construct either the fist valid regular 
-- triangulation or the first face of the convex hull. Functons are generics and
-- the algorithm works for 2D and 3D, but no attempt was made for higher dimensions.
getFirstEdge :: (PointND a)=> Plane a -> SetPoint a -> [PointPointer] -> [PointPointer] -> Maybe (PointPointer, PointPointer)
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
getMaxDistPoint :: (PointND a)=> Plane a -> SetPoint a -> [PointPointer] -> Maybe (PointPointer, Double)
getMaxDistPoint _ _ []         = Nothing
getMaxDistPoint divPlane sP (i:is) = Just $ foldl' func d1 is
  where
    d1   = (i, dist i)
    dist = norm . getProjOnPlane divPlane . (sP !.)
    func old x
      | d > (snd old) = (x, d)
      | otherwise     = old
      where d = dist x
              
-- | Finds the most distant point along the direction refdir.
getMaxDistPointOnDir ::(PointND a)=> a -> SetPoint a -> [PointPointer] -> Maybe (PointPointer, Double) 
getMaxDistPointOnDir _      _  []     = Nothing
getMaxDistPointOnDir refdir sP (i:is) = Just $ foldl' func d1 is
  where
    d1       = (i, dist i)
    dist x   = norm $ projAonB (sP !. x) refdir
    func old x
        | d > (snd old) = (x, d)
        | otherwise     = old
        where d = dist x              

-- | Given two initial points, it climbs outward the set of points using the touchPlane function (the plane which conteins
-- the two points and the direction perpenicular to refdir and conteined on the divison plane)
climber :: (PointND a)=> a -> Plane a -> SetPoint a -> PointPointer -> PointPointer
        -> [PointPointer] -> [PointPointer] -> Maybe (PointPointer, PointPointer)
climber refdir divPlane sP pi1 pi2 ps1 ps2 = goTop pi1 pi2
  where
    goTop p1 p2 = do
      facePlane <- touchPlane refdir divPlane (sP !. p1) (sP !. p2)
      let
        getPP     = pointSetPartition (whichSideOfPlane facePlane) sP . cleanedPS
        cleanedPS = filter (\x -> x /= p1 && x /= p2)
        pp1       = getPP ps1
        pp2       = getPP ps2
        
        onB1_1    = pointsOnB1 pp1
        onB2_1    = pointsOnB2 pp1
        onPlane_1 = pointsOnPlane pp1
        
        onB1_2    = pointsOnB1 pp2
        onB2_2    = pointsOnB2 pp2
        onPlane_2 = pointsOnPlane pp2
        
        -- move touchPlane futher out
        moveUp pa pb = let
          n1 = nextP pa ps1
          n2 = nextP pb ps2
          in goTop n1 n2
            
        -- get the closest point to p
        moveClosestTo p ps = let
          dist = powerDist (sP!p) . (sP!)
          np   = minimumBy (\a b -> dist a `compare` dist b) ps
          in return (p, np)
        
        -- get the closest point to divPlane
        getClosestToDiv ps = let
          projection p = sP!.p &. (normalize.planeNormal) divPlane
          dist x       = abs $ projection x - planeDist divPlane 
          in minimumBy (\a b -> dist a `compare` dist b) ps

        -- get the farest point in the outward direction
        nextP p ps
          | mvdir > 0 = getFarest p (filter (\x -> faceDist x > 0 && x /= p) ps)
          | otherwise = getFarest p (filter (\x -> faceDist x < 0 && x /= p) ps)

        -- find move direction of facePlane (outwards the set of points)
        mvdir = (planeNormal facePlane) &. (getProjOnPlane divPlane refdir)

        -- get the farest point from facePlane
        getFarest p [] = p
        getFarest _ ps = maximumBy (\a b -> faceDist a `compare` faceDist b) ps

        faceDist x =  (sP!.x &- sP!.p1) &. (planeNormal facePlane)

      case (onB1_1, onPlane_1, onB2_1, onB1_2, onPlane_2, onB2_2) of
      --(11, 1_, 12, 21, 2_, 22)  
        ([], [], [], [], [], [])     -> Nothing
        
        ([], [], _ , [], [], _ )     -> return (p1, p2)
        (_ , [], [], _ , [], [])     -> return (p1, p2)

        (_ , [] , [], _ , onP, [])   -> moveClosestTo p1 (p2:onP)
        (_ , onP, [], _ , [] , [])   -> moveClosestTo p2 (p1:onP)
        ([], [] , _ , [], onP, _ )   -> moveClosestTo p1 (p2:onP)
        ([], onP, _ , [], [] , _ )   -> moveClosestTo p2 (p1:onP)
        ([], onP1, _ , [], onP2, _ ) -> let p = getClosestToDiv (p1:onP1)
                                        in moveClosestTo p (p2:onP2)
        
        _                            -> moveUp p1 p2
          

