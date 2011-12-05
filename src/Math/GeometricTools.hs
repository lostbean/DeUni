{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}

module Math.GeometricTools where

import Prelude hiding (null, lookup)
import Data.Vec hiding (map, length, fromList, fold, get, Map)
import Control.Applicative ((<$>))
import Control.Monad.State.Lazy
import Data.Maybe
import System.Random
import Data.Array.Diff hiding (elems)


import Math.Types

pointSetPartition::(Point -> Position) -> SetPoint -> [PointPointer] -> PointPartition
pointSetPartition func sP ps =  convert $ splitInBox ([],[],[]) ps
    where
        convert (p1,p2,pA) = PointPartition p1 p2 pA
        splitInBox (p1,p2,pA) []     = (p1,p2,pA)
        splitInBox (p1,p2,pA) (x:xs) = case func (sP!x) of
            B1      -> splitInBox (x:p1,p2,pA) xs
            B2      -> splitInBox (p1,x:p2,pA) xs
            OnPlane -> splitInBox (p1,p2,x:pA) xs
            _       -> splitInBox (p1,p2,pA)   xs

whichBoxIsIt::BoxPair -> Point -> Position
whichBoxIsIt pairBox p
    | inbox1           = B1
    | inbox2           = B2
    | inbox1 && inbox2 = OnPlane
    | otherwise        = None
    where
        inbox1 = isInBox (halfBox1 pairBox) p
        inbox2 = isInBox (halfBox2 pairBox) p
        isInBox box (Vec3D x y z) =
            between (xMin box) (xMax box) x
         && between (yMin box) (yMax box) y
         && between (zMin box) (zMax box) z
        between min max x
            | min < max = (x >= min) && (max >= x)   -- >= will get points on the edge of the box and store if P1 those are on the commun face
            | min > max = (x <= min) && (max <= x)
            | otherwise = error ("Zero size box: " ++ show (pairBox))

whichSideOfPlane::Plane -> Point -> Position
whichSideOfPlane plane p = case compare projection dist of
    EQ -> OnPlane
    GT -> B1
    LT -> B2
    where
        projection = (planeNormal plane) `dot` p
        dist = planeDist plane

-- TODO verify correctness
genBox::(SubUnit a b, Ord a, Show a)=>Box -> StateMBC a (Plane, BoxPair)
genBox box@Box { xMax, xMin, yMax, yMin, zMax, zMin } = do
    seed <- randomSeed <$> get
    let (out, newSeed) = dirSelector seed
    modify (\x -> x { randomSeed=newSeed })
    return out
    where
        dirSelector seed
            | dir == 1  = ((Plane (Vec3D 1 0 0) newX, BoxPair box { xMax = newX } box { xMin = newX }), s4)
            | dir == 2  = ((Plane (Vec3D 0 0 1) newY, BoxPair box { zMax = newY } box { zMin = newY }), s4)
            | otherwise = ((Plane (Vec3D 0 1 0) newZ, BoxPair box { yMax = newZ } box { yMin = newZ }), s4)
            where
                (dir, s1) = randomR (1,3::Int) seed
                (a,   s2) = randomR (0,1) s1
                (b,   s3) = randomR (0,1) s2
                (c,   s4) = randomR (0,1) s3
                newX      = (xMax - xMin)*a + xMin
                newY      = (yMax + yMin)*b + yMin
                newZ      = (zMax + zMin)*c + zMin

isInBox::Box -> Point -> Bool
isInBox box (Vec3D x y z) = (xMax box > x && xMin box < x) &&
                            (yMax box > y && yMin box < y) &&
                            (zMax box > z && zMin box < z)

-- TODO verify correctness
genPlane::Box -> (Plane, BoxPair)
genPlane box@Box { xMax, xMin, yMax, yMin, zMax, zMin }
    | (deltaX >= (max deltaY deltaZ)) = (Plane (Vec3D 1 0 0) halfX, BoxPair box { xMax = halfX } box { xMin = halfX })
    | ((max deltaX deltaY) <= deltaZ) = (Plane (Vec3D 0 0 1) halfZ, BoxPair box { zMax = halfZ } box { zMin = halfZ })
    | otherwise                       = (Plane (Vec3D 0 1 0) halfY, BoxPair box { yMax = halfY } box { yMin = halfY })
    where
        deltaX = abs (xMax - xMin)
        deltaY = abs (yMax - yMin)
        deltaZ = abs (zMax - zMin)
        halfX = (xMax + xMin)/2
        halfY = (yMax + yMin)/2
        halfZ = (zMax + zMin)/2

-- | Project a vector\point on the plane that goes throw the oringe.
--   It discard the distance on Plane data. It assumes that the plane pass throw the oringe
getProjOnPlane::Plane -> Point -> Point
getProjOnPlane plane p = projOnPlane p
    where
        nd = planeNormal plane
        -- Projection A in B = B * (A°B)/(B°B)
        projNDPlane x = nd * (Vec3D k k k)
            where k = (dot x nd) / (dot nd nd)
        projOnPlane x = x - (projNDPlane x)


getPlane::(Point,Point,Point) -> Maybe Plane
getPlane (a,b,c)
    | nSize == 0 = Nothing
    | d >= 0     = Just $ Plane { planeNormal = normN
                                , planeDist   = d }
    | d < 0      = Just $ Plane { planeNormal = inv normN
                                , planeDist   = -d }
    where
        n     = pack $ (unpack (b - a)) `cross` (unpack (c - a))
        nSize = sqrt (dot n n)
        normN = normalize n
        d     = normN `dot` a
        inv n = (Vec3D (-1) (-1) (-1)) * n
