-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| First Face |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module Math.FirstSeed
       ( makeFirstFace ) where

import Prelude hiding (null, lookup)
import Data.Vec hiding (map, length, fromList, fold, get, Map)
import Data.List (map, foldl', filter, head, (\\), minimumBy, maximumBy)
import qualified Data.List as L
import Control.Applicative ((<$>))
import Data.Maybe
import Data.Array.Diff hiding (elems)

import Math.GeometricTools
import Math.Types

makeFirstFace::Plane -> SetPoint -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe Face
makeFirstFace alpha sP sideA sideB ps = do
    (pA, pB) <- getFirstEdge alpha sP sideA sideB
    (pC, nd) <- getThrirdPoint sP pA pB ps
    return $ Face (pA, pB, pC) nd


getFirstEdge::Plane -> SetPoint -> [PointPointer] -> [PointPointer] -> Maybe (PointPointer, PointPointer)
getFirstEdge divPlane sP ps1 ps2 = do
    (p1, d1) <- getFirstPoint divPlane sP ps1
    (p2, d2) <- getFirstPoint divPlane sP ps2
    if d1 > d2
        then climber B1 divPlane sP p1 p2 ps1 ps2
        else climber B2 divPlane sP p1 p2 ps1 ps2


getFirstPoint::Plane -> SetPoint -> [PointPointer] -> Maybe (PointPointer, Double)
getFirstPoint _ _ []         = Nothing
getFirstPoint alpha sP (x:xs) = Just $ foldl' func d1 xs
    where
    d1     = (x, dist x)
    dist x = norm $ getProjOnPlane alpha (sP!x)
    func old x
        | d > (snd old) = (x, d)
        | otherwise     = old
        where d = dist x

climber::Position -> Plane -> SetPoint -> PointPointer -> PointPointer -> [PointPointer] -> [PointPointer] -> Maybe (PointPointer, PointPointer)
climber mostDistPoint divPlane sP p1 p2 ps1 ps2 = goTop p1 p2
    where
        goTop p1 p2 = do
            ppPS1     <- getPP (sP!p1) (sP!p1 + projPlaneND) (sP!p2) ps1
            ppPS2     <- getPP (sP!p1) (sP!p1 + projPlaneND) (sP!p2) ps2
            facePlane <- getPlane (sP!p1,  (sP!p1 + projPlaneND), sP!p2)
            let
                --In a ideal world (p+projPlanaND) should also be execluded
                ps1B1 = (pointsOnB1 ppPS1) \\ [p1, p2]
                ps1B2 = (pointsOnB2 ppPS1) \\ [p1, p2]
                ps2B1 = (pointsOnB1 ppPS2) \\ [p1, p2]
                ps2B2 = (pointsOnB2 ppPS2) \\ [p1, p2]

                okExit = return (p1, p2)

                move
                    | moveDir   = goTop (selByDist p1 ps1B1) (selByDist p2 ps2B1)
                    | otherwise = goTop (selByDist p1 ps1B2) (selByDist p2 ps2B2)

                moveDir = 0 < (dot (planeNormal facePlane) projPOnDivPlane)

                selByDist x [] = x
                selByDist _ ps = maximumBy (\a b -> dist a `compare` dist b) ps

                dist x = abs $ dot (sP!x - sP!p1) (planeNormal facePlane)

            case (ps1B1, ps1B2, ps2B1, ps2B2) of
                ([], _ , [], _ ) -> okExit
                (_ , [], _ , []) -> okExit
                _                -> move

        getPP::Point -> Point -> Point -> [PointPointer] -> Maybe PointPartition
        getPP _  _  _  [] = Nothing
        getPP p1 p2 p3 ps = func <$> getPlane (p1, p2, p3)
                where func x = pointSetPartition (whichSideOfPlane x) sP ps 

        projPOnDivPlane
            | mostDistPoint == B1 = calcProjOnDivPlane (sP!p1)
            | mostDistPoint == B2 = calcProjOnDivPlane (sP!p2)

        projPlaneND = (pack.normalize.(cross (unpack $ planeNormal divPlane)).unpack) projPOnDivPlane
        calcProjOnDivPlane = getProjOnPlane divPlane


getThrirdPoint::SetPoint -> PointPointer -> PointPointer -> [PointPointer] -> Maybe (PointPointer, Point)
getThrirdPoint sP pA pB ps = scan ps
    where
    cleanList x = ps \\ [pA, pB, x]
    scan [] = Nothing
    scan (x:xs)
        | x == pA || x == pB = scan xs
        | otherwise = case getPlane (sP!pA, sP!pB, sP!x) of
            Just plane
                | (L.null.pointsOnB1) pp &&
                  (L.null.pointsOnB2) pp -> Nothing
                | (L.null.pointsOnB1) pp -> return (x, nd)
                | (L.null.pointsOnB2) pp -> return (x, ind)
                | otherwise              -> scan xs
                where pp = pointSetPartition (whichSideOfPlane plane) sP (cleanList x)
                      nd = planeNormal plane
                      ind = nd * (Vec3D (-1) (-1) (-1))
            Nothing    -> scan xs
