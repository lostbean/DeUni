-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| DeHull |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Math.DeHull where

import Prelude hiding (null, lookup)
import Data.Vec hiding (map, length, fromList, fold, get, Map)
import Data.List (map, foldl')
import Control.Applicative ((<$>))
import Data.Maybe
import Data.Array.Diff hiding (elems)

import Math.GeometricTools
import Math.Types
import Math.FirstSeed

instance SubUnit Edge Face where
    buildUnit      = makeFace
    build1stUnit   = makeFirstFace
    getAllSubUnits = extractAllFaceEdges
    subUnitPos     = edgePos


edgePos::BoxPair -> SetPoint -> Edge -> Position
edgePos pairBox sP (Edge a b) = case (findPos $ sP!a, findPos $ sP!b) of
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


extractAllFaceEdges::SetPoint -> Face -> [ActiveSubUnit Edge]
extractAllFaceEdges sP sigma = fsAll
    where
        (a,b,c) = facePoints sigma
        nd      = refND sigma
        fsAll   = [ ActiveUnit (Edge a b) c nd
                  , ActiveUnit (Edge a c) b nd
                  , ActiveUnit (Edge b c) a nd ]

makeFace::ActiveSubUnit Edge -> SetPoint -> [PointPointer] -> Maybe Face
makeFace _ _ [] = Nothing
makeFace e sP ps = do
    refPoint  <- get1stAng ps
    (pC, ang) <- findNext refPoint
    nd        <- planeNormal <$> getPlane (sP!pA, sP!pB, sP!pC)
    return $ buildFace pC (func (sP!pC) ang nd)
    where
        pA = (pointL.activeUnit) e
        pB = (pointR.activeUnit) e
        oldND = assocND e
        buildFace x nd = Face (pA, pB, x) nd

        get1stAng []    = Nothing
        get1stAng (p:ps)
            | isNaN ang = get1stAng ps
            | otherwise = Just $ (p, ang)
            where ang = calcAngBetweenSimplex e sP p

        scanMax old x
            | isNaN ang       = old
            | ang < (snd old) = (x, ang)
            | otherwise       = old
            where ang = calcAngBetweenSimplex e sP x

        findNext a1 = return $ foldl' scanMax a1 ps
        ind x     = ((*) (Vec3D (-1) (-1) (-1))) x
        disND p x   = calcAngBetweenSimplex $ fakeEdge pA pB (sP!p + x)
        fakeEdge a b c = ActiveUnit (Edge a b) a c
        func::Point -> Double -> Point -> Point
        func p ang nd
            | ang > 0   = if dot nd oldND > 0 then ind nd else nd
            | ang < 0   = if dot nd oldND > 0 then nd else ind nd
            | otherwise = nd --if (disND p nd) > ang then nd else ind nd


calcAngBetweenSimplex::ActiveSubUnit Edge -> SetPoint -> PointPointer -> Double
calcAngBetweenSimplex ae sP p
    | pA==p || pB==p || pC==p = 0/0
    | otherwise               = dot (normalToEdge vOnOldFace) (normalToEdge vOnNewFace)
    where
        pA = (pointL.activeUnit) ae
        pB = (pointR.activeUnit) ae
        pC = assocP ae
        pB'        = sP!pB
        edge       = sP!pA - pB'
        vOnOldFace = sP!pC - pB'
        vOnNewFace = sP!p  - pB'
        -- Projection A in B = B * (A°B)/(B°B)
        projToEdge x = edge * (Vec3D k k k)
            where k = (dot x edge) / (dot edge edge)
        normalToEdge x = normalize $ x - (projToEdge x)

 
