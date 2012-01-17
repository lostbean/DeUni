-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Hull3D |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE TypeFamilies #-}

module DeUni.Dim3.Hull2D where

import Prelude hiding (null, lookup)
import Data.List (map, foldl')
import Control.Applicative ((<$>))
import Data.Maybe
import Data.Array.Diff hiding (elems)

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim3.Base2D
import Math.Vector

instance Buildable S1 Point2D where
  type Sub S1    = S0
  buildUnit      = makeFace
  build1stUnit   = makeFirstFace
  getAllSubUnits = extractAllFaceEdges
  subUnitPos     = edge3DPos 

makeFirstFace::Plane Point2D -> SetPoint Point2D -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (S1 Point2D)
makeFirstFace alpha sP sideA sideB ps = do
  (pA, pB) <- getFirstEdge alpha sP sideA sideB
  nd       <- plane2DNormal <$> calcPlane sP (S1 pA pB _)
  return $ Edge2D pA pB nd
  
edge3DPos pairBox sP e = edgePos pairBox sP (edge3DL.activeUnit $ e) (edge3DR.activeUnit $ e)

extractAllFaceEdges::SetPoint Point3D -> S1 Point3D -> [ActiveSubUnit S1 Point3D]
extractAllFaceEdges sP sigma = fsAll
    where
        (a,b,c) = face3DPoints sigma
        nd      = face3DND sigma
        fsAll   = [ ActiveUnit (Edge3D a b) c nd
                  , ActiveUnit (Edge3D a c) b nd
                  , ActiveUnit (Edge3D b c) a nd ]

makeFace::ActiveSubUnit S1 Point3D -> SetPoint Point3D -> [PointPointer] -> Maybe (S1 Point3D)
makeFace _ _ [] = Nothing
makeFace e sP ps = do
    refPoint  <- get1stAng ps
    (pC, ang) <- findNext refPoint
    let face = Face3D { face3DPoints = (pA, pB, pC), face3DND = undefined }
    nd        <- planeNormal <$> calcPlane sP face
    return $ buildFace pC (func (sP!.pC) ang nd)
    where
        pA = (edge3DL.activeUnit) e
        pB = (edge3DR.activeUnit) e
        oldND = assocND e
        buildFace x nd = Face3D { face3DPoints = (pA, pB, x), face3DND = nd }

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
        disND p x   = calcAngBetweenSimplex $ fakeEdge pA pB (sP!.p &+ x)
        fakeEdge a b c = ActiveUnit (Edge3D a b) a c
        func::Point3D -> Double -> Point3D -> Point3D
        func p ang nd
            | ang > 0   = if nd &. oldND > 0 then neg nd else nd
            | ang < 0   = if nd &. oldND > 0 then nd else neg nd
            | otherwise = nd --if (disND p nd) > ang then nd else ind nd


calcAngBetweenSimplex::ActiveSubUnit S1 Point3D -> SetPoint Point3D -> PointPointer -> Double
calcAngBetweenSimplex ae sP p
    | pA==p || pB==p || pC==p = 0/0
    | otherwise               = (normalToEdge vOnOldFace) &. (normalToEdge vOnNewFace)
    where
        pA = (edge3DL.activeUnit) ae
        pB = (edge3DR.activeUnit) ae
        pC = assocP ae
        pB'        = sP!.pB
        edge       = sP!.pA &- pB'
        vOnOldFace = sP!.pC &- pB'
        vOnNewFace = sP!.p  &- pB'
        -- Projection A in B = B * (A°B)/(B°B)
        projToEdge x = edge &* k
            where k = (x &. edge) / (edge &. edge)
        normalToEdge x = normalize $ x &- (projToEdge x)

 
