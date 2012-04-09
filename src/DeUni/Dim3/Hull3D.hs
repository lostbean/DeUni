-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Hull3D |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE TypeFamilies #-}

module DeUni.Dim3.Hull3D where

import Prelude hiding (null, lookup)
import Data.List (map, foldl',null, (\\))
import Control.Applicative ((<$>))
import Data.Maybe

import Hammer.Math.Vector

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim3.Base3D

instance Buildable S1 Point3D where
  type Sub S1    = S0
  buildUnit      = makeFace
  build1stUnit   = makeFirstFace
  getAllSubUnits = extractAllFaceEdges
  subUnitPos     = edge3DPos 

makeFirstFace::Plane Point3D -> SetPoint Point3D -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (S1 Point3D)
makeFirstFace alpha sP sideA sideB ps = do
  (pA, pB) <- getFirstEdge alpha sP sideA sideB
  (pC, nd) <- getThrirdPoint sP pA pB ps
  return $ Face3D (pA, pB, pC)
  
edge3DPos pairBox sP e = edgePos pairBox sP (edge3DL.activeUnit $ e) (edge3DR.activeUnit $ e)

extractAllFaceEdges::SetPoint Point3D -> S1 Point3D -> [ActiveSubUnit S1 Point3D]
extractAllFaceEdges sP sigma = 
  let (a,b,c) = face3DPoints sigma
  in  [ ActiveUnit (Edge3D a b) c undefined
      , ActiveUnit (Edge3D b c) a undefined 
      , ActiveUnit (Edge3D c a) b undefined ]

makeFace::ActiveSubUnit S1 Point3D -> SetPoint Point3D -> [PointPointer] -> Maybe (S1 Point3D)
makeFace _ _ [] = Nothing
makeFace e sP ps = do
  (pC, ang) <- findNext <$> get1stAng ps
  return $ Face3D { face3DPoints = (pA, pB, pC) }
  where
    pA = (edge3DL.activeUnit) e
    pB = (edge3DR.activeUnit) e

    get1stAng []     = Nothing
    get1stAng (p:ps) = case calcAngBetweenSimplex e sP p of 
      Nothing  -> get1stAng ps
      Just ang -> Just (p, ang)

    scanMin old@(p0, ang0) p = case calcAngBetweenSimplex e sP p of
      Nothing -> old
      Just ang
        | p0 == p    -> old
        | truncation > abs (ang - ang0) -> let
          dist = getSignDist sP pA pB
          in if dist p0 > dist p
             then (p, ang)
             else old
        | ang < ang0 -> (p, ang)
        | otherwise  -> old

    findNext a1@(p, ang) = foldl' scanMin a1 ps


calcAngBetweenSimplex::ActiveSubUnit S1 Point3D -> SetPoint Point3D -> PointPointer -> Maybe Double
calcAngBetweenSimplex ae sP p
  | pA==p || pB==p || pC==p = Nothing
  | otherwise               = return $ (normalToEdge vOnOldFace) &. (normalToEdge vOnNewFace)
  where
    pA             = (edge3DL.activeUnit) ae
    pB             = (edge3DR.activeUnit) ae
    pC             = assocP ae
    pB'            = sP!.pB
    edge           = sP!.pA &- pB'
    vOnOldFace     = sP!.pC &- pB'
    vOnNewFace     = sP!.p  &- pB'
    normalToEdge x = normalofAtoB x edge

 
