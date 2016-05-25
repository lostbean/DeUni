{-# LANGUAGE
    FlexibleContexts
  , RecordWildCards
  , NamedFieldPuns
  , MultiParamTypeClasses
  , TypeSynonymInstances
  , TypeFamilies
  #-}
{-# OPTIONS_GHC -fno-warn-missing-signatures #-}

module DeUni.Dim3.Hull3D where

import Prelude hiding (null, lookup)
import Data.List (foldl')

import Hammer.Math.Algebra

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim3.Base3D

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Hull3D |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

instance Buildable S1 Point3D where
  type Sub S1    = S0
  buildUnit      = makeFace
  build1stUnit   = makeFirstFace
  getAllSubUnits = extractAllFaceEdges
  subUnitPos     = edge3DPos

makeFirstFace :: Plane Point3D -> SetPoint Point3D -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (S1 Point3D)
makeFirstFace alpha sP sideA sideB ps = do
  (pA, pB) <- getFirstEdge alpha sP sideA sideB
  (pC,  _) <- getThrirdPoint sP pA pB ps
  return $ Face3D (pA, pB, pC)

edge3DPos pairBox sP e = let
  l = edge3DL . activeUnit $ e
  r = edge3DR . activeUnit $ e
  in edgePos pairBox sP l r

extractAllFaceEdges :: SetPoint Point3D -> S1 Point3D -> [ActiveSubUnit S1 Point3D]
extractAllFaceEdges _ sigma =
  let (a,b,c) = face3DPoints sigma
  in  [ ActiveUnit (Edge3D a b) c undefined
      , ActiveUnit (Edge3D b c) a undefined
      , ActiveUnit (Edge3D c a) b undefined ]

makeFace :: ActiveSubUnit S1 Point3D -> SetPoint Point3D -> [PointPointer] -> Maybe (S1 Point3D)
makeFace _ _ [] = Nothing
makeFace e sP ps = do
  (pC, _) <- findNext <$> get1stAng ps
  return $ Face3D { face3DPoints = (pA, pB, pC) }
  where
    pA = (edge3DL.activeUnit) e
    pB = (edge3DR.activeUnit) e

    get1stAng []     = Nothing
    get1stAng (pp:pps) = case calcAngBetweenSimplex e sP pp of
      Nothing  -> get1stAng pps
      Just ang -> Just (pp, ang)

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

    findNext a1 = foldl' scanMin a1 ps


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
