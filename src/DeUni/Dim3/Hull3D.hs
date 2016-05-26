{-# LANGUAGE
    FlexibleContexts
  , MultiParamTypeClasses
  , TypeFamilies
  #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
module DeUni.Dim3.Hull3D where

import Prelude
import Data.List (foldl')

import Linear.Vect

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim3.Base3D

-- ===============================| Hull3D |=================================

instance Buildable S1 Vec3 where
  type Sub S1    = S0
  buildUnit      = makeFace
  build1stUnit   = makeFirstFace
  getAllSubUnits = extractAllFaceEdges
  subUnitPos     = edge3DPos

makeFirstFace :: Plane Vec3 -> SetPoint Vec3 -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (S1 Vec3)
makeFirstFace alpha sP sideA sideB ps = do
  (pA, pB) <- getFirstEdge alpha sP sideA sideB
  (pC,  _) <- getThrirdPoint sP pA pB ps
  return $ Face3D (pA, pB, pC)

edge3DPos :: (PointND a, Buildable simplex Vec3, Sub simplex ~ S0)
          => BoxPair a -> SetPoint a -> ActiveSubUnit simplex Vec3 -> Position
edge3DPos pairBox sP e = let
  l = edge3DL . activeUnit $ e
  r = edge3DR . activeUnit $ e
  in edgePos pairBox sP l r

extractAllFaceEdges :: SetPoint Vec3 -> S1 Vec3 -> [ActiveSubUnit S1 Vec3]
extractAllFaceEdges _ sigma =
  let (a,b,c) = face3DPoints sigma
  in  [ ActiveUnit (Edge3D a b) c undefined
      , ActiveUnit (Edge3D b c) a undefined
      , ActiveUnit (Edge3D c a) b undefined ]

makeFace :: ActiveSubUnit S1 Vec3 -> SetPoint Vec3 -> [PointPointer] -> Maybe (S1 Vec3)
makeFace _ _ [] = Nothing
makeFace e sP ps = do
  (pC, _) <- findNext <$> get1stAng ps
  return Face3D { face3DPoints = (pA, pB, pC) }
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
        | isMainlyZero (ang - ang0) -> let
          dist = getSignDist sP pA pB
          in if dist p0 > dist p
             then (p, ang)
             else old
        | ang < ang0 -> (p, ang)
        | otherwise  -> old

    findNext a1 = foldl' scanMin a1 ps


calcAngBetweenSimplex :: ActiveSubUnit S1 Vec3 -> SetPoint Vec3 -> PointPointer -> Maybe Double
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
