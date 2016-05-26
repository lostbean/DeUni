{-# LANGUAGE
    FlexibleContexts
  , MultiParamTypeClasses
  , TypeFamilies
  #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
module DeUni.Dim3.Delaunay3D where

import Prelude hiding (null, lookup)
import Data.List      (null, (\\))
import Data.Vector    ((!))

import Linear.Vect

import DeUni.Dim3.Hull3D ()

import DeUni.GeometricTools
import DeUni.Types
import DeUni.Dim3.Base3D
import DeUni.Dim3.ReTri3D

-- ============================| Delaunay3D |============================

instance Buildable S2 Vec3 where
  type Sub S2    = S1
  buildUnit      = makeSimplex
  build1stUnit   = makeFirstSimplex
  getAllSubUnits = extractAllSimplexFaces
  subUnitPos     = face3DPos

makeFirstSimplex :: Plane Vec3 -> SetPoint Vec3 -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (S2 Vec3)
makeFirstSimplex alpha sP sideA sideB ps = do
  face  <- build1stUnit alpha sP sideA sideB ps
  plane <- calcPlane sP face
  let
    newND
      | (null.pointsOnB1) pp = nd
      | otherwise            = neg nd
    (a,b,c) = face3DPoints face
    nd = plane3DNormal plane
    psClean = ps \\ [a, b, c]
    pp = pointSetPartition (whichSideOfPlane plane) sP psClean
    actFace = ActiveUnit { activeUnit = face, assocP = undefined, assocND = newND }
  makeSimplex actFace sP ps

face3DPos :: (PointND a, Buildable simplex Vec3, Sub simplex ~ S1)
          => BoxPair a -> SetPoint a -> ActiveSubUnit simplex Vec3 -> Position
face3DPos pairBox sP face = let (a, b, c) = (face3DPoints.activeUnit) face in facePos pairBox sP a b c

extractAllSimplexFaces :: SetPoint Vec3 -> S2 Vec3 -> [ActiveSubUnit S2 Vec3]
extractAllSimplexFaces sP sigma = map toSimplexFace fsAll
  where
    (a,b,c,d) = tetraPoints sigma
    fsAll     = [((a,b,d), c), ((a,d,c), b), ((d,b,c), a), ((a,b,c), d)]
    toSimplexFace fp@(f, x) = let
      nd = outterND fp
      in ActiveUnit { activeUnit = Face3D f, assocP = x, assocND = nd }
    outterND ((na, nb, nc), x) = let
      nd = normalize (sP!.nb &- sP!.na) &^ (sP!.nc &- sP!.na)
      in if (sP!.na &- sP!.x) &. nd > 0 then neg nd else nd

makeSimplex :: ActiveSubUnit S2 Vec3 -> SetPoint Vec3 -> [PointPointer] -> Maybe (S2 Vec3)
makeSimplex actFace sP ps = do
  minR  <- findMinRadius
  buildSimplexFace minR
  where
    buildSimplexFace (_, d) =
      let (rad, center) = getCircumSphere (sP!a) (sP!b) (sP!c) (sP!d)
      in return Tetrahedron { circumSphereCenter = center
                            , circumSphereRadius = rad
                            , tetraPoints  = (a,b,c,d) }
    -- | Remove points from face to avoid get 0.0 in findMin
    cleanP        = filter (\i -> isSideOk i && (i /= a) && (i /= b) && (i /= c)) ps
    findMinRadius = findMinimunButZero getFaceDist cleanP
    getFaceDist   = fst . getFaceDistCenter (sP!a) (sP!b) (sP!c) . (sP!)
    isSideOk i    = 0 < (sP!.a &- sP!.i) &. nd
    (a, b, c)     = (face3DPoints.activeUnit) actFace
    nd            = assocND actFace
