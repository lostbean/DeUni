-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Delaunay3D |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE TypeFamilies #-}


module DeUni.Dim3.Delaunay3D where

import Prelude hiding (null, lookup)
import Data.List (map, foldl', null, (\\))
import Control.Applicative ((<$>))
import Data.Maybe
import Data.Vector ((!))

import Hammer.Math.Vector

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim3.Base3D
import DeUni.Dim3.Hull3D
import DeUni.Dim3.ReTri3D


instance Buildable S2 Point3D where
  type Sub S2    = S1
  buildUnit      = makeSimplex
  build1stUnit   = makeFirstSimplex
  getAllSubUnits = extractAllSimplexFaces
  subUnitPos     = face3DPos

makeFirstSimplex::Plane Point3D -> SetPoint Point3D -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (S2 Point3D)
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

face3DPos pairBox sP face = let (a, b, c) = (face3DPoints.activeUnit) face in facePos pairBox sP a b c

extractAllSimplexFaces::SetPoint Point3D -> S2 Point3D -> [ActiveSubUnit S2 Point3D]
extractAllSimplexFaces sP sigma = map toSimplexFace fsAll
  where
    (a,b,c,d) = tetraPoints sigma
    fsAll     = [((a,b,d), c), ((a,d,c), b), ((d,b,c), a), ((a,b,c), d)]
    toSimplexFace fp@(f, x) = let nd = outterND fp
      in ActiveUnit { activeUnit = Face3D f, assocP = x, assocND = nd }
    outterND ((a,b,c), x) = let nd = normalize (sP!.b &- sP!.a) &^ (sP!.c &- sP!.a) 
      in if (sP!.a &- sP!.x) &. nd > 0 then neg nd else nd

makeSimplex::ActiveSubUnit S2 Point3D -> SetPoint Point3D -> [PointPointer] -> Maybe (S2 Point3D)
makeSimplex actFace sP ps = do
  minR  <- findMinRadius
  sigma <- buildSimplexFace minR
  return sigma
  where
    buildSimplexFace (_, d) =
      let (radius, center) = getCircumSphere (sP!a) (sP!b) (sP!c) (sP!d) 
      in return Tetrahedron { circumSphereCenter = center
                            , circumSphereRadius = radius
                            , tetraPoints  = (a,b,c,d) }
    -- | Remove points from face to avoid get 0.0 in findMin
    cleanP        = filter (\i -> (isSideOk i) && (i /= a) && (i /= b) && (i /= c)) ps
    findMinRadius = findMinimunButZero getFaceDist sP cleanP
    getFaceDist   = fst . getFaceDistCenter (sP!a) (sP!b) (sP!c) . (sP!)
    isSideOk i    = 0 < (sP!.a &- sP!.i) &. nd
    face@(a,b,c)  = (face3DPoints.activeUnit) actFace
    nd            = assocND actFace
 