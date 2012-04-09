-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Delaunay2D |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE TypeFamilies #-}


module DeUni.Dim2.Delaunay2D where

import Prelude hiding (null, lookup)
import Data.List (map, foldl', null, (\\))
import Control.Applicative ((<$>))
import Data.Maybe
import Data.Vector ((!))

import Hammer.Math.Vector

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim2.Base2D

import DeUni.Dim2.ReTri2D

instance Buildable S2 Point2D where
  type Sub S2    = S1
  buildUnit      = makeSimplex
  build1stUnit   = makeFirstSimplex
  getAllSubUnits = extractAllFaceEdges
  subUnitPos     = edge3DPos


makeFirstSimplex::Plane Point2D -> SetPoint Point2D -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (S2 Point2D)
makeFirstSimplex alpha sP sideA sideB ps = do
  (pA, pB) <- getFirstEdge alpha sP sideA sideB
  let edge = Edge2D pA pB
  plane    <- calcPlane sP edge
  let 
    newND
      | (null.pointsOnB1) pp = neg nd
      | otherwise            = nd
    nd = plane2DNormal plane                               
    psClean = ps \\ [pA, pB]
    pp = pointSetPartition (whichSideOfPlane plane) sP psClean
    actFace = ActiveUnit { activeUnit = edge, assocP = undefined, assocND = newND }
  makeSimplex actFace sP ps

  
edge3DPos pairBox sP e = edgePos pairBox sP (edge2DL.activeUnit $ e) (edge2DR.activeUnit $ e)

extractAllFaceEdges::SetPoint Point2D -> S2 Point2D -> [ActiveSubUnit S2 Point2D]
extractAllFaceEdges sP sigma = map toSimplexFace fsAll
  where
    (a,b,c) = face2DPoints sigma
    fsAll   = [((a,b), c), ((b,c), a), ((c,a), b)]
    toSimplexFace fp@((a,b), x) = ActiveUnit { activeUnit=(Edge2D a b), assocP=x, assocND=outterND fp }
    outterND ((a,b), x) = case plane2D (sP!.a) (sP!.b) of
      Just p -> let nd = plane2DNormal p
                in if (sP!.x &- sP!.a) &. nd > 0 then neg nd else nd
      _      -> error "Delaunay2D: Bad face!!!"
      
makeSimplex::ActiveSubUnit S2 Point2D -> SetPoint Point2D -> [PointPointer] -> Maybe (S2 Point2D)
makeSimplex actFace sP ps = do
  minR  <- findMinRadius
  sigma <- buildSimplexFace minR
  return sigma
  where
    buildSimplexFace (dist, i) = 
      let (radius, center) = getCircumCircle (sP!a) (sP!b) (sP!i)
      in return $ Face2D { circleCenter = center
                         , circleRadius = radius
                         , face2DPoints = (a,b,i) }
    -- | Remove points from face to avoid get 0.0 in findMin
    cleanP        = filter (\i -> (isSideOk i) && (i /= a) && (i /= b)) ps
    findMinRadius = findMinimunButZero getFaceDist sP cleanP
    getFaceDist   = fst . getFaceDistCenter (sP!a) (sP!b) . (sP!)
    isSideOk i    = (sP!.i &- sP!.a) &. nd > 0
    edge          = activeUnit actFace
    a             = edge2DR edge
    b             = edge2DL edge
    nd            = assocND actFace

    

