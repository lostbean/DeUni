{-# LANGUAGE
    FlexibleContexts
  , RecordWildCards
  , NamedFieldPuns
  , MultiParamTypeClasses
  , TypeSynonymInstances
  , TypeFamilies
  #-}
{-# OPTIONS_GHC -fno-warn-missing-signatures #-}

module DeUni.Dim2.Delaunay2D where

import Prelude hiding (null, lookup)
import Data.List      (null, (\\))
import Data.Vector    ((!))

import Hammer.Math.Algebra

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim2.Base2D

import DeUni.Dim2.ReTri2D

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Delaunay2D |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

instance Buildable S2 Point2D where
  type Sub S2    = S1
  buildUnit      = makeSimplex
  build1stUnit   = makeFirstSimplex
  getAllSubUnits = extractAllFaceEdges
  subUnitPos     = edge3DPos


makeFirstSimplex :: Plane Point2D -> SetPoint Point2D -> [PointPointer]
                 -> [PointPointer] -> [PointPointer] -> Maybe (S2 Point2D)
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

edge3DPos pairBox sP e = let
  l = edge2DL . activeUnit $ e
  r = edge2DR . activeUnit $ e
  in edgePos pairBox sP l r

extractAllFaceEdges :: SetPoint Point2D -> S2 Point2D -> [ActiveSubUnit S2 Point2D]
extractAllFaceEdges sP sigma = map toSimplexFace fsAll
  where
    (a,b,c) = face2DPoints sigma
    fsAll   = [((a,b), c), ((b,c), a), ((c,a), b)]
    toSimplexFace fp@((p1, p2), x) = ActiveUnit { activeUnit = (Edge2D p1 p2)
                                                , assocP = x
                                                , assocND = outterND fp }
    outterND ((p1, p2), x) = case plane2D (sP !. p1) (sP !. p2) of
      Just p
        | (sP!.x &- sP!.p1) &. nd > 0 -> neg nd
        | otherwise                   -> nd
        where nd = plane2DNormal p
      _ -> error "Delaunay2D: Bad face!!!"

makeSimplex :: ActiveSubUnit S2 Point2D -> SetPoint Point2D -> [PointPointer] -> Maybe (S2 Point2D)
makeSimplex actFace sP ps = do
  minR  <- findMinRadius
  sigma <- buildSimplexFace minR
  return sigma
  where
    buildSimplexFace (_, i) =
      let (rad, center) = getCircumCircle (sP!a) (sP!b) (sP!i)
      in return $ Face2D { circleCenter = center
                         , circleRadius = rad
                         , face2DPoints = (a, b, i) }
    -- | Remove points from face to avoid get 0.0 in findMin
    cleanP        = filter (\i -> (isSideOk i) && (i /= a) && (i /= b)) ps
    findMinRadius = findMinimunButZero getFaceDist cleanP
    getFaceDist   = fst . getFaceDistCenter (sP!a) (sP!b) . (sP!)
    isSideOk i    = (sP!.i &- sP!.a) &. nd > 0
    edge          = activeUnit actFace
    a             = edge2DR edge
    b             = edge2DL edge
    nd            = assocND actFace
