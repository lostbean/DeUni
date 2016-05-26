{-# LANGUAGE
    FlexibleContexts
  , MultiParamTypeClasses
  , TypeFamilies
  #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
module DeUni.Dim2.Delaunay2D where

import Prelude
import Data.List ((\\))
import Data.Vector ((!))

import Linear.Vect

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim2.Base2D

import DeUni.Dim2.ReTri2D

-- ==============================| Delaunay2D |===================================

instance Buildable S2 Vec2 where
  type Sub S2    = S1
  buildUnit      = makeSimplex
  build1stUnit   = makeFirstSimplex
  getAllSubUnits = extractAllFaceEdges
  subUnitPos     = edge3DPos

makeFirstSimplex :: Plane Vec2 -> SetPoint Vec2 -> [PointPointer]
                 -> [PointPointer] -> [PointPointer] -> Maybe (S2 Vec2)
makeFirstSimplex alpha sP sideA sideB ps = do
  (pA, pB) <- getFirstEdge alpha sP sideA sideB
  let edge = Edge2D pA pB
  plane    <- calcPlane sP edge
  let
    newND
      | (null . pointsOnB1) pp = neg nd
      | otherwise              = nd
    nd = plane2DNormal plane
    psClean = ps \\ [pA, pB]
    pp = pointSetPartition (whichSideOfPlane plane) sP psClean
    actFace = ActiveUnit { activeUnit = edge, assocP = undefined, assocND = newND }
  makeSimplex actFace sP ps

edge3DPos :: (PointND a, Buildable simplex Vec2, Sub simplex ~ S1)
          => BoxPair a -> SetPoint a -> ActiveSubUnit simplex Vec2 -> Position
edge3DPos pairBox sP e = let
  l = edge2DL . activeUnit $ e
  r = edge2DR . activeUnit $ e
  in edgePos pairBox sP l r

extractAllFaceEdges :: SetPoint Vec2 -> S2 Vec2 -> [ActiveSubUnit S2 Vec2]
extractAllFaceEdges sP sigma = map toSimplexFace fsAll
  where
    (a,b,c) = face2DPoints sigma
    fsAll   = [((a,b), c), ((b,c), a), ((c,a), b)]
    toSimplexFace fp@((p1, p2), x) = ActiveUnit { activeUnit = Edge2D p1 p2
                                                , assocP     = x
                                                , assocND    = outterND fp
                                                }
    outterND ((p1, p2), x) = case plane2D (sP !. p1) (sP !. p2) of
      Just p
        | (sP!.x &- sP!.p1) &. nd > 0 -> neg nd
        | otherwise                   -> nd
        where nd = plane2DNormal p
      _ -> error "Delaunay2D: Bad face!!!"

makeSimplex :: ActiveSubUnit S2 Vec2 -> SetPoint Vec2 -> [PointPointer] -> Maybe (S2 Vec2)
makeSimplex actFace sP ps = do
  minR  <- findMinRadius
  buildSimplexFace minR
  where
    buildSimplexFace (_, i) =
      let (rad, center) = getCircumCircle (sP!a) (sP!b) (sP!i)
      in return Face2D { circleCenter = center
                       , circleRadius = rad
                       , face2DPoints = (a, b, i)
                       }
    -- | Remove points from face to avoid get 0.0 in findMin
    cleanP        = filter (\i -> isSideOk i && (i /= a) && (i /= b)) ps
    findMinRadius = findMinimunButZero getFaceDist cleanP
    getFaceDist   = fst . getFaceDistCenter (sP!a) (sP!b) . (sP!)
    isSideOk i    = (sP!.i &- sP!.a) &. nd > 0
    edge          = activeUnit actFace
    a             = edge2DR edge
    b             = edge2DL edge
    nd            = assocND actFace
