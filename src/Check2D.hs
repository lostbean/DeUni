-----------------------------------------------------------------------------
--
-- Module      :  DeUniChecker
-- Copyright   :
-- License     :  AllRightsReserved
--
-- Maintainer  :
-- Stability   :
-- Portability :
--
-- |
--
-----------------------------------------------------------------------------


{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}

module Check2D where

import Test.QuickCheck
import Control.Applicative
import Control.Monad
import Data.Array.Diff
import qualified Data.IntMap as IM
import qualified Data.Set as S
import qualified Data.List as L

import Hammer.Math.Vector

import DeUni.DeWall
import DeUni.Types
import DeUni.FirstSeed
import DeUni.GeometricTools
import DeUni.Dim2.Base2D
import DeUni.Dim2.Delaunay2D
import DeUni.Dim2.ReTri2D

import VTKGenSimplex


runChecker =  do
  let myArgs = Args {replay = Nothing, maxSuccess = 1000, maxDiscard = 5000, maxSize = 1000, chatty = True}
  print "Testing 1st edge.."    
  quickCheckWith myArgs prop_1stEdge
  
  print "Testing CircumCircle.."    
  quickCheckWith myArgs prop_CircumCircle
  
  print "Testing 1st Face.."    
  quickCheckWith myArgs prop_1stSimplex
  
  print "Testing Delaunay.."    
  quickCheckWith myArgs prop_Delaunay


instance Arbitrary (Box Point2D) where
  arbitrary = liftM2 getBox max min
    where 
      getBox max min = Box2D {xMax2D=max, xMin2D=min, yMax2D=max, yMin2D=min}
      min = choose (-500,-200)
      max = choose (500,200)

instance Arbitrary Vec2 where
  arbitrary = liftM2 Vec2 p p
    where p = choose (-200,200)

instance (Arbitrary a) => Arbitrary (WPoint a) where
  arbitrary = liftM2 WPoint s arbitrary
    where s = choose (1, 8)
  
instance (Ix a, Integral a, Arbitrary b) => Arbitrary (DiffArray a b) where
  arbitrary   =
    (\x -> listArray (0,fromIntegral (length x - 1)) x) <$> arbitrary 


error_precisson = (10e-3)

msgFail text = printTestCase ("\x1b[7m Fail: " ++ show text ++ "! \x1b[0m")


prop_Delaunay::Box Point2D -> SetPoint Point2D -> Property
prop_Delaunay box sp = (length ps) > 4 ==> whenFail (writeVTKfile "Delaunay2D_err.vtu" sp wall) fulltest
  where
    fulltest   = testWall .&&. testHull .&&. testSize
    (wall, st) = runDelaunay2D box sp ixps
    testw x    = testProperFace sp x
    testh x    = testHullEdge sp x
    testWall   = testIM testw wall
    testHull   = testSet testh (S.map activeUnit $ externalFaces st)
    ixps       = indices sp
    ps         = elems sp
    testSize   = msgFail "gen obj /= num add" $ IM.size wall == count st
    

testIM test map
  | IM.null map = err
  | otherwise   = let (x, xs) = IM.deleteFindMin map
                  in IM.fold (\a b -> b .&&. test a) (test x) xs
  where
    err          = msgFail "empty output" False

testSet test set
  | S.null set = err
  | otherwise  = let (x, xs) = S.deleteFindMin set
                 in S.fold (\a b -> b .&&. test a) (test x) xs
  where
    err          = msgFail "empty output" False


prop_CircumCircle::WPoint Vec2 -> WPoint Vec2 -> WPoint Vec2 -> Property
prop_CircumCircle a b c = msgFail ("bad center", map (powerDist (WPoint r center)) [a,b,c]) test
  where
    test         = and $ map testcenter [a,b,c]
    testcenter i = error_precisson > (abs $ powerDist i (WPoint r center))
    (r, center)  = getCircumCircle a b c


prop_1stSimplex::Box Point2D -> SetPoint Point2D -> Property
prop_1stSimplex box sP = pretest ==> test
  where
    (plane, pairBox) = cutBox box []
    p                = indices sP 
    pp               = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1               = pointsOnB1 pp
    p2               = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest          = (length p) > 4 && p1 /= [] && p2 /= []
    test             = case makeFirstSimplex plane sP p1 p2 p of
        Just sigma -> testProperFace sP sigma
        _          -> msgFail "non-gen 1st Face3D" False


prop_1stEdge::Box Point2D -> SetPoint Point2D -> Property
prop_1stEdge box sP = pretest ==> test
  where
    (plane, pairBox) = cutBox box []
    p                = indices sP
    pp               = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1               = pointsOnB1 pp
    p2               = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest          = (length p) > 4 && p1 /= [] && p2 /= []
    test             = case getFirstEdge plane sP p1 p2 of
        Just (pA, pB) -> let edge = Edge2D pA pB in testHullEdge sP edge
        _         -> msgFail "non-gen 1st Face3D" False


testProperFace::SetPoint Point2D -> S2 Point2D -> Property
testProperFace sP sigma = msgFail ("non empty circle", [pA,pB,pC], map (powerDist wC.(sP!)) (indices sP)) isSphereOK
                     -- .&&. prop_CircumCircle (sP!a) (sP!b) (sP!c)
  where
    (pA,pB,pC)     = face2DPoints sigma
    center         = circleCenter sigma
    radius         = circleRadius sigma
    wC             = WPoint radius center
    cleanP         = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC)) (indices sP)
    -- Test if the CircumSphere is empty
    isSphereOK     = and $ map testEmptySph cleanP
    testEmptySph i = 0 < (powerDist (sP!i) wC)


testHullEdge::SetPoint Point2D -> S1 Point2D -> Property
testHullEdge sP edge = test
  where
    pA     = edge2DR edge
    pB     = edge2DL edge
    pp     = (\x -> pointSetPartition (whichSideOfPlane x) sP cleanP) <$> plane
    plane  = calcPlane sP edge
    cleanP = filter (\i -> (i /= pA) && (i /= pB)) $ indices sP
    test   = case pp of
      Nothing -> msgFail "no plane from face" False
      Just x  -> case (pointsOnB1 x, pointsOnB2 x, pointsOnPlane x) of
        ([],[],[]) -> msgFail "no points on partition" False
        ([],[],_)  -> msgFail "all points on plane" False
        ([],_,_)   -> label "face on B1" True
        (_,[],_)   -> label "face on B2" True
        _          -> msgFail ("non-Hull face", x, edge) False