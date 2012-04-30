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
import Data.Monoid ((<>))
import Data.Maybe (isJust)
import Data.IntMap (IntMap)
import Data.Set (Set)
import Data.Map (Map)
import qualified Data.IntMap as IM
import qualified Data.Map as Map
import qualified Data.Set as S
import qualified Data.List as L
import qualified Data.Vector as Vec
import Data.Vector (Vector, (!))

import Hammer.Math.Vector hiding (Vector)

import DeUni.DeWall
import DeUni.Types
import DeUni.FirstSeed
import DeUni.GeometricTools
import DeUni.Dim2.Base2D
import DeUni.Dim2.Delaunay2D
import DeUni.Dim2.ReTri2D

import RenderSVG


runChecker =  do
  let myArgs = Args {replay = Nothing, maxSuccess = 1000, maxDiscard = 500, maxSize = 1000, chatty = True}
  print "Testing 1st edge.."    
  quickCheckWith myArgs prop_1stEdge
  
  print "Testing CircumCircle.."    
  quickCheckWith myArgs prop_CircumCircle
  
  print "Testing 1st Face.."    
  quickCheckWith myArgs prop_1stSimplex
  
  print "Testing Delaunay.."    
  quickCheckWith myArgs prop_Delaunay
  

instance Arbitrary (Box Point2D) where
  arbitrary = let      
    getBoxSqr = liftM2 getBox max min
    getBoxRect = do
      max1 <- max
      max2 <- max
      min1 <- min
      min2 <- min
      return $ Box2D {xMax2D=max1, xMin2D=min1, yMax2D=max2, yMin2D=min2}
    getBox max min = Box2D {xMax2D=max, xMin2D=min, yMax2D=max, yMin2D=min}
    min = frequency [(1, elements [-550,-250]), (2, choose (-550,-250))]
    max = frequency [(1, elements [550,  250]), (2, choose (550,  250))]
    in  oneof [getBoxSqr, getBoxRect]

instance Arbitrary Vec2 where
  arbitrary = let p = choose (-200,200) in liftM2 Vec2 p p
       
instance (Arbitrary a) => Arbitrary (WPoint a) where
  arbitrary = let
    s = choose (1, 10)
    in liftM2 WPoint s arbitrary
  
instance Arbitrary (Vector (WPoint Point2D)) where
  arbitrary = let
    s = choose (1, 10)
    vecWall = let
      p1 = elements [-200,200]
      p2 = choose (-200,200)
      in oneof [ liftM2 Vec2 p1 p2, liftM2 Vec2 p2 p1 ]
    wpWall = frequency [ (20, arbitrary), (1, liftM2 WPoint s vecWall) ]
    in Vec.fromList <$> oneof [ listOf arbitrary
                              , listOf wpWall ]


error_precisson = (10e-3)

msgFail text = printTestCase ("\x1b[7m Fail: " ++ show text ++ "! \x1b[0m")


prop_Delaunay::Box Point2D -> SetPoint Point2D -> Property
prop_Delaunay box sp = (length ps) > 4 ==> whenFail (log) fulltest
  where
    fulltest   = testWall .&&. testHull .&&. testSize .&&. testClo
    (wall, st) = runDelaunay2D box sp ixps
    hull       = S.map activeUnit $ externalFaces st
    testw x    = testProperFace sp x
    testh x    = testHullEdge sp x
    testWall   = testIM testw wall
    testHull   = testSet testh hull
    (testClo, cloS2) = testClosure wall hull
    ixps       = let size = Vec.length sp in if size <= 0 then [] else [0 .. size - 1]
    ps         = Vec.toList sp
    testSize   = msgFail "gen obj /= num add" $ IM.size wall == count st
    log = let
      opens = map (testEdge.fst) $ Map.toList $ Map.filter (==Open) $ cloS2
      
      op  = closeUpOnBox box $
            renderBox2D box
         <> renderSetS2Triangle sp wall
         <> renderSetPoint2D sp
         <> renderSetPair sp opens
      
      dia2 = closeUpOnBox box $
             renderBox2D box
          <> renderSetS2Triangle sp wall
          <> renderSetPoint2D sp
          <> renderSetS1 sp hull
          
      dia3 = closeUpOnBox box $
             renderBox2D box
          <> renderSetS2Triangle sp wall
          <> renderSetS2Circle wall
          
      in do
        renderSVG "Delaunay+hull.svg" (sizeSpec (Just 500, Just 500)) dia2
        renderSVG "Delaunay+circle.svg" (sizeSpec (Just 500, Just 500)) dia3
        renderSVG "openEdges.svg" (sizeSpec (Just 500, Just 500)) op

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
    p   = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    pp               = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1               = pointsOnB1 pp
    p2               = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest          = (length p) > 4 && p1 /= [] && p2 /= []
    test             = case makeFirstSimplex plane sP p1 p2 p of
        Just sigma -> testProperFace sP sigma
        _          -> msgFail "non-gen 1st Face2D" False


prop_1stEdge::Box Point2D -> SetPoint Point2D -> Property
prop_1stEdge box sP = pretest ==> test
  where
    (plane, pairBox) = cutBox box []
    p   = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    pp               = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1               = pointsOnB1 pp
    p2               = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest          = (length p) > 4 && p1 /= [] && p2 /= []
    test             = case getFirstEdge plane sP p1 p2 of
        Just (pA, pB) -> let edge = Edge2D pA pB in testHullEdge sP edge
        _             -> msgFail "non-gen 1st Edge2D" False


testProperFace::SetPoint Point2D -> S2 Point2D -> Property
testProperFace sP sigma =  whenFail (log) $ msgFail ("non empty sphere", testAllPoints, sigma, sP!pA, sP!pB, sP!pC) isSphereOK
  where
    ps             = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    (pA,pB,pC)     = face2DPoints sigma
    center         = circleCenter sigma
    radius         = circleRadius sigma
    wC             = WPoint radius center
    cleanP         = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC)) ps
    -- Test if the CircumSphere is empty
    isSphereOK     = null testAllPoints
    testAllPoints  = filter isJust $ map testEmptySph cleanP
    testEmptySph i
      | 0 < powerDist (sP!i) wC = Nothing
      | otherwise = Just (i, powerDist (sP!i) wC)
    log = let
      box = Box2D {xMax2D = 500, xMin2D = -500, yMax2D = 500, yMin2D = -500}
      s   = IM.singleton 1 sigma
      dia = closeUpOnBox box
          $ renderSetPoint2D sP
         <> renderSetS2Triangle sP s
         <> renderSetS2Circle s
      dia2 = renderSetPoint2D sP
         <> renderSetS2Triangle sP s
         <> renderSetS2Circle s
      in do
        renderSVG "ProperFace.svg" (sizeSpec (Just 500, Just 500)) dia
        renderSVG "ProperFace2.svg" (sizeSpec (Just 500, Just 500)) dia2


testHullEdge::SetPoint Point2D -> S1 Point2D -> Property
testHullEdge sP edge = whenFail (log) test
  where
    ps     = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    pA     = edge2DR edge
    pB     = edge2DL edge
    plane  = calcPlane sP edge
    cleanP = filter (\i -> (i /= pA) && (i /= pB)) ps
    test   = case plane of
      Nothing -> msgFail "no plane from face" False
      Just x  -> let
        pp = pointSetPartition (whichSideOfPlane x) sP cleanP
        in case (pointsOnB1 pp, pointsOnB2 pp, pointsOnPlane pp) of
          ([],[],[]) -> msgFail "no points on partition" False
          ([],[],_ ) -> msgFail "all points on plane" False
          ([],_ ,_ ) -> label "face on B1" True
          (_,[] ,_ ) -> label "face on B2" True
          (b1,b2, _) -> let 
            fclean = filter (\p -> (10e-8) > ((sP!.pA &- sP!.p) &. (plane2DNormal x))) 
            in msgFail ("non-Hull face", x, edge) False
               --if null (fclean b1) || null (fclean b2)
               --then label "face on B1" True
               --else msgFail ("non-Hull face", x, edge) False
    log = let
      hull = S.singleton edge
      dia  = renderSetPoint2D sP
         <> renderSetS1 sP hull
      in do
        renderSVG "HullEdge.svg" (sizeSpec (Just 500, Just 500)) dia
        
        
data TestClosure = Open | Closed | Error deriving (Show, Eq)
newtype ClosureID = CloID {testEdge::(Int, Int)} deriving (Show, Eq)
instance Ord ClosureID where
  compare (CloID (a1,a2)) (CloID (b1,b2)) = compEdge a1 a2 b1 b2
  
testClosure::IntMap (S2 Point2D) -> Set (S1 Point2D) -> (Property, Map ClosureID TestClosure)
testClosure ss2 ss1 = (testError, cloS2)
  where
    testError = Map.foldl (\acc x -> acc .&&. ftest x) (property True) cloS2
    ftest x = (msgFail ("Missing face!!! :" ++ (show . Map.filter (==Open) $ cloS2)) $ x /= Open)
         .&&. (msgFail ("More than 2 faces per edge!!!") $ x /= Error)
    
    cloS1 = S.foldl funcS1 Map.empty ss1
    funcS1 acc s1 = let
      a           = edge2DR s1
      b           = edge2DL s1
      func _ Open = Closed
      func _ _    = Error
      in Map.insertWith func (CloID (a,b)) Open acc
      
    cloS2 = IM.foldl funcS2 cloS1 ss2
    funcS2 acc s2 = let
      (a, b, c)   = face2DPoints s2
      func _ Open = Closed
      func _ _    = Error
      add p1 p2 = Map.insertWith func (CloID (p1, p2)) Open
      in add a b $ add b c $ add c a acc