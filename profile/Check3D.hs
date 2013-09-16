{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}

module Check3D where

import qualified Data.IntMap as IM
import qualified Data.Map    as Map
import qualified Data.Set    as S
import qualified Data.List   as L
import qualified Data.Vector as Vec

import Test.QuickCheck
import Control.Applicative
import Control.Monad

import Data.Maybe  (isJust)
import Data.IntMap (IntMap)
import Data.Set    (Set)
import Data.Vector (Vector, (!))

import Hammer.Math.Algebra

import DeUni.DeWall
import DeUni.Types
import DeUni.FirstSeed
import DeUni.GeometricTools
import DeUni.Dim3.Base3D
import DeUni.Dim3.Delaunay3D
import DeUni.Dim3.Hull3D
import DeUni.Dim3.ReTri3D

import VTKRender

runChecker =  do
  let myArgs = Args { replay = Nothing
                    , maxSuccess = 1000
                    , maxDiscardRatio = 1
                    , maxSize = 1000
                    , chatty = True }
  
  print "Testing 1st face.."    
  quickCheckWith myArgs prop_1stFace
  
  print "Testing Circumsphere.."    
  quickCheckWith myArgs prop_CircumSphere
  
  print "Testing 1st tetrahedron.."    
  quickCheckWith myArgs prop_1stSimplex

  print "Testing Convex Hull.."    
  quickCheckWith myArgs prop_ConvHull
  
  print "Testing Delaunay.."    
  quickCheckWith myArgs prop_Delaunay


instance Arbitrary (Box Point3D) where
  arbitrary = let      
    getBoxSqr = liftM2 getBox max min
    getBoxRect = do
      max1 <- max
      max2 <- max
      max3 <- max
      min1 <- min
      min2 <- min
      min3 <- min
      return $ Box3D {xMax3D=max1, xMin3D=min1, yMax3D=max2, yMin3D=min2, zMax3D=max3, zMin3D=min3}
    getBox max min = Box3D {xMax3D=max, xMin3D=min, yMax3D=max, yMin3D=min, zMax3D=max, zMin3D=min}
    min = frequency [(1, elements [-550,-250]), (2, choose (-550,-250))]
    max = frequency [(1, elements [550,  250]), (2, choose (550,  250))]
    in  oneof [getBoxSqr, getBoxRect]
        
instance Arbitrary Vec3 where
  arbitrary = let p = choose (-200,200) in liftM3 Vec3 p p p
       
instance (Arbitrary a) => Arbitrary (WPoint a) where
  arbitrary = let
    s = choose (1, 10)
    in liftM2 WPoint s arbitrary
  
instance Arbitrary (Vector (WPoint Point3D)) where
  arbitrary = let
    s = choose (1, 10)
    vecWall = let
      p1 = elements [-200,200]
      p2 = choose (-200,200)
      in oneof [ liftM3 Vec3 p1 p2 p2
               , liftM3 Vec3 p2 p1 p2 
               , liftM3 Vec3 p2 p2 p1 ]
    wpWall = frequency [ (20, arbitrary), (1, liftM2 WPoint s vecWall) ]
    in Vec.fromList <$> oneof [ listOf arbitrary
                              , listOf wpWall ]


error_precisson = (10e-2)

msgFail text  = printTestCase ("\x1b[7m Fail: " ++ show text ++ "! \x1b[0m")

plotFail file sp obj = let
  ps = Vec.convert $ Vec.map point sp
  in whenFail (writeVTKfile file ps obj)

prop_ConvHull::Box Point3D -> SetPoint Point3D -> Property
prop_ConvHull box sp = (length ps) > 4 ==> plotFail "Hull3D_err.vtu" sp hull fulltest
  where
    fulltest    = testHull .&&. testClosure .&&. testSize
    (hull, st)  = runHull3D box sp ixps
    testh x     = testHullFace sp x
    testHull    = testIM testh hull
    testSize    = msgFail "gen obj /= num add" $ IM.size hull == count st
    testClosure = msgFail "open Hull" $ (S.null.externalFaces) st
    ixps        = let size = Vec.length sp in if size <= 0 then [] else [0 .. size - 1]
    ps          = Vec.toList sp


prop_Delaunay::Box Point3D -> SetPoint Point3D -> Property
prop_Delaunay box sp = (length ps) > 4 ==> plotFail "Delaunay3D_err.vtu" sp wall fulltest
  where
    fulltest   = testWall .&&. testHull .&&. testSize .&&. testClo
    (wall, st) = runDelaunay3D box sp ixps
    testw x    = testProperTetrahedron sp x
    testh x    = testHullFace sp x
    testWall   = testIM testw wall
    testHull   = testSet testh (S.map activeUnit $ externalFaces st)
    testClo    = testClosure wall (S.map activeUnit $ externalFaces st)
    ixps       = let size = Vec.length sp in if size <= 0 then [] else [0 .. size - 1]
    ps         = Vec.toList sp
    testSize   = msgFail "gen obj /= num add" $ IM.size wall == count st
    

testIM test map
  | IM.null map = err
  | otherwise  = let (x, xs) = IM.deleteFindMin map
                 in IM.fold (\a b -> b .&&. test a) (test $ snd x) xs
  where
    err          = msgFail "empty output" False

testSet test set
  | S.null set = err
  | otherwise  = let (x, xs) = S.deleteFindMin set
                 in S.fold (\a b -> b .&&. test a) (test x) xs
  where
    err          = msgFail "empty output" False


prop_CircumSphere::WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> Property
prop_CircumSphere a b c d = msgFail "bad center" test
  where
    test         = and $ map testcenter [a,b,c,d]
    testcenter i = error_precisson > (abs $ powerDist i (WPoint r center))
    (r, center)  = getCircumSphere a b c d


prop_1stSimplex::Box Point3D -> SetPoint Point3D -> Property
prop_1stSimplex box sP = pretest ==> test
  where
    (plane, pairBox) = cutBox box []
    p  = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1] 
    pp = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1 = pointsOnB1 pp
    p2 = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest = (length p) > 4 && p1 /= [] && p2 /= []
    test = case makeFirstSimplex plane sP p1 p2 p of
        Just sigma -> let
          wall = IM.singleton 1 sigma
          test = testProperTetrahedron sP sigma
          in plotFail "1stSimplex_err.vtu" sP wall test
        _          -> msgFail "non-gen 1st Tetra3D" False


prop_1stFace::Box Point3D -> SetPoint Point3D -> Property
prop_1stFace box sP = pretest ==> test
  where
    (plane, pairBox) = cutBox box []
    p  = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    pp = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1 = pointsOnB1 pp
    p2 = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest = (length p) > 4 && p1 /= [] && p2 /= []
    test = case makeFirstFace plane sP p1 p2 p of
        Just face -> testHullFace sP face
        _         -> msgFail "non-gen 1st Face3D" False


testProperTetrahedron::SetPoint Point3D -> S2 Point3D -> Property
testProperTetrahedron sP sigma = msgFail ("non empty sphere", testAllPoints, sigma) isSphereOK
  where
    ps             = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    (pA,pB,pC,pD)  = tetraPoints sigma
    center         = circumSphereCenter sigma
    radius         = circumRadius sigma
    wC             = WPoint radius center
    cleanP         = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC) && (i /= pD)) ps
    -- Test if the CircumSphere is empty
    isSphereOK     = null testAllPoints
    testAllPoints  = filter isJust $ map testEmptySph cleanP
    testEmptySph i
      | 0 < powerDist (sP!i) wC = Nothing
      | otherwise = Just (i, powerDist (sP!i) wC)

testHullFace::SetPoint Point3D -> S1 Point3D -> Property
testHullFace sP face = case calcPlane sP face of
  Nothing -> msgFail "no plane from face" False
  Just pl -> let
    ps     = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    (pA,pB,pC) = face3DPoints face
    cleanP = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC)) ps
    pp = pointSetPartition (whichSideOfPlane pl) sP cleanP
    
    testSigma = let
      newND
        | (null.pointsOnB1) pp = nd
        | otherwise            = neg nd
      nd = plane3DNormal pl                      
      actFace = ActiveUnit { activeUnit = face, assocP = undefined, assocND = newND }
      in case makeSimplex actFace sP ps of
        Just sigma -> testProperTetrahedron sP sigma
        _          -> msgFail "No possible tetraredron for hull face!" False
         
    testHull = case (pointsOnB1 pp, pointsOnB2 pp, pointsOnPlane pp) of
      ([],[],[]) -> msgFail "no points on partition" False
      ([],[],_)  -> msgFail "all points on plane" False
      ([],_,_)   -> label "face on B1" True
      (_,[],_)   -> label "face on B2" True
      _          -> msgFail ("non-Hull face", pp, face3DPoints face, map (\i -> (i, planeNormal pl &. (sP!.i &- sP!.pA))) ps) False
      
    in testHull .&&. testSigma        
        
data TestClosure = Open | Closed | Error deriving (Show, Eq)
newtype ClosureID = CloID (Int, Int, Int) deriving (Show, Eq)
instance Ord ClosureID where
  compare (CloID a) (CloID b) = compFace a b
  
testClosure::IntMap (S2 Point3D) -> Set (S1 Point3D) -> Property
testClosure ss2 ss1 = testError
  where
    testError = Map.foldl (\acc x -> acc .&&. ftest x) (property True) cloS2
    ftest x = (msgFail ("Missing face!!!: " ++ (show . Map.filter (==Open) $ cloS2)) $ x /= Open)
         .&&. (msgFail ("More than 2 faces per edge!!!") $ x /= Error)
    
    cloS1 = S.foldl funcS1 Map.empty ss1
    funcS1 acc s1 = let
      face        = CloID $ face3DPoints s1
      func _ Open = Closed
      func _ _    = Error
      in Map.insertWith func face Open acc
      
    cloS2 = IM.foldl funcS2 cloS1 ss2
    funcS2 acc s2  = let
      (a, b, c, d) = tetraPoints s2
      func _ Open  = Closed
      func _ _     = Error
      add p1 p2 p3 = Map.insertWith func (CloID (p1,p2,p3)) Open
      in add a b c $ add a b d $ add b c d $ add c a d acc
         
         
