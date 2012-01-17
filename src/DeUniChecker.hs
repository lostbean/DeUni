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

module DeUniChecker where

import Test.QuickCheck
import Control.Applicative
import Control.Monad
import Data.Array.Diff
import qualified Data.IntMap as IM
import qualified Data.Set as S
import qualified Data.List as L

import DeUni.DeWall
import Math.Vector
import DeUni.Types
import DeUni.GeometricTools
import DeUni.Dim3.Base3D
import DeUni.Dim3.Delaunay3D
import DeUni.Dim3.Hull3D

import VTKGenSimplex

instance Arbitrary (Box Point3D) where
  arbitrary = liftM2 getBox max min
    where 
      getBox max min = Box3D {xMax3D=max, xMin3D=min, yMax3D=max, yMin3D=min, zMax3D=max, zMin3D=min}
      min = choose (-200,-100)
      max = choose (200,100)

instance Arbitrary Vec3 where
  arbitrary = liftM3 Vec3 p p p
    where p = choose (-100,100)

instance (Arbitrary a) => Arbitrary (WPoint a) where
  arbitrary = WPoint 0 <$> arbitrary
  
instance (Ix a, Integral a, Arbitrary b) => Arbitrary (DiffArray a b) where
  arbitrary   =
    (\x -> listArray (0,fromIntegral (length x - 1)) x) <$> arbitrary 


error_precisson = (10e-8)

msgFail text = printTestCase ("\x1b[7m Fail: " ++ show text ++ "! \x1b[0m")

prop_ConvHull::Box Point3D -> SetPoint Point3D -> Property
prop_ConvHull box sp = (length ps) > 4 ==> whenFail (writeVTKfile "Hull3D_err.vtu" sp hull) fulltest
  where
    fulltest    = testHull .&&. testClosure .&&. testSize
    (hull, st)  = runHull3D box sp ixps
    testh x     = testHullFace sp x
    testHull    = testIM testh hull
    testSize    = msgFail "gen obj /= num add" $ IM.size hull == count st
    testClosure = msgFail "open Hull" $ (S.null.externalFaces) st
    ixps        = indices sp
    ps          = elems sp


prop_Delaunay::Box Point3D -> SetPoint Point3D -> Property
prop_Delaunay box sp = (length ps) > 4 ==> whenFail (writeVTKfile "Delaunay3D_err.vtu" sp wall) fulltest
  where
    fulltest   = testWall .&&. testHull .&&. testSize
    (wall, st) = runDelaunay3D box sp ixps
    testw x    = testProperTetrahedron sp x
    testh x    = testHullFace sp x
    testWall   = testIM testw wall
    testHull   = testSet testh (S.map activeUnit $ externalFaces st)
    ixps       = indices sp
    ps         = elems sp
    testSize   = msgFail "gen obj /= num add" $ IM.size wall == count st
    

testIM test map
  | IM.null map = err
  | otherwise  = let (x, xs) = IM.deleteFindMin map
                 in IM.fold (\a b -> b .&&. test a) (test x) xs
  where
    err          = msgFail "empty output" False

testSet test set
  | S.null set = err
  | otherwise  = let (x, xs) = S.deleteFindMin set
                 in S.fold (\a b -> b .&&. test a) (test x) xs
  where
    err          = msgFail "empty output" False


prop_Projection::Point3D -> Point3D -> Property
prop_Projection a b = msgFail "bad projection" $ c &. b < error_precisson
  where
    c = normalofAtoB a b 
    

prop_partition::Box Point3D -> SetPoint Point3D -> Property
prop_partition box sP = msgFail "bad points partition" (test1 || test2)
  where
    (plane, pairBox) = cutBox box []
    p   = indices sP 
    ppb = pointSetPartition (whichBoxIsIt pairBox) sP p
    ppp = pointSetPartition (whichSideOfPlane plane) sP p
    test1 = ((L.sort $ pointsOnB1 ppb) == (L.sort $ pointsOnB1 ppp))
         && ((L.sort $ pointsOnB2 ppb) == (L.sort $ pointsOnB2 ppp))
    test2 = ((L.sort $ pointsOnB1 ppb) == (L.sort $ pointsOnB2 ppp))
         && ((L.sort $ pointsOnB2 ppb) == (L.sort $ pointsOnB1 ppp))


prop_CircumSphere::(Vec3,Vec3,Vec3,Vec3) -> Property
prop_CircumSphere (a,b,c,d) = msgFail "CircumCenter" test
  where
    test   = and $ map testcenter [a,b,c,d]
    testcenter i = error_precisson > (abs $ (norm $ center &- i) - radius)
    (radius, center) = getCircumSphere (a, b ,c) d


prop_1stSimplex::Box Point3D -> SetPoint Point3D -> Property
prop_1stSimplex box sP = pretest ==> test
  where
    (plane, pairBox) = cutBox box []
    p  = indices sP 
    pp = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1 = pointsOnB1 pp
    p2 = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest = (length p) > 4 && p1 /= [] && p2 /= []
    test = case makeFirstSimplex plane sP p1 p2 p of
        Just sigma -> testProperTetrahedron sP sigma
        _          -> msgFail "non-gen 1st Tetra3D" False


prop_1stFace::Box Point3D -> SetPoint Point3D -> Property
prop_1stFace box sP = pretest ==> test
  where
    (plane, pairBox) = cutBox box []
    p  = indices sP
    pp = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1 = pointsOnB1 pp
    p2 = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest = (length p) > 4 && p1 /= [] && p2 /= []
    test = case makeFirstFace plane sP p1 p2 p of
        Just face -> testHullFace sP face
        _         -> msgFail "non-gen 1st Face3D" False


testProperTetrahedron::SetPoint Point3D -> S2 Point3D -> Property
testProperTetrahedron sP sigma = msgFail ("bad tetrahedron center ", map testC [pA,pB,pC,pD]) isCenterOK
                            .&&. msgFail ("non empty sphere", [pA,pB,pC,pD]) isSphereOK
  where
    (pA,pB,pC,pD)  = tetraPoints sigma
    center         = circumSphereCenter sigma
    radius         = norm $ sP!.pA &- center
    cleanP         = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC) && (i /= pD)) (indices sP)
    -- | Test if the it is the center of the simplex
    isCenterOK     = and $ map testCenter [pA,pB,pC,pD]
    testCenter i   = error_precisson > (abs $ (norm $ center &- sP!.i) - radius)
    testC i        = (abs $ (norm $ center &- sP!.i) - radius)
    -- | Test if the CircumSphere is empty
    isSphereOK     = and $ map testEmptySph cleanP
    testEmptySph i = radius < norm (sP!.i &- center)



testHullFace::SetPoint Point3D -> S1 Point3D -> Property
testHullFace sP face = test
  where
    (pA,pB,pC) = face3DPoints face
    pp     = (\x -> pointSetPartition (whichSideOfPlane x) sP cleanP) <$> plane
    plane  = calcPlane sP face
    cleanP = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC)) $ indices sP
    test   = case pp of
      Nothing -> msgFail "no plane from face" False
      Just x  -> case (pointsOnB1 x, pointsOnB2 x, pointsOnPlane x) of
        ([],[],[]) -> msgFail "no points on partition" False
        ([],[],_)  -> msgFail "all points on plane" False
        ([],_,_)   -> label "face on B1" True
        (_,[],_)   -> label "face on B2" True
        _          -> msgFail ("non-Hull face", x, face3DPoints face) False