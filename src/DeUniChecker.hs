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

module DeUniChecker where

import Test.QuickCheck.Arbitrary
import Test.QuickCheck.Property
import Test.QuickCheck.Gen
import Data.Vec hiding (map, length)
import Control.Monad
import Math.DeUni



import Debug.Trace
debug :: Show a => String -> a -> a
debug s x = trace (s ++ show x) x


instance Arbitrary Box where
  arbitrary = liftM2 getBox max min
    where getBox max min = Box {xMax=max, xMin=min, yMax=max, yMin=min, zMax=max, zMin=min}
          min = choose (-200,-100)
          max = choose (200,100)

instance Arbitrary Vec3D where
  arbitrary = liftM3 Vec3D p p p
    where p = choose (-100,100)

error_precisson = (10e-8)

prop_ConvHull::Box -> [Vec3D] -> Property
prop_ConvHull box ps = (length ps) > 4 ==> fullTest
    where
    fullTest     = (not $ null ps) && (and $ map (testFace ps) (runDeHull box ps))


prop_Delaunay::Box -> [Vec3D] -> Property
prop_Delaunay box ps = (length ps) > 4 ==> fullTest && sizeTest
    where
    fullTest      = and $ map (testProperTetrahedron ps) wall
    sizeTest      = length wall == genD
    (wall, genD)  = runDeWall box ps


prop_CircumSphere::(Vec3D,Vec3D,Vec3D,Vec3D) -> Bool
prop_CircumSphere (a,b,c,d) = and $ map test [a,b,c,d]
    where
    test i = error_precisson > (abs $ (norm $ center - i) - radius)
    (radius, center) = getCircumSphere (a, b ,c) d


prop_1stSimplex::Box -> [Vec3D] -> Property
prop_1stSimplex box p = (length p) > 4 && p1 /= [] && p2 /= [] ==> test
    where
    (plane, pairBox) = genPlane box
    pp = pointSetPartition (whichBoxIsIt pairBox) p
    p1 = pointsOnB1 pp
    p2 = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    test = case makeFirstSimplex plane p1 p2 p of
        Just sigma -> testProperTetrahedron p sigma
        _          -> False


prop_1stFace::Box -> [Vec3D] -> Property
prop_1stFace box p = (length p) > 4 && p1 /= [] && p2 /= [] ==> test
    where
    (plane, pairBox) = genPlane box
    pp = pointSetPartition (whichBoxIsIt pairBox) p
    p1 = pointsOnB1 pp
    p2 = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    test = case makeFirstFace plane p1 p2 p of
        Just face -> testFace p face
        _         -> False

testProperTetrahedron::[Vec3D] -> Simplex -> Bool
testProperTetrahedron ps sigma = isCenterOK && isSphereOK
    where
    (pA,pB,pC,pD)  = setCellID sigma
    center         = circumSphereCenter sigma
    radius         = norm $ pA - center
    cleanP         = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC) && (i /= pD)) ps
    -- | Test if the it is the center of the simplex
    isCenterOK     = and $ map testCenter [pA,pB,pC,pD]
    testCenter i   = error_precisson > (abs $ (norm $ center - i) - radius)
    -- | Test if the CircumSphere is empty
    isSphereOK     = and $ map testEmptySph cleanP
    testEmptySph i = radius < norm (i - center)


testFace::[Vec3D] -> Face -> Bool
testFace ps face = isNDOK && isHullOK
    where
    (pA,pB,pC)   = facePoints face
    nd           = refND face
    cleanP       = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC)) ps
    -- | Test Normal vector. The inner product must to be zero between ND and the vector in the face
    isNDOK       = and $ map testND [(pA,pB),(pB,pC),(pC,pA)]
    testND (a,b) = error_precisson > (abs $ dot nd (a - b))
    -- | Test if the face is a hull. All others points must lie on the opposite side defined by the ND
    isHullOK     = and $ map testHull cleanP
    testHull i   = (dot nd (i - pA)) < 0

