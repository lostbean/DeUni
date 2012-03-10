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

module CheckCommon where

import Test.QuickCheck
import Control.Applicative
import Control.Monad
import qualified Data.IntMap as IM
import qualified Data.Set as S
import qualified Data.List as L
import qualified Data.Vector as Vec
import Data.Vector (Vector)

import Hammer.Math.Vector hiding (Vector)

import DeUni.DeWall
import DeUni.Types
import DeUni.GeometricTools
import DeUni.Dim3.Base3D
import DeUni.Dim2.ReTri2D

import VTKRender

runChecker =  do
  let myArgs = Args {replay = Nothing, maxSuccess = 1000, maxDiscard = 5000, maxSize = 1000, chatty = True}
  print "Testing Projection.."    
  quickCheckWith myArgs prop_Projection
  
  print "Testing Parttition.."    
  quickCheckWith myArgs prop_partition
  
  print "Testing QR decomposition.."    
  quickCheckWith myArgs prop_QR
  

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
  
instance (Arbitrary a) => Arbitrary (Vector a) where
  arbitrary   = liftM2 Vec.generate arbitrary arbitrary


instance Arbitrary Vec2 where
  arbitrary = liftM2 Vec2 p p
    where p = choose (-100,100)
          
instance Arbitrary Mat2 where
  arbitrary = liftM2 Mat2 arbitrary arbitrary

error_precisson = (10e-4)

msgFail text = printTestCase ("\x1b[7m Fail: " ++ show text ++ "! \x1b[0m")
    

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

test_Mat2 a b = and bs
  where
    bs = map (\x -> abs x < error_precisson) [x11,x21,x12,x22]
    (Mat2 (Vec2 x11 x21) (Vec2 x12 x22)) = a &- b
    

prop_QR::Mat2 -> Property
prop_QR a = testQ .&&. testQR
  where
    testR  = msgFail "Bad R"  $ error_precisson > (abs $ det r - x11*x22)
      where (Mat2 (Vec2 x11 x21) (Vec2 x12 x22)) = transpose r.*.r
    testQ  = msgFail "Bad Q"  $ test_Mat2 idmtx (transpose q.*.q)
    testQR = msgFail "Bad QR" $ test_Mat2 a (r.*.q)
    (q,r)  = qrDecomp a


prop_Projection::Point3D -> Point3D -> Property
prop_Projection a b = msgFail "bad projection" $ c &. b < error_precisson
  where
    c = normalofAtoB a b 
    

prop_partition::Box Point3D -> SetPoint Point3D -> Property
prop_partition box sP = msgFail "bad points partition" (test1 || test2)
  where
    (plane, pairBox) = cutBox box []
    p   = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    ppb = pointSetPartition (whichBoxIsIt pairBox) sP p
    ppp = pointSetPartition (whichSideOfPlane plane) sP p
    test1 = ((L.sort $ pointsOnB1 ppb) == (L.sort $ pointsOnB1 ppp))
         && ((L.sort $ pointsOnB2 ppb) == (L.sort $ pointsOnB2 ppp))
    test2 = ((L.sort $ pointsOnB1 ppb) == (L.sort $ pointsOnB2 ppp))
         && ((L.sort $ pointsOnB2 ppb) == (L.sort $ pointsOnB1 ppp))

