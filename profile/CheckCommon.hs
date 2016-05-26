{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}

module CheckCommon where

import qualified Data.IntMap as IM
import qualified Data.Set    as S
import qualified Data.List   as L
import qualified Data.Vector as Vec

import Data.Vector (Vector)
  
import Test.QuickCheck
import Control.Applicative
import Control.Monad

import Linear.Vect

import DeUni.DeWall
import DeUni.Types
import DeUni.GeometricTools
import DeUni.Dim3.Base3D
import DeUni.Dim2.ReTri2D

import VTKRender

runChecker =  do
  let myArgs = Args { replay = Nothing
                    , maxSuccess = 1000
                    , maxDiscardRatio = 5
                    , maxSize = 1000
                    , chatty = True }

  print "Testing Projection.."    
  quickCheckWith myArgs prop_Projection
  
  print "Testing Parttition.."    
  quickCheckWith myArgs prop_partition
  

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
  arbitrary =Vec.fromList <$> arbitrary

instance Arbitrary Vec2 where
  arbitrary = liftM2 Vec2 p p
    where p = choose (-100,100)
          
instance Arbitrary Mat2 where
  arbitrary = liftM2 Mat2 arbitrary arbitrary

error_precisson = (10e-4)

msgFail text = printTestCase ("\x1b[7m Fail: " ++ show text ++ "! \x1b[0m")
    
testIM :: (a -> Gen Prop) -> IM.IntMap a -> Gen Prop
testIM test map
  | IM.null map = err
  | otherwise   = let (x, xs) = IM.deleteFindMin map
                  in IM.fold (\a acc -> acc .&&. test a) (test $ snd x) xs
  where
    err = msgFail "empty output" False

testSet :: (a -> Gen Prop) -> S.Set a -> Gen Prop
testSet test set
  | S.null set = err
  | otherwise  = let (x, xs) = S.deleteFindMin set
                 in S.fold (\a b -> b .&&. test a) (test x) xs
  where
    err          = msgFail "empty output" False

prop_Projection::Point3D -> Point3D -> Property
prop_Projection a b = msgFail "bad projection" $ c &. b < error_precisson
  where c = normalofAtoB a b 
    
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

