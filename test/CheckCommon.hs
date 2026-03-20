{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -Wno-orphans #-}

module CheckCommon where

import qualified Data.IntMap as IM
import qualified Data.List as L
import qualified Data.Set as S
import qualified Data.Vector as Vec

import Data.Vector (Vector)

import Control.Monad
import Test.QuickCheck

import Linear.Mat
import Linear.Vect

import DeUni.DeWall
import DeUni.FirstSeed

-- | Standard testing arguments
myArgs :: Args
myArgs =
    Args
        { replay = Nothing
        , maxSuccess = 1000
        , maxDiscardRatio = 10
        , maxSize = 1000
        , chatty = True
        , maxShrinks = 100
        }

-- Arbitrary Instances

instance Arbitrary (Box Point3D) where
    arbitrary = liftM2 getBox max_ val_min
      where
        getBox mx mn = Box3D{xMax3D = mx, xMin3D = mn, yMax3D = mx, yMin3D = mn, zMax3D = mx, zMin3D = mn}
        val_min = choose (-200, -100)
        max_ = choose (200, 100)

instance Arbitrary (Box Point2D) where
    arbitrary =
        let
            getBoxSqr = liftM2 getBox max_ mn_
            getBoxRect = do
                max1 <- max_
                max2 <- max_
                min1 <- mn_
                min2 <- mn_
                return $ Box2D{xMax2D = max1, xMin2D = min1, yMax2D = max2, yMin2D = min2}
            getBox mx mn = Box2D{xMax2D = mx, xMin2D = mn, yMax2D = mx, yMin2D = mn}
            mn_ = frequency [(1, elements [-550, -250]), (2, choose (-550, -250))]
            max_ = frequency [(1, elements [550, 250]), (2, choose (550, 250))]
         in
            oneof [getBoxSqr, getBoxRect]

instance Arbitrary Plane3D where
    arbitrary = do
        nd <- arbitrary
        dst <- choose (-200, 200)
        return $ Plane3D nd dst

instance Arbitrary Vec3D where
    arbitrary = liftM3 Vec3 p p p
      where
        p = choose (-100, 100)

instance Arbitrary Vec2D where
    arbitrary = liftM2 Vec2 p p
      where
        p = choose (-100, 100)

instance (Arbitrary (p Double)) => Arbitrary (WPoint p) where
    arbitrary = do
        w <- choose (0, 10)
        p <- arbitrary
        return $ WPoint w p

instance Arbitrary (Vector (WPoint Point3D)) where
    arbitrary = genPointVector

instance Arbitrary (Vector (WPoint Point2D)) where
    arbitrary = genPointVector

genPointVector :: (PointND p, Arbitrary (p Double), Norm Double p) => Gen (Vector (WPoint p))
genPointVector = do
    let
        wpNormal = arbitrary
        wpWall = do
            s <- choose (1, 10)
            -- This is simplified for the polymorphic version, but we can specialize if needed
            WPoint s <$> arbitrary

        collinear = do
            p1 <- arbitrary
            p2 <- arbitrary
            -- Using a simplified collinear generator that works for both 2D and 3D
            let diff = p2 &- p1
            if normsqr diff < 1e-3
                then collinear
                else do
                    let dir = normalize diff
                    ts <- replicateM 5 (choose (-100, 100))
                    let pts = WPoint 0 p1 : [WPoint 0 (p1 &+ (t *& dir)) | t <- ts]
                    return $ Vec.fromList pts

        coincident = do
            p <- arbitrary
            return $ Vec.fromList (replicate 6 (WPoint 0 p))

    frequency
        [ (70, Vec.fromList <$> (choose (5, 20) >>= \n -> replicateM n wpNormal))
        , (10, Vec.fromList <$> (choose (5, 20) >>= \n -> replicateM n wpWall))
        , (10, collinear)
        , (10, coincident)
        ]

instance Arbitrary (Mat2 Double) where
    arbitrary = liftM2 Mat2 arbitrary arbitrary

-- Helper functions

{- | Check if a point set has enough distinct points to be non-degenerate.
Returns True if there are at least n distinct points (by position).
-}
hasDistinctPoints :: (Eq (p Double), PointND p) => Int -> Vector (WPoint p) -> Bool
hasDistinctPoints n sp = length (L.nub positions) >= n
  where
    positions = map point (Vec.toList sp)

error_precisson :: Double
error_precisson = 10e-4

msgFail :: (Testable prop, Show a) => a -> prop -> Property
msgFail text = counterexample ("\x1b[7m Fail: " ++ show text ++ "! \x1b[0m")

testIM :: (a -> Property) -> IM.IntMap a -> Property
testIM test mp
    | IM.null mp = label "degenerate input (empty output)" True
    | otherwise =
        let (x, xs) = IM.deleteFindMin mp
         in IM.foldr (\a acc -> acc .&&. test a) (test $ snd x) xs

testSet :: (a -> Property) -> S.Set a -> Property
testSet test st
    | S.null st = label "degenerate input (empty output)" True
    | otherwise =
        let (x, xs) = S.deleteFindMin st
         in S.foldr (\a b -> b .&&. test a) (test x) xs

-- Properties

prop_Projection :: Vec3D -> Vec3D -> Property
prop_Projection a b = msgFail "bad projection" $ c &. b < error_precisson
  where
    c = normalofAtoB a b

prop_partition :: Box Point3D -> SetPoint Point3D -> Property
prop_partition box sP = msgFail "bad points partition" (test1 || test2)
  where
    (plane, pairBox) = cutBox box []
    p = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    ppb = pointSetPartition (whichBoxIsIt pairBox) sP p
    ppp = pointSetPartition (whichSideOfPlane plane) sP p
    test1 =
        ((L.sort $ pointsOnB1 ppb) == (L.sort $ pointsOnB1 ppp))
            && ((L.sort $ pointsOnB2 ppb) == (L.sort $ pointsOnB2 ppp))
    test2 =
        ((L.sort $ pointsOnB1 ppb) == (L.sort $ pointsOnB2 ppp))
            && ((L.sort $ pointsOnB2 ppb) == (L.sort $ pointsOnB1 ppp))

prop_projAonB_idempotent :: Vec3D -> Vec3D -> Property
prop_projAonB_idempotent a b =
    (normsqr b > error_precisson) ==>
        let p = projAonB a b
            pp = projAonB p b
         in msgFail "projAonB is not idempotent" $ norm (p &- pp) < error_precisson

prop_projAonB_orthogonal :: Vec3D -> Vec3D -> Property
prop_projAonB_orthogonal a b =
    (normsqr b > error_precisson) ==>
        let p = projAonB a b
            orth = a &- p
         in msgFail "projAonB is not orthogonal" $ abs (orth &. b) < error_precisson

prop_powerDist_symmetric :: WPoint Point3D -> WPoint Point3D -> Property
prop_powerDist_symmetric a b =
    msgFail "powerDist is not symmetric" $ abs (powerDist a b - powerDist b a) < error_precisson

prop_powerDist_euclidean :: Vec3D -> Vec3D -> Property
prop_powerDist_euclidean a b =
    let wa = WPoint 0 a
        wb = WPoint 0 b
        pd = powerDist wa wb
        ed = normsqr (a &- b)
     in msgFail "powerDist != distSqr when weights are zero" $ abs (pd - ed) < error_precisson

prop_whichSideOfPlane_consistency :: Plane Point3D -> Vec3D -> Property
prop_whichSideOfPlane_consistency plane p =
    let side = whichSideOfPlane plane p
        nd = normalize (planeNormal plane)
        dist = p &. nd - planeDist plane
     in case side of
            B1 -> msgFail "B1 but dist < 0" $ dist >= -error_precisson
            B2 -> msgFail "B2 but dist > 0" $ dist <= error_precisson
            OnPlane -> msgFail "OnPlane but dist /= 0" $ abs dist < error_precisson
            _ -> property True

prop_getMaxDistPoint :: Plane Point3D -> SetPoint Point3D -> Property
prop_getMaxDistPoint plane sP =
    (Vec.length sP > 0) ==>
        let p = [0 .. Vec.length sP - 1]
         in case getMaxDistPoint plane sP p of
                Nothing -> property False
                Just (_, dMax) ->
                    let dists = map (norm . getProjOnPlane plane . (sP !.)) p
                     in msgFail "Not the max dist point" $ dMax >= maximum dists - error_precisson

prop_getMaxDistPointOnDir :: Vec3D -> SetPoint Point3D -> Property
prop_getMaxDistPointOnDir dir sP =
    (Vec.length sP > 0 && normsqr dir > error_precisson) ==>
        let p = [0 .. Vec.length sP - 1]
            refdir = normalize dir
         in case getMaxDistPointOnDir refdir sP p of
                Nothing -> property False
                Just (_, dMax) ->
                    let dists = map (\x -> norm $ projAonB (sP !. x) refdir) p
                     in msgFail "Not the max dist point on dir" $ dMax >= maximum dists - error_precisson

prop_compEdge_consistent :: Int -> Int -> Int -> Int -> Property
prop_compEdge_consistent a b c d =
    let res = compEdge a b c d
        set1 = S.fromList [a, b]
        set2 = S.fromList [c, d]
        amax = max a b
        amin = min a b
        bmax = max c d
        bmin = min c d
     in case res of
            EQ -> msgFail "EQ but sets differ" $ set1 == set2
            LT -> msgFail "LT but set1 >= set2" $ set1 /= set2 && (amax < bmax || (amax == bmax && amin < bmin))
            GT -> msgFail "GT but set1 <= set2" $ set1 /= set2 && (amax > bmax || (amax == bmax && amin > bmin))

prop_compFace_consistent :: (Int, Int, Int) -> (Int, Int, Int) -> Property
prop_compFace_consistent (a1, a2, a3) (b1, b2, b3) =
    let res = compFace (a1, a2, a3) (b1, b2, b3)
        set1 = S.fromList [a1, a2, a3]
        set2 = S.fromList [b1, b2, b3]
        sort3 (x, y, z) = case L.sort [x, y, z] of
            [i, j, k] -> (k, j, i) -- fast3DSort is descending
            _ -> error "sort3: impossible"
     in case res of
            EQ -> msgFail "EQ but sets differ" $ set1 == set2
            LT -> msgFail "LT but a >= b" $ sort3 (a1, a2, a3) < sort3 (b1, b2, b3)
            GT -> msgFail "GT but a <= b" $ sort3 (a1, a2, a3) > sort3 (b1, b2, b3)

prop_findMinimunButZero :: [Double] -> Property
prop_findMinimunButZero ds =
    let ps = [0 .. length ds - 1]
        func i = ds !! i
        res = findMinimunButZero func ps
        nonZeros = filter (\x -> x /= 0) ds
     in case (res, null nonZeros) of
            (Nothing, True) -> property True
            (Just (val, _), False) -> msgFail "Not the minimum non-zero" $ val == minimum nonZeros
            _ -> property False
