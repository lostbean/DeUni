{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE UndecidableInstances #-}

module Check3D where

import qualified Data.IntMap as IM
import qualified Data.Map as Map
import qualified Data.Set as S
import qualified Data.Vector as Vec
import qualified Data.Vector.Unboxed as VU

import Test.QuickCheck

import Data.IntMap (IntMap)
import Data.Maybe (isJust)
import Data.Set (Set)
import Data.Vector ((!))

import Linear.Vect

import DeUni.DeWall
import DeUni.Dim3.Delaunay3D
import DeUni.Dim3.Hull3D
import DeUni.Dim3.ReTri3D

import CheckCommon
import VTKRender (writeVTKfile)

runChecker :: IO ()
runChecker = do
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

    print "Testing Delaunay with 4 points.."
    quickCheckWith myArgs prop_DelaunayMinPoints

plotFail :: (Testable prop) => FilePath -> SetPoint Point3D -> a -> prop -> Property
plotFail file sp obj test =
    let
        pts :: VU.Vector (Vec3 Double)
        pts = Vec.convert $ Vec.map point sp
     in
        whenFail (writeVTKfile file pts obj) test

prop_ConvHull :: Box Point3D -> SetPoint Point3D -> Property
prop_ConvHull box sp = (length ps) > 4 ==> plotFail "Hull3D_err.vtu" sp hull fulltest
  where
    fulltest = testHull .&&. testClo .&&. testSize
    (hull, st) = runHull3D box sp ixps
    testh x = testHullFace sp x
    testHull = testIM testh hull
    testSize = msgFail "gen obj /= num add" $ IM.size hull == count st
    testCloVal = msgFail "open Hull" $ (S.null . externalFaces) st
    testClo = testCloVal
    ixps = let size = Vec.length sp in if size <= 0 then [] else [0 .. size - 1]
    ps = Vec.toList sp

prop_Delaunay :: Box Point3D -> SetPoint Point3D -> Property
prop_Delaunay box sp = (length ps) > 4 ==> plotFail "Delaunay3D_err.vtu" sp wall fulltest
  where
    fulltest = testWall .&&. testHull .&&. testSize .&&. testClo
    (wall, st) = runDelaunay3D box sp ixps
    testw x = testProperTetrahedron sp x
    testh x = testHullFace sp x
    testWall = testIM testw wall
    testHull = testSet testh (S.map activeUnit $ externalFaces st)
    testClo = testClosure wall (S.map activeUnit $ externalFaces st)
    ixps = let size = Vec.length sp in if size <= 0 then [] else [0 .. size - 1]
    ps = Vec.toList sp
    testSize = msgFail "gen obj /= num add" $ IM.size wall == count st

prop_DelaunayMinPoints :: Box Point3D -> (WPoint Point3D, WPoint Point3D, WPoint Point3D, WPoint Point3D) -> Property
prop_DelaunayMinPoints box (p1, p2, p3, p4) =
    let sp = Vec.fromList [p1, p2, p3, p4]
     in (p1 /= p2 && p2 /= p3 && p3 /= p4 && p1 /= p4) ==>
            let (wall, _) = runDelaunay3D box sp [0, 1, 2, 3]
             in (IM.size wall == 1) .||. (IM.size wall == 0)

prop_CircumSphere :: WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> WPoint Point3D -> Property
prop_CircumSphere a b c d = msgFail "bad center" test
  where
    test = and $ map testcenter [a, b, c, d]
    testcenter i = error_precisson > (abs $ powerDist i (WPoint r center))
    (r, center) = getCircumSphere a b c d

prop_1stSimplex :: Box Point3D -> SetPoint Point3D -> Property
prop_1stSimplex box sP = pretest ==> test
  where
    (plane, pairBox) = cutBox box []
    p = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    pp = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1 = pointsOnB1 pp
    p2 = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest = (length p) > 4 && p1 /= [] && p2 /= []
    test = case makeFirstSimplex plane sP p1 p2 p of
        Just sigma -> testProperTetrahedron sP sigma
        _ -> msgFail "non-gen 1st Tetra3D" False

prop_1stFace :: Box Point3D -> SetPoint Point3D -> Property
prop_1stFace box sP = pretest ==> test
  where
    (plane, pairBox) = cutBox box []
    p = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    pp = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1 = pointsOnB1 pp
    p2 = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    pretest = (length p) > 4 && p1 /= [] && p2 /= []
    test = case makeFirstFace plane sP p1 p2 p of
        Just face -> testHullFace sP face
        _ -> msgFail "non-gen 1st Face3D" False

testProperTetrahedron :: SetPoint Point3D -> S2 Point3D -> Property
testProperTetrahedron sP sigma = msgFail ("non empty sphere", testAllPoints, sigma) isSphereOK
  where
    ps = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
    (pA, pB, pC, pD) = tetraPoints sigma
    center = circumSphereCenter sigma
    rad = circumRadius sigma
    wC = WPoint rad center
    cleanP = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC) && (i /= pD)) ps
    -- Test if the CircumSphere is empty
    isSphereOK = null testAllPoints
    testAllPoints = filter isJust $ map testEmptySph cleanP
    testEmptySph i
        | 0 < powerDist (sP ! i) wC = Nothing
        | otherwise = Just (i, powerDist (sP ! i) wC)

testHullFace :: SetPoint Point3D -> S1 Point3D -> Property
testHullFace sP face = case calcPlane sP face of
    Nothing -> msgFail "no plane from face" False
    Just pl ->
        let
            ps = let size = Vec.length sP in if size <= 0 then [] else [0 .. size - 1]
            (pA, pB, pC) = face3DPoints face
            cleanP = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC)) ps
            pp = pointSetPartition (whichSideOfPlane pl) sP cleanP

            testSigma =
                let
                    newND
                        | (null . pointsOnB1) pp = nd
                        | otherwise = neg nd
                    nd = plane3DNormal pl
                    actFace = ActiveUnit{activeUnit = face, assocP = undefined, assocND = newND}
                 in
                    case makeSimplex actFace sP ps of
                        Just sigma -> testProperTetrahedron sP sigma
                        _ -> msgFail "No possible tetraredron for hull face!" False

            testHull = case (pointsOnB1 pp, pointsOnB2 pp, pointsOnPlane pp) of
                ([], [], []) -> msgFail "no points on partition" False
                ([], [], _) -> msgFail "all points on plane" False
                ([], _, _) -> label "face on B1" True
                (_, [], _) -> label "face on B2" True
                _ -> msgFail ("non-Hull face", pp, face3DPoints face, map (\i -> (i, planeNormal pl &. (sP !. i &- sP !. pA))) ps) False
         in
            testHull .&&. testSigma

data TestClosure = Open | Closed | Error deriving (Show, Eq)
newtype ClosureID = CloID (Int, Int, Int) deriving (Show, Eq)
instance Ord ClosureID where
    compare (CloID a) (CloID b) = compFace a b

testClosure :: IntMap (S2 Point3D) -> Set (S1 Point3D) -> Property
testClosure ss2 ss1 = testError
  where
    testError = Map.foldl (\acc x -> acc .&&. ftest x) (property True) cloS2
    ftest x =
        (msgFail ("Missing face!!!: " ++ (show . Map.filter (== Open) $ cloS2)) $ x /= Open)
            .&&. (msgFail ("More than 2 faces per edge!!!") $ x /= Error)

    cloS1 = S.foldl funcS1 Map.empty ss1
    funcS1 acc s1 =
        let
            face = CloID $ face3DPoints s1
            func _ Open = Closed
            func _ _ = Error
         in
            Map.insertWith func face Open acc

    cloS2 = IM.foldr funcS2 cloS1 ss2
    funcS2 s2 acc =
        let
            (a, b, c, d) = tetraPoints s2
            func _ Open = Closed
            func _ _ = Error
            add p1 p2 p3 = Map.insertWith func (CloID (p1, p2, p3)) Open
         in
            add a b c $ add a b d $ add b c d $ add c a d acc
