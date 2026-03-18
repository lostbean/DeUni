module Main where

import Test.Hspec
import Test.QuickCheck

import qualified Check2D as C2D
import qualified Check3D as C3D
import qualified CheckCommon as C2C

main :: IO ()
main = hspec $ do
    describe "Common Geometric Properties" $ do
        it "projAonB is idempotent" $
            property C2C.prop_projAonB_idempotent
        it "projAonB is orthogonal" $
            property C2C.prop_projAonB_orthogonal
        it "powerDist is symmetric" $
            property C2C.prop_powerDist_symmetric
        it "powerDist is consistent with Euclidean distance for zero weights" $
            property C2C.prop_powerDist_euclidean
        it "whichSideOfPlane is consistent with normal and distance" $
            property C2C.prop_whichSideOfPlane_consistency
        it "getMaxDistPoint identifies the furthest point" $
            property C2C.prop_getMaxDistPoint
        it "getMaxDistPointOnDir identifies the furthest point along a direction" $
            property C2C.prop_getMaxDistPointOnDir
        it "compEdge provides a consistent total ordering" $
            property C2C.prop_compEdge_consistent
        it "compFace provides a consistent total ordering" $
            property C2C.prop_compFace_consistent
        it "findMinimunButZero identifies the smallest non-zero value" $
            property C2C.prop_findMinimunButZero

    describe "2D Delaunay and Geometric Tests" $ do
        it "Calculates correct circumcircles" $
            property C2D.prop_CircumCircle
        it "Finds a valid first edge" $
            property C2D.prop_1stEdge
        it "Generates a valid first simplex" $
            property C2D.prop_1stSimplex
        it "Delaunay triangulation handles minimum point set (3 points)" $
            property C2D.prop_DelaunayMinPoints
        -- prop_Delaunay is very sensitive to degeneracies, running it with fewer tests or as a separate spec if it's too slow
        it "Delaunay triangulation property (Empty Circumcircle)" $
            withMaxSuccess 100 C2D.prop_Delaunay

    describe "3D Delaunay and Convex Hull Tests" $ do
        it "Calculates correct circumspheres" $
            property C3D.prop_CircumSphere
        it "Generates a valid first face" $
            property C3D.prop_1stFace
        it "Generates a valid first tetrahedron" $
            property C3D.prop_1stSimplex
        it "Delaunay triangulation handles minimum point set (4 points)" $
            property C3D.prop_DelaunayMinPoints
        it "Convex Hull property" $
            withMaxSuccess 100 C3D.prop_ConvHull
        it "Delaunay triangulation property (Empty Circumsphere)" $
            withMaxSuccess 100 C3D.prop_Delaunay
