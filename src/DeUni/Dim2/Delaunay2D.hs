-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Delaunay3D |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE TypeFamilies #-}


module DeUni.Dim3.Delaunay2D where

import Prelude hiding (null, lookup)
import Data.List (map, foldl')
import Control.Applicative ((<$>))
import Data.Maybe
import Data.Array.Diff hiding (elems)

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim3.Base2D
import DeUni.Dim3.Hull2D
import Math.Vector



instance Buildable S2 Point2D where
  type Sub S2    = S1
  buildUnit      = makeSimplex
  build1stUnit   = makeFirstSimplex
  getAllSubUnits = extractAllSimplexFaces
  subUnitPos     = face3DPos

makeFirstSimplex::Plane Point3D -> SetPoint Point3D -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (S2 Point3D)
makeFirstSimplex alpha sP sideA sideB ps = do
    face <- build1stUnit alpha sP sideA sideB ps
    let actFace = ActiveUnit { activeUnit=face, assocP=(-1), assocND=face3DND face }
    makeSimplex actFace sP ps

face3DPos pairBox sP face = let (a, b, c) = (face3DPoints.activeUnit) face in facePos pairBox sP a b c

extractAllSimplexFaces::SetPoint Point3D -> S2 Point3D -> [ActiveSubUnit S2 Point3D]
extractAllSimplexFaces sP sigma = map toSimplexFace fsAll
    where
    (a,b,c,d) = tetraPoints sigma
    fsAll     = [((a,b,d), c), ((a,d,c), b), ((d,b,c), a), ((a,b,c), d)]
    toSimplexFace fp@(f, x) = ActiveUnit { activeUnit=(Face3D f nd), assocP=x, assocND=nd }
        where nd = outterND fp
    outterND ((a,b,c), x) = if (sP!.a &- sP!.x) &. nd > 0 then neg nd else nd
        where nd = normalize (sP!.b &- sP!.a) &^ (sP!.c &- sP!.a)

makeSimplex::ActiveSubUnit S2 Point3D -> SetPoint Point3D -> [PointPointer] -> Maybe (S2 Point3D)
makeSimplex actFace sP ps = do
    minR  <- findMinRadius
    sigma <- buildSimplexFace minR
    return sigma
    --if debug ("sigma: " ++ show sigma) $ testProperTetrahedron sP ps sigma then return sigma else (error "Fuck!!!!")
    where
        buildSimplexFace (_, d) = return Tetrahedron { circumSphereCenter = center
                                                     , circumRadius = radius
                                                     , tetraPoints  = (a,b,c,d) }
            where (radius, center) = getCircumSphere (sP!.a, sP!.b, sP!.c) (sP!.d)

        -- | Remove points from face to avoid get 0.0 in findMin
        cleanP        = filter (\i -> (isSideOk i) && (i /= a) && (i /= b) && (i /= c)) ps
        findMinRadius = findMinimunButZero (getRadius sP actFace) sP cleanP
        isSideOk i    = 0 < (sP!.a &- sP!.i) &. nd
        face@(a,b,c)  = (face3DPoints.activeUnit) actFace
        nd            = (face3DND.activeUnit) actFace


getRadius::SetPoint Point3D -> ActiveSubUnit S2 Point3D -> Point3D -> Double
getRadius sP actFace i
    | (getSide center)       && (getSide i)       = radius
    | (not $ getSide center) && (not $ getSide i) = radius
    | otherwise                                   = (-radius)
    where
        nd               = (face3DND.activeUnit) actFace
        getSide x        = 0 > nd &. (sP!.a &- x)
        face@(a,b,c)     = (face3DPoints.activeUnit) actFace
        (radius, center) = getCircumSphere (sP!.a, sP!.b, sP!.c) i


getCircumSphere::(Point3D, Point3D, Point3D) -> Point3D -> (Double, Point3D)
getCircumSphere (a, b, c) d = (radius, center)
    where
        radius   = abs $ (norm q)/div
        center   = a &+ (q &* (1/div))

        ref      = a
        deltaA   = b &- ref
        deltaB   = c &- ref
        deltaC   = d &- ref
        crossB_C = deltaB &^ deltaC
        crossC_A = deltaC &^ deltaA
        crossA_B = deltaA &^ deltaB
        x        = (normsqr deltaA) *& crossB_C
        w        = (normsqr deltaB) *& crossC_A
        t        = (normsqr deltaC) *& crossA_B
        div      = 2 * (deltaA &. crossB_C)
        q        = x &+ w &+ t


-- | Performance can be improve by removing the duplicate call to "func" in "dropZero" and the first "(func x, x)"
-- | OBS: Not the closest to zero. In that case
findClosestButZero::(Point3D -> Double) -> SetPoint Point3D -> [PointPointer] -> Maybe (Double, PointPointer)
findClosestButZero func = findMinimunButZero (abs.func)


-- | Performance can be improve by removing the duplicate call to "func" in "dropZero" and the first "(func x, x)"
-- | OBS: Not the closest to zero. In that case
findMinimunButZero::(Point3D -> Double) -> SetPoint Point3D -> [PointPointer] -> Maybe (Double, PointPointer)
findMinimunButZero func sP ps = case pStartWithNoZero of
    []     -> Nothing
    (x:xs) -> Just $ foldl' (\pair i -> foldMaybe pair (func' i, i)) (func' x, x) xs
    where
      func' = func.(sP!.)
      pStartWithNoZero = dropWhile dropZero ps
      dropZero = (flip$(==).func') 0
      foldMaybe new@(n, i) old@(nOld, iOld)
        | n == 0 = old
        | n > nOld = old
        | n < nOld = new
        | otherwise = error $ "Multiple points on circle or sphere! " ++ show new 

