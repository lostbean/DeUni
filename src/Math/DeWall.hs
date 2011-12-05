-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| DeWall |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Math.DeWall where

import Prelude hiding (null, lookup)
import Data.Vec hiding (map, length, fromList, fold, get, Map)
import Data.List (map, foldl')
import Control.Applicative ((<$>))
import Data.Maybe
import Data.Array.Diff hiding (elems)

import Math.GeometricTools
import Math.Types
import Math.FirstSeed


instance SubUnit Face Simplex where
    buildUnit      = makeSimplex
    build1stUnit   = makeFirstSimplex
    getAllSubUnits = extractAllSimplexFaces
    subUnitPos     = facePos

makeFirstSimplex::Plane -> SetPoint -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe Simplex
makeFirstSimplex alpha sP sideA sideB ps = do
    face <- makeFirstFace alpha sP sideA sideB ps
    let actFace = ActiveUnit { activeUnit=face, assocP=(-1), assocND=refND face }
    makeSimplex actFace sP ps

facePos::BoxPair -> SetPoint -> Face -> Position
facePos pairBox sP (Face (a,b,c) _ ) = case (findPos $ sP!a, findPos $ sP!b, findPos $ sP!c) of
    (B1,B1,B1)                -> B1
    (B2,B2,B2)                -> B2
    (B1,OnPlane,OnPlane)      -> B1
    (B2,OnPlane,OnPlane)      -> B2
    (OnPlane,B1,OnPlane)      -> B1
    (OnPlane,B2,OnPlane)      -> B2
    (OnPlane,OnPlane,B1)      -> B1
    (OnPlane,OnPlane,B2)      -> B2
    (B1,B1,OnPlane)           -> B1
    (B2,B2,OnPlane)           -> B2
    (OnPlane,B1,B1)           -> B1
    (OnPlane,B2,B2)           -> B2
    (B1,OnPlane,B1)           -> B1
    (B2,OnPlane,B2)           -> B2
    (OnPlane,OnPlane,OnPlane) -> None
    (None,_,_)                -> None
    (_,None,_)                -> None
    (_,_,None)                -> None
    _                         -> CrossPlane
    where findPos  = whichBoxIsIt pairBox


extractAllSimplexFaces::SetPoint -> Simplex -> [ActiveSubUnit Face]
extractAllSimplexFaces sP sigma = map toSimplexFace fsAll
    where
    (a,b,c,d) = setCellID sigma
    fsAll  = [((a,b,d), c), ((a,d,c), b), ((d,b,c), a), ((a,b,c), d)]
    toSimplexFace fp@(f, x) = ActiveUnit { activeUnit=(Face f nd), assocP=x, assocND=nd }
        where nd = outterND fp
    outterND ((a,b,c), x) = if (dot (sP!a - sP!x) nd) > 0 then inv nd else nd
        where
        nd    = pack $ normalize (unpack (sP!b - sP!a)) `cross` (unpack (sP!c - sP!a))
        inv v = (Vec3D (-1) (-1) (-1)) * v

makeSimplex::ActiveSubUnit Face -> SetPoint -> [PointPointer] -> Maybe Simplex
makeSimplex actFace sP ps = do
    minR <- findMinRadius
    sigma <- buildSimplexFace minR
    return sigma
    --if debug ("sigma: " ++ show sigma) $ testProperTetrahedron sP ps sigma then return sigma else (error "Fuck!!!!")
    where
        buildSimplexFace (_, d) = return Simplex { circumSphereCenter = center
                                                 , circumRadius = radius
                                                 , setCellID = (a,b,c,d) }
            where (radius, center) = getCircumSphere (sP!a, sP!b, sP!c) (sP!d)

        -- | Remove points from face to avoid get 0.0 in findMin
        cleanP        = filter (\i -> (isSideOk i) && (i /= a) && (i /= b) && (i /= c)) ps
        findMinRadius = findMinimunButZero (getRadius sP actFace) sP cleanP
        isSideOk i    = 0 < dot (sP!a - sP!i) nd
        face@(a,b,c)  = (facePoints.activeUnit) actFace
        nd            = (refND.activeUnit) actFace


getRadius::SetPoint -> ActiveSubUnit Face -> Point -> Double
getRadius sP actFace i
    | (getSide center)       && (getSide i)       = radius
    | (not $ getSide center) && (not $ getSide i) = radius
    | otherwise                                   = (-radius)
    where
        nd               = (refND.activeUnit) actFace
        getSide x        = 0 > (nd `dot` (sP!a - x))
        face@(a,b,c)     = (facePoints.activeUnit) actFace
        (radius, center) = getCircumSphere (sP!a, sP!b, sP!c) i


getCircumSphere::(Point, Point, Point) -> Point -> (Double, Point)
getCircumSphere (a, b, c) d = (radius, center)
    where
        radius = abs $ (norm q)/div
        center = a + (q/(pack $ vec div))

        ref = a
        deltaA = unpack (b - ref)
        deltaB = unpack (c - ref)
        deltaC = unpack (d - ref)
        crossB_C = (deltaB `cross` deltaC)
        crossC_A = (deltaC `cross` deltaA)
        crossA_B = (deltaA `cross` deltaB)
        x = ((norm2 deltaA) * crossB_C)
        w = ((norm2 deltaB) * crossC_A)
        t = ((norm2 deltaC) * crossA_B)
        norm2 x = vec n
            where n = dot x x
        div = 2 * (deltaA `dot` crossB_C)
        q = pack (x+w+t)


-- | Performance can be improve by removing the duplicate call to "func" in "dropZero" and the first "(func x, x)"
-- | OBS: Not the closest to zero. In that case
findClosestButZero::(Point -> Double) -> SetPoint -> [PointPointer] -> Maybe (Double, PointPointer)
findClosestButZero func = findMinimunButZero (abs.func)


-- | Performance can be improve by removing the duplicate call to "func" in "dropZero" and the first "(func x, x)"
-- | OBS: Not the closest to zero. In that case
findMinimunButZero::(Point -> Double) -> SetPoint -> [PointPointer] -> Maybe (Double, PointPointer)
findMinimunButZero func sP ps = case pStartWithNoZero of
    []     -> Nothing
    (x:xs) -> Just $ foldl' (\pair i -> foldMaybe pair (func' i, i)) (func' x, x) xs
    where
      func' = func.(sP!)
      pStartWithNoZero = dropWhile dropZero ps
      dropZero = (flip$(==).func') 0
      foldMaybe new@(n, i) old@(nOld, iOld)
        | n == 0 = old
        | n > nOld = old
        | n < nOld = new
        | otherwise = error $ "Multiple points on circle or sphere! " ++ show new 
