-----------------------------------------------------------------------------
--
-- Module      :  DeHull
-- Copyright   :
-- License     :  AllRightsReserved
--
-- Maintainer  :  Edgar Gomes
-- Stability   :  dev
-- Portability :
--
-- |
--
-----------------------------------------------------------------------------

{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE OverlappingInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE TypeFamilies #-}


module DeUni
( runDeHull
, SimplexHullFace (facePoints, outterND)
) where


import Prelude hiding (null, lookup)
import Data.Vec hiding (map, length, fromList, fold, get, Map)
import Data.List (map, foldl', filter, head, (\\), minimumBy, maximumBy)
import qualified Data.List as L
import Data.Set ( Set, deleteFindMax, member, empty, null, delete, insert, fromList, fold, elems, union, findMax )
import Data.Maybe
import Control.Applicative ((<$>))
import Control.Monad.State.Lazy




-- TODO Move Box to other module


import Debug.Trace
debug :: Show a => String -> a -> a
debug s x = trace (s ++ show x) x


-- >>>>>>>>>>> Type definitions <<<<<<<<<<<<<<<<<<<
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- | Define a point in 3D (x,y,z)
type Point = Vec3D

instance Ord Vec3D where
    compare p1@(Vec3D a1 b1 c1) p2@(Vec3D a2 b2 c2)
        | a1 > a2   = GT
        | a1 < a2   = LT
        | b1 > b2   = GT
        | b1 < b2   = LT
        | c1 > c2   = GT
        | c1 < c2   = LT
        | otherwise = EQ


-- | Create a structure for face (Triangles in the case of 3D DT) and store the orientation
--   of the face related to the previous generated simplex that create the face.

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data Simplex = Simplex
    { circumSphereCenter :: Point
    , setCellID          :: (Point, Point, Point, Point)
    } deriving (Show, Eq)

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Instances for edge
data Edge = Edge
    { pointL :: Point
    , pointR :: Point
    } deriving (Show)

instance Ord Edge where
    compare a@(Edge a1 a2) b@(Edge b1 b2)
        | amax > bmax = GT
        | amax < bmax = LT
        | amin > bmin = GT
        | amin < bmin = LT
        | otherwise   = EQ
        where
            amax = (eUp a1 a2)
            amin = (eBot a1 a2)
            bmax = (eUp b1 b2)
            bmin = (eBot b1 b2)
            eUp  x y = if compare x y == GT then x else y
            eBot x y = if compare x y == GT then y else x

instance Eq Edge where
    x == y = compare x y == EQ


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data Face = Face
   { facePoints :: (Point, Point, Point)
   , refND      :: Point
   } deriving Show

instance Ord Face where
    compare a b = compare a' b'
        where
        a' = (fast3DSort.facePoints) a
        b' = (fast3DSort.facePoints) b
        fast3DSort::(Point, Point, Point) -> (Point, Point, Point)
        fast3DSort face@(a, b, c)
           | (a >= b) && (b >= c) = face
           | otherwise = (a', b', c')
            where
                minab = min a b
                maxab = max a b
                a'    = max (maxab) c
                b'    = max (min (maxab) c) (minab)
                c'    = min (minab) c

instance Eq Face where
    x == y = compare x y == EQ


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- ActiveSubUnit
data (SubUnit a, Ord a)=>ActiveSubUnit a = ActiveUnit
    { activeUnit :: a
    , assocP     :: Point
    , assocND    :: Point
    } deriving (Show, Eq)

instance (SubUnit a, Ord a)=>Ord (ActiveSubUnit a) where
    compare a b = compare e1 e2
        where
            e1 = activeUnit a
            e2 = activeUnit b

type SetActiveSubUnits a = Set a


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- | Define the a plane that dissect the space in two half-space.
--   All the elements of DT will be constructed recursively over this plane.
--   Plane is represented by a vector from oring (0,0,0) and a scalar such that k*(a,b,c)
--   is the closest point in the plane to the oring
data Plane = Plane
    { planeNormal::Vec3D
    , planeDist  ::Double
    } deriving (Show, Eq)

data PointPartition = PointPartition
    { pointsOnB1   :: [Point]
    , pointsOnB2   :: [Point]
    , pointsOnPlane:: [Point]
    } deriving (Show)

-- | Define possible possitions of the elements for the 1st half-space (Box1=B1),
--   2nd (Box2=B2) and intersect by the plane (B1B2Plane).
data Position = B1
              | B2
              | OnPlane
              | CrossPlane
              | None
              deriving (Show, Eq)

data Box = Box
    { xMax::Double
    , xMin::Double
    , yMax::Double
    , yMin::Double
    , zMax::Double
    , zMin::Double
    } deriving (Show)

data BoxPair = BoxPair
    { halfBox1::Box
    , halfBox2::Box
    } deriving (Show)

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- | Group the data that must be update along the computation (State).
--   Use of state monad will make it clear and keep the purity of the code.
data (SubUnit a, Ord a)=>SubUnitsSets a = SubUnitsSets
    { aflAlpha, aflBox1, aflBox2 :: SetActiveSubUnits a
    , externalFaces              :: [a]
    } deriving (Show)

type StateMBC a = State (SubUnitsSets a)

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class SubUnit subUnit where
    type Unit subUnit :: *
    buildUnit      :: ActiveSubUnit subUnit -> [Point] -> Maybe (Unit subUnit)
    build1stUnit   :: Plane -> [Point] -> [Point] -> [Point] -> Maybe (Unit subUnit)
    getAllSubUnits :: Unit subUnit -> [ActiveSubUnit subUnit]
    subUnitPos     :: BoxPair -> subUnit -> Position

instance SubUnit Edge where
    type Unit Edge = Face
    buildUnit      = makeFace
    build1stUnit   = makeFirstFace
    getAllSubUnits = extractAllFaceEdges
    subUnitPos     = edgePos

instance SubUnit Face where
    type Unit Face = Simplex
    buildUnit      = makeSimplex
    build1stUnit   = undefined
    getAllSubUnits = extractAllSimplexFaces
    subUnitPos     = facePos
{-
class SubUnit subUnit unit | subUnit -> unit, unit -> subUnit where
    buildUnit      :: ActiveSubUnit subUnit -> [Point] -> Maybe unit
    build1stUnit   :: Plane -> [Point] -> [Point] -> [Point] -> Maybe unit
    getAllSubUnits :: unit -> [ActiveSubUnit subUnit]
    subUnitPos     :: BoxPair -> subUnit -> Position

instance SubUnit Edge Face where
    buildUnit      = makeFace
    build1stUnit   = makeFirstFace
    getAllSubUnits = extractAllFaceEdges
    subUnitPos     = edgePos

instance SubUnit Face Simplex where
    buildUnit      = makeSimplex
    build1stUnit   = undefined
    getAllSubUnits = extractAllSimplexFaces
    subUnitPos     = facePos
-}
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%| Exposed functions |%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runDeHull::Box -> [Point] -> [Face]
runDeHull box ps = fst $ runState (mbc ps (empty::SetActiveSubUnits Edge) box) initState


initState = SubUnitsSets
    { aflAlpha      = empty
    , aflBox1       = empty
    , aflBox2       = empty
    , externalFaces = []
    }


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%| Non-Exposed functions |%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-- | Marriage Before Conquer
mbc::(SubUnit a)=>[Point] -> SetActiveSubUnits a -> Box -> StateMBC a [Unit a]
mbc p afl box = do
    cleanAFLs
    if (null afl)
        then do
            (unit, afl) <- case build1stUnit plane p1 p2 p of
                Just un  -> return (([un], fromList $ getAllSubUnits un)::(SubUnit a)=>([Unit a],SetActiveSubUnits a))
                _        -> return ([]   , empty)
            analyze1stUnit unit afl
        else do
            mapM_ (splitAF pairBox) (elems afl)
            units <- getUnitsOnPlane p pairBox plane
            analyzeUnit pairBox p1 p2 units
    where
        (plane, pairBox) = genPlane box
        pPartition       = pointSetPartition (whichBoxIsIt pairBox) p
        p1               = pointsOnB1 pPartition
        p2               = (pointsOnB2 pPartition) ++ (pointsOnPlane pPartition)
        cleanAFLs        = modify (\x -> x { aflAlpha=empty, aflBox1=empty, aflBox2=empty })

        analyze1stUnit::(SubUnit a)=>[Unit a] -> SetActiveSubUnits a -> StateMBC a [Unit a]
        analyze1stUnit unit afl
            | (not.null) afl       = do
                us  <- mbc p afl box
                return (us ++ unit)

            | p1 == [] && p2 /= [] = do
                us2 <- mbc p2 afl (halfBox2 pairBox)
                return $ us2 ++ unit

            | p1 /= [] && p2 == [] = do
                us1 <- mbc p1 afl (halfBox1 pairBox)
                return $ us1 ++ unit

            | otherwise = return []


        --analyzeUnit::BoxPair -> [Point] -> [Point] -> [b] -> StateMBC a [b]
        analyzeUnit pairBox p1 p2 units = get >>= redirect
            where
            redirect st
                | null afl1 && null afl2 = return units
                | (null) afl1 = do
                    us <- mbc p2 afl2 box2
                    return (us ++ units)
                | (null) afl2 = do
                    us <- mbc p1 afl1 box1
                    return (us ++ units)
                | otherwise   = do
                    -- The state is independant and can discarted as it will be
                    -- ereased at the bigein of the next recursive func call
                    return (us1 ++ us2 ++ units)
                where
                !afl1 = aflBox1 st
                !afl2 = aflBox2 st
                box1 = halfBox1 pairBox
                box2 = halfBox2 pairBox
                us1 = evalState (mbc p1 afl1 box1) st
                us2 = evalState (mbc p2 afl2 box2) st


-- Simplex Wall Construction
getUnitsOnPlane::(SubUnit a)=>[Point] -> BoxPair -> Plane -> StateMBC a [Unit a]
getUnitsOnPlane p pairBox plane = do
    st <- get
    if null (aflAlpha st)
        then do
            return []
        else do
            subUnit <- getOneSubUnit st
            recursion (buildUnit subUnit p) subUnit
    where
        getOneSubUnit       = return.findMax.aflAlpha
        getOthersSubUnits x = return.(L.delete x).getAllSubUnits
        removeSubUnit su    = modify (\x -> x { aflAlpha = delete su (aflAlpha x) })
        recursion t subUnit = case t of
            Just sig -> do
                getOthersSubUnits subUnit sig >>= mapM_ (splitAF pairBox)
                removeSubUnit subUnit
                s <- getUnitsOnPlane p pairBox plane
                return (sig:s)
            _ -> do
                modify (\x -> x { externalFaces = subUnit : (externalFaces x) })
                removeSubUnit subUnit
                getUnitsOnPlane p pairBox plane


splitAF::(SubUnit a)=>BoxPair -> a -> StateMBC a ()
splitAF pairBox e = case subUnitPos pairBox (activeUnit e) of
    B1         -> upP1
    B2         -> upP2
    CrossPlane -> upAlpha
    _          -> return ()
    where
        upP1     = updateSubUnit e modP1
        upP2     = updateSubUnit e modP2
        upAlpha  = updateSubUnit e modAlpha
        modP1    = (aflBox1 , \mod x -> x { aflBox1  = mod })
        modP2    = (aflBox2 , \mod x -> x { aflBox2  = mod })
        modAlpha = (aflAlpha, \mod x -> x { aflAlpha = mod })

        updateSubUnit edge (func, mod) = do
            set <- func <$> get
            case member edge set of
                False -> do
                    modify (mod $ insert edge set)
                    return ()
                True  -> do
                    modify (mod $ delete edge set)
                    return ()


pointSetPartition::(Point -> Position) -> [Point] -> PointPartition
pointSetPartition func setPoint =  convert $ splitInBox ([],[],[]) setPoint
    where
        convert (p1,p2,pA) = PointPartition p1 p2 pA
        splitInBox (p1,p2,pA) []     = (p1,p2,pA)
        splitInBox (p1,p2,pA) (x:xs) = case func x of
            B1      -> splitInBox (x:p1,p2,pA) xs
            B2      -> splitInBox (p1,x:p2,pA) xs
            OnPlane -> splitInBox (p1,p2,x:pA) xs
            _       -> splitInBox (p1,p2,pA)   xs

whichBoxIsIt::BoxPair -> Point -> Position
whichBoxIsIt pairBox p
    | inbox1           = B1
    | inbox2           = B2
    | inbox1 && inbox2 = OnPlane
    | otherwise        = None
    where
        inbox1 = isInBox (halfBox1 pairBox) p
        inbox2 = isInBox (halfBox2 pairBox) p
        isInBox box (Vec3D x y z) =
            between (xMin box) (xMax box) x
         && between (yMin box) (yMax box) y
         && between (zMin box) (zMax box) z
        between min max x
            | min < max = (x >= min) && (max >= x)   -- >= will get points on the edge of the box and store if P1 those are on the commun face
            | min > max = (x <= min) && (max <= x)
            | otherwise = error ("Zero size box: " ++ show (pairBox))

whichSideOfPlane::Plane -> Point -> Position
whichSideOfPlane plane p = case compare projection dist of
    EQ -> OnPlane
    GT -> B1
    LT -> B2
    where
        projection = (planeNormal plane) `dot` p
        dist = planeDist plane


-- TODO verify correctness
genPlane::Box -> (Plane, BoxPair)
genPlane box@Box { xMax, xMin, yMax, yMin, zMax, zMin }
    | (deltaX >= (max deltaY deltaZ)) = (Plane (Vec3D 1 0 0) halfX, BoxPair box { xMax = halfX } box { xMin = halfX })
    | ((max deltaX deltaY) <= deltaZ) = (Plane (Vec3D 0 0 1) halfZ, BoxPair box { zMax = halfZ } box { zMin = halfZ })
    | otherwise                       = (Plane (Vec3D 0 1 0) halfY, BoxPair box { yMax = halfY } box { yMin = halfY })
    where
        deltaX = abs (xMax - xMin)
        deltaY = abs (yMax - yMin)
        deltaZ = abs (zMax - zMin)
        halfX = (xMax + xMin)/2
        halfY = (yMax + yMin)/2
        halfZ = (zMax + zMin)/2

-- | Project a vector\point on the plane that goes throw the oringe.
--   It discard the distance on Plane data. It assumes that the plane pass throw the oringe
getProjOnPlane::Plane -> Point -> Point
getProjOnPlane plane p = projOnPlane p
    where
        nd = planeNormal plane
        -- Projection A in B = B * (A°B)/(B°B)
        projNDPlane x = nd * (Vec3D k k k)
            where k = (dot x nd) / (dot nd nd)
        projOnPlane x = x - (projNDPlane x)


getPlane::(Point,Point,Point) -> Maybe Plane
getPlane (a,b,c)
    | nSize == 0 = Nothing
    | d >= 0     = Just $ Plane { planeNormal = normN
                                , planeDist   = d }
    | d < 0      = Just $ Plane { planeNormal = inv normN
                                , planeDist   = -d }
    where
        n     = pack $ (unpack (b - a)) `cross` (unpack (c - a))
        nSize = sqrt (dot n n)
        normN = normalize n
        d     = normN `dot` a
        inv n = (Vec3D (-1) (-1) (-1)) * n


-- %%%%%%%%%%%%%%%%%%%%%% First Face %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makeFirstFace::Plane -> [Point] -> [Point] -> [Point] -> Maybe Face
makeFirstFace alpha sideA sideB ps = do
    (pA, pB) <- getFirstEdge alpha sideA sideB
    (pC, nd) <- getThrirdPoint pA pB ps
    return $ Face (pA, pB, pC) nd


getFirstEdge::Plane -> [Point] -> [Point] -> Maybe (Point, Point)
getFirstEdge divPlane ps1 ps2 = do
    (p1, d1) <- getFirstPoint divPlane ps1
    (p2, d2) <- getFirstPoint divPlane ps2
    if d1 > d2
        then climber B1 divPlane p1 p2 ps1 ps2
        else climber B2 divPlane p1 p2 ps1 ps2


getFirstPoint::Plane -> [Point] -> Maybe (Point, Double)
getFirstPoint _ []         = Nothing
getFirstPoint alpha (x:xs) = Just $ foldl' func d1 xs
    where
    d1     = (x, dist x)
    dist x = norm $ getProjOnPlane alpha x
    func old x
        | d > (snd old) = (x, d)
        | otherwise     = old
        where d = dist x

climber::Position -> Plane -> Point -> Point -> [Point] -> [Point] -> Maybe (Point, Point)
climber mostDistPoint divPlane p1 p2 ps1 ps2 = goTop p1 p2
    where
        goTop p1 p2 = do
            ppPS1     <- getPP p1 (p1 + projPlaneND) p2 ps1
            ppPS2     <- getPP p1 (p1 + projPlaneND) p2 ps2
            facePlane <- getPlane (p1,  (p1 + projPlaneND), p2)
            let
                ps1B1 = pointsOnB1 ppPS1
                ps1B2 = pointsOnB2 ppPS1
                ps2B1 = pointsOnB1 ppPS2
                ps2B2 = pointsOnB2 ppPS2

                okExit = return (p1, p2)

                move
                    | moveDir   = goTop (selByDist p1 ps1B1) (selByDist p2 ps2B1)
                    | otherwise = goTop (selByDist p1 ps1B2) (selByDist p2 ps2B2)

                moveDir = 0 < (dot (planeNormal facePlane) projPOnDivPlane)

                selByDist x [] = x
                selByDist _ ps = maximumBy (\a b -> dist a `compare` dist b) ps

                dist x = abs $ dot (x-p1) (planeNormal facePlane)

            case (ps1B1, ps1B2, ps2B1, ps2B2) of
                ([], _ , [], _ ) -> okExit
                (_ , [], _ , []) -> okExit
                _                -> move

        getPP::Point -> Point -> Point -> [Point] -> Maybe PointPartition
        getPP _  _  _  [] = Nothing
        getPP p1 p2 p3 ps = func <$> getPlane (p1, p2, p3)
                where func x = pointSetPartition (whichSideOfPlane x) (ps \\ [p1, p2, p3])

        projPOnDivPlane
            | mostDistPoint == B1 = calcProjOnDivPlane p1
            | mostDistPoint == B2 = calcProjOnDivPlane p2

        projPlaneND = (pack.normalize.(cross (unpack $ planeNormal divPlane)).unpack) projPOnDivPlane
        calcProjOnDivPlane = getProjOnPlane divPlane


getThrirdPoint::Point -> Point -> [Point] -> Maybe (Point, Point)
getThrirdPoint pA pB setPoint = scan setPoint
    where
    cleanList x = setPoint \\ [pA, pB, x]
    scan [] = Nothing
    scan (p:ps)
        | p == pA || p == pB = scan ps
        | otherwise = case getPlane (pA, pB, p) of
            Just plane
                | (L.null.pointsOnB1) pp &&
                  (L.null.pointsOnB2) pp -> Nothing
                | (L.null.pointsOnB1) pp -> return (p, nd)
                | (L.null.pointsOnB2) pp -> return (p, ind)
                | otherwise              -> scan ps
                where pp = pointSetPartition (whichSideOfPlane plane) (cleanList p)
                      nd = planeNormal plane
                      ind = nd * (Vec3D (-1) (-1) (-1))
            Nothing    -> scan ps



-- %%%%%%%%%%%%%%%%%%%%%% DeHull %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edgePos::BoxPair -> Edge -> Position
edgePos pairBox (Edge a b) = case (findPos a, findPos b) of
        (B1,B1)      -> B1
        (B2,B2)      -> B2
        (B1,OnPlane) -> B1
        (B2,OnPlane) -> B2
        (OnPlane,B1) -> B1
        (OnPlane,B2) -> B2
        (None,_)     -> None
        (_,None)     -> None
        _            -> CrossPlane
        where findPos  = whichBoxIsIt pairBox


extractAllFaceEdges::Face -> [ActiveSubUnit Edge]
extractAllFaceEdges sigma = fsAll
    where
        (a,b,c) = facePoints sigma
        nd      = refND sigma
        fsAll   = [ ActiveUnit (Edge a b) c nd
                  , ActiveUnit (Edge a c) b nd
                  , ActiveUnit (Edge b c) a nd ]

makeFace::ActiveSubUnit Edge -> [Point] -> Maybe Face
makeFace _ []     = Nothing
makeFace e ps = do
    refPoint  <- get1stAng ps
    (pC, ang) <- findNext refPoint
    nd        <- planeNormal <$> getPlane (pA, pB, pC)
    return $ buildSimplex pC (func pC ang nd)
    where
        pA = (pointL.activeUnit) e
        pB = (pointR.activeUnit) e
        oldND = assocND e
        buildSimplex x nd = Face (pA, pB, x) nd

        get1stAng []    = Nothing
        get1stAng (p:ps)
            | isNaN ang = get1stAng ps
            | otherwise = Just $ (p, ang)
            where ang = calcAngBetweenSimplex e p

        scanMax old x
            | isNaN ang       = old
            | ang < (snd old) = (x, ang)
            | otherwise       = old
            where ang = calcAngBetweenSimplex e x

        findNext a1 = return $ foldl' scanMax a1 ps
        ind x     = ((*) (Vec3D (-1) (-1) (-1))) x
        disND p x   = calcAngBetweenSimplex $ fakeEdge pA pB (p + x)
        fakeEdge a b c = ActiveUnit (Edge a b) c c
        func::Point -> Double -> Point -> Point
        func p ang nd
            | ang > 0   = if dot nd oldND > 0 then ind nd else nd
            | ang < 0   = if dot nd oldND > 0 then nd else ind nd
            | otherwise = nd --if (disND p nd) > ang then nd else ind nd


calcAngBetweenSimplex::ActiveSubUnit Edge -> Point -> Double
calcAngBetweenSimplex ae p
    | pA==p || pB==p || pC==p = 0/0
    | otherwise               = dot (normalToEdge vOnOldFace) (normalToEdge vOnNewFace)
    where
        pA = (pointL.activeUnit) ae
        pB = (pointR.activeUnit) ae
        pC = assocP ae
        edge = pA - pB
        vOnOldFace = pC - pB
        vOnNewFace = p  - pB
        -- Projection A in B = B * (A°B)/(B°B)
        projToEdge x = edge * (Vec3D k k k)
            where k = (dot x edge) / (dot edge edge)
        normalToEdge x = normalize $ x - (projToEdge x)



-- %%%%%%%%%%%%%%%%%%%%%% DeWall %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

facePos::BoxPair -> Face -> Position
facePos pairBox (Face (a,b,c) _ ) = case (findPos a, findPos b, findPos c) of
    (B1,B1,B1)           -> B1
    (B2,B2,B2)           -> B2
    (B1,OnPlane,OnPlane) -> B1
    (B2,OnPlane,OnPlane) -> B2
    (OnPlane,B1,OnPlane) -> B1
    (OnPlane,B2,OnPlane) -> B2
    (OnPlane,OnPlane,B1) -> B1
    (OnPlane,OnPlane,B2) -> B2
    (None,_,_)           -> None
    (_,None,_)           -> None
    (_,_,None)           -> None
    _                    -> CrossPlane
    where findPos  = whichBoxIsIt pairBox


extractAllSimplexFaces::Simplex -> [ActiveSubUnit Face]
extractAllSimplexFaces sigma = map toSimplexFace fsAll
    where
        (a,b,c,d) = setCellID sigma
        fsAll  = [((a,b,d), c), ((a,d,c), b), ((d,b,c), a), ((a,b,c), d)]
        toSimplexFace (f, x) = ActiveUnit { activeUnit=(Face f innerND), assocP=x, assocND=innerND }
            where innerND = undefined --getPlane f whichSideOfPlane x


makeSimplex::ActiveSubUnit Face -> [Point] -> Maybe Simplex
makeSimplex actFace p = do
    minR <- findMinRadius
    buildSimplexFace minR
    where
        buildSimplexFace (_, d) = return Simplex {circumSphereCenter = center, setCellID=(a,b,c,d)}
            where (_, center) = getCircumSphere (a, b, c) d

        -- | Remove points from face to avoid get 0.0 in findMin
        cleanP        = filter (\i -> (i /= a) && (i /= b) && (i /= c)) p
        findMinRadius = findMinimunButZero (getRadius actFace) (filter isSideOk cleanP)
        isSideOk i    = undefined
        face@(a,b,c)  = (facePoints.activeUnit) actFace
        nd            = assocND actFace


getRadius::ActiveSubUnit Face -> Point -> Double
getRadius actFace i
    | (getSide center) && (getSide i) = radius
    | otherwise                       = (-radius)
    where
        nd               = assocND actFace
        getSide          = (0 >) . (dot nd) . (\x-> a - x)
        face@(a,b,c)     = (facePoints.activeUnit) actFace
        (radius, center) = getCircumSphere face i


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
findClosestButZero::(Point -> Double) -> [Point] -> Maybe (Double, Point)
findClosestButZero func = findMinimunButZero (abs.func)


-- | Performance can be improve by removing the duplicate call to "func" in "dropZero" and the first "(func x, x)"
-- | OBS: Not the closest to zero. In that case
findMinimunButZero::(Point -> Double) -> [Point] -> Maybe (Double, Point)
findMinimunButZero func p = case pStartWithNoZero of
    []     -> Nothing
    (x:xs) -> Just $ foldl' (\pair i -> foldMaybe pair (func i, i)) (func x, x) xs
    where
        pStartWithNoZero = dropWhile dropZero p
        dropZero = (flip$(==).func) 0
        foldMaybe new@(n, i) old@(nOld, iOld)
                | n == 0 = old
                | n > nOld = old
                | n < nOld = new
                | otherwise = error $ "Multiple points on circle or sphere! " ++ show new











