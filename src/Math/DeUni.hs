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
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}


module Math.DeUni
( runDeHull
, runDeWall
, Face      (facePoints, refND)
, Simplex   (circumSphereCenter, setCellID)
, getCircumSphere
, Box       (..)


-- Remove after test
, genPlane
, pointSetPartition
, whichBoxIsIt
, pointsOnB1
, pointsOnB2
, pointsOnPlane
, makeFirstSimplex
, makeFirstFace
, makeSimplex
, ActiveSubUnit (..)
, whichSideOfPlane
, getPlane
, isInBox
) where


import Prelude hiding (null, lookup)
import Data.Vec hiding (map, length, fromList, fold, get, Map)
import Data.List (map, foldl', filter, head, (\\), minimumBy, maximumBy)
import qualified Data.List as L
import Data.Set ( Set, deleteFindMax, member, empty, null, delete, insert, fromList, fold, elems, union, findMax )
import qualified Data.IntMap as IM
import Data.IntMap (IntMap)
import Data.Maybe
import Control.Applicative ((<$>))
import Control.Monad.State.Lazy
import System.Random
import Data.Array.Diff hiding (elems)


import Debug.Trace
debug :: Show a => String -> a -> a
debug s x = x --trace (s ++ show x) x


-- >>>>>>>>>>> Type definitions <<<<<<<<<<<<<<<<<<<
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- | Define a point in 3D (x,y,z)
type Point           = Vec3D
type SetPoint        = DiffArray PointPointer Point
type SetSimplex      = IntMap Simplex
type SetFace         = IntMap Face
newtype PointPointer = PointPointer Int deriving (Eq, Ord, Num, Ix, Real, Enum, Integral, Show)

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
    , setCellID          :: (PointPointer, PointPointer, PointPointer, PointPointer)
    } deriving (Show, Eq)

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Instances for edge
data Edge = Edge
    { pointL :: PointPointer
    , pointR :: PointPointer
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
   { facePoints :: (PointPointer, PointPointer, PointPointer)
   , refND      :: Point
   } deriving Show

instance Ord Face where
    compare a b = compare a' b'
        where
        a' = (fast3DSort.facePoints) a
        b' = (fast3DSort.facePoints) b
        fast3DSort::(PointPointer, PointPointer, PointPointer) -> (PointPointer, PointPointer, PointPointer)
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
data ActiveSubUnit a = ActiveUnit
    { activeUnit :: a
    , assocP     :: PointPointer
    , assocND    :: Point
    } deriving (Show, Eq)

instance (Ord a)=>Ord (ActiveSubUnit a) where
    compare a b = compare e1 e2
        where
            e1 = activeUnit a
            e2 = activeUnit b

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
    { pointsOnB1   :: [PointPointer]
    , pointsOnB2   :: [PointPointer]
    , pointsOnPlane:: [PointPointer]
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
data SubUnitsSets a = SubUnitsSets
    { aflAlpha, aflBox1, aflBox2 :: Set (ActiveSubUnit a)
    , externalFaces              :: [ActiveSubUnit a]
    , randomSeed                 :: StdGen
    , count                      :: Int
    , setPoint                   :: SetPoint
    } deriving (Show)

type SetActiveSubUnits a = Set (ActiveSubUnit a)
type StateMBC a          = State (SubUnitsSets a)

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class SubUnit subUnit unit | subUnit -> unit, unit -> subUnit where
    buildUnit      :: ActiveSubUnit subUnit -> SetPoint -> [PointPointer] -> Maybe unit
    build1stUnit   :: Plane -> SetPoint -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe unit
    getAllSubUnits :: SetPoint -> unit -> [ActiveSubUnit subUnit]
    subUnitPos     :: BoxPair -> SetPoint -> subUnit -> Position

instance SubUnit Edge Face where
    buildUnit      = makeFace
    build1stUnit   = makeFirstFace
    getAllSubUnits = extractAllFaceEdges
    subUnitPos     = edgePos

instance SubUnit Face Simplex where
    buildUnit      = makeSimplex
    build1stUnit   = makeFirstSimplex
    getAllSubUnits = extractAllSimplexFaces
    subUnitPos     = facePos
--}
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Exposed functions |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runDeHull::Box -> SetPoint -> [PointPointer] -> (SetFace, SubUnitsSets Edge)
runDeHull box sP ps = runState (mbc ps (empty::SetActiveSubUnits Edge) box) init
  where
    init = initState sP


runDeWall::Box -> SetPoint -> [PointPointer] -> (SetSimplex, SubUnitsSets Face)
runDeWall box sP ps = runState (mbc ps (empty::SetActiveSubUnits Face) box) init
  where
    init = initState sP

initState sP = SubUnitsSets
    { aflAlpha      = empty
    , aflBox1       = empty
    , aflBox2       = empty
    , externalFaces = []
    , randomSeed    = mkStdGen 240779
    , count         = 0
    , setPoint      = sP
    }


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Non-Exposed functions |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


-- | Marriage Before Conquer
mbc::(SubUnit a b, Ord a, Show a)=>[PointPointer] -> SetActiveSubUnits a -> Box -> StateMBC a (IntMap b)
mbc p afl box = do
    cleanAFLs
    (plane, pairBox) <- genBox box
    sP <- liftM setPoint get
    let
        pp = pointSetPartition (whichBoxIsIt pairBox) sP p
        p1 = pointsOnB1 pp
        p2 = (pointsOnB2 pp) ++ (pointsOnPlane pp)
    if (debug "outter: " $ null afl)
        then do
            us <- case build1stUnit plane sP p1 p2 p of
                Just unit  -> do
                    mapM_ (splitAF pairBox) (getAllSubUnits sP unit)
                    units <- getUnitsOnPlane p pairBox plane
                    cnt   <- liftM count get
                    modify (\x -> x { count = cnt + 1 })
                    return $ IM.insert cnt unit units
                _        -> return IM.empty
            analyzeUnit pairBox p1 p2 us
        else do
            mapM_ (splitAF pairBox) (elems afl)
            units <- getUnitsOnPlane p pairBox plane
            analyzeUnit pairBox p1 p2 units
    where
        cleanAFLs = modify (\x -> x { aflAlpha=empty, aflBox1=empty, aflBox2=empty })
        
        analyzeUnit pairBox p1 p2 units = get >>= redirect
            where
            redirect st = case debug "split: " (IM.null units, null afl1, null afl2, L.null p1, L.null p2) of
                -- The state is independant and can discarted as it will be
                -- ereased at the bigein of the next recursive func call
              --(units, afl1,  afl2,  p1,    p2   )
                (True,  True,  True,  False, True ) -> mbc p1 afl1 box1
                (True,  True,  True,  True,  False) -> mbc p2 afl2 box2
                (True,  True,  True,  _ ,    _    ) -> return IM.empty
                (True,  False, True,  False, True ) -> mbc p1 afl1 box1
                (True,  True,  False, True,  False) -> mbc p2 afl2 box2
                (True,  _,     _,     False, False) -> mbc p (union afl1 afl2) box
                (False, True,  True,  _ ,    _    ) -> return units
                (False, True,  False, _ ,    _    ) -> (IM.union units) <$> mbc p2 afl2 box2
                (False, False, True,  _ ,    _    ) -> (IM.union units) <$> mbc p1 afl1 box1
                (False, False, False, _ ,    _    ) -> do
                  let c = count st
                      (us1, st1) = runState (mbc p1 afl1 box1) st
                      (us2, st2) = runState (mbc p2 afl2 box2) st
                  modify (\x -> x { count = (count st1) - c + (count st2) })
                  return (us1 `IM.union` us2 `IM.union` units)
                where
                afl1 = aflBox1 st
                afl2 = aflBox2 st
                box1 = halfBox1 pairBox
                box2 = halfBox2 pairBox
                --us1  = evalState (mbc p1 afl1 box1) st
                --us2  = evalState (mbc p2 afl2 box2) st


-- Simplex Wall Construction
getUnitsOnPlane::(SubUnit a b, Ord a, Show a)=>[PointPointer] -> BoxPair -> Plane -> StateMBC a (IntMap b)
getUnitsOnPlane p pairBox plane = do
    st <- get
    sP <- liftM setPoint get
    let
      getOneActSubUnit    = return.findMax.aflAlpha
      getOthersSubUnits x = return.(L.delete x).(getAllSubUnits sP)
      removeSubUnit su    = modify (\x -> x { aflAlpha = delete su (aflAlpha x) })
      recursion t actSubUnit = case t of
        Just sig -> do
          getOthersSubUnits actSubUnit sig >>= mapM_ (splitAF pairBox)
          removeSubUnit actSubUnit
          s <- getUnitsOnPlane p pairBox plane
          cnt <- liftM count get
          modify (\x -> x { count = cnt + 1 })
          return $ IM.insert cnt sig s
        _ -> do
          modify (\x -> x { externalFaces = actSubUnit : (externalFaces x) })
          removeSubUnit actSubUnit
          getUnitsOnPlane p pairBox plane
    if null (debug "inner: " $ aflAlpha st)
        then do
            return IM.empty
        else do
            actSubUnit <- getOneActSubUnit st
            recursion (buildUnit actSubUnit sP p) actSubUnit


splitAF::(SubUnit a b, Ord a)=>BoxPair -> ActiveSubUnit a -> StateMBC a ()
splitAF pairBox e = do
  sP <- liftM setPoint get
  case subUnitPos pairBox sP (activeUnit e) of
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


pointSetPartition::(Point -> Position) -> SetPoint -> [PointPointer] -> PointPartition
pointSetPartition func sP ps =  convert $ splitInBox ([],[],[]) ps
    where
        convert (p1,p2,pA) = PointPartition p1 p2 pA
        splitInBox (p1,p2,pA) []     = (p1,p2,pA)
        splitInBox (p1,p2,pA) (x:xs) = case func (sP!x) of
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
genBox::(SubUnit a b, Ord a, Show a)=>Box -> StateMBC a (Plane, BoxPair)
genBox box@Box { xMax, xMin, yMax, yMin, zMax, zMin } = do
    seed <- randomSeed <$> get
    let (out, newSeed) = dirSelector seed
    modify (\x -> x { randomSeed=newSeed })
    return out
    where
        dirSelector seed
            | dir == 1  = ((Plane (Vec3D 1 0 0) newX, BoxPair box { xMax = newX } box { xMin = newX }), s4)
            | dir == 2  = ((Plane (Vec3D 0 0 1) newY, BoxPair box { zMax = newY } box { zMin = newY }), s4)
            | otherwise = ((Plane (Vec3D 0 1 0) newZ, BoxPair box { yMax = newZ } box { yMin = newZ }), s4)
            where
                (dir, s1) = randomR (1,3::Int) seed
                (a,   s2) = randomR (0,1) s1
                (b,   s3) = randomR (0,1) s2
                (c,   s4) = randomR (0,1) s3
                newX      = (xMax - xMin)*a + xMin
                newY      = (yMax + yMin)*b + yMin
                newZ      = (zMax + zMin)*c + zMin

isInBox::Box -> Point -> Bool
isInBox box (Vec3D x y z) = (xMax box > x && xMin box < x) &&
                            (yMax box > y && yMin box < y) &&
                            (zMax box > z && zMin box < z)

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

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| First Face |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

makeFirstFace::Plane -> SetPoint -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe Face
makeFirstFace alpha sP sideA sideB ps = do
    (pA, pB) <- getFirstEdge alpha sP sideA sideB
    (pC, nd) <- getThrirdPoint sP pA pB ps
    return $ Face (pA, pB, pC) nd


getFirstEdge::Plane -> SetPoint -> [PointPointer] -> [PointPointer] -> Maybe (PointPointer, PointPointer)
getFirstEdge divPlane sP ps1 ps2 = do
    (p1, d1) <- getFirstPoint divPlane sP ps1
    (p2, d2) <- getFirstPoint divPlane sP ps2
    if d1 > d2
        then climber B1 divPlane sP p1 p2 ps1 ps2
        else climber B2 divPlane sP p1 p2 ps1 ps2


getFirstPoint::Plane -> SetPoint -> [PointPointer] -> Maybe (PointPointer, Double)
getFirstPoint _ _ []         = Nothing
getFirstPoint alpha sP (x:xs) = Just $ foldl' func d1 xs
    where
    d1     = (x, dist x)
    dist x = norm $ getProjOnPlane alpha (sP!x)
    func old x
        | d > (snd old) = (x, d)
        | otherwise     = old
        where d = dist x

climber::Position -> Plane -> SetPoint -> PointPointer -> PointPointer -> [PointPointer] -> [PointPointer] -> Maybe (PointPointer, PointPointer)
climber mostDistPoint divPlane sP p1 p2 ps1 ps2 = goTop p1 p2
    where
        goTop p1 p2 = do
            ppPS1     <- getPP (sP!p1) (sP!p1 + projPlaneND) (sP!p2) ps1
            ppPS2     <- getPP (sP!p1) (sP!p1 + projPlaneND) (sP!p2) ps2
            facePlane <- getPlane (sP!p1,  (sP!p1 + projPlaneND), sP!p2)
            let
                --In a ideal world (p+projPlanaND) should also be execluded
                ps1B1 = (pointsOnB1 ppPS1) \\ [p1, p2]
                ps1B2 = (pointsOnB2 ppPS1) \\ [p1, p2]
                ps2B1 = (pointsOnB1 ppPS2) \\ [p1, p2]
                ps2B2 = (pointsOnB2 ppPS2) \\ [p1, p2]

                okExit = return (p1, p2)

                move
                    | moveDir   = goTop (selByDist p1 ps1B1) (selByDist p2 ps2B1)
                    | otherwise = goTop (selByDist p1 ps1B2) (selByDist p2 ps2B2)

                moveDir = 0 < (dot (planeNormal facePlane) projPOnDivPlane)

                selByDist x [] = x
                selByDist _ ps = maximumBy (\a b -> dist a `compare` dist b) ps

                dist x = abs $ dot (sP!x - sP!p1) (planeNormal facePlane)

            case (ps1B1, ps1B2, ps2B1, ps2B2) of
                ([], _ , [], _ ) -> okExit
                (_ , [], _ , []) -> okExit
                _                -> move

        getPP::Point -> Point -> Point -> [PointPointer] -> Maybe PointPartition
        getPP _  _  _  [] = Nothing
        getPP p1 p2 p3 ps = func <$> getPlane (p1, p2, p3)
                where func x = pointSetPartition (whichSideOfPlane x) sP ps 

        projPOnDivPlane
            | mostDistPoint == B1 = calcProjOnDivPlane (sP!p1)
            | mostDistPoint == B2 = calcProjOnDivPlane (sP!p2)

        projPlaneND = (pack.normalize.(cross (unpack $ planeNormal divPlane)).unpack) projPOnDivPlane
        calcProjOnDivPlane = getProjOnPlane divPlane


getThrirdPoint::SetPoint -> PointPointer -> PointPointer -> [PointPointer] -> Maybe (PointPointer, Point)
getThrirdPoint sP pA pB ps = scan ps
    where
    cleanList x = ps \\ [pA, pB, x]
    scan [] = Nothing
    scan (x:xs)
        | x == pA || x == pB = scan xs
        | otherwise = case getPlane (sP!pA, sP!pB, sP!x) of
            Just plane
                | (L.null.pointsOnB1) pp &&
                  (L.null.pointsOnB2) pp -> Nothing
                | (L.null.pointsOnB1) pp -> return (x, nd)
                | (L.null.pointsOnB2) pp -> return (x, ind)
                | otherwise              -> scan xs
                where pp = pointSetPartition (whichSideOfPlane plane) sP (cleanList x)
                      nd = planeNormal plane
                      ind = nd * (Vec3D (-1) (-1) (-1))
            Nothing    -> scan xs


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| DeHull |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edgePos::BoxPair -> SetPoint -> Edge -> Position
edgePos pairBox sP (Edge a b) = case (findPos $ sP!a, findPos $ sP!b) of
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


extractAllFaceEdges::SetPoint -> Face -> [ActiveSubUnit Edge]
extractAllFaceEdges sP sigma = fsAll
    where
        (a,b,c) = facePoints sigma
        nd      = refND sigma
        fsAll   = [ ActiveUnit (Edge a b) c nd
                  , ActiveUnit (Edge a c) b nd
                  , ActiveUnit (Edge b c) a nd ]

makeFace::ActiveSubUnit Edge -> SetPoint -> [PointPointer] -> Maybe Face
makeFace _ _ [] = Nothing
makeFace e sP ps = do
    refPoint  <- get1stAng ps
    (pC, ang) <- findNext refPoint
    nd        <- planeNormal <$> getPlane (sP!pA, sP!pB, sP!pC)
    return $ buildFace pC (func (sP!pC) ang nd)
    where
        pA = (pointL.activeUnit) e
        pB = (pointR.activeUnit) e
        oldND = assocND e
        buildFace x nd = Face (pA, pB, x) nd

        get1stAng []    = Nothing
        get1stAng (p:ps)
            | isNaN ang = get1stAng ps
            | otherwise = Just $ (p, ang)
            where ang = calcAngBetweenSimplex e sP p

        scanMax old x
            | isNaN ang       = old
            | ang < (snd old) = (x, ang)
            | otherwise       = old
            where ang = calcAngBetweenSimplex e sP x

        findNext a1 = return $ foldl' scanMax a1 ps
        ind x     = ((*) (Vec3D (-1) (-1) (-1))) x
        disND p x   = calcAngBetweenSimplex $ fakeEdge pA pB (sP!p + x)
        fakeEdge a b c = ActiveUnit (Edge a b) a c
        func::Point -> Double -> Point -> Point
        func p ang nd
            | ang > 0   = if dot nd oldND > 0 then ind nd else nd
            | ang < 0   = if dot nd oldND > 0 then nd else ind nd
            | otherwise = nd --if (disND p nd) > ang then nd else ind nd


calcAngBetweenSimplex::ActiveSubUnit Edge -> SetPoint -> PointPointer -> Double
calcAngBetweenSimplex ae sP p
    | pA==p || pB==p || pC==p = 0/0
    | otherwise               = dot (normalToEdge vOnOldFace) (normalToEdge vOnNewFace)
    where
        pA = (pointL.activeUnit) ae
        pB = (pointR.activeUnit) ae
        pC = assocP ae
        pB'        = sP!pB
        edge       = sP!pA - pB'
        vOnOldFace = sP!pC - pB'
        vOnNewFace = sP!p  - pB'
        -- Projection A in B = B * (A°B)/(B°B)
        projToEdge x = edge * (Vec3D k k k)
            where k = (dot x edge) / (dot edge edge)
        normalToEdge x = normalize $ x - (projToEdge x)


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| DeWall |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if debug ("sigma: " ++ show sigma) $ testProperTetrahedron sP ps sigma then return sigma else (error "Fuck!!!!")
    where
        buildSimplexFace (_, d) = return Simplex {circumSphereCenter = center, setCellID=(a,b,c,d)}
            where (_, center) = getCircumSphere (sP!a, sP!b, sP!c) (sP!d)

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


-- ^^^^^^^^^^^^^^^ GARBAGE
testProperTetrahedron::SetPoint -> [PointPointer] -> Simplex -> Bool
testProperTetrahedron sP ps sigma = isCenterOK && isSphereOK
    where
    error_precisson = (10e-8)
    (pA,pB,pC,pD)  = setCellID sigma
    center         = circumSphereCenter sigma
    radius         = norm $ sP!pA - center
    cleanP         = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC) && (i /= pD)) ps
    -- | Test if the it is the center of the simplex
    isCenterOK     = and $ map testCenter [pA,pB,pC,pD]
    testCenter i   = error_precisson > (abs $ (norm $ center - sP!i) - radius)
    -- | Test if the CircumSphere is empty
    isSphereOK     = and $ map testEmptySph cleanP
    testEmptySph i = radius < norm (sP!i - center)








