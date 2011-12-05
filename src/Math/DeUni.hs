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
, reRunDeWall
, Face         (facePoints, refND)
, Simplex      (circumSphereCenter, circumRadius, setCellID)
, PointPointer (..)
, Point
, SetPoint
, SetSimplex
, SetFace
, Box          (..)
, StateVarsMBC (..)
, SetActiveSubUnits
, ActiveSubUnit (..)
, extractAllSimplexFaces

-- Remove after test
, getCircumSphere
, genPlane
, pointSetPartition
, whichBoxIsIt
, pointsOnB1
, pointsOnB2
, pointsOnPlane
, makeFirstSimplex
, makeSimplex
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

import Math.WeightedPoint as WP
import Math.DeWall
import Math.DeHull
import Math.Types
import Math.GeometricTools


import Debug.Trace
debug :: Show a => String -> a -> a
debug s x = x --trace (s ++ show x) x


-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Exposed functions |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reRunDeWall::StateVarsMBC Face -> Box -> SetPoint -> [PointPointer] -> SetActiveSubUnits Face -> (SetSimplex, StateVarsMBC Face)
reRunDeWall st box sP ps faces = runState (mbc ps faces box) init
  where
    init = st { aflAlpha=empty, aflBox1=empty, aflBox2=empty, setPoint=sP }

runDeHull::Box -> SetPoint -> [PointPointer] -> (SetFace, StateVarsMBC Edge)
runDeHull box sP ps = runState (mbc ps (empty::SetActiveSubUnits Edge) box) init
  where
    init = initState sP


runDeWall::Box -> SetPoint -> [PointPointer] -> (SetSimplex, StateVarsMBC Face)
runDeWall box sP ps = runState (mbc ps (empty::SetActiveSubUnits Face) box) init
  where
    init = initState sP

initState sP = StateVarsMBC
    { aflAlpha      = empty
    , aflBox1       = empty
    , aflBox2       = empty
    , externalFaces = empty
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
                      (us2, st2) = runState (mbc p2 afl2 box2) st1
                  --modify (\x -> x { count = (count st1) - c + (count st2) })
                  modify (\x -> x {externalFaces=(externalFaces st1) `union` (externalFaces st2)})
                  modify (\x -> x { count = count st2 })
                  return (us1 `IM.union` us2 `IM.union` units)
                where
                afl1 = aflBox1 st
                afl2 = aflBox2 st
                box1 = halfBox1 pairBox
                box2 = halfBox2 pairBox


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
          modify (\x -> x { externalFaces = insert actSubUnit (externalFaces x) })
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


-- ^^^^^^^^^^^^^^^ Internal (in module) Test
testProperTetrahedron::SetPoint -> [PointPointer] -> Simplex -> Bool
testProperTetrahedron sP ps sigma = if isRadiusOK && isCenterOK && isSphereOK
                                    then True
                                    else error $ "Puta Merda " ++ show (isCenterOK, isSphereOK)
    where
    error_precisson = (10e-8)
    (pA,pB,pC,pD)  = setCellID sigma
    center         = circumSphereCenter sigma
    radius         = circumRadius sigma
    isRadiusOK     = error_precisson > (abs $ radius - (norm $ sP!pA - center))
    cleanP         = filter (\i -> (i /= pA) && (i /= pB) && (i /= pC) && (i /= pD)) ps
    -- | Test if the it is the center of the simplex
    isCenterOK     = and $ map testCenter [pA,pB,pC,pD]
    testCenter i   = error_precisson > (abs $ (norm $ center - sP!i) - radius)
    -- | Test if the CircumSphere is empty
    isSphereOK     = and $ map testEmptySph cleanP
    testEmptySph i = if radius < norm (sP!i - center)
                     then True
                     else error $ "Puta Merda " ++ show (i, sP!i, sigma)

