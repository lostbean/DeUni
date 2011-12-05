-- >>>>>>>>>>> Type definitions <<<<<<<<<<<<<<<<<<<
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


module Math.Types where

import Prelude
import Data.Vec
import Data.Set (Set)
import Data.IntMap (IntMap)
import Data.Maybe (Maybe)
import Data.Array.Diff (DiffArray, Ix)
import Control.Monad.State.Lazy (State)
import System.Random (StdGen)

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
    , circumRadius       :: Double  
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
data StateVarsMBC a = StateVarsMBC
    { aflAlpha, aflBox1, aflBox2 :: Set (ActiveSubUnit a)
    , externalFaces              :: Set (ActiveSubUnit a)
    , randomSeed                 :: StdGen
    , count                      :: Int
    , setPoint                   :: SetPoint
    } deriving (Show)

type SetActiveSubUnits a = Set (ActiveSubUnit a)
type StateMBC a          = State (StateVarsMBC a)

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class SubUnit subUnit unit | subUnit -> unit, unit -> subUnit where
    buildUnit      :: ActiveSubUnit subUnit -> SetPoint -> [PointPointer] -> Maybe unit
    build1stUnit   :: Plane -> SetPoint -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe unit
    getAllSubUnits :: SetPoint -> unit -> [ActiveSubUnit subUnit]
    subUnitPos     :: BoxPair -> SetPoint -> subUnit -> Position
