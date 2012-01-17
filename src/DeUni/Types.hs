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
{-# LANGUAGE UndecidableInstances #-}


module DeUni.Types where

import Prelude
import Data.Set (Set)
import Data.IntMap (IntMap)
import Data.Maybe (Maybe)
import Data.Array.Diff (DiffArray, Ix, (!))
import Control.Monad.State.Lazy (State)

import Math.Vector


import Debug.Trace
debug :: Show a => String -> a -> a
debug s x = trace (s ++ show x) x


-- | Define a point in 2D and 3D
type Point3D         = Vec3
type Point2D         = Vec2
newtype PointPointer = PointPointer Int deriving (Eq, Ord, Num, Ix, Real, Enum, Integral, Show)

data WPoint p = WPoint
              { weigth :: Double
              , point  :: p
              } deriving (Show,Eq)
                       
-- | Create a structure for face (Triangles in the case of 3D DT) and store the orientation
--   of the face related to the previous generated simplex that create the face.

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                
instance (PointND a) => Ord (S0 a) where
  compare = compS0
instance (PointND a) => Eq (S0 a) where
  x == y = compS0 x y == EQ

instance (PointND a) => Ord (S1 a) where
  compare = compS1
instance (PointND a) => Eq (S1 a) where
  x == y = compS1 x y == EQ

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- Instances for edge
compEdge a1 a2 b1 b2
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

compFace a b = compare a' b'
  where
    a' = fast3DSort a
    b' = fast3DSort b
    fast3DSort face@(a, b, c)
      | (a >= b) && (b >= c) = face
      | otherwise = (a', b', c')
        where
          minab = min a b
          maxab = max a b
          a'    = max (maxab) c
          b'    = max (min (maxab) c) (minab)
          c'    = min (minab) c

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- ActiveSubUnit
data (Buildable simplex dim) => ActiveSubUnit simplex dim = ActiveUnit
    { activeUnit :: (Sub simplex) dim
    , assocP     :: PointPointer
    , assocND    :: Point3D
    }

instance (Ord (Sub simplex dim), Buildable simplex dim) => Ord (ActiveSubUnit simplex dim) where
  compare a b = compare e1 e2
    where
      e1 = activeUnit a
      e2 = activeUnit b

instance (Ord (Sub simplex dim), Buildable simplex dim) => Eq (ActiveSubUnit simplex dim) where
  a == b = compare e1 e2 == EQ
    where
      e1 = activeUnit a
      e2 = activeUnit b

-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- | Define the a plane that dissect the space in two half-space.
--   All the elements of DT will be constructed recursively over this plane.
--   Plane is represented by a vector from oring (0,0,0) and a scalar such that k*(a,b,c)
--   is the closest point in the plane to the oring

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

data BoxPair a = BoxPair
    { halfBox1::Box a
    , halfBox2::Box a
    }

type SetPoint a = DiffArray PointPointer (WPoint a)

(!.)::SetPoint a -> PointPointer -> a
sP !. ix = point $ sP ! ix 


-- | Group the data that must be update along the computation (State).
--   Use of state monad will make it clear and keep the purity of the code.
data StateVarsMBC simplex dim = StateVarsMBC
    { aflAlpha, aflBox1, aflBox2 :: Set (ActiveSubUnit simplex dim)
    , externalFaces              :: Set (ActiveSubUnit simplex dim)
    , count                      :: Int
    , setPoint                   :: SetPoint dim
    }

type SetActiveSubUnits simplex dim = Set (ActiveSubUnit simplex dim)
type StateMBC simplex dim          = State (StateVarsMBC simplex dim)


class (PointND dim) => Buildable simplex dim where
  type Sub simplex  :: * -> *
  buildUnit         :: ActiveSubUnit simplex dim -> SetPoint dim -> [PointPointer] -> Maybe (simplex dim)
  build1stUnit      :: Plane dim -> SetPoint dim -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (simplex dim)
  getAllSubUnits    :: Maybe (ActiveSubUnit simplex dim) -> SetPoint dim -> (simplex dim) -> [ActiveSubUnit simplex dim]
  subUnitPos        :: BoxPair dim -> SetPoint dim -> ActiveSubUnit simplex dim -> Position


class (Vector p, DotProd p, CrossProd p, Show p) => PointND p where
  data Box p   :: *
  data Plane p :: *
  data S0 p    :: *
  data S1 p    :: *
  data S2 p    :: *
  compS0       :: S0 p -> S0 p -> Ordering
  compS1       :: S1 p -> S1 p -> Ordering
  isInBox      :: Box p -> p -> Bool
  calcPlane    :: SetPoint p -> S1 p -> Maybe (Plane p)
  planeNormal  :: Plane p -> p
  planeDist    :: Plane p -> Double
  makePlane    :: p -> Double -> Plane p
  touchPlane   :: p -> Plane p -> p -> p -> Maybe (Plane p)
  cutBox       :: Box p -> [Position] -> (Plane p, BoxPair p)
