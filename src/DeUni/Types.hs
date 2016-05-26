{-# LANGUAGE
    FlexibleContexts
  , TypeFamilies
  , MultiParamTypeClasses
  , UndecidableInstances
  #-}
module DeUni.Types where

import Prelude
import Data.Set                 (Set)
import Data.Vector              (Vector, (!))
import Control.Monad.State.Lazy (State)

import Linear.Vect

-- =================================================================================================

-- | Define a point in 2D and 3D
type PointPointer = Int

-- | Create a structure for face (Triangles in the case of 3D DT) and store the orientation
--   of the face related to the previous generated simplex that create the face.
data WPoint p
  = WPoint
  { weight :: Double
  , point  :: p Double
  }

instance (Show (v Double)) => Show (WPoint v) where
  show (WPoint w p) = "(" ++ show w ++ ")" ++ show p

instance (Eq (v Double)) => Eq (WPoint v) where
  (WPoint w1 p1) == (WPoint w2 p2) = p1 == p2 && w1 == w2

-- =================================================================================================

instance (PointND a) => Ord (S0 a) where
  compare = compS0
instance (PointND a) => Eq (S0 a) where
  x == y = compS0 x y == EQ

instance (PointND a) => Ord (S1 a) where
  compare = compS1
instance (PointND a) => Eq (S1 a) where
  x == y = compS1 x y == EQ

-- =================================================================================================

compEdge :: (Ord a)=> a -> a -> a -> a -> Ordering
compEdge a1 a2 b1 b2
  | amax > bmax = GT
  | amax < bmax = LT
  | amin > bmin = GT
  | amin < bmin = LT
  | otherwise   = EQ
  where
    amax = eUp  a1 a2
    amin = eBot a1 a2
    bmax = eUp  b1 b2
    bmin = eBot b1 b2
    eUp  x y = if compare x y == GT then x else y
    eBot x y = if compare x y == GT then y else x

compFace :: (Ord a)=> (a, a, a) -> (a, a, a) -> Ordering
compFace a b = compare a' b'
  where
    a' = fast3DSort a
    b' = fast3DSort b

fast3DSort :: (Ord a)=> (a, a, a) -> (a, a, a)
fast3DSort face@(a, b, c)
  | (a >= b) && (b >= c) = face
  | otherwise            = (a', b', c')
  where
    minab = min a b
    maxab = max a b
    a'    = max maxab c
    b'    = max (min maxab c) minab
    c'    = min minab c

-- =================================================================================================
-- ActiveSubUnit
data ActiveSubUnit simplex v = ActiveUnit
    { activeUnit :: (Sub simplex) v
    , assocP     :: PointPointer
    , assocND    :: v Double
    }

instance (Ord (Sub simplex v), Buildable simplex v) => Ord (ActiveSubUnit simplex v) where
  compare a b = compare e1 e2
    where
      e1 = activeUnit a
      e2 = activeUnit b

instance (Ord (Sub simplex v), Buildable simplex v) => Eq (ActiveSubUnit simplex v) where
  a == b = e1 == e2
    where
      e1 = activeUnit a
      e2 = activeUnit b

-- =================================================================================================
-- | Define the a plane that dissect the space in two half-space.
--   All the elements of DT will be constructed recursively over this plane.
--   Plane is represented by a vector from origin (0,0,0) and a scalar such that k*(a,b,c)
--   is the closest point in the plane to the origin

data PointPartition = PointPartition
    { pointsOnB1   :: [PointPointer]
    , pointsOnB2   :: [PointPointer]
    , pointsOnPlane:: [PointPointer]
    } deriving (Show)

-- | Define possible possitions of the elements for the 1st half-space (Box1=B1),
--   2nd (Box2=B2) and intersect by the plane (B1B2Plane).
data Position
  = B1
  | B2
  | OnPlane
  | CrossPlane
  | None
  deriving (Show, Eq)

data BoxPair a
  = BoxPair
  { halfBox1 :: Box a
  , halfBox2 :: Box a
  }

type SetPoint v = Vector (WPoint v)

(!.)::SetPoint a -> PointPointer -> a Double
sP !. ix = point $ sP ! ix


-- | Group the data that must be update along the computation (State).
--   Use of state monad will make it clear and keep the purity of the code.
data StateVarsMBC simplex v
  = StateVarsMBC
  { aflAlpha, aflBox1, aflBox2 :: Set (ActiveSubUnit simplex v)
  , externalFaces              :: Set (ActiveSubUnit simplex v)
  , count                      :: Int
  , setPoint                   :: SetPoint v
  }

type SetActiveSubUnits simplex v = Set   (ActiveSubUnit simplex v)
type StateMBC          simplex v = State (StateVarsMBC  simplex v)


class (PointND v) => Buildable simplex v where
  type Sub simplex :: (* -> *) -> *
  buildUnit        :: ActiveSubUnit simplex v -> SetPoint v -> [PointPointer] -> Maybe (simplex v)
  build1stUnit     :: Plane v -> SetPoint v -> [PointPointer] -> [PointPointer] -> [PointPointer] -> Maybe (simplex v)
  getAllSubUnits   :: SetPoint v -> simplex v -> [ActiveSubUnit simplex v]
  subUnitPos       :: BoxPair v -> SetPoint v -> ActiveSubUnit simplex v -> Position


class (LinearMap Double p, DotProd Double p, Show (p Double)) => PointND p where
  data Box p   :: *
  data Plane p :: *
  data S0 p    :: *
  data S1 p    :: *
  data S2 p    :: *
  compS0       :: S0 p -> S0 p -> Ordering
  compS1       :: S1 p -> S1 p -> Ordering
  circumOrigin :: S2 p -> p Double
  circumRadius :: S2 p -> Double
  isInBox      :: Box p -> p Double -> Bool
  calcPlane    :: SetPoint p -> S1 p -> Maybe (Plane p)
  planeNormal  :: Plane p -> p Double
  planeDist    :: Plane p -> Double
  makePlane    :: p Double -> Double -> Plane p
  touchPlane   :: p Double -> Plane p -> p Double -> p Double -> Maybe (Plane p)
  cutBox       :: Box p -> [Position] -> (Plane p, BoxPair p)
