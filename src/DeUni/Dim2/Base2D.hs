{-# LANGUAGE
    FlexibleInstances
  , RecordWildCards
  , TypeFamilies
  #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
module DeUni.Dim2.Base2D where

import Linear.Vect

import DeUni.Types

instance PointND Vec2 where
  data Box Vec2 = Box2D
    { xMax2D :: Double
    , xMin2D :: Double
    , yMax2D :: Double
    , yMin2D :: Double
    } deriving (Show)

  data Plane Vec2 = Plane2D
    { plane2DNormal :: Vec2D
    , plane2DDist   :: Double
    } deriving (Show)

  data S0 Vec2 = P2D
    { point2D :: PointPointer
    } deriving (Show)

  data S1 Vec2 = Edge2D
    { edge2DL  :: PointPointer
    , edge2DR  :: PointPointer
    } deriving (Show)

  data S2 Vec2 = Face2D
    { face2DPoints :: (PointPointer, PointPointer, PointPointer)
    , circleRadius :: Double
    , circleCenter :: Vec2D
    } deriving (Show)

  compS0 _ _ = error "[DeUni] 2D S0 comparison not defined."

  compS1 a b = compEdge (edge2DL a) (edge2DR a) (edge2DL b) (edge2DR b)

  circumOrigin = circleCenter

  circumRadius = circleRadius

  isInBox box (Vec2 x y) = let
    between minV maxV v
      -- will get points on the edge of the box and store if P1
      -- those are on the commun face
      | minV < maxV = (v >= minV) && (maxV >= v)
      | minV > maxV = (v <= minV) && (maxV <= v)
      | otherwise = error ("Zero size box: " ++ show box)
    in between (xMin2D box) (xMax2D box) x &&
       between (yMin2D box) (yMax2D box) y

  planeNormal = plane2DNormal
  planeDist   = plane2DDist
  makePlane   = Plane2D

  calcPlane sp edge = let
    a = edge2DL edge
    b = edge2DR edge
    in plane2D (sp!.a) (sp!.b)

  touchPlane _ _ = plane2D

  cutBox box subB
    | null subB = smartBox box box
    | otherwise = func box subB
    where
      func sub [] = smartBox box sub
      func sub (p:ps) = case p of
        B1 -> func (halfBox1.snd $ smartBox sub sub) ps
        B2 -> func (halfBox2.snd $ smartBox sub sub) ps
        _  -> func sub ps

instance Show (BoxPair Vec2) where
  show b = show (halfBox1 b, halfBox2 b)

smartBox :: Box Vec2 -> Box Vec2 -> (Plane Vec2, BoxPair Vec2)
smartBox box Box2D{..}
  | deltaX >= deltaY = (Plane2D (Vec2 1 0) halfX, cutX)
  | otherwise        = (Plane2D (Vec2 0 1) halfY, cutY)
  where
    cutX = BoxPair box { xMax2D = halfX } box { xMin2D = halfX }
    cutY = BoxPair box { yMax2D = halfY } box { yMin2D = halfY }
    deltaX = abs (xMax2D - xMin2D)
    deltaY = abs (yMax2D - yMin2D)
    halfX = (xMax2D + xMin2D)/2
    halfY = (yMax2D + yMin2D)/2


plane2D :: Vec2D -> Vec2D -> Maybe (Plane Vec2)
plane2D a b
  | nSize == 0 = Nothing
  | d >= 0     = Just $ makePlane normN d
  | otherwise  = Just $ makePlane (neg normN) (-d)
  where
    (Vec2 x y) = b &- a
    n          = Vec2 (-y) x
    nSize = vlen n
    -- Double normalization to avoid floating point operations
    -- errors in some computers. Critical in case of multiple
    -- points algined in a plane e.g. on a face of the box
    normN = (normalize . normalize) n
    d     = normN &. a
