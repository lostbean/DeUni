{-# LANGUAGE
    FlexibleContexts
  , RecordWildCards
  , TypeFamilies
  #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
module DeUni.Dim3.Base3D where

import qualified Data.List as L

import Data.List ((\\))
import Data.Vector ((!))

import Linear.Vect

import DeUni.GeometricTools
import DeUni.Types
import DeUni.Dim2.ReTri2D

instance PointND Vec3 where
  data Box Vec3 =  Box3D
    { xMax3D :: Double
    , xMin3D :: Double
    , yMax3D :: Double
    , yMin3D :: Double
    , zMax3D :: Double
    , zMin3D :: Double
    } deriving (Show, Eq)

  data Plane Vec3 = Plane3D
    { plane3DNormal :: Vec3D
    , plane3DDist   :: Double
    } deriving (Show, Eq)

  data S0 Vec3 = Edge3D
    { edge3DL :: PointPointer
    , edge3DR :: PointPointer
    } deriving (Show)

  data S1 Vec3 = Face3D
   { face3DPoints :: (PointPointer, PointPointer, PointPointer)
   } deriving (Show)

  data S2 Vec3 = Tetrahedron
   { circumSphereCenter :: Vec3D
   , circumSphereRadius :: Double
   , tetraPoints        :: (PointPointer, PointPointer, PointPointer, PointPointer)
   } deriving (Show, Eq)

  compS0 a b = compEdge (edge3DL a) (edge3DR a) (edge3DL b) (edge3DR b)

  compS1 a b = compFace (face3DPoints a) (face3DPoints b)

  circumOrigin = circumSphereCenter

  circumRadius = circumSphereRadius

  isInBox box (Vec3 x y z) = let
    between minV maxV v
      --will get points on the edge of the box and store if P1 those are on the commun face
      | minV < maxV = (v >= minV) && (maxV >= v)
      | minV > maxV = (v <= minV) && (maxV <= v)
      | otherwise = error ("Zero size box: " ++ show box)
    in between (xMin3D box) (xMax3D box) x
    && between (yMin3D box) (yMax3D box) y
    && between (zMin3D box) (zMax3D box) z


  planeNormal = plane3DNormal
  planeDist   = plane3DDist
  makePlane   = Plane3D

  calcPlane sp face
    | nSize == 0 = Nothing
    | d >= 0     = Just $ makePlane normN d
    | otherwise  = Just $ makePlane (inv normN) (-d)
    where
      (a,b,c) = face3DPoints face
      n       = (sp!.b &- sp!.a) &^ (sp!.c &- sp!.a)
      nSize   = vlen n
      -- Double normalization to avoid floating point operations errors in some computers
      -- Critical in case of multiple points algined in a plane e.g. on a face of the box
      normN   = (normalize . normalize) n
      d       = normN &. (sp!.a)
      inv     = ((-1) *&)

  touchPlane refdir divPlane a b
    | nSize == 0 = Nothing
    | d >= 0     = Just $ makePlane normND  d
    | otherwise  = Just $ makePlane (neg normND) (-d)
    where
      refOnDivPlane = getProjOnPlane divPlane refdir
      fixRotDir     = normalize $ (planeNormal divPlane) &^ refOnDivPlane
      nd            = (b &- a) &^ fixRotDir
      nSize         = norm nd
      normND        = normalize nd
      d             = normND &. a

  cutBox inBox subB
    | null subB = smartBox inBox inBox
    | otherwise = func     inBox subB
    where
      func sub [] = smartBox inBox sub
      func sub (p:ps) = case p of
        B1 -> func (halfBox1.snd $ smartBox sub sub) ps
        B2 -> func (halfBox2.snd $ smartBox sub sub) ps
        _  -> func sub ps
      smartBox box Box3D{..}
        | deltaX >= max deltaY deltaZ = (Plane3D (Vec3 1 0 0) halfX, cutX)
        | deltaY >= max deltaX deltaZ = (Plane3D (Vec3 0 1 0) halfY, cutY)
        | otherwise                   = (Plane3D (Vec3 0 0 1) halfZ, cutZ)
        where
          cutX = BoxPair box { xMax3D = halfX } box { xMin3D = halfX }
          cutY = BoxPair box { yMax3D = halfY } box { yMin3D = halfY }
          cutZ = BoxPair box { zMax3D = halfZ } box { zMin3D = halfZ }
          deltaX = abs (xMax3D - xMin3D)
          deltaY = abs (yMax3D - yMin3D)
          deltaZ = abs (zMax3D - zMin3D)
          halfX = (xMax3D + xMin3D)/2
          halfY = (yMax3D + yMin3D)/2
          halfZ = (zMax3D + zMin3D)/2

getThrirdPoint :: (PointND Vec3) => SetPoint Vec3 -> PointPointer -> PointPointer -> [PointPointer] -> Maybe (PointPointer, Vec3D)
getThrirdPoint sP pA pB ps = do
  (_, i) <- findThird
  getND i
  where
    cleanList x = ps \\ [pA, pB, x]
    justHull    = filter isHull ps

    findThird = let
      dist = getSignDist sP pA pB
      in findMinimunButZero' dist justHull

    isHull x
      | x == pA || x == pB = False
      | otherwise = let
        face = Face3D { face3DPoints = (pA, pB, x) }
        in case calcPlane sP face of
          Just plane
            | (L.null.pointsOnB1) pp -> True
            | (L.null.pointsOnB2) pp -> True
            | otherwise              -> False
            where pp = pointSetPartition (whichSideOfPlane plane) sP (cleanList x)
          Nothing                    -> False

    getND x = let
      face = Face3D { face3DPoints = (pA, pB, x) }
      func plane
        | (L.null.pointsOnB1) pp = return (x, nd)
        | otherwise              = return (x, neg nd)
        where pp = pointSetPartition (whichSideOfPlane plane) sP (cleanList x)
              nd = planeNormal plane
      in calcPlane sP face >>= func

-- | Rotate points from a plane with normal nd to a plane with normal = (0,0,1)
-- using axi-angle and Rodriges' equation.
rotate :: Vec3D -> Vec3D -> Vec3D
rotate nd x = let
  v    = Vec3 0 0 1
  w    = nd &^ v
  cosV = nd &. v
  sinV = sqrt (1 - cosV * cosV)
  in x &* cosV &+ (w &^ x) &* sinV &+ w &* ((w &. x)*(1 - cosV))

-- | Get the signed distance from c to the center of the edge <a,b> which are circumscribed
-- by a circle.
getSignDist :: SetPoint Vec3 -> PointPointer -> PointPointer -> PointPointer -> Maybe Double
getSignDist sp a b c = let
  face = Face3D { face3DPoints = (a, b, c) }
  getIn2D plane = let
    rot = rotate (planeNormal plane)
    -- rotate the face <a,b,c>, where a,b,c = R3, to a face with normal // (0,0,1)
    (Vec3 x1 y1 _) = rot (sp!.a)
    (Vec3 x2 y2 _) = rot (sp!.b)
    (Vec3 x3 y3 _) = rot (sp!.c)
    -- and cast the points to R2 by excluding their z values
    pA = (sp!a) {point = Vec2 x1 y1}
    pB = (sp!b) {point = Vec2 x2 y2}
    pC = (sp!c) {point = Vec2 x3 y3}
    -- find the circumcircle and the signed distance
    in fst $ getFaceDistCenter pA pB pC

  in getIn2D <$> calcPlane sp face
