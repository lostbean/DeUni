-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Hull3D |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE TypeFamilies #-}

module DeUni.Dim3.Hull2D where

import Prelude hiding (null, lookup)
import Data.List (map, foldl')
import Control.Applicative ((<$>))
import Data.Maybe
import Data.Array.Diff hiding (elems)

import Hammer.Math.Vector

import DeUni.GeometricTools
import DeUni.Types
import DeUni.FirstSeed
import DeUni.Dim3.Base2D


-- | It is not possible to implement convex hull 2D due the subUints
-- are point and therefor can't cross the alpha plane.

{--
instance Buildable S1 Point2D where
  type Sub S1    = S0
  buildUnit      = makeFace
  build1stUnit   = makeFirstFace
  getAllSubUnits = extractAllFaceEdges
  subUnitPos     = edge3DPos 
--}
 
