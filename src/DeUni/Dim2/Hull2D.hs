{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE NamedFieldPuns #-}
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%| Hull3D |%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeSynonymInstances #-}

module DeUni.Dim3.Hull2D where

import Control.Applicative ((<$>))
import Data.Array.Diff hiding (elems)
import Data.List (foldl', map)
import Data.Maybe
import Prelude hiding (lookup, null)

import Linear.Vect

import DeUni.Dim3.Base2D
import DeUni.FirstSeed
import DeUni.GeometricTools
import DeUni.Types

{- | It is not possible to implement convex hull 2D due the subUints
are point and therefor can't cross the alpha plane.
-}

{--
instance Buildable S1 Point2D where
  type Sub S1    = S0
  buildUnit      = makeFace
  build1stUnit   = makeFirstFace
  getAllSubUnits = extractAllFaceEdges
  subUnitPos     = edge3DPos
--}
