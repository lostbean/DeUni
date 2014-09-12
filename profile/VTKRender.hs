

{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}

module VTKRender
       ( writeVTKfile
       ) where

import qualified Data.Vector.Unboxed as VU

import Hammer.Math.Algebra
import Hammer.VTK

import DeUni.Types
import DeUni.Dim3.Base3D
import DeUni.Dim2.Base2D
import DeUni.Dim2.Delaunay2D
import DeUni.Dim3.Delaunay3D
import DeUni.Dim3.Hull3D

writeVTKfile file ps cs = writeUniVTKfile file $
                          mkUGVTK "RegularTriangulation" ps cs

type PointPointer = Int

instance RenderCell (S1 Point3D) where
  makeCell cell = let (a,b,c) = face3DPoints cell in VU.fromList [a,b,c]
  getType _ = VTK_TRIANGLE

instance RenderCell (S2 Point3D) where
  makeCell cell = let (a,b,c,d) = tetraPoints cell in VU.fromList [a,b,c,d]
  getType _ = VTK_TETRA
              
instance RenderCell (S2 Point2D) where
  makeCell cell = let (a,b,c) = face2DPoints cell in VU.fromList [a,b,c]
  getType _ = VTK_TRIANGLE

--instance (RenderElemVTK a) => RenderElemVTK (WPoint a) where
--  renderPoint = renderPoint . point






