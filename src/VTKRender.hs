

{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}

module VTKRender
       ( writeVTKfile
       ) where

import qualified Data.Vector as Vec

import Hammer.Math.Vector
import Hammer.Render.VTK.VTKRender
import Hammer.Render.VTK.Base

import DeUni.Types
import DeUni.Dim3.Base3D
import DeUni.Dim2.Base2D
import DeUni.Dim2.Delaunay2D
import DeUni.Dim3.Delaunay3D
import DeUni.Dim3.Hull3D

writeVTKfile file ps cs = writeUniVTKfile (text2Path file) $ mkUGVTK "RegularTriangulation" ps cs

type PointPointer = Int

instance RenderCell (S1 Point3D) where
  makeCell cell = let (a,b,c) = face3DPoints cell in Vec.fromList [a,b,c]
  getType _ = VTK_TRIANGLE

instance RenderCell (S2 Point3D) where
  makeCell cell = let (a,b,c,d) = tetraPoints cell in Vec.fromList [a,b,c,d]
  getType _ = VTK_TETRA
              
instance RenderCell (S2 Point2D) where
  makeCell cell = let (a,b,c) = face2DPoints cell in Vec.fromList [a,b,c]
  getType _ = VTK_TRIANGLE

instance (RenderPoint a) => RenderPoint (WPoint a) where
  renderPoint = renderPoint . point






