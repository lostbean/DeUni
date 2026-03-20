{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# OPTIONS_GHC -Wno-orphans #-}

module VTKRender (
    writeVTKfile,
) where

import qualified Data.Vector.Unboxed as VU

import Hammer.VTK

import DeUni.Dim2.Base2D
import DeUni.Dim3.Base3D
import DeUni.Types

writeVTKfile :: (Monad m) => a -> b -> c -> m ()
writeVTKfile _ _ _ = return ()

instance RenderCell (S1 Point3D) where
    makeCell cell = let (a, b, c) = face3DPoints cell in VU.fromList [a, b, c]
    getType _ = VTK_TRIANGLE

instance RenderCell (S2 Point3D) where
    makeCell cell = let (a, b, c, d) = tetraPoints cell in VU.fromList [a, b, c, d]
    getType _ = VTK_TETRA

instance RenderCell (S2 Point2D) where
    makeCell cell = let (a, b, c) = face2DPoints cell in VU.fromList [a, b, c]
    getType _ = VTK_TRIANGLE
