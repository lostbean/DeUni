-----------------------------------------------------------------------------
--
-- Module      :  Main
-- Copyright   :
-- License     :  AllRightsReserved
--
-- Maintainer  :
-- Stability   :
-- Portability :
--
-- |
--
-----------------------------------------------------------------------------

module Main ( main ) where

import Check3D as C3D
import Check2D as C2D
import CheckCommon as C2C

import Test.QuickCheck
import Data.Array.Diff

import DeUni.DeWall
import Math.Vector
import DeUni.Types
import DeUni.Dim3.Base3D
import DeUni.GeometricTools


main = do
  C2C.runChecker
  print "=========== Common test finished =================="
  C2D.runChecker
  print "=========== 2D test finished ======================"
  C3D.runChecker
  print "=========== 3D test finished ======================"
