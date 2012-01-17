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

import DeUniChecker

import Test.QuickCheck
import Data.Array.Diff

import DeUni.DeWall
import Math.Vector
import DeUni.Types
import DeUni.Dim3.Base3D
import DeUni.GeometricTools


main = do
  let myArgs = Args {replay = Nothing, maxSuccess = 1000, maxDiscard = 5000, maxSize = 1000, chatty = True}
  
  print "Testing Projection.."    
  quickCheckWith myArgs prop_Projection
  
  print "Testing Parttition.."    
  quickCheckWith myArgs prop_partition
  
  print "Testing 1st face.."    
  quickCheckWith myArgs prop_1stFace
  
  print "Testing Circumsphere.."    
  quickCheckWith myArgs prop_CircumSphere
  
  print "Testing 1st tetrahedron.."    
  quickCheckWith myArgs prop_1stSimplex

  print "Testing Convex Hull.."    
  quickCheckWith myArgs prop_ConvHull
  
  print "Testing Delaunay.."    
  quickCheckWith myArgs prop_Delaunay
  


  --let (sim,i) = runDeWall box ps
  --print $ (length sim,i)

box = Box3D { xMax3D = 167.8502429742366, xMin3D = -117.74551048874855
                  , yMax3D = 167.8502429742366, yMin3D = -117.74551048874855
                  , zMax3D = 167.8502429742366, zMin3D = -117.74551048874855 }
vs = [ Vec3 (97.39372693002224) (73.40617212466896) (2.6973414700478315)
           , Vec3 (58.095817361027) (-98.87312105856836) (82.42319310083985)
           , Vec3 (-94.2528698593378) (-99.21665950678289) (-59.45456698536873)
           , Vec3 (12.748646875843406) (66.25062758103013) (67.8034599404782)
           , Vec3 (75.50557223148644) (-9.201711229979992) (-23.980601178482175) ]
sP = (\x -> listArray (0,fromIntegral (length x - 1)) x) vs :: DiffArray PointPointer Vec3