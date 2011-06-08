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

import Math.DeUni
import DeUniChecker
import Test.QuickCheck
import Data.Vec hiding (length)


main = do
  let myArgs = Args {replay = Nothing, maxSuccess = 1000, maxDiscard = 5000, maxSize = 1000, chatty = True}
  --quickCheckWith myArgs prop_1stFace
  --quickCheckWith myArgs prop_ConvHull

  --quickCheckWith myArgs prop_CircumSphere
  --quickCheckWith myArgs prop_1stSimplex
  let box = Box {xMax = 167.8502429742366, xMin = -117.74551048874855, yMax = 167.8502429742366, yMin = -117.74551048874855, zMax = 167.8502429742366, zMin = -117.74551048874855}
  let ps = [Vec3D (97.39372693002224) (73.40617212466896) (2.6973414700478315) , Vec3D (58.095817361027) (-98.87312105856836) (82.42319310083985) , Vec3D (-94.2528698593378) (-99.21665950678289) (-59.45456698536873) , Vec3D (12.748646875843406) (66.25062758103013) (67.8034599404782) , Vec3D (75.50557223148644) (-9.201711229979992) (-23.980601178482175)]
  let (sim,i) = runDeWall box ps
  print $ (length sim,i)
  quickCheckWith myArgs prop_Delaunay