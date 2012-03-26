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


import DeUni.DeWall
import Criterion.Main
import RandomGenerator
import qualified Data.Vector as Vec
import Data.Random

main = let
  vol  = 20
  var  = 1
  anis = (1,1,1)
  dist n gen = genFullRandomGrainDistribution gen n vol (Data.Random.normal 15 1) anis

  test (DistributedPoints box arr) = let
    len  = Vec.length arr
    pset = [0..len-1]
    in runDelaunay3D box arr pset
  
  in do
     gen <- getRandomGen (Just 2908)
     dist1 <- dist 1000 gen
     -- dist2 <- dist 5000 gen
     defaultMain [
       bgroup "fib" [ bench "1000" $ whnf test dist1
                    -- , bench "5000" $ whnf test dist2
                    ]
       ]

