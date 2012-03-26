
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE RecordWildCards #-}

module RandomGenerator
( genFullRandomGrainDistribution
, genPoint
, getRandomGen
, DistributedPoints(..)
, calcBoxSize
) where

import Control.Monad (replicateM, liftM, foldM)
import Control.Monad.Loops (iterateUntil)
import Control.Monad.ST (runST)
import Data.Vector (Vector, (!), fromList)
import Data.IORef
import Data.Random
import Data.Random.RVar
import Data.Random.Source.StdGen
import System.Random.Mersenne.Pure64

import qualified Hammer.Math.Vector as AlgLin
import Hammer.Math.Vector hiding (Vector)
import DeUni.Types
import DeUni.Dim3.Base3D
import DeUni.Dim2.Base2D


data DistributedPoints a =
  DistributedPoints { box      ::Box a
                    , setPoint ::SetPoint a }

class (AlgLin.Vector v, AlgLin.Pointwise v)=> GenRandom v where
  type Ratio v                   :: *
  calcBoxSize                    :: Double -> Ratio v -> Box v
  boxDim                         :: Box v -> v
  genPoint                       :: IORef PureMT -> RVar Double -> IO v
       

instance GenRandom Point2D where
  type Ratio Point2D = (Double, Double)
  
  calcBoxSize avgVolume (xRatio, yRatio) = let
    refRatio = xRatio
    k1       = yRatio/refRatio
    a        = (avgVolume/k1)**(1/2)
    in Box2D {xMax2D = a, xMin2D = 0, yMax2D = a*k1, yMin2D = 0}

  boxDim Box2D{..} = let
    dx = (xMax2D - xMin2D)
    dy = (yMax2D - yMin2D)
    in Vec2 dx dy

  genPoint gen f = do
    -- To avoid repetition on the pseudo-random generator, use one external gen
    -- wrapped in an StateMonad. Or for internal gen use : "gen <- getRandomGen"
    a <- sampleFrom gen f
    b <- sampleFrom gen f
    return $ Vec2 a b

instance GenRandom Point3D where
  type Ratio Point3D = (Double, Double, Double)
  
  calcBoxSize avgVolume (xRatio, yRatio, zRatio) =
    let
    refRatio = xRatio
    k1       = yRatio/refRatio
    k2       = zRatio/refRatio
    a        = (avgVolume/(k1*k2))**(1/3)
    in Box3D {xMax3D=a, xMin3D=0, yMax3D=a*k1, yMin3D=0, zMax3D=a*k2, zMin3D=0}

  boxDim Box3D{..} = let
    dx = (xMax3D - xMin3D)
    dy = (yMax3D - yMin3D)
    dz = (zMax3D - zMin3D)
    in Vec3 dx dy dz

  genPoint gen f = do
    -- To avoid repetition on the pseudo-random generator, use one external gen
    -- wrapped in an StateMonad. Or for internal gen use : "gen <- getRandomGen"
    a <- sampleFrom gen f
    b <- sampleFrom gen f
    c <- sampleFrom gen f
    return $ Vec3 a b c


genFullRandomGrainDistribution::(GenRandom v)=> IORef PureMT -> Int -> Double -> RVar Double -> Ratio v -> IO (DistributedPoints v)
genFullRandomGrainDistribution gen targetNGrains avgVolume dist ratio = 
  let
    totalVolume = avgVolume*(fromIntegral targetNGrains)
    box         = calcBoxSize totalVolume ratio
    delta       = boxDim box
    getPoint    = do
      p  <- genPoint gen stdUniform
      gs <- genGrainSize gen dist
      return $ WPoint gs (delta &! p) 
    in replicateM targetNGrains getPoint >>= return.(DistributedPoints box).fromList


genGrainSize::IORef PureMT -> RVar Double -> IO Double
genGrainSize gen f = sampleFrom gen f

getRandomGen::Maybe Int -> IO (IORef PureMT)
getRandomGen x = case x of
    Nothing -> newPureMT >>= newIORef
    (Just seed) ->  (return $ pureMT (fromIntegral seed)) >>= newIORef
