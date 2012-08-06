{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RecordWildCards #-}

module RenderSVG where

import qualified Data.ByteString.Lazy as BS
import qualified Blaze.ByteString.Builder as B
import qualified Data.Vector as Vec
import qualified Data.List as L
import Data.Vector (Vector)
import Data.IntMap (IntMap)
import qualified Data.IntMap as IM
import Data.Map (Map)
import qualified Data.Map as Map
import Data.Set (Set)
import qualified Data.Set as S
import Data.Colour (AlphaColour,withOpacity)

import Diagrams.Prelude hiding (width, height, interval)
import Diagrams.Backend.SVG

import DeUni.Types
import DeUni.Dim2.Base2D
import Hammer.Math.Vector hiding (Vector)


sizeSpec (width, height) = case (width, height) of
  (Nothing, Nothing) -> Absolute
  (Just w, Nothing)  -> Width (fromIntegral w)
  (Nothing, Just h)  -> Height (fromIntegral h)
  (Just w, Just h)   -> Dims (fromIntegral w) (fromIntegral h)

renderSVG :: String -> SizeSpec2D -> Diagram SVG R2 -> IO ()
renderSVG fileName sizeSpec dia = let
  build = renderDia SVG (SVGOptions fileName sizeSpec) dia
  in BS.writeFile fileName (B.toLazyByteString build)
                      
closeUpOnBox::Box Point2D -> Diagram SVG R2 -> Diagram SVG R2
closeUpOnBox Box2D{..} = let
  r = r2 (xMax2D - xMin2D, yMax2D - yMin2D)
  p = p2 (xMin2D, yMin2D)
  in view p r

renderBox2D::Box Point2D -> Diagram SVG R2
renderBox2D Box2D{..} = let
  dx = abs (xMax2D - xMin2D)
  dy = abs (yMax2D - yMin2D)
  boxSize = r2 (xMin2D + dx/2, yMin2D + dy/2)
  in rect dx dy
     # translate boxSize
     # lc blue
     # lw 0.05

renderSetS1::Vector (WPoint Point2D) -> Set ( S1 Point2D) -> Diagram SVG R2
renderSetS1 ps ss = let
  func acc s1 = let
    a = edge2DR s1
    b = edge2DL s1
    in acc <> seg (ps!.a) (ps!.b)
  seg a b = let
    seg = v2p a ~~ v2p b
    in strokeT seg
       # lw 0.3
       # lc red
       # translate (v2r a)
  in S.foldl' func mempty ss
     
renderSetPair::Vector (WPoint Point2D) -> [(Int, Int)] -> Diagram SVG R2
renderSetPair ps ss = let
  func acc (a, b) = acc <> seg (ps!.a) (ps!.b)
  seg a b = let
    seg = v2p a ~~ v2p b
    in strokeT seg
       # lw 0.6
       # lc green
       # translate (v2r a)
  in L.foldl' func mempty ss
     
renderSetS2::Vector (WPoint Point2D) -> IntMap (S2 Point2D) -> Diagram SVG R2
renderSetS2 ps ss = let
  func acc x = let
    (a,b,c) = face2DPoints x
    in acc <> renderTri (ps!.a) (ps!.b) (ps!.c)
  in IM.foldl func mempty ss
     
renderSetS2Triangle::Vector (WPoint Point2D) -> IntMap (S2 Point2D) -> Diagram SVG R2
renderSetS2Triangle ps ss = let
  func acc x = let
    (a,b,c) = face2DPoints x
    in acc <> renderTri (ps!.a) (ps!.b) (ps!.c)  
  in IM.foldl func mempty ss

renderSetS2Circle::IntMap (S2 Point2D) -> Diagram SVG R2
renderSetS2Circle ss = let
  func acc x = acc <> renderCircle (circleCenter x) (circleRadius x ** 0.5) (green `withOpacity` 0.15)
                    # fcA (green `withOpacity` 0.15)
  in IM.foldl func mempty ss

renderTri::Point2D -> Point2D -> Point2D -> Diagram SVG R2
renderTri a b c = let
  ab = v2p a ~~ v2p b
  bc = v2p b ~~ v2p c
  ca = v2p c ~~ v2p a
  tri = ab <> bc <> ca
  in strokeT tri
     # fcA (yellow `withOpacity` 0.30)
     # lw 0.1
     # lc orange
     # translate (v2r a)

renderSetPoint2D::Vector (WPoint Point2D) -> Diagram SVG R2
renderSetPoint2D ps = Vec.foldl' (\acc x -> acc <> renderCircle (point x) (weight x) (red `withOpacity` 0.15)) mempty ps

v2r (Vec2 x y) = r2 (x,y)

v2p (Vec2 x y) = p2 (x,y)


renderCircle::Point2D -> Double -> AlphaColour Double -> Diagram SVG R2
renderCircle point diameter color = let
  v = v2r point                      
  in circle diameter
     # lc black
     # fcA color
     # translate v 
  <> circle 0.2
     # fc black
     # translate v
     
