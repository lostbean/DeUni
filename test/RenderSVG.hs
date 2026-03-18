{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

module RenderSVG where

import qualified Data.ByteString.Lazy as BS
import qualified Data.IntMap as IM
import qualified Data.List as L
import qualified Data.Map as Map
import qualified Data.Set as S
import qualified Data.Text as T
import qualified Data.Text.Lazy.Encoding as TL
import qualified Data.Vector as Vec

import Data.Colour (AlphaColour, withOpacity)
import Data.IntMap (IntMap)
import Data.Map (Map)
import Data.Set (Set)
import Data.Vector (Vector)

import Diagrams.Backend.SVG
import Diagrams.Prelude hiding (Box, Vector, height, interval, view, width)
import Graphics.Svg.Core (renderText)

import DeUni.DeWall
import Linear.Vect

sizeSpec (width, height) = case (width, height) of
    (Nothing, Nothing) -> mkSizeSpec2D Nothing Nothing
    (Just w, Nothing) -> mkSizeSpec2D (Just $ fromIntegral w) Nothing
    (Nothing, Just h) -> mkSizeSpec2D Nothing (Just $ fromIntegral h)
    (Just w, Just h) -> mkSizeSpec2D (Just $ fromIntegral w) (Just $ fromIntegral h)

renderSVG :: String -> SizeSpec V2 Double -> Diagram SVG -> IO ()
renderSVG fileName szSpec dia =
    let
        build = renderDia SVG (SVGOptions szSpec Nothing "" [] True) dia
     in
        BS.writeFile fileName (TL.encodeUtf8 $ renderText build)

closeUpOnBox :: Box Point2D -> Diagram SVG -> Diagram SVG
closeUpOnBox Box2D{..} dia =
    let
        dx = xMax2D - xMin2D
        dy = yMax2D - yMin2D
        boxCenter = p2 ((xMax2D + xMin2D) / 2, (yMax2D + yMin2D) / 2)
        clipRect = rect dx dy # moveTo boxCenter
     in
        dia # clipBy clipRect

renderBox2D :: Box Point2D -> Diagram SVG
renderBox2D Box2D{..} =
    let
        dx = abs (xMax2D - xMin2D)
        dy = abs (yMax2D - yMin2D)
        boxSize = r2 (xMin2D + dx / 2, yMin2D + dy / 2)
     in
        rect dx dy
            # translate boxSize
            # lc blue
            # lw 0.05

renderSetS1 :: Vector (WPoint Point2D) -> Set (S1 Point2D) -> Diagram SVG
renderSetS1 ps ss =
    let
        func acc s1 =
            let
                a = edge2DR s1
                b = edge2DL s1
             in
                acc <> seg (ps !. a) (ps !. b)
        seg a b =
            let
                sg = v2p a ~~ v2p b
             in
                strokeT sg
                    # lw 0.3
                    # lc red
                    # translate (v2r a)
     in
        S.foldl' func mempty ss

renderSetPair :: Vector (WPoint Point2D) -> [(Int, Int)] -> Diagram SVG
renderSetPair ps ss =
    let
        func acc (a, b) = acc <> seg (ps !. a) (ps !. b)
        seg a b =
            let
                sg = v2p a ~~ v2p b
             in
                strokeT sg
                    # lw 0.6
                    # lc green
                    # translate (v2r a)
     in
        L.foldl' func mempty ss

renderSetS2 :: Vector (WPoint Point2D) -> IntMap (S2 Point2D) -> Diagram SVG
renderSetS2 ps ss =
    let
        func acc x =
            let
                (a, b, c) = face2DPoints x
             in
                acc <> renderTri (ps !. a) (ps !. b) (ps !. c)
     in
        IM.foldl func mempty ss

renderSetS2Triangle :: Vector (WPoint Point2D) -> IntMap (S2 Point2D) -> Diagram SVG
renderSetS2Triangle ps ss =
    let
        func acc x =
            let
                (a, b, c) = face2DPoints x
             in
                acc <> renderTri (ps !. a) (ps !. b) (ps !. c)
     in
        IM.foldl func mempty ss

renderSetS2Circle :: IntMap (S2 Point2D) -> Diagram SVG
renderSetS2Circle ss =
    let
        func acc x =
            acc
                <> renderCircle (circleCenter x) (circleRadius x ** 0.5) (green `withOpacity` 0.15)
                    # fcA (green `withOpacity` 0.15)
     in
        IM.foldl func mempty ss

renderTri :: Vec2D -> Vec2D -> Vec2D -> Diagram SVG
renderTri a b c =
    let
        ab = v2p a ~~ v2p b
        bc = v2p b ~~ v2p c
        ca = v2p c ~~ v2p a
        tri = ab <> bc <> ca
     in
        strokeT tri
            # fcA (yellow `withOpacity` 0.30)
            # lw 0.1
            # lc orange
            # translate (v2r a)

renderSetPoint2D :: Vector (WPoint Point2D) -> Diagram SVG
renderSetPoint2D ps = Vec.foldl' (\acc x -> acc <> renderCircle (point x) (weight x) (red `withOpacity` 0.15)) mempty ps

v2r (Vec2 x y) = r2 (x, y)

v2p (Vec2 x y) = p2 (x, y)

renderCircle :: Vec2D -> Double -> AlphaColour Double -> Diagram SVG
renderCircle pt diameter color =
    let
        v = v2r pt
     in
        circle diameter
            # lc black
            # fcA color
            # translate v
            <> circle 0.2
                # fc black
                # translate v
