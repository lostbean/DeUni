

{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE MultiParamTypeClasses #-}



module VTKGenSimplex
       ( renderVTK
       , writeVTKfile
       ) where

import Control.Monad.State.Lazy
import Data.Array.Diff (listArray, (!), elems, bounds, (//), DiffArray)

import Data.IntMap (IntMap)
import qualified Data.IntMap as IM

import Text.XML
import qualified Text.XML as X

import Data.Text (Text, pack, append)
import qualified Data.Text as T

import Math.Vector
import DeUni.Types
import DeUni.Dim3.Base3D
import DeUni.Dim3.Delaunay3D
import DeUni.Dim3.Hull3D
import TemplateVTKXML


type Render = State VTKRender
type SetSimplex = IntMap (S2 Point3D)

class (Buildable slx dim, RenderPoints dim) => RenderVTK slx dim where
  render::(slx dim) -> Render ()
  
class RenderPoints dim where
  point2text::Maybe dim -> Text

data VTKRender = VTKRender
                 { cellConn     :: [Int]
                 , cellOffset   :: [Int]
                 , cellPointer  :: Int
                 , cellType     :: [Int] }
                 deriving (Show)
                       
initState = VTKRender
                 { cellConn     = []
                 , cellOffset   = []
                 , cellPointer  = 0
                 , cellType     = [] }
  

writeVTKfile::(RenderVTK slx dim)=> FilePath -> SetPoint dim -> IntMap (slx dim) -> IO ()
writeVTKfile name arr gs = X.writeFile X.def name (renderVTK arr gs)

renderVTK::(RenderVTK slx dim)=> SetPoint dim -> IntMap (slx dim) -> Document
renderVTK arr x = renderVTKDoc np nc dataPoint dataCell p c
  where
    ss        = IM.toList x
    st        = execState (mapM_ (render.snd) ss) (initState)
    (lb,hb)   = bounds arr
    ps        = [lb..hb]
    np        = length ps
    nc        = IM.size x
    dataPoint = (renderScalarPointData "Pscalar" $ map fromIntegral [0..hb]):[]
    dataCell  = (renderScalarCellData "ID" $ map (fromIntegral.fst) ss):[]
    
    p = renderPoints (map (\i->if i >= lb && i <= hb
                               then Just $ arr!.i
                               else Nothing ) [0..hb])
    
    c = renderCells (cellConn st) (cellOffset st) (cellType st) [] []

renderPoints::(RenderPoints dim)=> [Maybe dim] -> Element
renderPoints points = 
  Element {
    elementName = Name {nameLocalName = "Points", nameNamespace = Nothing, namePrefix = Nothing},
    elementAttributes = [],
    elementNodes = [
      NodeContent "\n",
      NodeElement (
        Element {
           elementName = Name {nameLocalName = "DataArray", nameNamespace = Nothing, namePrefix = Nothing},
           elementAttributes = [
             (Name {nameLocalName = "type", nameNamespace = Nothing, namePrefix = Nothing}, "Float32"),
             (Name {nameLocalName = "Name", nameNamespace = Nothing, namePrefix = Nothing}, "Points"),
             (Name {nameLocalName = "NumberOfComponents", nameNamespace = Nothing, namePrefix = Nothing}, "3"),
             (Name {nameLocalName = "format", nameNamespace = Nothing, namePrefix = Nothing}, "ascii")],
           elementNodes = [
             NodeContent (T.unwords $ map point2text points)] }),
      NodeContent "\n"]}


instance RenderPoints Point3D where
  point2text x = case x of
    Just (Vec3 x y z) -> T.unwords [toTxt x, toTxt y, toTxt z]
    Nothing           -> pack "0.0 0.0 0.0"      
      
instance RenderVTK S1 Point3D where
  render face = do
    st <- get
    let 
      (a,b,c) = face3DPoints face
      ps = map fromIntegral [a,b,c]
    modify (\x -> x { cellConn     = ps ++ (cellConn x) })
    modify (\x -> x { cellOffset   = (cellOffset x) ++ [ (cellPointer x) + 3 ] })
    modify (\x -> x { cellPointer  = cellPointer x + 3 })
    modify (\x -> x { cellType     = 5:(cellType x) })

instance RenderVTK S2 Point3D where
  render sigma = do
    st <- get
    let 
      (a,b,c,d) = tetraPoints sigma
      ps = map fromIntegral [a,b,c,d]
    modify (\x -> x { cellConn     = ps ++ (cellConn x) })
    modify (\x -> x { cellOffset   = (cellOffset x) ++ [ (cellPointer x) + 4 ] })
    modify (\x -> x { cellPointer  = cellPointer x + 4 })
    modify (\x -> x { cellType     = 10:(cellType x) })






