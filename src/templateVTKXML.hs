{-# LANGUAGE OverloadedStrings #-}
module TemplateVTKXML
       ( renderVTKDoc
       , renderScalarPointData
       , renderScalarCellData 
       , renderCells
       , toTxt ) where

import Data.List (intersperse)
import Data.Text (Text, pack, append)
import qualified Data.Text as T
import Text.XML


toTxt::(Show a)=>a -> Text
toTxt = pack.show

renderVTKDoc::Int -> Int -> [Element] -> [Element] -> Element -> Element -> Document
renderVTKDoc nPoints nCells pointData cellData points cells =
  Document {
    documentPrologue = Prologue {prologueBefore = [], prologueDoctype = Nothing, prologueAfter = []}, 
    documentRoot = Element {
      elementName = Name {nameLocalName = "VTKFile", nameNamespace = Nothing, namePrefix = Nothing}, 
      elementAttributes = [(Name {nameLocalName = "type", nameNamespace = Nothing, namePrefix = Nothing}, "UnstructuredGrid")], 
      elementNodes = [ 
        NodeContent "\n",
        NodeElement (
          Element {
             elementName = Name {nameLocalName = "UnstructuredGrid", nameNamespace = Nothing, namePrefix = Nothing},
             elementAttributes = [], 
             elementNodes = [
               NodeContent "\n",
               NodeElement (
                 Element {
                    elementName = Name {nameLocalName = "Piece", nameNamespace = Nothing, namePrefix = Nothing},
                    elementAttributes = [
                      -- Number of points and cells
                      (Name {nameLocalName = "NumberOfPoints", nameNamespace = Nothing, namePrefix = Nothing}, toTxt nPoints),
                      (Name {nameLocalName = "NumberOfCells", nameNamespace = Nothing, namePrefix = Nothing},toTxt nCells)], 
                    elementNodes =
                      -- Insert Point data
                      map NodeElement pointData
                      ++
                      -- Insert Cell Data
                      map NodeElement cellData
                      ++
                      -- Insert points
                      [NodeElement ( points ),
                      -- Insert cells
                      NodeElement ( cells )] }),
               NodeContent "\n\n"]}),
        NodeContent "\n"																
        ]}, documentEpilogue = []}


renderScalarPointData::String -> [Double] -> Element
renderScalarPointData name scalarData =
  Element {
    elementName = Name {nameLocalName = "PointData", nameNamespace = Nothing, namePrefix = Nothing},
    elementAttributes = [(Name {nameLocalName = "Scalars", nameNamespace = Nothing, namePrefix = Nothing},pack name)],
    elementNodes = [
      NodeContent "\n",
      NodeElement (
        Element {
           elementName = Name {nameLocalName = "DataArray", nameNamespace = Nothing, namePrefix = Nothing},
           elementAttributes = [
             (Name {nameLocalName = "type", nameNamespace = Nothing, namePrefix = Nothing}, "Float32"),
             (Name {nameLocalName = "Name", nameNamespace = Nothing, namePrefix = Nothing}, "my_scalars"),
             (Name {nameLocalName = "format", nameNamespace = Nothing, namePrefix = Nothing}, "ascii")],
           elementNodes = [
             NodeContent "\n\t\t",
             NodeContent (T.unwords $ map (pack.show) scalarData)]}),
      NodeContent "\n"]}
  
renderScalarCellData::String -> [Double] -> Element
renderScalarCellData name scalarData = 
  Element {
    elementName = Name {nameLocalName = "CellData", nameNamespace = Nothing, namePrefix = Nothing},
    elementAttributes = [(Name {nameLocalName = "Scalars", nameNamespace = Nothing, namePrefix = Nothing}, pack name)],
    elementNodes = [
      NodeContent "\n",
      NodeElement (
        Element {
           elementName = Name {nameLocalName = "DataArray", nameNamespace = Nothing, namePrefix = Nothing},
           elementAttributes = [
             (Name {nameLocalName = "type", nameNamespace = Nothing, namePrefix = Nothing}, "Float32"),
             (Name {nameLocalName = "Name", nameNamespace = Nothing, namePrefix = Nothing}, "scalars"),
             (Name {nameLocalName = "format", nameNamespace = Nothing, namePrefix = Nothing}, "ascii")],
           elementNodes = [
             NodeContent "\n\t\t",
             NodeContent (T.unwords $ map (pack.show) scalarData) ]}),
      NodeContent "\n"]}
  
renderCells::[Int] -> [Int] -> [Int] -> [Int] -> [Int] -> Element
renderCells cellConn cellOffsets cellTypes faces faceOffsets = Element {
    elementName = Name {nameLocalName = "Cells", nameNamespace = Nothing, namePrefix = Nothing},
    elementAttributes = [],
    elementNodes = intersperse (NodeContent "\n") full }
  where
    full = (fconn cellConn ++ foof cellOffsets ++ ftype cellTypes ++ fface faces ++ ffaceoff faceOffsets) 
    fconn [] = [] 
    fconn l = [NodeElement (
        Element {
           elementName = Name {nameLocalName = "DataArray", nameNamespace = Nothing, namePrefix = Nothing},
           elementAttributes = [
             (Name {nameLocalName = "type", nameNamespace = Nothing, namePrefix = Nothing}, "Int64"),
             (Name {nameLocalName = "Name", nameNamespace = Nothing, namePrefix = Nothing}, "connectivity"),
             (Name {nameLocalName = "format", nameNamespace = Nothing, namePrefix = Nothing}, "ascii")],
           elementNodes = [
             NodeContent (T.unwords $ map toTxt l)] })]
              
    foff [] = []
    foof l = [NodeElement (
        Element {
           elementName = Name {nameLocalName = "DataArray", nameNamespace = Nothing, namePrefix = Nothing},
           elementAttributes = [
             (Name {nameLocalName = "type", nameNamespace = Nothing, namePrefix = Nothing}, "Int64"),
             (Name {nameLocalName = "Name", nameNamespace = Nothing, namePrefix = Nothing}, "offsets"),
             (Name {nameLocalName = "format", nameNamespace = Nothing, namePrefix = Nothing}, "ascii")],
           elementNodes = [
             NodeContent (T.unwords $ map toTxt l)] })]

    ftype [] = []
    ftype l = [NodeElement (
        Element {
           elementName = Name {nameLocalName = "DataArray", nameNamespace = Nothing, namePrefix = Nothing},
           elementAttributes = [
             (Name {nameLocalName = "type", nameNamespace = Nothing, namePrefix = Nothing}, "UInt8"),
             (Name {nameLocalName = "Name", nameNamespace = Nothing, namePrefix = Nothing}, "types"),
             (Name {nameLocalName = "format", nameNamespace = Nothing, namePrefix = Nothing}, "ascii")],
           elementNodes = [
             NodeContent (T.unwords $ map toTxt l)] })]
      
    fface [] = []
    fface l = [NodeElement (
        Element {
           elementName = Name {nameLocalName = "DataArray", nameNamespace = Nothing, namePrefix = Nothing},
           elementAttributes = [
             (Name {nameLocalName = "type", nameNamespace = Nothing, namePrefix = Nothing}, "Int64"),
             (Name {nameLocalName = "Name", nameNamespace = Nothing, namePrefix = Nothing}, "faces"),
             (Name {nameLocalName = "format", nameNamespace = Nothing, namePrefix = Nothing}, "ascii")],
           elementNodes = [
             NodeContent (T.unwords $ map toTxt l)] })]
      
    ffaceoff [] = []
    ffaceoff l = [NodeElement (
        Element {
           elementName = Name {nameLocalName = "DataArray", nameNamespace = Nothing, namePrefix = Nothing},
           elementAttributes = [
             (Name {nameLocalName = "type", nameNamespace = Nothing, namePrefix = Nothing}, "Int64"),
             (Name {nameLocalName = "Name", nameNamespace = Nothing, namePrefix = Nothing}, "faceoffsets"),
             (Name {nameLocalName = "format", nameNamespace = Nothing, namePrefix = Nothing}, "ascii")],
           elementNodes = [
             NodeContent (T.unwords $ map toTxt l)]})]