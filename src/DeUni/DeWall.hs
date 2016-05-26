{-# LANGUAGE FlexibleContexts #-}
module DeUni.DeWall
  ( runHull3D
  , runDelaunay
  , runDelaunay3D
  , runDelaunay2D
  , reRun
  , SetSimplex2D
  , SetSimplex3D
  , SetFace3D
  , module DeUni.Types
  , module DeUni.Dim2.Base2D
  , module DeUni.Dim3.Base3D
  , module DeUni.GeometricTools
  ) where

import Prelude
import Control.Monad.State.Lazy
import Data.IntMap (IntMap)
import Linear.Vect
import qualified Data.List   as L
import qualified Data.Set    as S
import qualified Data.IntMap as IM

import DeUni.Types
import DeUni.GeometricTools
import DeUni.Dim3.Base3D
import DeUni.Dim2.Base2D
import DeUni.Dim2.Delaunay2D ()
import DeUni.Dim3.Delaunay3D ()
import DeUni.Dim3.Hull3D     ()

-- =================================================================================================

type SetSimplex2D = IntMap (S2 Vec2)
type SetSimplex3D = IntMap (S2 Vec3)
type SetFace3D    = IntMap (S1 Vec3)

-- ===================================| Exposed functions |=========================================

reRun :: (Buildable s p, Ord (Sub s p))
      => StateVarsMBC s p -> Box p -> SetPoint p -> [PointPointer] -> SetActiveSubUnits s p -> (IntMap(s p), StateVarsMBC s p)
reRun st box sP ps faces = runState (mbc ps faces box []) initi
  where
    initi = st { aflAlpha=S.empty, aflBox1=S.empty, aflBox2=S.empty, setPoint=sP }

runHull3D :: Box Vec3 -> SetPoint Vec3 -> [PointPointer] -> (SetFace3D, StateVarsMBC S1 Vec3)
runHull3D box sP ps = let
  initi = initState sP
  in runState (mbc ps (S.empty::SetActiveSubUnits S1 Vec3) box []) initi

runDelaunay3D :: Box Vec3 -> SetPoint Vec3 -> [PointPointer] -> (SetSimplex3D, StateVarsMBC S2 Vec3)
runDelaunay3D box sP ps = let
  initi = initState sP
  in runState (mbc ps (S.empty::SetActiveSubUnits S2 Vec3) box []) initi

runDelaunay2D :: Box Vec2 -> SetPoint Vec2 -> [PointPointer] -> (SetSimplex2D, StateVarsMBC S2 Vec2)
runDelaunay2D box sP ps = let
  initi = initState sP
  in runState (mbc ps (S.empty::SetActiveSubUnits S2 Vec2) box []) initi

runDelaunay :: (Buildable S2 a)=> Box a -> SetPoint a -> [PointPointer] -> (IntMap (S2 a), StateVarsMBC S2 a)
runDelaunay box sP ps = let
  initi = initState sP
  in runState (mbc ps S.empty box []) initi

initState :: SetPoint v -> StateVarsMBC simplex v
initState sP = StateVarsMBC
  { aflAlpha      = S.empty
  , aflBox1       = S.empty
  , aflBox2       = S.empty
  , externalFaces = S.empty
  , count         = 0
  , setPoint      = sP
  }

-- ====================================| Non-Exposed functions |====================================

-- | Marriage Before Conquer
mbc :: (Buildable slx dim, Ord (Sub slx dim))=> [PointPointer] -> SetActiveSubUnits slx dim -> Box dim -> [Position] -> StateMBC slx dim (IntMap (slx dim))
mbc p afl box subbox = do
  cleanAFLs
  sP <- fmap setPoint get
  let
    (plane, pairBox) = cutBox box subbox
    pp               = pointSetPartition (whichBoxIsIt pairBox) sP p
    p1               = pointsOnB1 pp
    p2               = pointsOnB2 pp ++ pointsOnPlane pp
  mapM_ (splitAF pairBox) (S.elems afl)
  units <- if S.null afl  -- same as null aflAlpha && null aflBox1 && null aflBox2
    then get1stUnits pairBox plane p1 p2 p
    else getUnitsOnPlane pairBox plane p
  analyzeUnit pairBox p1 p2 units
    where
      analyzeUnit pairBox p1 p2 units = do
        st <- get
        let
          afl1 = aflBox1 st
          afl2 = aflBox2 st
          box1 = halfBox1 pairBox
          box2 = halfBox2 pairBox
        case (IM.null units, S.null afl1, S.null afl2, L.null p1, L.null p2) of
          -- The state is independant and can discarted as it will be
          -- ereased at the bigein of the next recursive func call
          --(units, afl1,  afl2,  p1,    p2   )
          (True,  _,     _,     False, True ) -> mbc p1 afl1 box1 []
          (True,  _,     _,     True,  False) -> mbc p2 afl2 box2 []
          (True,  False, True,  False, False) -> mbc p  afl1 box  (subbox ++ [B1])
          (True,  True,  False, False, False) -> mbc p  afl2 box  (subbox ++ [B2])
          (True,  False, False, False, False) -> do
            us1 <- mbc p afl1 box (subbox ++ [B1])
            us2 <- mbc p afl2 box (subbox ++ [B2])
            return (us1 `IM.union` us2 `IM.union` units)
          (True,  _,     _,     _ ,    _    ) -> return IM.empty
          (False, True,  True,  _ ,    _    ) -> return units
          (False, True,  False, _ ,    _    ) -> IM.union units <$> mbc p2 afl2 box2 []
          (False, False, True,  _ ,    _    ) -> IM.union units <$> mbc p1 afl1 box1 []
          (False, False, False, _ ,    _    ) -> do
            us1 <- mbc p1 afl1 box1 []
            us2 <- mbc p2 afl2 box2 []
            return (us1 `IM.union` us2 `IM.union` units)

cleanAFLs :: (Buildable slx dim, Ord (Sub slx dim)) => StateMBC slx dim ()
cleanAFLs = modify (\x -> x { aflAlpha=S.empty, aflBox1=S.empty, aflBox2=S.empty })


get1stUnits :: (Buildable slx dim, Ord (Sub slx dim)) => BoxPair dim -> Plane dim -> [PointPointer] ->[PointPointer] -> [PointPointer] -> StateMBC slx dim (IntMap (slx dim))
get1stUnits pairBox plane p1 p2 p = do
  sP <- fmap setPoint get
  case build1stUnit plane sP p1 p2 p of
    Just unit -> do
      mapM_ (splitAF pairBox) (getAllSubUnits sP unit)
      units <- getUnitsOnPlane pairBox plane p
      cnt   <- fmap count get
      modify (\x -> x { count = cnt + 1 })
      return $ IM.insert cnt unit units
    _         -> return IM.empty

-- Simplex Wall Construction
getUnitsOnPlane :: (Buildable slx dim, Ord (Sub slx dim)) => BoxPair dim -> Plane dim -> [PointPointer] -> StateMBC slx dim (IntMap (slx dim))
getUnitsOnPlane pairBox plane p = do
  st <- get
  sP <- fmap setPoint get
  if S.null (aflAlpha st)
    then return IM.empty
    else let
      actSubUnit       = S.findMax . aflAlpha $ st
      newUnit          = buildUnit actSubUnit sP p
      removeSubUnit su = modify (\x -> x { aflAlpha = S.delete su (aflAlpha x) })
      in case newUnit of
        Just sig -> do
          mapM_ (splitAF pairBox) (getAllSubUnits sP sig)
          sigs <- getUnitsOnPlane pairBox plane p
          cnt  <- fmap count get
          modify (\x -> x { count = cnt + 1 })
          return $ IM.insert cnt sig sigs
        _ -> do
          modify (\x -> x { externalFaces = S.insert actSubUnit (externalFaces x) })
          removeSubUnit actSubUnit
          getUnitsOnPlane pairBox plane p

splitAF :: (Buildable slx dim, Ord (Sub slx dim)) => BoxPair dim -> ActiveSubUnit slx dim -> StateMBC slx dim ()
splitAF pairBox e = do
  sP <- fmap setPoint get
  case subUnitPos pairBox sP e of
    B1         -> upP1
    B2         -> upP2
    CrossPlane -> upAlpha
    _          -> return ()
    where
      upP1     = updateSubUnit e modP1
      upP2     = updateSubUnit e modP2
      upAlpha  = updateSubUnit e modAlpha
      modP1    = (aflBox1 , \modf x -> x { aflBox1  = modf })
      modP2    = (aflBox2 , \modf x -> x { aflBox2  = modf })
      modAlpha = (aflAlpha, \modf x -> x { aflAlpha = modf })

      updateSubUnit edge (func, modf) = do
        set <- func <$> get
        if S.member edge set
          then modify (modf $ S.delete edge set)
          else modify (modf $ S.insert edge set)
        return ()
