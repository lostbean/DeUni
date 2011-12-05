module Math.WeightedPoint where

import Data.Vec hiding (reverse, foldl, map)
import qualified Data.Vec as V 

type Matrix = Mat33 Double


data WPoint = WPoint
            { weigth :: Double
            , point  :: Vec3 Double
            } deriving (Show)


a = WPoint 0 (1:.3:.0)
b = WPoint 0 (-3:.0:.0)
c = WPoint 0 (11:.0:.20)
d = WPoint 0 (0:.1:.0)
m1 :: Matrix
m1 = matFromLists [[12.0,6.0,-4.0],[-51.0,167.0,24.0],[4.0,-68.0,-41.0]]

getM :: WPoint -> WPoint -> WPoint -> WPoint -> Matrix
getM a b c d = (point b - point a) :. (point c - point a) :. (point d - point a) :. ()

getAlpha :: WPoint -> WPoint -> WPoint -> WPoint -> Vec3 Double
getAlpha a b c d = (fun b a) :. (fun c a) :. (fun d a) :. ()
  where fun x y = (normSq.point) x - (normSq.point) y - weigth x + weigth y
        
solveMu a r = (mux, muy, muz)
  where
    (ax, ay, az)    = unvec3 a
    (r1, r2, r3)    = unvec3 r
    (r11, r21, r31) = unvec3 r1
    (r12, r22, r32) = unvec3 r2
    (r13, r23, r33) = unvec3 r3
    mux = ax / r11
    muy = (ay - mux*r12) / r22
    muz = (az - mux*r13 - muy* r23) / r33

unvec3 (a:.b:.c:.()) = (a,b,c)        


center a b c d = (-0.5) * ((vec mux)*q1 + (vec muy)*q2 + (vec muz)*q3)
  where
    m               = getM a b c d
    (q,r)           = qrGramSchmidt m
    (mux, muy, muz) = solveMu (vec (-1) * (getAlpha a b c d)) r
    (q1, q2, q3)    = unvec3 q
    
       
delaunyDist a b c d
  | signDet r = dist   
  | otherwise = (-1) * dist
  where  
    dist            = muz*0.5 + (q3 `dot` (point a))
    m               = getM a b c d
    (q,r)           = qrGramSchmidt m
    (mux, muy, muz) = solveMu (vec (-1) * (getAlpha a b c d)) r
    (q1, q2, q3)    = unvec3 q
    
signDet r = r11*r22*r33 >= 0
  where
    (r1, r2, r3)    = unvec3 r
    (r11, r21, r31) = unvec3 r1
    (r12, r22, r32) = unvec3 r2
    (r13, r23, r33) = unvec3 r3
    


dsList a b c ds = map func ds
  where
    func         = delaunyDistPart qr upAlpha upR a
    -- Calc QR using a dump 'd' since it will be updated afterwards
    qr           = qrGramSchmidt $ getM a b c (WPoint 0 (0:.0:.0:.()))
    upAlpha x    = vec (-1) * (getAlpha a b c x)
    upR (q,r) d' = (r1:.r2:.r3:.())
      where
        (q1, q2, q3) = unvec3 q
        (r1, r2, _ ) = unvec3 r
        r3           = (fun q1:.fun q2:.fun q3:.())
        fun x        = dot x (point d' - point a)

delaunyDistPart qr@(q,_) upA upR a d
  | signDet r' = dist   
  | otherwise  = (-1) * dist
  where
    dist            = muz*0.5 + (q3 `dot` (point a))
    (mux, muy, muz) = solveMu (upA d) r'
    (q1, q2, q3)    = unvec3 q
    r'              = upR qr d



-- | Gram-Schmidt is a very simple process to find an orthonormal basis for a matrix.
-- It takes the columns of the matrix and subtracts the part of each vector that is
-- parallel to the previous columns. After performing this process, we are left with
-- the Q matrix of our decomposition.
toColumns   = toList
fromColumns = fromList

orthoBasis :: Matrix -> Matrix
orthoBasis = fromColumns . reverse . (foldl orthonormal []) . toColumns
  where 
    subtractParallel :: Vec3 Double -> Vec3 Double -> Vec3 Double
    subtractParallel v u = v - u * (vec (u `dot` v))
    --orthonormal :: Matrix -> Vec3 Double -> Matrix
    orthonormal vs v = (normalize $ foldl subtractParallel v vs) : vs

qrGramSchmidt :: Matrix -> (Matrix, Matrix)
qrGramSchmidt m = (q,r)
  where q = orthoBasis m
        r = m `multmm` (transpose q)




