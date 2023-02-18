{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}
module MarchingCubes.Internal
  ( module MarchingCubes.Internal
  ) where
import           Data.Array.Unboxed             ( IArray(..)
                                                , UArray
                                                )
import qualified Data.Array.Unboxed            as A
import           Data.Bits                      ( shiftL )
import qualified Data.Foldable                 as F
import           Data.List                      ( zipWith4 )
import           Data.Matrix                    ( Matrix(..)
                                                , elementwiseUnsafe
                                                , fromLists
                                                , getElem
                                                , getRow
                                                , matrix
                                                , minorMatrix
                                                , scaleMatrix
                                                )
import qualified Data.Matrix                   as M
import           Data.Sequence                  ( Seq
                                                , (|>)
                                                )
import qualified Data.Sequence                 as S
import qualified Data.Vector                   as V
import           Data.Vector.Unboxed            ( (!)
                                                , Unbox
                                                , Vector
                                                )
import qualified Data.Vector.Unboxed           as UV
import           MarchingCubes.Tables           ( crf
                                                , facePoints
                                                , indexArray
                                                )
import           MarchingCubes.Utils            ( arrayToMatrix
                                                , kro1
                                                , kro2
                                                , levelMatrix
                                                , whichIndicesAndItems
                                                )

facesNo7
  :: (Num a, Ord a, Unbox a)
  => Vector Int
  -> Vector Int
  -> Vector a
  -> Int
  -> Int
  -> [Int]
facesNo7 faces p1 values l j = map fun [0 .. l - 1]
 where
  fun i = if temp == 1 then shiftL 1 (j - 1) else 0
   where
    f  = abs (faces ! i) - 1
    e  = facePoints V.! f
    e1 = e ! 1
    e2 = e ! 2
    e3 = e ! 3
    e4 = e ! 4
    p  = p1 ! i - 2
    a  = values ! (p + e1)
    b  = values ! (p + e2)
    c  = values ! (p + e3)
    d  = values ! (p + e4)
    temp =
      (if faces ! i > 0 then 1 :: Int else -1)
        * (if a * b - c * d > 0 then 1 else -1)

faces7
  :: (RealFloat a, Unbox a)
  => Vector Int
  -> Vector Int
  -> Vector a
  -> Int
  -> Int
  -> [Int]
faces7 faces p1 values l j = map fun [0 .. l - 1]
 where
  fun i = if temp == 1 then shiftL 1 (j - 1) else 0
   where
    p         = (p1 ! i) - 1
    a0        = values ! p
    b0        = values ! (p + 3)
    c0        = values ! (p + 2)
    d0        = values ! (p + 1)
    a1        = values ! (p + 4)
    b1        = values ! (p + 7)
    c1        = values ! (p + 6)
    d1        = values ! (p + 5)
    a         = (a1 - a0) * (c1 - c0) - (b1 - b0) * (d1 - d0)
    b = c0 * (a1 - a0) + a0 * (c1 - c0) - d0 * (b1 - b0) - b0 * (d1 - d0)
    c         = a0 * c0 - b0 * d0
    tmax      = -b / 2 / a
    mxmm      = a * tmax * tmax + b * tmax + c
    mxmm'     = if isNaN mxmm then -1 else mxmm
    cond1     = a < 0
    cond2     = tmax > 0
    cond3     = tmax < 1
    cond4     = mxmm' > 0
    totalcond = cond1 && cond2 && cond3 && cond4
    temp =
      (if faces ! i > 0 then 1 :: Int else -1) * (if totalcond then 1 else -1)

faceType :: Real a => Matrix a -> a -> a -> Matrix Int
faceType mtrx level mx = elementwiseUnsafe (+) sum1 sum2
 where
  lm         = levelMatrix mtrx level (level < mx)
  m          = nrows mtrx
  n          = ncols mtrx
  minorMat   = minorMatrix m n lm
  sminorMat2 = scaleMatrix 2 (minorMatrix 1 n lm)
  sminorMat4 = scaleMatrix 4 (minorMatrix 1 1 lm)
  sminorMat8 = scaleMatrix 8 (minorMatrix m 1 lm)
  sum1       = elementwiseUnsafe (+) minorMat sminorMat2
  sum2       = elementwiseUnsafe (+) sminorMat4 sminorMat8

levCells
  :: (Real a, IArray UArray a)
  => UArray (Int, Int, Int) a
  -> a
  -> a
  -> Matrix Int
levCells a level mx = out
 where
  bottomTypes             = faceType (arrayToMatrix a 0) level mx
  (_, (nx', ny', nz'))    = bounds a
  -- nx = nx' + 1
  -- ny = ny' + 1
  -- nz = nz' + 1
  (lengths, cells, types) = go 0 S.empty bottomTypes S.empty S.empty
  go
    :: Int
    -> Seq Int
    -> Matrix Int
    -> Seq (Seq Int)
    -> Seq (Seq Int)
    -> (Seq Int, Seq (Seq Int), Seq (Seq Int))
  go k !lngths !bTypes !clls !tps
    | k == nz'  = (lngths, clls, tps)
    | otherwise = go (k + 1) (lngths |> l) tTypes (clls |> cll) (tps |> tp)
   where
    tTypes    = faceType (arrayToMatrix a (k + 1)) level mx
    cellTypes = elementwiseUnsafe (+) bTypes (scaleMatrix 16 tTypes)
    goodcells = whichIndicesAndItems cellTypes
    l         = S.length goodcells
    cll       = fmap (\(i, _) -> i + nx' * ny' * k + 1) goodcells
    tp        = fmap snd goodcells
  out = M.transpose (fromLists (concatMap f [0 .. nz' - 1]))
  f k = map (g k) [0 .. S.index lengths k - 1]
  g k l =
    [ c `mod` nx' + 1
    , (c `div` nx') `mod` ny' + 1
    , c `div` (nx' * ny') + 1
    , S.index (S.index types k) l
    ]
    where c = S.index (S.index cells k) l - 1

getBasic1 :: Vector Int -> Matrix Int -> Matrix Int
getBasic1 r vivjvk = elementwiseUnsafe (+) k1 k2
 where
  nR    = UV.length r
  cube1 = matrix nR 3 (\(i, j) -> getElem (r ! (i - 1) + 1) j vivjvk)
  k1    = kro1 indexArray nR
  k2    = kro2 cube1 8

getBasic2
  :: (Num a, Unbox a, IArray UArray a)
  => UArray (Int, Int, Int) a
  -> a
  -> Matrix Int
  -> Vector a
getBasic2 a level cubeco = UV.fromList values
 where
  f i j = getElem i j cubeco - 1
  values =
    [ a A.! (f i 1, f i 2, f i 3) - level | i <- [1 .. nrows cubeco - 1] ]
    ++ [0]

getTcase :: V.Vector Int -> Vector Int
getTcase types =
  UV.fromList [ crf ! (types V.! i) - 1 | i <- [0 .. V.length types - 1] ]

getR :: Vector Int -> Vector Int
getR tcase = UV.fromList $ F.toList $ go 0 S.empty
 where
  n = UV.length tcase
  go :: Int -> Seq Int -> Seq Int
  go i !out
    | i == n
    = out
    | otherwise
    = if tc
         == 1
         || tc
         == 2
         || tc
         == 5
         || tc
         == 8
         || tc
         == 9
         || tc
         == 11
         || tc
         == 14
      then
        go (i + 1) (out |> i)
      else
        go (i + 1) out
    where tc = tcase ! i

lambdaMu :: RealFrac a => [a] -> ([a], [a])
lambdaMu x1 = (lambda, mu)
 where
  lambda = map (\x -> fromInteger $ floor (x / 9)) x1
  mu     = map (1 -) lambda

average :: Num a => ([a], [a]) -> [a] -> [a] -> [a]
average (lambda, mu) = zipWith4 (\a b c d -> b * c + a * d) lambda mu

average7 :: Num a => ([a], [a]) -> [a] -> [a]
average7 (lambda, mu) = zipWith3 (\a b c -> b * c + a) lambda mu

average8 :: Num a => ([a], [a]) -> [a] -> [a]
average8 (lambda, mu) = zipWith3 (\a b c -> b * c - a) lambda mu

getPoints
  :: (Unbox a, RealFrac a)
  => Matrix Int
  -> Vector a
  -> [Int]
  -> [Int]
  -> [Int]
  -> Matrix a
getPoints cubeco values p1 x1 x2 = fromLists
  [out0, out1, out2, out3, out4, out5, out6, out7]
 where
  p1x1     = zipWith (+) p1 x1
  p1x2     = zipWith (+) p1 x2
  xx1      = map fromIntegral x1
  lambdamu = lambdaMu xx1
  v1       = map (\j -> fromIntegral $ getElem (j - 1) 1 cubeco) p1x1
  w1       = map (\j -> fromIntegral $ getElem j 1 cubeco) p1
  v2       = map (\j -> fromIntegral $ getElem (j - 1) 1 cubeco) p1x2
  w2       = map (\j -> fromIntegral $ getElem (j + 1) 1 cubeco) p1
  v3       = map (\j -> fromIntegral $ getElem (j - 1) 2 cubeco) p1x1
  w3       = map (\j -> fromIntegral $ getElem (j + 1) 2 cubeco) p1
  v4       = map (\j -> fromIntegral $ getElem (j - 1) 2 cubeco) p1x2
  w4       = map (\j -> fromIntegral $ getElem (j + 2) 2 cubeco) p1
  v5       = map (\j -> fromIntegral $ getElem (j - 1) 3 cubeco) p1x1
  w5       = map (\j -> fromIntegral $ getElem (j + 1) 3 cubeco) p1
  v6       = map (\j -> fromIntegral $ getElem (j - 1) 3 cubeco) p1x2
  w6       = map (\j -> fromIntegral $ getElem (j + 5) 3 cubeco) p1
  v7       = map (\j -> values ! (j - 2)) p1x1
  v8       = map (\j -> values ! (j - 2)) p1x2
  out0     = average lambdamu v1 w1
  out1     = average lambdamu v2 w2
  out2     = average lambdamu v3 w3
  out3     = average lambdamu v4 w4
  out4     = average lambdamu v5 w5
  out5     = average lambdamu v6 w6
  out6     = average7 lambdamu v7
  out7     = average8 lambdamu v8

calPoints :: Fractional a => Matrix a -> Matrix a
calPoints points = M.transpose $ fromLists [x, y, z]
 where
  x1 = getRow 1 points
  x2 = getRow 2 points
  y1 = getRow 3 points
  y2 = getRow 4 points
  z1 = getRow 5 points
  z2 = getRow 6 points
  v1 = getRow 7 points
  v2 = getRow 8 points
  s  = V.zipWith (\a b -> a / (a - b)) v1 v2
  x  = V.toList $ V.zipWith3 (\a b c -> a + c * (b - a)) x1 x2 s
  y  = V.toList $ V.zipWith3 (\a b c -> a + c * (b - a)) y1 y2 s
  z  = V.toList $ V.zipWith3 (\a b c -> a + c * (b - a)) z1 z2 s