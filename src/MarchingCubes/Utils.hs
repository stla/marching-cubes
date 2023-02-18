{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}
module MarchingCubes.Utils
  ( module MarchingCubes.Utils
  ) where
import           Data.Array.Unboxed             ( (!)
                                                , IArray(..)
                                                , UArray
                                                )
import           Data.List                      ( transpose )
import           Data.Matrix                    ( (<|>)
                                                , Matrix(..)
                                                , colVector
                                                , getElem
                                                , mapPos
                                                , matrix
                                                , toLists
                                                )
import qualified Data.Matrix                   as M
import           Data.Sequence                  ( (><)
                                                , Seq
                                                , (|>)
                                                )
import qualified Data.Sequence                 as S
import qualified Data.Vector                   as V
import           Data.Vector.Unboxed            ( Vector )
import qualified Data.Vector.Unboxed           as UV

levelMatrix :: Real a => Matrix a -> a -> Bool -> Matrix Int
levelMatrix mtrx level strict = mapPos (\_ x -> if lt x then 1 else 0) mtrx
  where lt = if strict then (<) level else (<=) level

arrayToMatrix :: IArray UArray a => UArray (Int, Int, Int) a -> Int -> Matrix a
arrayToMatrix arr k = matrix (nx + 1)
                             (ny + 1)
                             (\(i, j) -> arr ! (i - 1, j - 1, k))
  where (_, (nx, ny, _)) = bounds arr

whichIndicesAndItems :: Matrix Int -> Seq (Int, Int)
whichIndicesAndItems mtrx = go 1 S.empty
 where
  m = nrows mtrx + 1
  n = ncols mtrx + 1
  go :: Int -> Seq (Int, Int) -> Seq (Int, Int)
  go j !out | j == n    = out
            | otherwise = go (j + 1) (inner 1 j out)
  inner i !j !out | i == m    = out
                  | otherwise = inner (i + 1) j out'
   where
    out' =
      let x = getElem i j mtrx
      in  if x > 0 && x < 255
            then out |> ((j - 1) * (m - 1) + i - 1, x)
            else out

kro1 :: Matrix Int -> Int -> Matrix Int
kro1 mtrx n = matrix p ny f
 where
  nx = nrows mtrx
  ny = ncols mtrx
  p  = nx * n + 1
  f (i, j) = if i < p then getElem ((i - 1) `mod` nx + 1) j mtrx else 0

kro2 :: Matrix Int -> Int -> Matrix Int
kro2 mtrx n = matrix p ny f
 where
  nx         = nrows mtrx
  ny         = ncols mtrx
  p          = nx * n + 1
  replicates = foldr ((><) . S.replicate n) S.empty [1 .. nx] :: Seq Int
  f (i, j) = if i < p then getElem (S.index replicates (i - 1)) j mtrx else 0

replicateEach :: [a] -> [Int] -> Seq a
replicateEach list counts = foldr (><) S.empty sequences
  where sequences = zipWith S.replicate counts list

replicateEach' :: [a] -> Int -> Seq a
replicateEach' list n = foldr ((><) . S.replicate n) S.empty list

vector2matrix :: [a] -> Int -> Matrix a
vector2matrix list ny = M.fromList (length list `div` ny) ny list

cbind :: Matrix a -> [a] -> [a] -> Matrix a
cbind mtrx col1 col2 = mtrx <|> col1' <|> col2'
 where
  col1' = colVector (V.fromList col1)
  col2' = colVector (V.fromList col2)

subMatrix :: Matrix a -> [Int] -> Vector Int -> Matrix a
subMatrix mtrx rows cols = matrix
  (length rows)
  (UV.length cols)
  (\(i, j) -> getElem (rows !! (i - 1) + 1) (cols UV.! (j - 1) + 1) mtrx)

matrix2listMinusFirstColumn :: Matrix a -> [a]
matrix2listMinusFirstColumn mtrx =
  concat $ transpose $ tail (toLists $ M.transpose mtrx)

jthColumn :: Vector Int -> Int -> Int -> Vector Int
jthColumn vec ncol j = UV.map (\i -> vec UV.! (i * ncol + j))
                              (UV.enumFromN 0 nrow)
  where nrow = UV.length vec `div` ncol