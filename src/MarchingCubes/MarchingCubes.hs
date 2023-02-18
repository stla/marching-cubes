{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}
module MarchingCubes.MarchingCubes
  ( Voxel
  , XYZ
  , marchingCubes
  , makeVoxel
  ) where
import           Data.Array.Unboxed             ( IArray
                                                , UArray
                                                , listArray
                                                )
import qualified Data.Foldable                 as F
import           Data.List                      ( elemIndices
                                                , findIndices
                                                )
import           Data.Matrix                    ( (<->)
                                                , Matrix(ncols)
                                                , getCol
                                                , getRow
                                                , mapCol
                                                , submatrix
                                                )
import qualified Data.Matrix                   as M
import           Data.Maybe                     ( catMaybes
                                                , mapMaybe
                                                )
import qualified Data.Vector                   as V
import           Data.Vector.Unboxed            ( (!)
                                                , Unbox
                                                )
import qualified Data.Vector.Unboxed           as UV
import           MarchingCubes.Internal         ( calPoints
                                                , faces7
                                                , facesNo7
                                                , getBasic1
                                                , getBasic2
                                                , getPoints
                                                , getR
                                                , getTcase
                                                , levCells
                                                )
import           MarchingCubes.Tables           ( edgePoints
                                                , edgesLengths
                                                , edgesTable
                                                , edgesTable2
                                                , facesTable
                                                , specialInd
                                                , specialName
                                                , specialNedge
                                                , specialNface
                                                , specialPos
                                                )
import           MarchingCubes.Utils            ( cbind
                                                , jthColumn
                                                , matrix2listMinusFirstColumn
                                                , replicateEach
                                                , replicateEach'
                                                , subMatrix
                                                , vector2matrix
                                                )

type Bounds a = ((a, a), (a, a), (a, a))
type Dims = (Int, Int, Int)
type Voxel a = ((UArray Dims a, a), (Bounds a, Dims))
type XYZ a = (a, a, a)

-- | Make the voxel. 
makeVoxel
  :: (RealFloat a, IArray UArray a)
  => (XYZ a -> a) -- ^ the function defining the isosurface
  -> Bounds a     -- ^ bounds of the grid
  -> Dims         -- ^ numbers of subdivisions of the grid
  -> Voxel a
makeVoxel fun bds@((xm, xM), (ym, yM), (zm, zM)) dims@(nx, ny, nz) =
  ((listArray ((0, 0, 0), (nx - 1, ny - 1, nz - 1)) values, mxmm), (bds, dims))
 where
  x_ = [ xm + (xM - xm) * fracx i | i <- [0 .. nx - 1] ]
  fracx p = realToFrac p / (realToFrac nx - 1)
  y_ = [ ym + (yM - ym) * fracy i | i <- [0 .. ny - 1] ]
  fracy p = realToFrac p / (realToFrac ny - 1)
  z_ = [ zm + (zM - zm) * fracz i | i <- [0 .. nz - 1] ]
  fracz p = realToFrac p / (realToFrac nz - 1)
  values = [ fun (x, y, z) | x <- x_, y <- y_, z <- z_ ]
  mxmm   = maximum (filter (not . isNaN) values)

rescale :: Fractional a => (a, a) -> Int -> a -> a
rescale (minmm, maxmm) n w = minmm + (maxmm - minmm) * w / fromIntegral (n - 1)

rescaleMatrix :: Fractional a => Matrix a -> Bounds a -> Dims -> Matrix a
rescaleMatrix mtrx (xbds, ybds, zbds) (nx, ny, nz) = mtrx'''
 where
  mtrx'   = mapCol (\_ w -> rescale xbds nx (w - 1)) 1 mtrx
  mtrx''  = mapCol (\_ w -> rescale ybds ny (w - 1)) 2 mtrx'
  mtrx''' = mapCol (\_ w -> rescale zbds nz (w - 1)) 3 mtrx''

marchingCubes
  :: (RealFloat a, Unbox a, IArray UArray a) => Voxel a -> a -> Matrix a
marchingCubes ((voxel, mx), (bds, dims)) level = rescaleMatrix
  (maybe triangles1 (triangles1 <->) triangles2)
  bds
  dims
 where
  ijkt         = levCells voxel level mx
  vt           = getRow 4 ijkt
  tcase        = getTcase vt
  r            = getR tcase
  nR           = UV.length r
  ijk          = submatrix 1 3 1 (ncols ijkt) ijkt
  vivjvk       = M.transpose ijk
  cubeco       = getBasic1 r vivjvk
  values       = getBasic2 voxel level cubeco
  p1           = [ 8 * i + 1 | i <- [0 .. nR - 1] ]
  cases        = UV.map (\j -> vt V.! j - 1) r
  edgeslengths = UV.map (edgesLengths UV.!) cases
  p1rep        = F.toList $ replicateEach p1 (UV.toList edgeslengths)
  edges =
    [ (edgesTable V.! (cases ! i)) ! j
    | i <- [0 .. nR - 1]
    , j <- [0 .. edgeslengths ! i - 1]
    ]
  epoints        = map (\j -> edgePoints V.! (j - 1)) edges
  x1             = map (! 1) epoints
  x2             = map (! 2) epoints
  points         = getPoints cubeco values p1rep x1 x2
  triangles1     = calPoints points
  -- ~~## special cases ##~~
  specialCases   = map (`UV.elemIndices` tcase) specialName
  r3s            = filter (not . UV.null) specialCases
  cs             = findIndices (not . UV.null) specialCases
  setOfTriangles = catMaybes (zipWith special cs r3s)
  special c r3 = if null newtriangles
    then Nothing
    else Just $ foldl1 (<->) newtriangles
   where
    nR3     = UV.length r3
    cubeco3 = getBasic1 r3 vivjvk
    values3 = getBasic2 voxel level cubeco3
    p13     = UV.map (\i -> 8 * i + 1) (UV.enumFromN 0 nR3)
    cases3  = [ vt V.! (r3 ! i) - 1 | i <- [0 .. nR3 - 1] ]
    nedge   = specialNedge ! c
    faces3  = UV.concat $ map ((V.!) facesTable) cases3
    index3  = case c of
      0 -> facesNo7 faces3 p13 values3 nR3 1
      1 -> faces7 faces3 p13 values3 nR3 1
      _ -> index3'
     where
      nface    = specialNface ! c
      facelast = jthColumn faces3 nface (nface - 1)
      index3'  = loop 0 (faces7 facelast p13 values3 nR3 nface)
      loop j !idx
        | j == nface - 1
        = idx
        | otherwise
        = let facej = jthColumn faces3 nface j
          in  let temp = facesNo7 facej p13 values3 nR3 (j + 1)
              in  loop (j + 1) (zipWith (+) idx temp)
    edges3'      = UV.toList $ UV.concat $ map ((V.!) edgesTable2) cases3
    edges3       = vector2matrix edges3' nedge
    edgesp1index = cbind edges3 (UV.toList p13) index3
    ind3         = specialInd V.! c
    newtriangles = mapMaybe f [0 .. UV.length ind3 - 1]
    f j = triangles3
     where
      wrows      = elemIndices (ind3 ! j) index3
      triangles3 = if null wrows then Nothing else Just $ calPoints points3
       where
        wcols     = UV.cons nedge ((specialPos V.! c) V.! j)
        ed        = subMatrix edgesp1index wrows wcols
        col0ed    = V.toList $ getCol 1 ed
        col0edrep = F.toList $ replicateEach' col0ed (UV.length wcols - 1)
        edge1     = matrix2listMinusFirstColumn ed
        epoints'  = map (\k -> edgePoints V.! (k - 1)) edge1
        x1'       = map (! 1) epoints'
        x2'       = map (! 2) epoints'
        points3   = getPoints cubeco3 values3 col0edrep x1' x2'
  triangles2 = if null setOfTriangles
    then Nothing
    else Just $ foldl1 (<->) setOfTriangles
