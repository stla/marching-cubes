{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TupleSections    #-}
module MarchingCubes.Mesh
  ( makeMesh
  , Mesh (..)
  ) where
import           Data.Array.Unboxed             ( IArray
                                                , UArray
                                                )
import           Data.Bifunctor                 ( bimap )
import           Data.Foldable                  ( foldl' )
import           Data.IntMap.Strict             ( IntMap
                                                , elems
                                                , unionWith
                                                )
import qualified Data.IntMap.Strict            as IM
import           Data.List.Split                ( chunksOf )
import           Data.Map.Strict                ( insertLookupWithKey )
import qualified Data.Map.Strict               as M
import           Data.Matrix                    ( Matrix
                                                , nrows
                                                , toLists
                                                )
import           Data.Tuple.Extra               ( (***)
                                                , curry3
                                                , second
                                                , thd3
                                                )
import           Data.Vector.Unboxed            ( (!)
                                                , Unbox
                                                , Vector
                                                , fromList
                                                )
import           Linear                         ( V3(..)
                                                , (^+^)
                                                , (^-^)
                                                , cross
                                                , signorm
                                                )
import           MarchingCubes.MarchingCubes    ( Voxel
                                                , marchingCubes
                                                )

-- type Mesh a = ((Vector (a, a, a), [[Int]]), Vector (a, a, a))

data Mesh a = Mesh 
  { 
    _vertices :: Vector (a, a, a)
  , _faces :: [(Int, Int, Int)] 
  , _normals :: Vector (a, a, a)
  } 
  deriving Show

matrix2listOfTriples :: Matrix a -> [(a, a, a)]
matrix2listOfTriples mtrx =
  map (\row -> (row !! 0, row !! 1, row !! 2)) (toLists mtrx)

uniqueWithIndices :: (Unbox a, Ord a) => [a] -> (Vector a, Vector Int)
uniqueWithIndices list = (fromList *** fromList) $ go 0 M.empty list
 where
  go _ _ []       = ([], [])
  go i m (x : xs) = case insertLookupWithKey (curry3 thd3) x i m of
    (Nothing, m') -> bimap (x :) (i :) (go (i + 1) m' xs)
    (Just j , m') -> second (j :) (go i m' xs)

undupMesh
  :: (Unbox a, Ord a) => ([(a, a, a)], [[Int]]) -> (Vector (a, a, a), [[Int]])
undupMesh (vertices, faces) = (vertices', faces')
 where
  (vertices', idx) = uniqueWithIndices vertices
  faces'           = map (map (idx !)) faces

normal :: Floating a => (a, a, a) -> (a, a, a) -> (a, a, a) -> V3 a
normal v1 v2 v3 = signorm $ cross (v3' ^-^ v1') (v2' ^-^ v1') -- (NaN,NaN,NaN) if degenerate face
 where
  toV3 (x, y, z) = V3 x y z
  v1' = toV3 v1
  v2' = toV3 v2
  v3' = toV3 v3

faceNormals
  :: (Floating a, Unbox a) => Vector (a, a, a) -> [Int] -> IntMap (V3 a)
faceNormals vs face = IM.fromList nrmls
 where
  nrml  = normal (vs ! (face !! 0)) (vs ! (face !! 1)) (vs ! (face !! 2))
  nrmls = map (, nrml) face

normalsV3
  :: (Floating a, Unbox a) => (Vector (a, a, a), [[Int]]) -> IntMap (V3 a)
normalsV3 (vs, faces) = IM.map
  signorm
  (foldl' (unionWith (^+^)) IM.empty (map (faceNormals vs) faces))

normals
  :: (Floating a, Unbox a) => (Vector (a, a, a), [[Int]]) -> Vector (a, a, a)
normals mesh = fromList $ map fromV3 (elems $ normalsV3 mesh)
  where fromV3 (V3 x y z) = (x, y, z)

degenerateFace :: (Unbox a, Eq a) => Vector (a, a, a) -> [Int] -> Bool
degenerateFace vertices face = v1 == v2 || v1 == v3 || v2 == v3
 where
  v1 = vertices ! (face !! 0)
  v2 = vertices ! (face !! 1)
  v3 = vertices ! (face !! 2)

-- | Make mesh of isosurface.
makeMesh :: 
  (Unbox a, RealFloat a, IArray UArray a) 
  => Voxel a -- ^ voxel obtained with `makeVoxel`
  -> a       -- ^ isovalue
  -> Mesh a
makeMesh voxel level = Mesh 
  {
    _vertices = vertices', _faces = faces''', _normals = normals mesh
  }
 where
  mtrx                = marchingCubes voxel level
  vertices            = matrix2listOfTriples mtrx
  faces               = chunksOf 3 [0 .. nrows mtrx - 1]
  (vertices', faces') = undupMesh (vertices, faces)
  faces''             = filter (not . degenerateFace vertices') faces'
  faces'''            = map toTriplet faces''
  mesh                = (vertices', faces'')
  toTriplet ikj       = (ikj !! 0, ikj !! 2, ikj !! 1)