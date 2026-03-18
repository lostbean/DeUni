# DeUni

**Delaunay Triangulation and Convex Hull via Marriage Before Conquer**

DeUni is a Haskell library for computing 2D and 3D Delaunay triangulations (including regular/weighted triangulations) and 3D convex hulls. It implements the Marriage Before Conquer (MBC) variant of the divide-and-conquer algorithm, where the "marriage" step -- constructing simplices that straddle the dividing plane -- is performed before recursing into each half-space.

## Features

- 2D Delaunay triangulation (triangles from a planar point set)
- 3D Delaunay triangulation (tetrahedralization from a spatial point set)
- 3D convex hull (triangular surface mesh)
- Support for weighted (power) distances, enabling regular triangulations
- Dimension-generic core algorithm with concrete 2D and 3D instantiations
- Pure Haskell; computation runs in the `State` monad with no IO required

## Algorithm: Marriage Before Conquer (MBC)

The MBC algorithm is a divide-and-conquer strategy for constructing Delaunay triangulations and convex hulls. It differs from the classical approach by performing the merge ("marriage") step _before_ the recursive conquering of sub-problems.

### How it works

1. **Divide:** A cutting plane bisects the bounding box along its longest axis, partitioning all points into two half-spaces (B1 and B2).

1. **Marry:** Simplices that cross the cutting plane (the "alpha" set) are constructed first. Starting from an initial seed simplex found via `FirstSeed`, the algorithm builds all simplices whose sub-units (edges in 2D, faces in 3D) straddle the plane. Each new simplex's sub-units are classified into three Active Face Lists (AFLs):

   - `aflAlpha` -- sub-units crossing the cutting plane (processed next)
   - `aflBox1` -- sub-units entirely in half-space B1
   - `aflBox2` -- sub-units entirely in half-space B2

   A sub-unit appearing twice in an AFL is removed (it is shared by two simplices and therefore closed).

1. **Conquer:** Once all alpha simplices are built, the algorithm recurses into B1 with `aflBox1` and into B2 with `aflBox2`. Each recursive call re-bisects its sub-box and repeats.

1. **Termination:** Recursion ends when the AFL is empty and no further simplices can be constructed.

### Circumsphere / Circumcircle computation

DeUni computes circumspheres (3D) and circumcircles (2D) using QR decomposition (Gram-Schmidt orthogonalization) of the matrix formed by the simplex edge vectors. The `power distance` is used instead of Euclidean distance, supporting weighted point sets and regular triangulations.

## Module Reference

### Core modules

#### `DeUni.DeWall`

The main entry point. Re-exports all other modules and provides the top-level API functions:

| Function | Type | Description |
|---|---|---|
| `runDelaunay2D` | `Box Vec2 -> SetPoint Vec2 -> [PointPointer] -> (SetSimplex2D, StateVarsMBC S2 Vec2)` | Run 2D Delaunay triangulation |
| `runDelaunay3D` | `Box Vec3 -> SetPoint Vec3 -> [PointPointer] -> (SetSimplex3D, StateVarsMBC S2 Vec3)` | Run 3D Delaunay triangulation |
| `runDelaunay` | `(Buildable S2 a) => Box a -> SetPoint a -> [PointPointer] -> (IntMap (S2 a), StateVarsMBC S2 a)` | Generic Delaunay (works for any dimension with a `Buildable S2` instance) |
| `runHull3D` | `Box Vec3 -> SetPoint Vec3 -> [PointPointer] -> (SetFace3D, StateVarsMBC S1 Vec3)` | Run 3D convex hull |
| `reRun` | -- | Re-run MBC with a pre-existing state and a new set of active faces |

Type aliases:

- `SetSimplex2D = IntMap (S2 Vec2)` -- a map of 2D triangles
- `SetSimplex3D = IntMap (S2 Vec3)` -- a map of 3D tetrahedra
- `SetFace3D = IntMap (S1 Vec3)` -- a map of 3D triangular faces (hull)

#### `DeUni.Types`

Defines the core data structures and type classes.

**Key types:**

- `PointPointer = Int` -- index into the point vector
- `WPoint p` -- a weighted point with fields `weight :: Double` and `point :: p Double`
- `SetPoint v = Vector (WPoint v)` -- the input point set
- `ActiveSubUnit simplex v` -- an active sub-unit (edge or face) with its associated opposite point and outward normal
- `StateVarsMBC simplex v` -- the MBC state: three AFLs, external faces set, simplex counter, and the point set

**Key type classes:**

- `PointND p` -- defines associated data types for a dimension: `Box p`, `Plane p`, `S0 p` (vertex), `S1 p` (edge/face), `S2 p` (triangle/tetrahedron). Methods include `isInBox`, `calcPlane`, `cutBox`, `touchPlane`, `circumOrigin`, `circumRadius`, etc.
- `Buildable simplex v` -- defines how to build simplices: `buildUnit`, `build1stUnit`, `getAllSubUnits`, `subUnitPos`.

#### `DeUni.FirstSeed`

Finds the first valid edge spanning the division plane, which is the seed for the first simplex. The algorithm finds the most distant points from the cutting plane in each half-space and "climbs" outward until the edge lies on the convex hull boundary.

#### `DeUni.GeometricTools`

Geometric utility functions: `projAonB`, `normalofAtoB`, `whichSideOfPlane`, `pointSetPartition`, `powerDist`, etc.

### 2D modules

- **`DeUni.Dim2.Base2D`** -- `PointND Vec2` instance defining `Box2D`, `Plane2D`, `Edge2D`, `Face2D` (triangle with circumcircle)
- **`DeUni.Dim2.Delaunay2D`** -- `Buildable S2 Vec2` instance for 2D Delaunay triangulation
- **`DeUni.Dim2.ReTri2D`** -- Circumcircle computation using QR decomposition

### 3D modules

- **`DeUni.Dim3.Base3D`** -- `PointND Vec3` instance defining `Box3D`, `Plane3D`, `Edge3D`, `Face3D`, `Tetrahedron` (with circumsphere)
- **`DeUni.Dim3.Delaunay3D`** -- `Buildable S2 Vec3` instance for 3D Delaunay tetrahedralization
- **`DeUni.Dim3.Hull3D`** -- `Buildable S1 Vec3` instance for 3D convex hull construction
- **`DeUni.Dim3.ReTri3D`** -- Circumsphere computation using QR decomposition

## Key Data Types Summary

| Type | Dimension | Description |
|---|---|---|
| `WPoint v` | Any | Weighted point: `{ weight :: Double, point :: v Double }` |
| `Box2D` | 2D | Axis-aligned bounding box |
| `Box3D` | 3D | Axis-aligned bounding box |
| `Edge2D` | 2D | Edge (S1 in 2D): pair of point indices |
| `Face2D` | 2D | Triangle (S2 in 2D) with circumcircle |
| `Edge3D` | 3D | Edge (S0 in 3D): pair of point indices |
| `Face3D` | 3D | Triangular face (S1 in 3D): triple of point indices |
| `Tetrahedron` | 3D | Tetrahedron (S2 in 3D) with circumsphere |

## Usage Examples

### 2D Delaunay Triangulation

```haskell
import qualified Data.Vector as Vec
import DeUni.DeWall

box :: Box Vec2
box = Box2D { xMax2D = 100, xMin2D = -100, yMax2D = 100, yMin2D = -100 }

points :: SetPoint Vec2
points = Vec.fromList
  [ WPoint 0 (Vec2 10 20)
  , WPoint 0 (Vec2 (-30) 40)
  , WPoint 0 (Vec2 50 (-10))
  , WPoint 0 (Vec2 (-20) (-50))
  , WPoint 0 (Vec2 70 60)
  ]

indices :: [PointPointer]
indices = [0 .. Vec.length points - 1]

(triangles, finalState) = runDelaunay2D box points indices
-- triangles :: IntMap (S2 Vec2)
-- Each Face2D contains face2DPoints, circleRadius, circleCenter
```

### 3D Delaunay Triangulation

```haskell
import qualified Data.Vector as Vec
import DeUni.DeWall

box3d :: Box Vec3
box3d = Box3D
  { xMax3D = 100, xMin3D = -100
  , yMax3D = 100, yMin3D = -100
  , zMax3D = 100, zMin3D = -100
  }

points3d :: SetPoint Vec3
points3d = Vec.fromList
  [ WPoint 0 (Vec3 10 20 30)
  , WPoint 0 (Vec3 (-30) 40 10)
  , WPoint 0 (Vec3 50 (-10) 25)
  , WPoint 0 (Vec3 (-20) (-50) 15)
  , WPoint 0 (Vec3 70 60 (-40))
  ]

(tetrahedra, state) = runDelaunay3D box3d points3d [0..4]
```

### 3D Convex Hull

```haskell
import DeUni.DeWall

(hullFaces, hullState) = runHull3D box3d points3d [0..4]
-- hullFaces :: IntMap (S1 Vec3)
-- Each Face3D contains face3DPoints :: (PointPointer, PointPointer, PointPointer)
```

### Regular (Weighted) Triangulation

To compute a regular triangulation instead of standard Delaunay, set non-zero weights on the points:

```haskell
weightedPoints = Vec.fromList
  [ WPoint 5.0 (Vec3 10 20 30)   -- sphere of radius sqrt(5)
  , WPoint 2.0 (Vec3 (-30) 40 10)
  , WPoint 8.0 (Vec3 50 (-10) 25)
  ]
```

## Architecture Notes

The simplex hierarchy is encoded via associated type families in the `PointND` class:

- **S0** -- 0-simplex relative to the target: a point in 2D, an edge in 3D
- **S1** -- 1-simplex relative to the target: an edge in 2D, a triangular face in 3D
- **S2** -- 2-simplex (the target): a triangle in 2D, a tetrahedron in 3D

The `Buildable` class connects these levels: `Sub S2 = S1`. This design allows the core `mbc` function in `DeWall` to be completely generic -- it works identically for 2D triangulation, 3D triangulation, and 3D hull construction.

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| `base` | >= 4, < 5 | Standard library |
| `containers` | >= 0.4 | `Data.Set`, `Data.IntMap` for AFLs and simplex storage |
| `linear-vect` | >= 0.2, < 0.3 | Vector/matrix algebra (`Vec2`, `Vec3`, QR decomposition) |
| `mtl` | >= 2.0 | `State` monad for MBC algorithm state |
| `random` | >= 1.0 | Random number generation |
| `vector` | >= 0.10 | Efficient point set storage |

## Building

### With Cabal

```bash
cabal build
```

### With Nix (development shell)

```bash
nix develop
cabal build
```

### Building tests and benchmarks

Tests and benchmarks are gated behind the `test` flag:

```bash
cabal build -f test
cabal run DeUni-check    # QuickCheck test suite
cabal run DeUni-benchmark  # Criterion benchmarks
```

## Author

Edgar Gomes de Araujo (<talktoedgar@gmail.com>)

## License

MIT -- see [LICENSE](./LICENSE).
