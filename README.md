# DeUni

Delaunay triangulation and convex hull computation via the Marriage Before Conquer algorithm.

## What is this?

A Haskell library for computing:

- **2D Delaunay triangulations** -- given a set of points in a plane, produce a triangulation where no point lies inside any triangle's circumcircle
- **3D Delaunay triangulations** (tetrahedralization) -- the 3D analogue: partition space into tetrahedra where no point lies inside any tetrahedron's circumsphere
- **3D convex hulls** -- the smallest convex surface enclosing a 3D point set

It also supports **weighted (power) distances**, which produce regular triangulations -- a generalization of Delaunay where each point has an associated weight that influences the tessellation. This is useful when points represent spheres or grains of different sizes.

## Background: what is a Delaunay triangulation?

Given a set of points, there are many possible triangulations. The Delaunay triangulation is special: it maximizes the minimum angle across all triangles, producing "well-shaped" triangles that avoid long thin slivers. This property makes it the preferred triangulation for finite element methods, mesh generation, and constructing Voronoi diagrams (the Voronoi diagram is the geometric dual of the Delaunay triangulation).

The defining property: for every triangle (or tetrahedron in 3D), the circumscribed circle (or sphere) contains no other points from the set. This is called the "empty circumsphere" property.

## The algorithm: Marriage Before Conquer (MBC)

DeUni implements a divide-and-conquer strategy called Marriage Before Conquer, a variant of the DeWall algorithm. The name comes from the unusual ordering: the "marriage" (merge) step happens *before* the recursive "conquer" step.

### How it works

1. **Divide:** A cutting plane bisects the bounding box along its longest axis, splitting points into two half-spaces (B1 and B2).

2. **Marry (merge first):** Before recursing, the algorithm builds all simplices (triangles in 2D, tetrahedra in 3D) that straddle the cutting plane. Starting from an initial seed simplex, it grows the set of crossing simplices by processing an Active Face List (AFL). Each new simplex's sub-faces are classified:
   - Faces crossing the cutting plane go into the "alpha" AFL (processed next)
   - Faces entirely in B1 or B2 go into their respective AFLs
   - Duplicate faces (shared by two simplices, thus closed) are removed

3. **Conquer:** Recurse into each half-space with its accumulated AFL of open faces.

4. **Termination:** Recursion ends when the AFL is empty.

### Circumsphere computation

To determine whether a candidate point should form a new simplex, the algorithm computes circumcircles (2D) or circumspheres (3D) using QR decomposition (Gram-Schmidt orthogonalization) of the edge vector matrix. This approach naturally supports weighted points for regular triangulations.

## Example

```haskell
import qualified Data.Vector as Vec
import DeUni.DeWall

let box = Box2D { xMax2D = 100, xMin2D = -100, yMax2D = 100, yMin2D = -100 }

    points = Vec.fromList
      [ WPoint 0 (Vec2 10 20)
      , WPoint 0 (Vec2 (-30) 40)
      , WPoint 0 (Vec2 50 (-10))
      , WPoint 0 (Vec2 (-20) (-50))
      ]

    (triangles, _) = runDelaunay2D box points [0 .. Vec.length points - 1]
    -- triangles :: IntMap Face2D
```

## Where is it used?

- **VirMat** -- the virtual microstructure generator uses DeUni to compute the Delaunay triangulation of grain seed points, then extracts the Voronoi dual to form grain boundaries. The triangulation is periodically recomputed during sphere packing iterations to maintain neighbor relationships.

## How to build

```bash
# With Nix (recommended)
nix develop
cabal build --allow-newer

# With Cabal
cabal build

# Tests and benchmarks (behind the test flag)
cabal build -f test
cabal run DeUni-check        # QuickCheck test suite
cabal run DeUni-benchmark    # Criterion benchmarks
```

## License

MIT -- see [LICENSE](./LICENSE).
