#pragma once
// selective_refinement.h
// Header-only "Selective Refinement" stub that currently delegates to Dijkstra.
// The goal is to keep the call shape identical to Dijkstra so the renderer/UI
// can switch/overlay methods without code changes later.
//
// Future work inside this header:
//  - Add Steiner points on ORIGINAL edges around an initial path
//  - Optionally add ON-FACE edges (diagonals) inside incident faces
//  - Iterate until convergence (|L_i - L_{i-1}| < eps) or itmax
//
// NOTE: Comments are in English for public release.

#include <vector>
#include "dijkstra.h"  // reuse shortestPath(...) for the initial stub

struct SRParams {
  // Initial spacing relative to model scale (e.g., bbox diagonal ratio).
  float gamma = 0.04f;
  // Number of Steiner points per ORIGINAL edge around the current path.
  int   m     = 1;
  // Convergence threshold on path length change.
  float eps   = 1e-3f;
  // Maximum number of refinement iterations.
  int   itmax = 10;
  // Whether to add ON-FACE edges (diagonals) inside incident faces.
  bool  addOnFace = true;
};

// Forward declarations to avoid heavy includes here.
// geodesic.cpp defines these types and includes us after full definitions.
struct Mesh;
struct Graph;

// Public API: same shape as Dijkstra for easy swapping/overlay in the UI.
// For now, this just calls Dijkstra on the original 1-skeleton graph.
// Later, this function will build a refined graph (G_ref) around the path.
inline bool selectiveRefinementPath(const Mesh& M,
                                    const Graph& G0,
                                    int s, int t,
                                    std::vector<int>& outPath,
                                    const SRParams& P = {}) {
  (void)M; (void)P; // unused in the stub
  // Delegation: use the baseline Dijkstra result for now.
  return shortestPath(G0, s, t, outPath);
}