#pragma once
// selective_refinement.h
// Header-only Selective Refinement (phase 1):
// - Build initial path on the original 1-skeleton via Dijkstra
// - Around that path, subdivide ORIGINAL edges into m equal segments (Steiner points)
// - Build a refined graph Gref that replaces only those target edges with chains
// - Run Dijkstra on Gref
// - For compatibility with the current renderer, collapse the refined path back to
//   ORIGINAL vertex indices (steiner nodes are skipped in the returned path).
//
// Notes:
// - This header keeps comments in English for public release.
// - Minimal-change philosophy: geodesic.cpp remains unmodified at this step.
// - Next step for visualization: expose the steiner nodes/segments to draw them.

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <limits>
#include <algorithm>

#include "dijkstra.h"       // shortestPath(...)
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>

// Params for Selective Refinement
struct SRParams {
  float gamma = 0.04f;  // not used in phase 1 (kept for future spacing logic)
  int   m     = 1;      // number of Steiner points per target edge (m >= 1)
  float eps   = 1e-3f;  // not used in phase 1
  int   itmax = 10;     // not used in phase 1
  bool  addOnFace = true; // not used in phase 1
};

// Forward declarations; the full definitions are provided in geodesic.cpp
struct Mesh;
struct Graph;

// Local edge key (undirected) for hashing/uniqueness inside this header.
struct SR_EdgeKey {
  int a, b;
  SR_EdgeKey() : a(-1), b(-1) {}
  SR_EdgeKey(int i, int j) {
    if (i < j) { a = i; b = j; }
    else       { a = j; b = i; }
  }
  bool operator==(const SR_EdgeKey& o) const { return a == o.a && b == o.b; }
};
struct SR_EdgeKeyHash {
  std::size_t operator()(const SR_EdgeKey& k) const noexcept {
    // A simple pair-hash. Good enough for our use.
    return (std::hash<int>()(k.a) * 73856093u) ^ (std::hash<int>()(k.b) * 19349663u);
  }
};

// Build a refined graph around an initial path by subdividing ORIGINAL edges
// incident to the path's vertices into chains with m steiner points.
inline void sr_build_refined_graph_edge_subdivision(const Mesh& M,
                                                    const Graph& G0,
                                                    const std::vector<int>& path0,
                                                    int m,
                                                    Graph& Gref,
                                                    int& n_orig_out,
                                                    int& n_ref_out)
{
  const int n = static_cast<int>(G0.adj.size());
  n_orig_out = n;

  // 1) Collect the set of vertices that lie on the initial path.
  std::vector<char> on_path(n, 0);
  for (int v : path0) if (0 <= v && v < n) on_path[v] = 1;

  // 2) Collect UNIQUE edges from the original graph and mark which to subdivide.
  std::unordered_set<SR_EdgeKey, SR_EdgeKeyHash> all_edges;
  std::unordered_set<SR_EdgeKey, SR_EdgeKeyHash> sub_edges; // edges to subdivide
  all_edges.reserve(n * 3);
  sub_edges.reserve(n);

  // Helper to read an edge weight (i->j). We assume an undirected graph mirrored in adj.
  auto find_weight = [&](int i, int j) -> float {
    for (const auto& p : G0.adj[i]) if (p.first == j) return p.second;
    // Fallback: compute Euclidean length from mesh (should not happen if adj is consistent).
    return glm::length(M.V[i] - M.V[j]);
  };

  for (int u = 0; u < n; ++u) {
    for (const auto& kv : G0.adj[u]) {
      const int v = kv.first;
      if (v < 0 || v >= n || v == u) continue;
      SR_EdgeKey ek(u, v);
      if (all_edges.insert(ek).second) {
        // inserted as a unique undirected edge
        if (on_path[u] || on_path[v]) {
          sub_edges.insert(ek); // mark for subdivision (incident to the path)
        }
      }
    }
  }

  // 3) Prepare refined adjacency with extra nodes.
  const int m_per_edge = std::max(0, m);
  const int num_sub = static_cast<int>(sub_edges.size());
  const int num_new = m_per_edge * num_sub;
  const int n_ref = n + num_new;
  n_ref_out = n_ref;

  Gref.adj.clear();
  Gref.adj.resize(n_ref);

  // Map each sub-edge -> starting index of its steiner nodes (contiguous block of size m)
  std::unordered_map<SR_EdgeKey, int, SR_EdgeKeyHash> steiner_base;
  steiner_base.reserve(num_sub);

  // 4) Assign indices for Steiner nodes.
  int next_idx = n; // first new vertex index
  if (m_per_edge > 0) {
    for (const auto& ek : sub_edges) {
      steiner_base[ek] = next_idx;
      next_idx += m_per_edge;
    }
  }

  // 5) Add edges:
  //    - For edges NOT marked for subdivision: copy as-is into Gref.
  //    - For edges marked: replace with a chain i - s1 - s2 - ... - sm - j,
  //      with equal segment weights (original length / (m+1)).
  auto add_undirected = [&](int a, int b, float w) {
    Gref.adj[a].push_back({b, w});
    Gref.adj[b].push_back({a, w});
  };

  for (const auto& ek : all_edges) {
    const int i = ek.a, j = ek.b;
    const float w_ij = find_weight(i, j);

    auto it = sub_edges.find(ek);
    if (it == sub_edges.end() || m_per_edge == 0) {
      // Not subdivided: keep original edge.
      add_undirected(i, j, w_ij);
    } else {
      // Subdivide into (m+1) equal segments.
      const float seg_w = w_ij / float(m_per_edge + 1);
      int prev = i;
      int base = steiner_base[ek];
      // Create a chain of m nodes
      for (int k = 0; k < m_per_edge; ++k) {
        int cur = base + k;
        add_undirected(prev, cur, seg_w);
        prev = cur;
      }
      // Connect last steiner to j
      add_undirected(prev, j, seg_w);
    }
  }
}

// Collapse a refined path (with potential steiner nodes) back to original vertex IDs
// so that the current renderer (which only knows original V) can draw without changes.
// We simply drop any vertex >= n_orig. If consecutive kept vertices are equal after
// dropping, we coalesce.
inline void sr_collapse_to_original_vertices(const std::vector<int>& path_ref,
                                             int n_orig,
                                             std::vector<int>& path_orig_only)
{
  path_orig_only.clear();
  int last = -1;
  for (int v : path_ref) {
    if (v >= 0 && v < n_orig) {
      if (path_orig_only.empty() || v != last) {
        path_orig_only.push_back(v);
        last = v;
      }
    }
  }
  // Ensure at least endpoints appear if possible.
  if (path_orig_only.size() < 2 && !path_ref.empty()) {
    // nothing to do here; upstream will re-sample endpoints if unreachable
  }
}

// Public API (phase 1): subdivide ORIGINAL edges around the initial path and run Dijkstra
// on the refined graph. For compatibility, we return a path that contains only ORIGINAL
// vertex IDs (steiner nodes are not exposed yet).
inline bool selectiveRefinementPath(const Mesh& M,
                                    const Graph& G0,
                                    int s, int t,
                                    std::vector<int>& outPath,
                                    const SRParams& P = {})
{
  outPath.clear();
  const int n = static_cast<int>(G0.adj.size());
  if (n == 0 || s < 0 || s >= n || t < 0 || t >= n) return false;

  // Step 1: initial path on the original 1-skeleton.
  std::vector<int> path0;
  if (!shortestPath(G0, s, t, path0)) return false;

  // Step 2-3-4-5: build refined graph by subdividing edges incident to the initial path.
  Graph Gref;
  int n_orig = 0, n_ref = 0;
  sr_build_refined_graph_edge_subdivision(M, G0, path0, std::max(1, P.m), Gref, n_orig, n_ref);

  // Step 6: run Dijkstra on the refined graph.
  std::vector<int> path_ref;
  if (!shortestPath(Gref, s, t, path_ref)) {
    // Fallback to original path (should rarely happen)
    outPath = path0;
    return true;
  }

  // Step 7: collapse refined path back to ORIGINAL vertex IDs for current renderer.
  sr_collapse_to_original_vertices(path_ref, n_orig, outPath);
  if (outPath.size() < 2) {
    // If everything collapses away (e.g., path entirely through steiner nodes),
    // return the original-path endpoints as a safe fallback.
    outPath = path0;
  }
  return true;
}
