#pragma once
// selective_refinement.h
// Header-only Selective Refinement (phase 1: ORIGINAL-edge subdivision).
// This version adds optional debug outputs to visualize:
//  - Steiner points (positions)
//  - Refined edges (as segments)
//  - The SR path polyline (including Steiner nodes)
// Notes:
//  - Comments are in English for public release.
//  - The app can keep using original-vertex-only paths for drawing,
//    while debug overlays visualize Steiner nodes & refined chains.

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <limits>
#include <algorithm>
#include <tuple>

#include "dijkstra.h"
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>

struct SRParams {
  float gamma = 0.04f; // not used in phase 1
  int   m     = 1;     // number of Steiner points per target edge (>=1)
  float eps   = 1e-3f; // not used in phase 1
  int   itmax = 10;    // not used in phase 1
  bool  addOnFace = true; // not used in phase 1
};

// Forward declarations for interop with the app's types.
struct Mesh;
struct Graph;

// Undirected edge key for hashing/uniqueness.
struct SR_EdgeKey {
  int a, b;
  SR_EdgeKey() : a(-1), b(-1) {}
  SR_EdgeKey(int i, int j) { if (i < j) { a = i; b = j; } else { a = j; b = i; } }
  bool operator==(const SR_EdgeKey& o) const { return a == o.a && b == o.b; }
};
struct SR_EdgeKeyHash {
  std::size_t operator()(const SR_EdgeKey& k) const noexcept {
    return (std::hash<int>()(k.a) * 73856093u) ^ (std::hash<int>()(k.b) * 19349663u);
  }
};

// Internal LUT entry: map a steiner node index to its ORIGINAL edge and segment id.
struct SR_SteinerLUT {
  int i = -1; // original endpoint
  int j = -1; // original endpoint
  int k = -1; // 0..m-1 (which segment node along the chain)
  int m = 0;  // total steiner per edge
};

// Helper to compute the position of a refined vertex index.
inline glm::vec3 sr_vertex_position(const Mesh& M,
                                    int v,
                                    int n_orig,
                                    const std::vector<SR_SteinerLUT>* lut)
{
  if (v >= 0 && v < n_orig) {
    return M.V[v];
  }
  if (!lut) return glm::vec3(0.0f); // should not happen if caller passes lut when needed
  const SR_SteinerLUT& e = (*lut)[v - n_orig];
  const float t = float(e.k + 1) / float(e.m + 1);
  return (1.0f - t) * M.V[e.i] + t * M.V[e.j];
}

// Build refined graph by subdividing ORIGINAL edges incident to the initial path.
// Also fill 'lut' for steiner-node -> (i,j,k,m) reverse mapping.
inline void sr_build_refined_graph_edge_subdivision(const Mesh& M,
                                                    const Graph& G0,
                                                    const std::vector<int>& path0,
                                                    int m,
                                                    Graph& Gref,
                                                    int& n_orig_out,
                                                    int& n_ref_out,
                                                    std::vector<SR_SteinerLUT>* lut /*optional*/)
{
  const int n = (int)G0.adj.size();
  n_orig_out = n;

  // Mark vertices on the initial path.
  std::vector<char> on_path(n, 0);
  for (int v : path0) if (0 <= v && v < n) on_path[v] = 1;

  // Collect unique undirected edges and mark those incident to the path.
  std::unordered_set<SR_EdgeKey, SR_EdgeKeyHash> all_edges;
  std::unordered_set<SR_EdgeKey, SR_EdgeKeyHash> sub_edges;
  all_edges.reserve(n * 3);
  sub_edges.reserve(n);

  auto find_weight = [&](int i, int j) -> float {
    for (const auto& p : G0.adj[i]) if (p.first == j) return p.second;
    return glm::length(M.V[i] - M.V[j]); // fallback (should not happen)
  };

  for (int u = 0; u < n; ++u) {
    for (const auto& kv : G0.adj[u]) {
      const int v = kv.first;
      if (v < 0 || v >= n || v == u) continue;
      SR_EdgeKey ek(u, v);
      if (all_edges.insert(ek).second) {
        if (on_path[u] || on_path[v]) sub_edges.insert(ek);
      }
    }
  }

  // Decide refined graph size.
  const int m_per_edge = std::max(0, m);
  const int num_sub = (int)sub_edges.size();
  const int num_new = m_per_edge * num_sub;
  const int n_ref = n + num_new;
  n_ref_out = n_ref;

  Gref.adj.clear();
  Gref.adj.resize(n_ref);

  // Map sub-edge -> block base of its steiner nodes.
  std::unordered_map<SR_EdgeKey, int, SR_EdgeKeyHash> steiner_base;
  steiner_base.reserve(num_sub);

  if (lut) lut->clear();
  if (lut && m_per_edge > 0) lut->resize(num_new);

  int next_idx = n; // first steiner index
  if (m_per_edge > 0) {
    for (const auto& ek : sub_edges) {
      steiner_base[ek] = next_idx;
      // Fill LUT entries for this edge
      if (lut) {
        for (int k = 0; k < m_per_edge; ++k) {
          const int global_idx = next_idx + k;
          const int lut_idx    = global_idx - n;
          (*lut)[lut_idx] = SR_SteinerLUT{ ek.a, ek.b, k, m_per_edge };
        }
      }
      next_idx += m_per_edge;
    }
  }

  auto add_undirected = [&](int a, int b, float w) {
    Gref.adj[a].push_back({b, w});
    Gref.adj[b].push_back({a, w});
  };

  // Fill edges into Gref.
  for (const auto& ek : all_edges) {
    const int i = ek.a, j = ek.b;
    const float w_ij = find_weight(i, j);

    auto it = sub_edges.find(ek);
    if (it == sub_edges.end() || m_per_edge == 0) {
      // Keep original edge.
      add_undirected(i, j, w_ij);
    } else {
      // Replace by chain i - s1 - ... - sm - j
      const float seg_w = w_ij / float(m_per_edge + 1);
      int prev = i;
      const int base = steiner_base[ek];
      for (int k = 0; k < m_per_edge; ++k) {
        const int cur = base + k;
        add_undirected(prev, cur, seg_w);
        prev = cur;
      }
      add_undirected(prev, j, seg_w);
    }
  }
}

// Collapse a refined path back to ORIGINAL vertex IDs (drop steiner nodes).
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
}

// Public API (phase 1):
// - Build refined graph by subdividing ORIGINAL edges incident to the initial path.
// - Run Dijkstra on the refined graph.
// - Return original-only path for compatibility,
// - Optionally output debug overlays (steiner points, refined edges, SR polyline).
inline bool selectiveRefinementPath(const Mesh& M,
                                    const Graph& G0,
                                    int s, int t,
                                    std::vector<int>& outPath,
                                    const SRParams& P = {},
                                    // Optional debug outputs (pass nullptr to ignore)
                                    std::vector<glm::vec3>* dbgSteinerPoints = nullptr,
                                    std::vector<std::pair<glm::vec3, glm::vec3>>* dbgRefinedEdges = nullptr,
                                    std::vector<glm::vec3>* dbgSRPathPolyline = nullptr)
{
  outPath.clear();
  const int n = (int)G0.adj.size();
  if (n == 0 || s < 0 || s >= n || t < 0 || t >= n) return false;

  // 1) Initial path on original graph
  std::vector<int> path0;
  if (!shortestPath(G0, s, t, path0)) return false;

  // 2) Refined graph + LUT for steiner nodes
  Graph Gref;
  int n_orig = 0, n_ref = 0;
  std::vector<SR_SteinerLUT> lut; // size = num_new steiner nodes
  sr_build_refined_graph_edge_subdivision(M, G0, path0, std::max(1, P.m), Gref, n_orig, n_ref,
                                          (dbgSteinerPoints || dbgRefinedEdges || dbgSRPathPolyline) ? &lut : nullptr);

  // 3) Run Dijkstra on refined graph
  std::vector<int> path_ref;
  if (!shortestPath(Gref, s, t, path_ref)) {
    outPath = path0; // fallback
    return true;
  }

  // 4) Collapse to original-only indices for current renderer
  sr_collapse_to_original_vertices(path_ref, n_orig, outPath);
  if (outPath.size() < 2) outPath = path0;

  // 5) Optional debug overlays ------------------------------------------------
  if (dbgSteinerPoints || dbgRefinedEdges || dbgSRPathPolyline) {
    // Build a quick reverse map for edges to know which are subdivided
    // (we can infer from LUT as well).
    // 5.1 Fill dbgSteinerPoints
    if (dbgSteinerPoints) {
      dbgSteinerPoints->clear();
      dbgSteinerPoints->reserve(lut.size());
      for (size_t idx = 0; idx < lut.size(); ++idx) {
        const glm::vec3 p = sr_vertex_position(M, (int)(n_orig + idx), n_orig, &lut);
        dbgSteinerPoints->push_back(p);
      }
    }
    // 5.2 Fill dbgRefinedEdges: for each subdivided edge (i,j), emit chain segments
    if (dbgRefinedEdges) {
      dbgRefinedEdges->clear();
      // Rebuild the same edge sets to know which were subdivided, with m = P.m
      // (We could also deduce from LUT by grouping blocks of size m).
      std::unordered_set<SR_EdgeKey, SR_EdgeKeyHash> all_edges, sub_edges;
      std::vector<char> on_path(n_orig, 0);
      for (int v : path0) if (0 <= v && v < n_orig) on_path[v] = 1;
      for (int u = 0; u < n_orig; ++u) {
        for (const auto& kv : G0.adj[u]) {
          int v = kv.first;
          if (v < 0 || v >= n_orig || v == u) continue;
          SR_EdgeKey ek(u, v);
          if (all_edges.insert(ek).second) {
            if (on_path[u] || on_path[v]) sub_edges.insert(ek);
          }
        }
      }
      const int m_per_edge = std::max(1, P.m);
      // We need a way to fetch the assigned base index for each edge.
      // Recompute the same assignment order as in sr_build_refined_graph_edge_subdivision:
      std::unordered_map<SR_EdgeKey, int, SR_EdgeKeyHash> steiner_base;
      int next_idx = n_orig;
      for (const auto& ek : sub_edges) {
        steiner_base[ek] = next_idx;
        next_idx += m_per_edge;
      }
      // Emit segments i - s1, s1 - s2, ..., sm - j for each subdivided edge.
      for (const auto& ek : sub_edges) {
        const int i = ek.a, j = ek.b;
        const int base = steiner_base[ek];
        // positions
        glm::vec3 p_i = M.V[i];
        glm::vec3 p_j = M.V[j];
        // chain positions
        auto pos_of = [&](int idx) { return sr_vertex_position(M, idx, n_orig, &lut); };
        int prev = i;
        for (int k = 0; k < m_per_edge; ++k) {
          int cur = base + k;
          dbgRefinedEdges->push_back({ pos_of(prev), pos_of(cur) });
          prev = cur;
        }
        dbgRefinedEdges->push_back({ pos_of(prev), pos_of(j) });
      }
    }
    // 5.3 Fill dbgSRPathPolyline: map refined path indices to 3D positions
    if (dbgSRPathPolyline) {
      dbgSRPathPolyline->clear();
      dbgSRPathPolyline->reserve(path_ref.size());
      for (int v : path_ref) {
        dbgSRPathPolyline->push_back(sr_vertex_position(M, v, n_orig, &lut));
      }
    }
  }

  return true;
}
