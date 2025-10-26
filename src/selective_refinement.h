#pragma once
// selective_refinement.h
// Header-only Selective Refinement (iterative):
// - Iterate: (1) build a refined graph around the current original-only path
//            (2) run Dijkstra on the refined graph
//            (3) collapse to original-only path for the next iteration
// - Stop when relative path-length change < eps or iteration reaches itmax.
// - Supports:
//    * ORIGINAL-edge subdivision (m steiner per target edge)
//    * ON-FACE rungs inside incident triangles
// - Returns an original-only path for compatibility with the current renderer.
// - Optional debug overlays from the **final iteration**:
//    * dbgSteinerPoints (cyan), dbgRefinedEdges (cyan)
//    * dbgOnFaceEdges (magenta), dbgSRPathPolyline (yellow)
//
// Notes:
// - Comments are in English for public release.

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <limits>
#include <algorithm>
#include <tuple>
#include <cmath>

#include "dijkstra.h"
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>

struct SRParams {
  float gamma   = 0.04f;  // not used in this step
  int   m       = 2;      // steiner points per subdivided edge (>=1 recommended)
  float eps     = 1e-3f;  // relative convergence threshold on path length
  int   itmax   = 10;     // maximum refinement iterations
  bool  addOnFace = true; // add ON-FACE rung edges
};

// Forward declarations matching app's types.
struct Mesh;
struct Graph;

// Undirected edge key
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

// Reverse map entry for steiner id -> original edge & segment id.
struct SR_SteinerLUT {
  int i = -1; // original endpoint
  int j = -1; // original endpoint
  int k = -1; // 0..m-1 (which subdivision)
  int m = 0;  // total steiner per edge
};

// Position of a refined vertex id (original or steiner).
inline glm::vec3 sr_vertex_position(const Mesh& M,
                                    int v,
                                    int n_orig,
                                    const std::vector<SR_SteinerLUT>* lut)
{
  if (v >= 0 && v < n_orig) return M.V[v];
  if (!lut) return glm::vec3(0.0f);
  const SR_SteinerLUT& e = (*lut)[v - n_orig];
  const float t = float(e.k + 1) / float(e.m + 1);
  return (1.0f - t) * M.V[e.i] + t * M.V[e.j];
}

// Build refined graph around an original-only path: subdivide incident ORIGINAL edges,
// and optionally add ON-FACE rungs inside incident triangles.
// Optionally fills: 'lut' (steiner reverse map), 'dbgOnFaceSegs' (world-space rungs).
inline void sr_build_refined_graph_edge_subdivision(const Mesh& M,
                                                    const Graph& G0,
                                                    const std::vector<int>& path0,
                                                    int m,
                                                    Graph& Gref,
                                                    int& n_orig_out,
                                                    int& n_ref_out,
                                                    std::vector<SR_SteinerLUT>* lut /*optional*/,
                                                    std::vector<std::pair<glm::vec3, glm::vec3>>* dbgOnFaceSegs /*optional*/,
                                                    bool addOnFace /*toggle*/)
{
  const int n = (int)G0.adj.size();
  n_orig_out = n;

  // Mark vertices on the current original-only path.
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
  if (dbgOnFaceSegs) dbgOnFaceSegs->clear();

  int next_idx = n; // first steiner index
  if (m_per_edge > 0) {
    for (const auto& ek : sub_edges) {
      steiner_base[ek] = next_idx;
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
  auto pos_of = [&](int idx)->glm::vec3 {
    return sr_vertex_position(M, idx, n, lut);
  };

  // Fill edges into Gref: original or subdivided chains.
  for (const auto& ek : all_edges) {
    const int i = ek.a, j = ek.b;
    const float w_ij = find_weight(i, j);

    auto it = sub_edges.find(ek);
    if (it == sub_edges.end() || m_per_edge == 0) {
      add_undirected(i, j, w_ij);
    } else {
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

  // ON-FACE rungs between aligned subdivision nodes (optional).
  if (addOnFace && m_per_edge > 0) {
    auto add_onface_between = [&](const SR_EdgeKey& e1, const SR_EdgeKey& e2){
      auto it1 = steiner_base.find(e1);
      auto it2 = steiner_base.find(e2);
      if (it1 == steiner_base.end() || it2 == steiner_base.end()) return;
      const int b1 = it1->second;
      const int b2 = it2->second;
      for (int k = 0; k < m_per_edge; ++k) {
        const int u = b1 + k;
        const int v = b2 + k;
        const glm::vec3 pu = pos_of(u);
        const glm::vec3 pv = pos_of(v);
        const float w = glm::length(pu - pv);
        add_undirected(u, v, w);
        if (dbgOnFaceSegs) dbgOnFaceSegs->push_back({pu, pv});
      }
    };

    for (const auto& tri : M.F) {
      const int a = tri.x, b = tri.y, c = tri.z;
      SR_EdgeKey ab(a,b), bc(b,c), ca(c,a);
      add_onface_between(ab, bc);
      add_onface_between(bc, ca);
      add_onface_between(ca, ab);
    }
  }
}

// Collapse refined path back to ORIGINAL vertex IDs.
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

// Compute polyline length in 3D for a refined path (may include steiner).
inline float sr_polyline_length(const Mesh& M,
                                const std::vector<int>& path_ref,
                                int n_orig,
                                const std::vector<SR_SteinerLUT>* lut)
{
  if (path_ref.size() < 2) return 0.0f;
  float L = 0.0f;
  glm::vec3 prev = sr_vertex_position(M, path_ref[0], n_orig, lut);
  for (size_t i = 1; i < path_ref.size(); ++i) {
    glm::vec3 cur = sr_vertex_position(M, path_ref[i], n_orig, lut);
    L += glm::length(cur - prev);
    prev = cur;
  }
  return L;
}

// Public API (iterative SR with eps/itmax):
inline bool selectiveRefinementPath(const Mesh& M,
                                    const Graph& G0,
                                    int s, int t,
                                    std::vector<int>& outPath,
                                    const SRParams& P = {},
                                    // Optional debug outputs (final iteration)
                                    std::vector<glm::vec3>* dbgSteinerPoints = nullptr,
                                    std::vector<std::pair<glm::vec3, glm::vec3>>* dbgRefinedEdges = nullptr,
                                    std::vector<glm::vec3>* dbgSRPathPolyline = nullptr,
                                    std::vector<std::pair<glm::vec3, glm::vec3>>* dbgOnFaceEdges = nullptr)
{
  outPath.clear();
  const int n = (int)G0.adj.size();
  if (n == 0 || s < 0 || s >= n || t < 0 || t >= n) return false;

  // Initial path on the original graph
  std::vector<int> path_orig;
  if (!shortestPath(G0, s, t, path_orig)) return false;

  // Iteration state
  float lastL = std::numeric_limits<float>::infinity();
  std::vector<int> best_path_orig = path_orig;

  // Buffers for the last iteration (for debug export)
  std::vector<SR_SteinerLUT> last_lut;
  std::vector<glm::vec3> last_dbgSteiner;
  std::vector<std::pair<glm::vec3, glm::vec3>> last_dbgRefined;
  std::vector<std::pair<glm::vec3, glm::vec3>> last_dbgOnFace;
  std::vector<glm::vec3> last_dbgSRpoly;

  const int itmax = std::max(1, P.itmax);
  const int m_use = std::max(1, P.m);
  const float eps = std::max(0.0f, P.eps);

  for (int it = 0; it < itmax; ++it) {
    // 1) Build refined graph around current original-only path
    Graph Gref;
    int n_orig = 0, n_ref = 0;
    std::vector<SR_SteinerLUT> lut;
    std::vector<std::pair<glm::vec3, glm::vec3>> onface_debug;

    sr_build_refined_graph_edge_subdivision(
      M, G0, path_orig, m_use,
      Gref, n_orig, n_ref,
      (dbgSteinerPoints || dbgRefinedEdges || dbgSRPathPolyline || dbgOnFaceEdges) ? &lut : nullptr,
      (dbgOnFaceEdges ? &onface_debug : nullptr),
      P.addOnFace
    );

    // 2) Run Dijkstra on refined graph
    std::vector<int> path_ref;
    if (!shortestPath(Gref, s, t, path_ref)) {
      // If refined graph fails, keep the last good result
      break;
    }

    // 3) Measure refined path length (3D)
    const float curL = sr_polyline_length(M, path_ref, n_orig,
                          (dbgSteinerPoints || dbgRefinedEdges || dbgSRPathPolyline || dbgOnFaceEdges) ? &lut : nullptr);

    // Keep debug buffers from this iteration (overwrite each round)
    if (dbgSteinerPoints || dbgRefinedEdges || dbgSRPathPolyline || dbgOnFaceEdges) {
      // steiner points
      last_dbgSteiner.clear();
      last_dbgSteiner.reserve(lut.size());
      for (size_t idx = 0; idx < lut.size(); ++idx) {
        last_dbgSteiner.push_back(sr_vertex_position(M, (int)(n_orig + idx), n_orig, &lut));
      }
      // refined chains: rebuild assignment order consistent with builder
      last_dbgRefined.clear();
      {
        std::unordered_set<SR_EdgeKey, SR_EdgeKeyHash> all_edges, sub_edges;
        std::vector<char> on_path(n_orig, 0);
        for (int v : path_orig) if (0 <= v && v < n_orig) on_path[v] = 1;
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
        std::unordered_map<SR_EdgeKey, int, SR_EdgeKeyHash> steiner_base;
        int next_idx = n_orig;
        for (const auto& ek : sub_edges) {
          steiner_base[ek] = next_idx;
          next_idx += m_use;
        }
        auto pos_of = [&](int idx)->glm::vec3 { return sr_vertex_position(M, idx, n_orig, &lut); };
        for (const auto& ek : sub_edges) {
          const int i = ek.a, j = ek.b;
          const int base = steiner_base[ek];
          int prev = i;
          for (int k = 0; k < m_use; ++k) {
            int cur = base + k;
            last_dbgRefined.push_back({ pos_of(prev), pos_of(cur) });
            prev = cur;
          }
          last_dbgRefined.push_back({ pos_of(prev), pos_of(j) });
        }
      }
      // on-face
      last_dbgOnFace = std::move(onface_debug);
      // SR polyline
      last_dbgSRpoly.clear();
      last_dbgSRpoly.reserve(path_ref.size());
      for (int v : path_ref) last_dbgSRpoly.push_back(sr_vertex_position(M, v, n_orig, &lut));
      // store lut for potential external use (not exported directly)
      last_lut = std::move(lut);
    }

    // 4) Collapse refined path to original-only for the next iteration
    std::vector<int> next_path_orig;
    sr_collapse_to_original_vertices(path_ref, n_orig, next_path_orig);
    if (next_path_orig.size() < 2) {
      // degenerate; stop with current best
      break;
    }

    // 5) Convergence test (relative change)
    if (std::isfinite(lastL)) {
      const float rel = std::fabs(curL - lastL) / std::max(curL, 1e-6f);
      if (rel < eps) {
        best_path_orig = next_path_orig;
        break;
      }
    }
    // prepare next round
    lastL = curL;
    best_path_orig = next_path_orig;
    path_orig.swap(next_path_orig);
  }

  // Output the best original-only path
  outPath = best_path_orig;

  // Export debug overlays from the final iteration
  if (dbgSteinerPoints)   *dbgSteinerPoints   = std::move(last_dbgSteiner);
  if (dbgRefinedEdges)    *dbgRefinedEdges    = std::move(last_dbgRefined);
  if (dbgOnFaceEdges)     *dbgOnFaceEdges     = std::move(last_dbgOnFace);
  if (dbgSRPathPolyline)  *dbgSRPathPolyline  = std::move(last_dbgSRpoly);

  return (outPath.size() >= 2);
}
