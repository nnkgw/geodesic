#pragma once
// dijkstra.h
// Header-only Dijkstra shortest path on a graph's 1-skeleton.
// Keeps the exact signature used in geodesic.cpp for minimal changes.
//
// Graph requirements (as defined in geodesic.cpp):
//   struct Graph {
//     std::vector<std::vector<std::pair<int,float>>> adj;
//   };
//
// Usage:
//   #include "dijkstra.h"   // after Graph is defined
//   std::vector<int> path;
//   bool ok = shortestPath(G, src, dst, path);
//
// Notes:
// - Indices are assumed to be 0-based.
// - Edge weights are non-negative.
// - This is a standard Dijkstra with a binary heap (priority_queue).

#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <algorithm>

struct Graph; // forward declaration; the full definition must be visible before including this header

// Core solver: identical signature to your existing function for a zero-diff call site.
// Returns:
//   - true  : path found; 'outPath' is filled from s -> t
//   - false : t is unreachable
// Complexity: O((V + E) log V) with a binary heap.
inline bool shortestPath(const Graph& G, int s, int t, std::vector<int>& outPath)
{
    const int n = static_cast<int>(G.adj.size());
    if (s < 0 || s >= n || t < 0 || t >= n || n == 0) {
        outPath.clear();
        return false;
    }

    const float INF = std::numeric_limits<float>::infinity();
    std::vector<float> dist(n, INF);
    std::vector<int>   prev(n, -1);

    using Node = std::pair<float,int>; // (distance, vertex)
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

    dist[s] = 0.0f;
    pq.push({0.0f, s});

    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u]) continue;     // stale entry
        if (u == t) break;             // early exit when target is popped

        // relax all outgoing edges
        for (const auto& nbr : G.adj[u]) {
            const int v = nbr.first;
            const float w = nbr.second; // non-negative edge weight assumed
            const float nd = d + w;
            if (nd < dist[v]) {
                dist[v] = nd;
                prev[v] = u;
                pq.push({nd, v});
            }
        }
    }

    if (!std::isfinite(dist[t])) {
        outPath.clear();
        return false; // unreachable
    }

    // Reconstruct path t -> s, then reverse
    outPath.clear();
    for (int cur = t; cur != -1; cur = prev[cur]) outPath.push_back(cur);
    std::reverse(outPath.begin(), outPath.end());
    return true;
}
