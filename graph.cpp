#include "graph.h"

#define forEachEdge(u, e) for (auto e = _first[(u)]; e < _first[(u) + 1]; ++e)

void Graph::buildFromBipartite() {
    int a, b;
    num_left = 0; num_right = 0;
    while (read_int(&a) == 1 && read_int(&b) == 1) {
        --a; --b;
        _edges.push_back(make_pair(a, b));
        num_left = max(num_left, a + 1);
        num_right = max(num_right, b + 1);
    }

    // Add source and sink
    int num_edges = _edges.size();
    int n = num_left + num_right + 2;
    src = n - 2; sink = n - 1;
    printf("n - 2 = %'i num_left = %'i num_right = %'i m = %'i\n", n - 2, num_left, num_right, num_edges);

    // Relabel right vertices
    for (int i = 0; i < num_edges; ++i) _edges[i].second += num_left;
    for (int i = 0; i < num_edges; ++i) { assert( _edges[i].first >= 0 && _edges[i].first < num_left ); assert( _edges[i].second >= 0 && _edges[i].second < n); }

    // Add edges from source and to sink
    for (int u = 0; u < num_left; ++u) _edges.push_back(pii(src, u));
    for (int u = num_left; u < src; ++u) _edges.push_back(pii(u, sink));
}

// Initializes _first, capacities, and edges; reserves labeling space
void Graph::prepare() {
    int num_edges = _edges.size();

    printf("Count-sorting the edges...\n");
    int n = num_left + num_right + 2;
    vector<int> cnt(n);
    for (int i = 0; i < num_edges; ++i) {
        ++cnt[ _edges[i].first ];       // one direct from u
        ++cnt[ _edges[i].second ];      // one reverse from v
    }
    for (int u = 1; u < n; ++u) cnt[u] += cnt[u - 1];
    assert( cnt[n - 1] == 2 * num_edges );

    printf("Reserving...\n");
    edges.resize(2 * num_edges);
    capacities.resize(2 * num_edges);
    _first.clear();
    _first.push_back(&edges[0]);
    for (int c:cnt) _first.push_back(&edges[c]);

    printf("Sorting...\n");
    for (int i = num_edges - 1; i >= 0; --i) {
        if ((num_edges - i) % 5000000 == 0) { printf("  ... %'i", num_edges - i); fflush(stdout); }
        int u = _edges[i].first, v = _edges[i].second;

        int pos = --cnt[u], pos_rev = --cnt[v];
        edges[pos]     = Edge({ v, &edges[pos_rev], 1 });
        edges[pos_rev] = Edge({ u, &edges[pos    ], 0 });

        capacities[pos] = 1;
    }
    puts("");

    _edges.clear(); _edges.shrink_to_fit();
    label.resize(n);
    start = std::chrono::steady_clock::now();
}

void Graph::print_edge(const Edge* e) const {
    printf("u = %i v = %i cap = %lf flow = %lf\n", e->rev->v, e->v, (double)capacity(e), (double)edge_flow(e));
    e = e->rev;
    printf("u = %i v = %i cap = %lf flow = %lf\n", e->rev->v, e->v, (double)capacity(e), (double)edge_flow(e));
}
