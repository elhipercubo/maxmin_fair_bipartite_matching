void Graph::recompute_heights() {
#ifdef BUCKETS
//        printf("recomputing heights\n");
    int n = num_vertices();
    buckets.assign(n, n, &label[0], &excess[0]);
//        printf("done recomputing heights\n");
#endif

#ifdef CURRENT
    current.resize(num_vertices());
    for (int u = 0; u < num_vertices(); ++u) current[u] = _first[u];
#endif
}

void Graph::init_flow(bool use_fifo) {
#if !defined(BUCKETS) || defined(FIFO)
    use_fifo = true;
#endif
    using_fifo = use_fifo;

    int n = num_vertices();

    label.assign(n, 0);
    excess.assign(n, 0);
    label[src] = 1;                   // temporary so that push doesn't complain
    recompute_heights();
    forEachEdge(src, e) if (is_forward_edge(e)) {
        // TODO: we're ignoring vertices with out degree 0
        if (deg(e->v) && e->res_cap > 0) {
            push(src, e, e->res_cap);
        }    
    }
    label[src] = n;

    for (int u = 0; u < num_left; ++u) label[u] = 2;
    for (int u = num_left; u < src; ++u) label[u] = 1;
    label[src] = n; label[sink] = 0;
    recompute_heights();
    initialized = true;
    phase = 0;
}

bool Graph::push(int u, Edge* e, Captype amount) {
    assert(amount > 0);
    assert(e->res_cap > 0);
    assert(u == src || label[e->v] == label[e->rev->v] - 1);

    ++pushes;

    e->res_cap -= amount;
    e->rev->res_cap += amount;

    if ( ( excess[e->v] += amount) == amount && e->v != sink ) {
        activate(e->v);
    }
    if ( ( excess[u] -= amount ) == 0 ) {
        ++sat_pushes;
        // Saturating push. DO NOT deactivate v; has already been done by max_flow()!
        return true;
    }
    return false;
}

void Graph::relabel(int u, int d) {
    if (d > num_vertices()) d = num_vertices();
    assert(label[u] < d);
    ++relabels;
#ifdef BUCKETS
    int cnt = buckets.relabel(u, d);
    gap_node_cnt += cnt;
    if (cnt > 0) {
        assert(label[u] == num_vertices());
        ++num_gaps;
    }
#else
    label[u] = d;
#endif
}

void Graph::discharge(int u) {
    while (true) {
//        printf("%i\n", buckets.num_included);
        assert(buckets.num_included >= 0);
//        if (label[u] >= buckets.num_included) { if (label[u] < num_vertices()) relabel(u, num_vertices()); ++skipped; return; }
        if (label[u] >= num_vertices()/* - skipped*/) { ++skipped; return; }

        int d = num_vertices();

        ++work;
        ++work_since_gap;
#ifdef CURRENT
        Edge *l = current[u];
        for (; current[u] < _first[u + 1]; ++current[u]) {
            auto e = current[u];
            ++work;
            ++work_since_gap;
            ++edges_looked;
            if (is_admissible(u, e)) {
                if (push(u, e, min(e->res_cap, excess[u]) )) return;
            } else {
                if (e->res_cap > 0) d = min(d, label[e->v]);
            }
        }
        // Now d contains the minimum label of a neighbour after or at _first[u].
        // Update d to contain the previous neighbours as well. 
        for (Edge* e = _first[u]; e < l; ++e) {
            if (e->res_cap > 0) d = min(d, label[e->v]);
            ++edges_looked;
        }
        current[u] = _first[u];
#else
        forEachEdge(u, e) {
            ++work;
            ++work_since_gap;
            ++edges_looked;
            if (is_admissible(u, e)) {
                if (push(u, e, min(e->res_cap, excess[u]) )) return;
            } else {
                if (e->res_cap > 0) d = min(d, label[e->v]);
            }
        }
#endif

        // Now d contains the minimum label of a neighbour
        relabel(u, d + 1);
    }
}

void Graph::print_stats() {
    if (using_fifo) printf("====== Using FIFO ========\n"); else            printf("====== Using MAXD ========\n");
    set<int> h;
    for (int u = 0; u < num_vertices(); ++u) h.insert(label[u]);
    printf("Statistics:\n");
    double t = ((std::chrono::steady_clock::now() - start)).count() * 1e-9;
	cout << "  time elapsed: " <<  t << " seconds" << endl;
    printf("  flow achieved = %lf %lli\n", (double)excess[sink], excess[sink]);
    printf("  flow I'd like to achieve = %lf\n\n", (double)-excess[src]);
    printf("  flow achieved per time = %lf\n", (double)excess[sink] / t);
    printf("  work = %lli\n", work);

    printf("  saturating pushes = %'i\n", sat_pushes);
    printf("  non-saturating pushes = %'i\n", pushes - sat_pushes);
    printf("  pushes = %'i\n", pushes);
    printf("  relabels = %'i\n", relabels);
    printf("  activations = %'i\n", activations);
    printf("  skipped activations = %'i (%lf)\n", skipped, (double)skipped / activations);
    printf("  edges looked at = %'lli (%lf)\n", edges_looked, double(edges_looked) / edges.size());
    printf("  different heights = %'i\n", (int)h.size());
    printf("  num gaps = %'i\n", num_gaps);
    printf("  gaps node count = %'i\n", gap_node_cnt);
    printf("  num updates = %'i\n", num_updates);

    assert(-excess[src] >= excess[sink]);
}

void Graph::remove_from_queue(int u) {
    if (using_fifo)
        active.pop_front();
}

void Graph::global_update() {
    ++num_updates;

    last_relabels = relabels;
    work = work_since_gap = 0;

    int n = num_vertices();
    int infinity = num_vertices();
    fill(&label[0], &label[n], infinity);

    int* q = (int*) malloc(n * sizeof(int));
    auto qstart = q, qend = q;

    int v = phase == 0 ? sink : src;
    label[v] = 0;
    *qend++ = v;
    if (phase == 0) {
        while (qstart != qend) {
            int u = *qstart++, d = label[u] + 1;
            forEachEdge(u, e)
                if (e->rev->res_cap > 0) {
                    ++edges_looked;
                    int v = e->v;
                    if (label[v] > d) {
                        *qend++ = v;
                        label[v] = d;
                    }
                }
        }
        recompute_heights();
    } else {
        while (qstart != qend) {
            int u = *qstart++, d = label[u] + 1;
            forEachEdge(u, e)
                if (e->res_cap > 0) {
                    ++edges_looked;
                    int v = e->v;
                    if (label[v] > d) {
                        *qend++ = v;
                        label[v] = d;
                    }
                }
        }
    }

    recompute_heights();
    free(q);
}

Captype Graph::max_flow(int* c1, int* c2, vector<int>& dead) {
    if (!initialized) { puts("not initialized!\n"); exit(1); }

    work = work_since_gap = skipped = 0;
    int n = num_vertices();

    bool do_updates = false;
#if defined(GLOBAL_HEURISTIC) || !defined(BUCKETS)
    do_updates = true;
#endif
    while (some_active()) {
        int u = select_active();
#ifdef BUCKETS
        assert(buckets.aMax == label[u]);
#endif
        remove_from_queue(u);
        if (label[u] < num_vertices()) {
            assert(excess[u] > 0);
#ifdef BUCKETS
            assert(label[u] <= buckets.aMax);
            buckets.deactivate(u);
#endif
            discharge(u);
            assert(excess[u] == 0 || label[u] == num_vertices());
            if (do_updates && work > 3.3 * (long long)edges.size() && relabels > last_relabels + n / 10) global_update();
        }
        if (!(excess[u] == 0 || label[u] == num_vertices())) {
            buckets.print();
            printf("%i %lf %i\n", u, (double)excess[u] , label[u]);
            exit(1);
        }

        assert(excess[u] == 0 || label[u] == num_vertices());
    }

    Captype flow = excess[sink];

    // Compute distances from src
    global_update();

    int s_left = 0, s_right = 0, t_left = 0, t_right = 0;
    int d = n - 1;
    for (int u = 0; u < n; ++u) if (u != src && u != sink) {
        if      (label[u] > d && u < num_left) ++s_left;
        else if (label[u] > d && u >= num_left) ++s_right;
        else if (label[u] < d && u < num_left) ++t_left;
        else if (label[u] < d && u >= num_left) ++t_right;
    }
//    printf("   Left: source %i sink %i\n", s_left, t_left);
//    printf("   Right: source %i sink %i\n", s_right, t_right);
//    printf("Min vertex cover: sink_left + source_right = %i\n", t_left + s_right);  // all neighbours of s_left are in s_right
    if (s_left != 0) {
        *c1 = s_left; *c2 = s_right;
//        printf("   lambda <= %i / %i = %lf\n", *c2, *c1, (double)(*c2) / *c1);

        dead.clear();
        for (int u = 0; u < n; ++u) if (u != src && u != sink) {
//            if (is_dead[u]) continue;

            if ( (label[u] > d && u < num_left) || (label[u] > d && u >= num_left) )
                dead.push_back(u);
        }
    }

    return flow;
}




/*
void Graph::clearFlow() {
    // Restore capacities
    for (int u = 0; u < num_vertices(); ++u)
        forEachEdge(u, e)
            e->res_cap = capacity(e);
}
*/

