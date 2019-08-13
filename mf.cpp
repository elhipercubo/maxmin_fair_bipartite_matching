#include "util.h"
using namespace std;


// ****************** Compiling options *****************
// FULL_OUTPUT, BUCKETS, FIFO/MAX_D, CURRENT, GLOBAL_HEURISTIC

#define BUCKETS
//#define FIFO

#define CURRENT


#ifdef CURRENT
    #warning "CURRENT!"
#endif

#ifdef FIFO
    #warning "FIFO!"

    #ifndef BUCKETS
        #define GLOBAL_HEURISTIC
    #endif
#else
    #warning "MAXD!"
    #define BUCKETS
#endif

#ifdef GLOBAL_HEURISTIC
    #warning "GLOBAL_HEURISTIC!"
#endif
#ifdef BUCKETS
    #include "BucketList.cpp"
    #warning "BUCKETS!"
#else
    #warning "NO BUCKETS!"
#endif


#include "graph.cpp"
#include "maxflow.cpp"

Graph G;

void Graph::fair_decomposition() {
    int n = num_vertices();
    int c1, c2;
    vector<int> dead;
    // By Mendelhson-Dulmage, we can just keep the right vertices on a maximum matching
    for (int u = 0; u < n; ++u)
        forEachEdge(u, e) if (is_forward_edge(e)) {
            capacities[edgeIndex(e)] = num_left; // for the case lambda = 1
            e->res_cap = capacity(e);
        } else e->res_cap = 0;

    vector<int> block_id(n, 0);
    vector<Rational> lb, ub, next;
    lb.push_back(0);
    ub.push_back(1);
    next.push_back(-1);

    vector<bool> complete(1);
    for (int it = 0; ; ++it) {
        printf("===== iteration %i\n", it);
        int num_blocks = lb.size();

        vector<int> left_size(num_blocks), right_size(num_blocks);
        vector<Captype> p_lambda(num_blocks);

        for (int u = 0; u < num_left; ++u)
            ++left_size[block_id[u]];
        for (int u = num_left; u < src; ++u)
            ++right_size[block_id[u]];

        vector<int> numer(num_blocks), denom(num_blocks);
        for (int i = 0; i < num_blocks; ++i) {
            int g = __gcd(left_size[i], right_size[i]);
            numer[i] = right_size[i] / g;
            denom[i] = left_size[i] / g;
        }

        for (int i = 0; i < num_blocks; ++i)  {
            if (lb[i] < ub[i]) {
                next[i] = Rational(min(numer[i], denom[i])) / denom[i];

                if (num_blocks == 1 && it == 0/* && right_size[i] >= left_size[i]*/) {
                    next[i] = 1 - Rational(1) / denom[i]; // hack for first time
                }

                cout << "i = " << i << " left = " << left_size[i] << " right = " << right_size[i] << " lb = " << lb[i] << " next = " << next[i] << " ub = " << ub[i] << endl;

                assert(next[i] >= lb[i] && next[i] <= ub[i]);
                p_lambda[i] = ceil(next[i] * denom[i]);
            }
        }


        // Erase edges between different blocks
        for (int u = 0; u < n; ++u)
            forEachEdge(u, e) if (u != src && e->v != sink && block_id[u] != block_id[e->v]) {
//                    printf("removing %i, %i\n", u, e->v);
                capacities[edgeIndex(e)] = 0;
            }

        vector<bool> just_completed(num_blocks);
        for (int i = 0; i < num_blocks; ++i)
            if (lb[i] == ub[i] && !complete[i]) {
                printf("%i just completed!\n", i);
                complete[i] = true;
                just_completed[i] = true;
            }

        // Just in case
        /*
        for (int u = 0; u < src; ++u) if (!complete[block_id[u]])
            forEachEdge(u, e)
                if (capacity(e) > 0 && e->v != sink && block_id[u] == block_id[e->v])
                    capacities[edgeIndex(e)] = 1000000000;//left_size[ block_id[u] ];
                    */

        // Restore capacities for incomplete blocks
        forEachEdge(src, e) 
            if (!complete[block_id[e->v]])
                e->res_cap = capacities[edgeIndex(e)] = p_lambda[ block_id[e->v] ];

        forEachEdge(sink, e) 
            if (!complete[block_id[e->v]])
                e->rev->res_cap = capacities[edgeIndex(e->rev)] = denom[ block_id[e->v] ];

        for (int u = 0; u < src; ++u) if (!complete[block_id[u]]) 
            forEachEdge(u, e) {
                e->res_cap = capacity(e);
                if (e->v < src && block_id[e->v] != block_id[u]) {
                    // Remove normal edge between two different blocks
                    e->res_cap = 0;
                }
            }

        // New capacities for just completed blocks
        forEachEdge(src, e) 
            if (just_completed[block_id[e->v]])
                e->res_cap = 0;

        forEachEdge(sink, e) 
            if (just_completed[block_id[e->v]] && e->rev->res_cap > 0)
                e->rev->res_cap *= -1;  // for the reverse BFS from the sink, correct distances

        printf("before:\n");
        for (int i = 0; i < num_blocks; ++i) {
            printf("  block %i %s left = %i right = %i ",
                    i, lb[i] == ub[i] ? "COMPLETE!" : "",
                    left_size[i], right_size[i]
                    );
            cout << "     lb = " << lb[i] << " next = " << next[i] << " ub = " << ub[i] << endl;
                printf("  [%lg, %lg] %lg\n", 
                        (double)lb[i].numerator() / lb[i].denominator(),
                        (double)next[i].numerator() / next[i].denominator(),
                        (double)ub[i].numerator() / ub[i].denominator()
                        );
        }

        printf("Computing max flow...\n");
        initialized = false;
        init_flow();
        dead.empty();
        double total = max_flow(&c1, &c2, dead);
        print_stats();
        printf("total flow = %lf c1 = %i c2 = %i dead = %i\n", (double)total, c1, c2, (int)dead.size());

        vector<Captype> flow(num_blocks);
        forEachEdge(sink, e) {
            int id = block_id[e->v];
            flow[id] += edge_flow(e->rev);
        }

        // Now if the cut has something on the left, dead are those with small neighborhoods (A)
        // still (moved) are the rest (bar A), with high neihborhoods
        vector<vector<int>> moved(num_blocks), moved_left(num_blocks);
        vector<bool> still(src, true);
        for (int u:dead) still[u] = false;

        for (int u = 0; u < src; ++u)
            if (still[u] && !complete[block_id[u]]) {
                moved[block_id[u]].push_back(u);
                if (u < num_left)
                    moved_left[block_id[u]].push_back(u);
//                    printf("    moving %i\n", u);
            }

        bool finished = true;
        for (int i = 0; i < num_blocks; ++i) {
            if (complete[i]) continue;

            printf("flow[%i] = %lf left = %i right = %i %lli %lli\n", i, (double)flow[i], left_size[i], right_size[i], left_size[i] * numer[i], right_size[i] * denom[i]);
            assert(flow[i] <= (long long)left_size[i] * denom[i]);
//                printf("%lg %lg %lli\n", (double)flow[i], left_size[i] * ub[i] * precision, flow[i] - floor(left_size[i] * ub[i] * precision));
            if (flow[i] == left_size[i] * numer[i]) {
                // optimal
                printf("block %i complete! lambda = %lf\n", i, (double)p_lambda[i] / denom[i]);
                lb[i] = ub[i] = next[i];
            } else {
                finished = false;
                ostringstream stream1, stream2;
                assert(!moved_left[i].empty());
//                    assert(left_size[i] > moved_left[i].size());
                if ((int)moved[i].size() == left_size[i] + right_size[i]) {
                    //only when 1?
                    // saturated left
//                        assert(false);
                    lb[i] = next[i];
                } else if (!moved[i].empty() && lb[i] < 1) {
                    for (int u:moved[i])
                        block_id[u] = lb.size();

                    lb.push_back(next[i]);
                    ub.push_back(ub[i]);
                    next.push_back(ub[i]);
                    complete.push_back(false);
                }
            }
        }
        fflush(stdout);

        if (finished) {
            // Restore edges from src
            forEachEdge(src, e) {
                capacities[edgeIndex(e)] = p_lambda[ block_id[e->v] ];
                e->res_cap = 0;
            }

            forEachEdge(sink, e) 
                if (e->rev->res_cap < 0) {
//                        printf("reversing");
                    e->rev->res_cap *= -1;
                }

            int giant = -1, no_blocks = 0;
            for (int i = 0; i < num_blocks; ++i)
                if (right_size[i] >= left_size[i]) {
                    if (giant < 0) giant = i;
                    else {
                        printf("giant %i\n", i);
                        ++no_blocks;
                        left_size[giant] += left_size[i];
                        right_size[giant] += right_size[i];

                        left_size[ i ] = n;
                        right_size[ i ] = 2 * n;
                    }
                }
            printf("giant = %i\n", giant);
            if (giant >= 0)
                for (int u = 0; u < n; ++u)
                    if (right_size[ block_id[u] ] >= left_size[ block_id[u] ])
                        block_id[u] = giant;

            vector<int> _a(num_blocks), _b(num_blocks);
            for (int i = 0; i < num_blocks; ++i) {
                if (right_size[i] > n) continue;
//                    if (right_size[i] > left_size[i]) right_size[i] = left_size[i];
//                    printf("%i %i\n", left_size[i], right_size[i]);
                int g = __gcd(left_size[i], right_size[i]);
                _a[i] = left_size[i] / g;
                _b[i] = right_size[i] / g;
            }

            // Compute a real flow
            phase = 1;
            global_update();
            recompute_heights();
            double flow = max_flow(&c1, &c2, dead);
            printf("flow = %lf expected = %lf\n", flow , total);

            printf("num_blocks = %i\n", num_blocks);
            vector<int> order, new_id(num_blocks);
            for (int i = 0; i < num_blocks; ++i)
                order.push_back(i);

            sort(order.begin(), order.end(),
                    [right_size, left_size](int a, int b) -> bool {
                        return (long long)right_size[a] * left_size[b] <
                               (long long)left_size[a] * right_size[b];
                    });

            // Remove unmatched right vertices
            int cnt = 0;
            forEachEdge(sink, e) {
                if (edge_flow(e->rev) == 0) {
                    block_id[e->v] = -1;
                    ++cnt;
                }
            }
            printf("cnt = %i\n", cnt);
            if (giant >= 0) {
                right_size[giant] -= cnt;
                _a[giant] = _b[giant] = 1;
                if (left_size[giant] != right_size[giant]) {
                    printf("%i %i %i\n", giant, left_size[giant], right_size[giant]);
                    puts("no\n");
                    exit(1);
                }
            }
            printf("=========== removed %i; %i %i\n", cnt, num_right, num_right - cnt);

            vector<vector<int>> in_block(num_blocks);
            for (int u = 0; u < src; ++u) if (block_id[u] >= 0)
                in_block[block_id[u]].push_back(u);
            vector<int> num_edges(num_blocks);
            vector<int> mult_edges(num_blocks);

#ifdef FULL_OUTPUT
            printf("# %i %i %i\n", num_left, num_right, num_blocks - no_blocks);
#endif
            for (int i: order) {
                if (right_size[i] > n) continue;

                long long min_flow = 1;
                for (int u:in_block[i]) if (u < num_left)
                    forEachEdge(u, e) if (edge_flow(e) >= min_flow)
                        ++num_edges[i];
#ifdef FULL_OUTPUT
                printf("# %i %i\n", left_size[i], right_size[i]);
                // Print matchings
                printf("#"); for (int u:in_block[i]) if (u < num_left) printf(" %i", u); printf("\n");
                printf("#"); for (int u:in_block[i]) if (u >= num_left) printf(" %i", u); printf("\n");
                printf("# %i %i\n", num_edges[i], _a[i]);
#endif
                for (int u:in_block[i]) if (u < num_left)  {
                    forEachEdge(u, e) if (e->v != sink && edge_flow(e) >= min_flow) {
                        long long times = edge_flow(e);
                        if (i == giant) {
                            times = 1;
                        } else {
                            int g = __gcd(left_size[i], right_size[i]);
//                                printf("%i %i %lf times %lli flow %lli\n", u, e->v,(double) e->res_cap, times, edge_flow(e));

                            /*
                            if (times % g!= 0) {
                                puts("warra\n");
                                exit(1);
                            }
                            */
//                                times /= g;
                        }

//                            printf("%i %i %lf %lf %lf\n", u, e->v, (double) edge_flow(e) / precision, (double)edge_flow(e), (double) min_flow);
                        assert(times >= 1);
#ifdef FULL_OUTPUT
                        printf("# %i %i %i\n", u, e->v, times);
#endif
                        capacities[edgeIndex(e)] = times;
                        mult_edges[i] += times;
            //            printf("# %i %i %lf\n", u, e->v, (double)edge_flow(e) / precision);
                    } else capacities[edgeIndex(e)] = 0;
                }
            }
            int single_edges = 0, total_edges = 0, num_matchings = 1, max_mb = 0;
            for (int j = 0, pos = 0; j < num_blocks; ++j) {
                int i = order[j];
                if (right_size[i] > n) continue;

                ++pos;

                ostringstream stream;
                stream << Rational(right_size[i]) / left_size[i];
                printf("lambda_%-6i = %-8lf\t = %12s\t", pos, (double)right_size[i] / left_size[i], stream.str().c_str());
                printf("\tleft = %5i\tright = %5i\tedges = %5i\tmult_edges = %5i\n", left_size[i], right_size[i], num_edges[i], mult_edges[i]);
                assert(mult_edges[i] = left_size[i] * _b[i]);
                single_edges += num_edges[i];
                total_edges += mult_edges[i];
                num_matchings += _a[i] - 1;
                max_mb = max(max_mb, _a[i]);
            }
            printf("single_edges = %'i\ntotal_edges = %'i\nnum_matchings = %'i max_mb = %'i\n", single_edges, total_edges, num_matchings, max_mb);
            return;
        }
    }
}

int main(int argc, char* argv[]) {
    setlocale(LC_NUMERIC, "");
    auto start = std::chrono::steady_clock::now();

    G.buildFromBipartite();
    cout << "Reading graph: " << ((std::chrono::steady_clock::now() - start)).count() * 1e-9 << " seconds" << endl; start = std::chrono::steady_clock::now();

    G.prepare();
	cout << "Preparing: " << ((std::chrono::steady_clock::now() - start)).count() * 1e-9 << " seconds" << endl; start = std::chrono::steady_clock::now();

    G.fair_decomposition();
}
