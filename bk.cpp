// g++ -O3 bk.cpp -o bk && time ./bk bibsonomy-decomp 
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <numeric>
#include <locale.h>
#include <map>
#include <iostream>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>
#include <unistd.h>
using namespace std;

#define gc(f) getc_unlocked(f)
int isSpaceChar(char c) {
            return c == ' ' || c == '\n' || c == '\r' || c == '\t'  || c == '/';
        }
inline int read_int(FILE *f, int* ret)
{
    char ch;
    int val=0;
    ch=gc(f);
    while(isSpaceChar(ch))
            ch=gc(f);
    if (feof(f)) return 0;
    while(!isSpaceChar(ch))
    {
        val=(val*10)+(ch-'0');
        ch=gc(f);
    }
    *ret = val;
    return 1;
}

typedef pair<int, int> pii;
typedef long long ll;

typedef int Captype;

struct Edge {
  int dest;
  Captype res_cap;
  Edge(int a, Captype b) : dest(a), res_cap(b) {}
};

struct Block {
    int left, right, g, n, id, num_blocks, real_right;
    vector<vector<Edge>> edges;

    const vector<int> &block_left;
    const vector<int> &block_right;

    vector<int> match_l;
    vector<int> match_r;
    vector<int> index_r;
    vector<bool> visited_r;

    Block(const vector<int>& _block_left, const vector<int>& _block_right, int _g, int _id, int _nb) : block_left(_block_left), block_right(_block_right), id(_id), num_blocks(_nb) {
        left = _block_left.size();
        right = _block_right.size();
        real_right = right;
        g = _g;
        assert(left >= right);
        n = left;
    }

    bool alt_path(int v) {
        // w is left, v is right
        visited_r[v] = true;
        for (int i = 0; i < edges[v].size(); ++i) {
            int w = edges[v][i].dest;
            if (match_l[w] < 0 || (!visited_r[match_l[w]] && alt_path(match_l[w]) )) {
                /*
                assert(v >= 0);
                assert(v < right);
                assert(w >= 0);
                assert(w < left);
                */

                match_l[w] = v;
                match_r[v] = w;
                index_r[v] = i;
                return true;
            }
        }
        return false;
    }

    vector<vector<pii>> ret;
    void output_matching(int m, int iter) {
        vector<pii> v;
        for (int r = 0; r < m; ++r)
            if (match_r[r] >= 0)
                v.push_back(pii(block_left[ match_r[r] ], block_right[ r ]));
        ret.push_back(v);
    }

    void print_all_matchings() {
        for (int iter = 0; iter < ret.size(); ++iter) 
            for (auto p:ret[iter])
                printf("# %i %i %i %i\n", id, iter, p.first, p.second);
    }

    void process() {
//        process_works(); print_all_matchings(); return;
//        {        auto e = edges; process_works(); edges.swap(e); right = real_right; ret.clear(); process_works(); print_all_matchings(); return;} 

        auto e = edges;
        if (!process_doesnt_work()) {
            printf("doing again!\n");
            edges.swap(e);
            right = real_right;
            ret.clear();
            process_works();
        } else
            printf("it worked!\n");
        print_all_matchings();
    }

    bool process_doesnt_work() {
        if (right == 1) {
            // Trivial lottery case
            match_r.assign(right, -1);
            for (int i = 0; i < right; ++i) {
                match_r[i] = i;
                output_matching(right, i);
                match_r[i] = -1;
            }
            return true;
        }

        // Add phantom edges
        int phantom_r = left - right;
        edges.resize(n);
        /*
        for (int u = 0; u < left; ++u) // left
            for (int j = 0; j < phantom_r / g; ++j) {
                int v = u % g + right + j * g;  // right
                assert(v >= 0 && v < n);
                edges[v].push_back(u);
            }
            */

        // Randomize edge order
        for (int u = 0; u < n; ++u) random_shuffle(edges[u].begin(), edges[u].end());

    	auto start = std::chrono::steady_clock::now();
//        right = n; /// LAZY HACK THAT'S GOING TO BITE ME LATER!!!!

        for (int iter = 0; iter < left / g; ++iter) {
    	    auto start2 = std::chrono::steady_clock::now();

            // Sanity check
            /*
            vector<int> deg_l(left), deg_r(right);
            for (int u = 0; u < right; ++u)
                for (int v: edges[u]) {
                    assert(0 <= v && v < left);
                    assert(0 <= u && u < right);
                    ++deg_l[v];
                    ++deg_r[u];
                }
            for (int u = 0; u < right; ++u)
                assert(deg_r[u] == left / g - iter);
            for (int u = 0; u < left; ++u) {
                if (iter == 0)
                    assert(deg_l[u] == right / g);
                assert(deg_l[u] >= right / g - iter);
            }
            */
            /*
            for (int u = 0; u < right; ++u) {
                printf("edges right %i: ", u);
                for (int w: edges[u])
                    printf("%i ", w);
                printf("\n");
            }
            */

            match_l.assign(left, -1);
            match_r.assign(right, -1);
            index_r.assign(right, -1);        // to save time this could be initialized outside

            int u = 0;
            visited_r.assign(right, false);

    	    auto start3 = std::chrono::steady_clock::now();
            vector<int> to_try, next;
            for (int u = 0; u < right; ++u) to_try.push_back(u);
            while (true) {
                bool progress = false;
                for (int u:to_try) {
                    if (match_r[u] < 0) {
                        if (!visited_r[u] && alt_path(u)) 
                            progress = true;
                        else 
                            next.push_back(u);
                    } 
                }
                if (!progress) {
                    printf("something is wrong! no progress!\n");
                    return false;
                }
                if (next.empty()) break;
                visited_r.assign(right, false);
                to_try.swap(next);
                next.clear();
            }


            /*
            printf("found matching nr %i/%i (block %i/%i: left = %i right = %i g = %i)!\n", iter + 1, left / g, id, num_blocks, left, block_right.size(), g);
            float t = ((std::chrono::steady_clock::now() - start)).count() * 1e-9;
            printf("%lf seconds per matching\n", t / (iter + 1));
            t = ((std::chrono::steady_clock::now() - start2)).count() * 1e-9;
            printf("since iteration: %lf seconds\n", t);
            t = ((std::chrono::steady_clock::now() - start3)).count() * 1e-9;
            printf("since matching: %lf seconds\n", t);
            */
//            t = ((std::chrono::steady_clock::now() - main)).count() * 1e-9; printf("main time = %lf\n", t);
               
            output_matching(n - phantom_r, iter);

            for (int r = 0; r < n - phantom_r; ++r) {
                assert(match_r[r] >= 0 && match_r[r] <= left);
                assert(match_l[ match_r[r] ] == r);
//                if (match_r[r] >= 0) printf("#%i %i %i %i\n", id, iter, block_left[ match_r[r] ], block_right[ r ]);
            }
            /*
            for (int r = 0; r < right; ++r)
                if (match_r[r] >= 0) printf("left %i => right %i\n", match_r[r], r);
            printf("\n");
            */

            for (int u = 0; u < right; ++u) {
                // Remove this edge
                int i = index_r[u];
                assert(i >= 0 && i < edges[u].size());
                assert(match_r[u] == edges[u][i].dest);
                int w = match_r[u];

                /*
                    swap(edges[u][i], edges[u][edges[u].size() - 1]);
                    assert(edges[u].back().dest == w);
                if (--edges[u].back().res_cap == 0)
                    edges[u].pop_back();
                    */
                if (--edges[u][i].res_cap == 0) {
                    swap(edges[u][i], edges[u][edges[u].size() - 1]);
                    assert(edges[u].back().dest == w);
                    edges[u].pop_back();
                }
            }
        }
        return true;
    }

    void process_works() {
        if (right == 1) {
            // Trivial lottery case
            match_r.assign(right, -1);
            for (int i = 0; i < right; ++i) {
                match_r[i] = i;
                output_matching(right, i);
                match_r[i] = -1;
            }
            return;
        }

        // Add phantom edges
        int phantom_r = left - right;
        edges.resize(n);
        for (int u = 0; u < left; ++u) // left
            for (int j = 0; j < phantom_r / g; ++j) {
                int v = u % g + right + j * g;  // right
                assert(v >= 0 && v < n);
                edges[v].push_back(Edge(u, 1));
            }

        // Randomize edge order
        for (int u = 0; u < n; ++u) random_shuffle(edges[u].begin(), edges[u].end());

    	auto start = std::chrono::steady_clock::now();
        right = n; /// LAZY HACK THAT'S GOING TO BITE ME LATER!!!!

        for (int iter = 0; iter < left / g; ++iter) {
    	    auto start2 = std::chrono::steady_clock::now();

            match_l.assign(left, -1);
            match_r.assign(right, -1);
            index_r.assign(right, -1);        // to save time this could be initialized outside

            int u = 0;
            visited_r.assign(right, false);

    	    auto start3 = std::chrono::steady_clock::now();
            vector<int> to_try, next;

            for (int u = 0; u < right; ++u) to_try.push_back(u);
            while (true) {
                bool progress = false;
                for (int u:to_try) {
                    if (match_r[u] < 0) {
                        if (!visited_r[u] && alt_path(u)) 
                            progress = true;
                        else 
                            next.push_back(u);
                    } 
                }
                if (!progress) {
                    printf("something is wrong! no progress!\n");
                    exit(1);
                }
                if (next.empty()) break;
                visited_r.assign(right, false);
                to_try.swap(next);
                next.clear();
            }

            /*
            printf("found matching nr %i/%i (block %i/%i: left = %i right = %i g = %i)!\n", iter + 1, left / g, id, num_blocks, left, block_right.size(), g);
            float t = ((std::chrono::steady_clock::now() - start)).count() * 1e-9;
            printf("%lf seconds per matching\n", t / (iter + 1));
            t = ((std::chrono::steady_clock::now() - start2)).count() * 1e-9;
            printf("since iteration: %lf seconds\n", t);
            t = ((std::chrono::steady_clock::now() - start3)).count() * 1e-9;
            printf("since matching: %lf seconds\n", t);
            */
//            t = ((std::chrono::steady_clock::now() - main)).count() * 1e-9; printf("main time = %lf\n", t);
               
            output_matching(n - phantom_r, iter);

            for (int r = 0; r < n - phantom_r; ++r) {
                assert(match_r[r] >= 0 && match_r[r] <= left);
                assert(match_l[ match_r[r] ] == r);
//                if (match_r[r] >= 0) printf("#%i %i %i %i\n", id, iter, block_left[ match_r[r] ], block_right[ r ]);
            }
            /*
            for (int r = 0; r < right; ++r)
                if (match_r[r] >= 0) printf("left %i => right %i\n", match_r[r], r);
            printf("\n");
            */

            for (int u = 0; u < right; ++u) {
                // Remove this edge
                int i = index_r[u];
                assert(i >= 0 && i < edges[u].size());
                assert(match_r[u] == edges[u][i].dest);
                int w = match_r[u];

                    swap(edges[u][i], edges[u][edges[u].size() - 1]);
                    assert(edges[u].back().dest == w);
                if (--edges[u].back().res_cap == 0)
                    edges[u].pop_back();
                /*
                if (--edges[u][i].res_cap == 0) {
                    swap(edges[u][i], edges[u][edges[u].size() - 1]);
                    assert(edges[u].back().dest == w);
                    edges[u].pop_back();
                }
                */
            }
        }
    }
};


void read_decomp(FILE* file) {
    int num_left, num_right;        // number of left and right vertices (excluding src, sink)
    int num_blocks;
    int n = 0;

    read_int(file, &num_left);
    read_int(file, &num_right);
    read_int(file, &num_blocks);
    printf("num_left = %i num_right = %i num_blocks = %i\n", num_left, num_right, num_blocks);
    printf("# %i\n", num_blocks);

    n = num_left + num_right;
    int id = 0, a, b;
    double last = -1;
    ll output_size = 0;
    while (read_int(file, &a) == 1 && read_int(file, &b) == 1) {
        unordered_map<int, int> pos_l, pos_r;
        vector<int> block_left, block_right;

        int g = __gcd(a, b);
        output_size += a / g * b;
        printf("# %i %i\n", a / g, b);

        printf("\na = %i b = %i g = %i\n", a, b, g);
//        printf("id = %-5i a = %-7i b = %-20i _a = %-7i _b = %-7i\n", id, a, b, a / g, b / g);
        double lambda = (double)b / a;
        assert(last <= lambda && last <= 1);
        lambda = last;

        int label = 0;

        int u, v, times;
        printf("reading vertices of block %i...\n", id);
        for (int i = 0; i < a; ++i) {
            int ret = read_int(file, &u);
            assert(ret > 0);
            assert(u < num_left);

            block_left.push_back(u + 1);
            pos_l[u] = i;
        }
        for (int i = 0; i < b; ++i) {
            int ret = read_int(file, &u);
            assert(ret > 0);
            assert(num_left <= u && u < n);

            block_right.push_back(u + 1 - num_left);
            pos_r[u] = i;
        }
        printf("done!\n");
        int block_edges = 0, denom;
        read_int(file, &block_edges);
        read_int(file, &denom);
        assert(block_left.size() % denom == 0);

        g = block_left.size() / denom;
        assert(block_right.size() % g == 0);

        vector<vector<Edge>> edges(b);
        long long total_edges = 0;
        for (int i = 0; i < block_edges; ++i) {
            read_int(file, &u);
            read_int(file, &v);
            read_int(file, &times);

            assert(u != v);

            total_edges += times;
//            edges[pos_r[v]].push_back( Edge(pos_l[u], times) );
            for (int i = 0; i < times; ++i)
                edges[pos_r[v]].push_back( Edge(pos_l[u], 1) );
        }
        printf("unique_edges = %i multi_edges = %lli\n", block_edges, total_edges);

        assert(a == block_left.size());
        assert(b == block_right.size());
        Block block(block_left, block_right, g, id, num_blocks);
        block.edges.swap(edges);
        block.process();

        ++id;
    }
    printf("\noutput size = %lli\n\n", output_size);
}

#include "sys/resource.h"

int main(int argc, char* argv[]) {
    setlocale(LC_NUMERIC, "");
    struct rlimit rlim = { RLIM_INFINITY, RLIM_INFINITY };

    if ( setrlimit(RLIMIT_STACK, &rlim) == -1 ) {
        perror("setrlimit error");
        exit(1);
    }

    assert(argc >= 2);
    printf("Decomposition file %s\n", argv[1]);
	auto start = std::chrono::steady_clock::now();
    read_decomp(fopen(argv[1], "r"));
	cout << "Reading: " << ((std::chrono::steady_clock::now() - start)).count() * 1e-9 << " seconds" << endl; start = std::chrono::steady_clock::now();
//ProfilerStart("fair.log");
//    G.find_matchings();
//ProfilerStop();
	cout << "Matchings: " << ((std::chrono::steady_clock::now() - start)).count() * 1e-9 << " seconds" << endl; start = std::chrono::steady_clock::now();
}
