// Bucket lists of active nodes for use in the push-relabel maximum flow algorithm.
struct BucketList {
    int n,              // total number of nodes (excluding sentinel)
        infinity,       // height at which nodes are ignored (except src and sink)
        num_included,   // number of non-ignored nodes
        aMax,           // maximum height of active nodes
        dMax ;          // maximum height of included nodes (with height < infinity);
    int* level;
    long long* excess;

    // Each bucket is maintained as two doubly linked list containing the active nodes and the inactive nodes at that level.
    // The sink is always active at level 0.
    typedef struct _Node {
        struct _Node* next;
        struct _Node* prev;
    } Node;

    struct Bucket {
        Node* firstActive;
        Node* firstInactive;
    };

    vector<Bucket> buckets;     // list of buckets
    vector<Node> nodes;         // all nodes, active or not, including src, sink and a sentinel node
    Node* sentinelNode;         // sentinel node to mark the end of a list

    BucketList() : n(-1) {}

    int node_id(const Node* node) const { return node - &nodes[0]; }

    // Initialize the bucket list:
    //   -n: total number of nodes (including src = n - 2 and sink = n - 1)
    //   infinity_level: nodes with label >= this will be ignored
    //   labels: array with node labels (heights). Must be >= 1
    //   _excess: array with excess incoming flow.  This is a quick ugly hack, it is only used to
    //            identify active nodes as those with positive excess
    void assign(int _n, int infinity_level, int* labels, long long* _excess) {
        n = _n;
        infinity = infinity_level;
        aMax = -1;
        dMax = -1;
        level = labels;
        excess =_excess;
        buckets.resize(infinity);
        nodes.resize(n + 1);
        sentinelNode = &nodes[n];
        sentinelNode->next = sentinelNode->prev = sentinelNode;
        for (int d = 0; d < infinity; ++d)
            buckets[d].firstActive = buckets[d].firstInactive = sentinelNode;

        iAdd(0, n - 1);     // the sink is always active at level 0
        num_included = 2;
        for (int u = 0; u < n - 2; ++u)
            if (level[u] < infinity) {
                ++num_included;
                if (level[u] > dMax) dMax = level[u];

                if (excess[u] > 0)
                    aAdd(level[u], u);
                else
                    iAdd(level[u], u);
            }
    }

    // Check internal data consistency
    void check() const {
        for (int d = 0; d <= dMax; ++d){
            if(!(buckets[d].firstActive != sentinelNode || buckets[d].firstInactive != sentinelNode)) printf("  violation %i\n", d);;
            assert(buckets[d].firstActive != sentinelNode || buckets[d].firstInactive != sentinelNode);
        }

        for (const Bucket* l = &buckets[0]; l != &buckets[dMax + 1]; ++l) {
            int d = l - &buckets[0];
            assert(l->firstActive != sentinelNode || l->firstInactive != sentinelNode);
            printf("checking d = %i active\n", d);
            for (Node* i = l->firstActive; i != sentinelNode; i = i->next) {
         //       printf("%p %i\n", i, i - &nodes[0]);
                assert(level[node_id(i)] == d);
            }
            printf("checking d = %i inactive\n", d);
            for (Node* i = l->firstInactive; i != sentinelNode; i = i->next) {
                assert(level[node_id(i)] == d);
            }
        }
    }

    // Print the bucket list
    void print() const {
        printf("  n = %i infinity = %i aMax = %i dMax = %i\n", n, infinity, aMax, dMax);
        printf("    active:\n");
        for (const Bucket* l = &buckets[0]; l != &buckets[dMax + 1]; ++l) if (l->firstActive != sentinelNode) {
            printf("      d = %li:", l - &buckets[0]);
            int r = 0;
            for (Node* i = l->firstActive; i != sentinelNode && r < 10; i = i->next, ++r)
                printf(" %i", node_id(i));
            r = 0;
            for (Node* i = l->firstActive; i != sentinelNode; i = i->next, ++r);
            printf(" (total %i)\n", r);
        }
        printf("\n");

        printf("    inactive:\n");
        for (const Bucket* l = &buckets[0]; l != &buckets[dMax + 1]; ++l) if (l->firstInactive != sentinelNode) {
            printf("      d = %li:", l - &buckets[0]);
            int r = 0;
            for (Node* i = l->firstInactive; i != sentinelNode && r < 10; i = i->next, ++r)
                printf(" %i", node_id(i));
            r = 0;
            for (Node* i = l->firstInactive; i != sentinelNode; i = i->next, ++r);
            printf(" (total %i)\n", r);
        }
        printf("\n");
    }

    /*
       active means excess > 0 and label < infinity
       inactive means excess == 0 and label < infinity
       before discharging, we need to remove from active because the level may change
    */
    void add_to_front(Node* &first, Node* i) {
        i->next = first; i->prev = sentinelNode;
        first->prev = i;
        first = i;
    }

    void remove(Node* &first, Node* i) {
        Node* next = i->next;
        if (first == i) {
            i->next->prev = sentinelNode;
            first = next;
        } else {
            Node* prev = i->prev;
            prev->next = i->next;
            next->prev = prev;
        }
    }

    void aAdd(int d, int u) {
        assert(d < infinity);
//        printf("aAdd %i %i\n", d, u);
        add_to_front(buckets[d].firstActive, &nodes[u]);
        if (d > aMax) aMax = d;
        if (dMax < aMax) dMax = aMax;
    }

    void iAdd(int d, int u) {
        assert(d < infinity);
//        printf("iAdd %i %i\n", d, u);
        add_to_front(buckets[d].firstInactive, &nodes[u]);
    }

    void aRemove(int d, int u) {
        assert(d < infinity);
//        printf("aRemove %i %i\n", d, u);
        remove(buckets[d].firstActive, &nodes[u]);
    }

    void iDelete(int d, int u) {
        assert(d < infinity);
//        printf("iDelete %i %i\n", d, u);
        remove(buckets[d].firstInactive, &nodes[u]);
    }

    void activate(int u) {
//        printf("activating %u\n", u);
        iDelete(level[u], u);
        aAdd(level[u], u);
    }

    int max_active() const {
        assert(buckets[aMax].firstActive != sentinelNode);
        return node_id( buckets[aMax].firstActive );
    }

    void deactivate(int u) {
 //       printf("deactivating %u\n", u);
        aRemove(level[u], u);
        while (aMax >= 0 && buckets[aMax].firstActive == sentinelNode) --aMax;
        iAdd(level[u], u);
    }

    bool some_active() const { return aMax > 0; }

    // u must not be in inactive list (temporarily)
    int relabel(int u, int d) {
        assert(level[u] < infinity);
        assert(excess[u] > 0);

        int gap = level[u];
        iDelete(gap, u);
        level[u] = d;

        int cnt = 0;
        if (buckets[gap].firstActive == sentinelNode && buckets[gap].firstInactive == sentinelNode) {
//            printf("gap at height %i\n ==== relabeling %i from %i to %i; infinity = %i\n", gap, u, gap, d, infinity);
            level[u] = infinity;
            cnt = 1;

            // Gap heuristic
            for (d = gap + 1; d <= dMax; ++d) {
                Bucket* l = &buckets[d];
                // This does nothing for highest level selection
                for (Node* i = l->firstActive; i != sentinelNode; i = i->next) {
                    level[node_id(i)] = infinity;
                    ++cnt;
                }

                for (Node* i = l->firstInactive; i != sentinelNode; i = i->next) {
                    level[node_id(i)] = infinity;
                    ++cnt;
                }
                buckets[d].firstActive = buckets[d].firstInactive = sentinelNode;
            }

            dMax = gap - 1;
            if (aMax >= dMax) {
                aMax = dMax;
                while (aMax >= 0 && buckets[aMax].firstActive == sentinelNode) --aMax;
            }
//            printf("removed %i nodes\n", cnt);
            assert(cnt > 0);
            num_included -= cnt;
//            printf("new label = %i\n", level[u]);
        }

        // put it back to inactive list if no gap
        if (level[u] < infinity) {
            iAdd(d, u);
            if (dMax < d) dMax = d;
        }

        return cnt;
    }
};
