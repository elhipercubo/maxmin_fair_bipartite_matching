typedef long long Captype;

// 24 bytes (int=4, pointer=8, Captype=8)
struct Edge {
    int v;              // edge from v to e->rev->v
    Edge* rev;          // pointer to edge from e->rev->v to v
    Captype res_cap;    // residual capacity ( = capacity[e] - current flow ) for forward edges; 0 for backward edges
};

typedef pair<int, int> pii;
class Graph {
  public:
    vector<pii> _edges;             // original edges, only used before prepare()
    
    vector<Edge> edges;             // edges with residual capacities
    vector<Captype> capacities;     // original capacities

    vector<Edge*> _first;           // pointer to first edge of each vertex
    vector<int> label;              // label of each vertex (for push-relabel)
    vector<Captype> excess;         // excess incoming flow of each vertex

    bool initialized;               // has there been a call to init_flow?

    chrono::steady_clock::time_point start;     // for timing statistics

    int src, sink;                  // we NEED src, sink > all other vertices
    int num_left, num_right;        // number of left and right vertices (excluding src, sink)

    // Push-relabel statistics
    int pushes, sat_pushes, relabels, activations, num_gaps, gap_node_cnt, num_updates, skipped;

    long long edges_looked, work, last_relabels, work_since_gap;

    // Phase of push-relabel algorithm (0 or 1).
    // To find min-cuts, phase 0 suffices; to find flows, phase 0 needs to be followed by phase 1.
    int phase = 0;

#ifdef CURRENT
    // Next edge to look at for each vertex
    vector<Edge*> current;
#endif

    // Queue of active nodes. Only used if  FIFO is defined.
    deque<int> active;
    bool using_fifo = false;

    // Bucket lists containing active nodes at each possible height (label)
    BucketList buckets;

    Graph() { pushes = relabels = activations = num_gaps = gap_node_cnt = num_updates = skipped = last_relabels = 0; initialized = false; }

    void print_edge(const Edge* e) const;

    // Read a bipartite graph from stdin (where ids of left and right vertices collide) and do the following:
    //   - set num_left, num_right
    //   - add src, sink and edge pairs from src to L and from R to sink
    //   - relabel the right vertices so that their ids come after those of the left vertices
    void buildFromBipartite();

    // Take _edges (containing the list of forward edges as pairs of integers, including those from src and sink) and does the following:
    //   - compute e->rev for each edge
    //   - sort the array of edges so that edges from a single vertex are consecutive
    //   - set capacities of all edges to 1
    //   - compute _first (first edge from each vertex)
    //   - allocate space for 'label'
    //   - set timer ('start')
    //   - free space for _edges (not needed anymore)
    // All other functions assume prepare() has been called after reading the input.
    void prepare();

    // Total number of vertices (including src and sink)
    int num_vertices() const { return _first.size() - 1; }

    // Total number of forward edges in the graph (that is, not including the backward edges using for max flow computations)
    int num_edges() const { return edges.size() / 2; }

    // Degree of a vertex
    int deg(int v) const { return _first[v + 1] - _first[v]; }

    // Transform between edge pointers and indices in the 'edges' array 
    int edgeIndex(const Edge* e) const { return e - &edges[0]; }

    // Original capacity of an edge
    Captype capacity(const Edge* e) const { return capacities[edgeIndex(e)]; }

    // Is e a forward edge?
    bool is_forward_edge(const Edge* e) const { return capacity(e) > 0; }

    // Current pre-flow being routed along this edge. Returns 0 for backward edges
    Captype edge_flow(const Edge* e) const { return is_forward_edge(e) ? capacity(e) - e->res_cap : 0; }

    // Can I push flow along edge e = (u, e->v) where the label decreases by exactly one?
    bool is_admissible(int u, const Edge* e) const { return e->res_cap > 0 && label[e->v] == label[u] - 1; }

    // Initialize bucket list with all nodes of heights < n. Activate the sink and all nodes of positive excess.
    // Also resets current[] to the first edge from each vertex.
    void recompute_heights();

    // Are there active nodes?
    bool some_active() {
        if (using_fifo)
            return !active.empty();
        else
            return buckets.some_active();
    }

    // Select the next active node to process
    int select_active() {
        if (using_fifo)
            return active.front();
        else
            return buckets.max_active();
    }

    // Activate a node != src, sink
    void activate(int v) {
        assert(excess[v] > 0);
        assert(label[v] < num_vertices());

        ++activations;
#ifdef BUCKETS
        buckets.activate(v);
#endif
        if (using_fifo) active.push_back(v);
    }

    // Remove a node from the active FIFO queue (but not from bucket lists). Only relevant if using_fifo = TRUE.
    void remove_from_queue(int);

    // Initializes a max flow computation from src to sink. We assume there are no vertices with degree 0.
    //   - reserve space for label, excess
    //   - call recompute_heights
    //   - greedily push flow from src to each forward edge from it
    //   - set label[src] = n, label[sink] = 0, label[left vertices] = 2 and label[right vertices] = 1
    //   - set initialized = true and phase = 0
    //   - set using_fifo 
    void init_flow(bool use_fifo = false);

    // Push a positive amount of flow along an admissible edge e = u->v (v != src) or an edge src->v
    // Active v if v != sink was inactive. (Does not perform any deactivation.)
    bool push(int u, Edge* e, Captype amount);

    // Relabels u != src, sink with 0 < label[u] < n to d. If BUCKETS is defined, use the gap heuristic.
    void relabel(int u, int d);

    // Discharge a node with positive excess by pushing flow to its neighbours through admissible edges and relabelling if necessary, 
    // until the excess becomes 0 or its label becomes n
    // u must have been deactivated from the bucket list
    void discharge(int u);

    // Computes distances from src (if phase = 0) or sink (if phase = 1)
    // in the residual graph. Updates distances (labels) accordingly.
    void global_update();

    Captype max_flow(int* c1, int* c2, vector<int>& dead);

    // Accumulated max flow statistics
    void print_stats();

    // Find the fair decomposition of the graph
    void fair_decomposition();
};
