#include<iostream>
#include<vector>
#include<memory>
#include<cmath>
#include<algorithm>
#include<queue>
#include<numeric>

// /home/zoujiu/Downloads/0907gz/god3/ disp

struct pointSingle {
    double x, y, z;
};

class KD_Tree {
public:
    struct Node {
        int left = -1, right = -1;
        int pos;
        char dim = 0;
        Node() {
            left = -1;
            right = -1;
        };
        Node(int pos_, char dim_) {
            this->left = -1;
            this->right = -1;
            this->pos = pos_;
            this->dim = dim_;
        }
    };
    std::vector<Node> nodes_;
    std::vector<int>  order_;
    const std::vector<pointSingle>& cloud_;
    int root_ = -1;
    double twice_tolerance;
    KD_Tree(const std::vector<pointSingle>& points_, double ClusterTolerance_) :
                cloud_(points_) {
        int number = (int)points_.size();
        order_.resize(number);
        std::iota(order_.begin(), order_.end(), 0);   // [0, 1, 2, ..., n - 1]

        this->twice_tolerance = ClusterTolerance_ * ClusterTolerance_;
        nodes_.reserve(number);
        root_ = this->KD_Tree_build(0, (int)order_.size(), 0);
    }
    ~KD_Tree() {
    }
    inline double distance_sq(const pointSingle& a, const pointSingle& c) {
        double dx = a.x - c.x;
        double dy = a.y - c.y;
        double dz = a.z - c.z;  
        return dx * dx + dy * dy + dz * dz;
    }
    int KD_Tree_build(int l, int r, int depth) {
        if( l >= r ) return -1;
        int dim = depth % 3;
        int mid = l + (r - l) / 2;
        std::function<bool(int, int)> compare = [&](int a, int c) {
            if( dim == 0 ) {
                return cloud_[a].x < cloud_[c].x;
            } else if( dim == 1 ) {
                return cloud_[a].y < cloud_[c].y;
            } else {
                return cloud_[a].z < cloud_[c].z;
            }
        };
        std::nth_element(order_.begin() + l, order_.begin() + mid,
                            order_.begin() + r, compare);
        Node ptr = Node(order_[mid], dim);
        int n = (int)nodes_.size();
        nodes_.push_back(ptr);
        nodes_[n].left = KD_Tree_build(l, mid, depth + 1);
        nodes_[n].right = KD_Tree_build(mid + 1, r, depth + 1);
        return n;
    }

    void extract(int nodeId, const pointSingle& query,
                 std::vector<int>& out_indices) {
        if(nodeId < 0) return;
        const Node& node = nodes_[nodeId];
        const pointSingle& p = cloud_[node.pos];
        double dist = distance_sq(query, p);
        if(dist <= this->twice_tolerance) {
            out_indices.push_back(node.pos);
        }
        int dim = node.dim;
        double diff;
        if(dim == 0) {
            diff = query.x - p.x;
        } else if(dim == 1) {
            diff = query.y - p.y;
        } else {
            diff = query.z - p.z;
        }
        // diff * diff <= ClusterTolerance * ClusterTolerance
        // diff * diff + unknow * unknow = Euclidean_distance * Euclidean_distance < ClusterTolerance * ClusterTolerance
        if(diff * diff <= this->twice_tolerance) {
            extract(nodes_[nodeId].left, query, out_indices);
            extract(nodes_[nodeId].right, query, out_indices);
        } else if(diff < 0) { // close to left, in left part 
            extract(nodes_[nodeId].left, query, out_indices);
        } else { // close to right, in right part 
            extract(nodes_[nodeId].right, query, out_indices);
        }
    }
};

class EuclideanCluster {
public:
    std::unique_ptr<KD_Tree> kdtree = nullptr;
    float ClusterTolerance_ = 0.2;
    int MinClusterSize_ = 2;
    int MaxClusterSize_ = 5000;
    EuclideanCluster(float ClusterTolerance, 
                        int MinClusterSize,
                        int MaxClusterSize) : 
                        ClusterTolerance_(ClusterTolerance),
                        MinClusterSize_(MinClusterSize),
                        MaxClusterSize_(MaxClusterSize) {}
    ~EuclideanCluster() {
    }
    std::vector<std::vector<int>> cluster(const std::vector<pointSingle>& allpoint) {
        int n = allpoint.size();
        if(n == 0) return {};
        this->kdtree = std::make_unique<KD_Tree>(allpoint, ClusterTolerance_);
        std::vector<std::vector<int>> ret;
        std::vector<char> status(n, 0);
        int size, now, index;
        std::vector<int> out_indices_all, out_indices;
        for(int i = 0; i < n; i++) {
            if(status[i] == 0) {
                out_indices_all = {i};
                status[i] = 1;
                std::queue<int> qe;
                qe.push(i);
                while(!qe.empty()) {
                    now = qe.front();
                    qe.pop();
                    out_indices.clear();
                    kdtree->extract(kdtree->root_, allpoint[now], out_indices);
                    size = (int)out_indices.size();
                    for(int ind = 0; ind < size; ind++) {
                        index = out_indices[ind];
                        if(status[index] == 0) {
                            qe.push(index);
                            status[index] = 1;
                            out_indices_all.push_back(index);
                            if((int) out_indices_all.size() >= MaxClusterSize_) break; // PCL continue not break
                        }
                    }
                    if((int) out_indices_all.size() >= MaxClusterSize_) break; // PCL continue not break
                }
                size = (int)out_indices_all.size();
                if( size >= MinClusterSize_ && size <= MaxClusterSize_) {
                    // std::sort(out_indices.begin(), out_indices.end());
                    ret.push_back(std::move(out_indices_all));
                }
            }
        }
        return ret;
    }
};

// ---------------- verification helpers (brute-force reference) ----------------
static double ref_dist2(const pointSingle& a, const pointSingle& b) {
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

// O(n) full scan: all points within tol of q (sorted)
static std::vector<int> brute_radius(const std::vector<pointSingle>& cloud,
                                     const pointSingle& q, double tol) {
    std::vector<int> out;
    double t2 = tol * tol;
    for(int i = 0; i < (int)cloud.size(); ++i)
        if(ref_dist2(cloud[i], q) <= t2) out.push_back(i);
    std::sort(out.begin(), out.end());
    return out;
}

// O(n^2) union-find: exact connected components of the radius graph,
// then filtered by [minSize, maxSize] (full-growth / PCL semantics)
static std::vector<std::vector<int>> brute_cluster(const std::vector<pointSingle>& cloud,
                                                   double tol, int minSize, int maxSize) {
    int n = (int)cloud.size();
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x) {
        while(parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    auto unite = [&](int a, int b) {
        a = find(a); b = find(b);
        if(a != b) parent[a] = b;
    };
    double t2 = tol * tol;
    for(int i = 0; i < n; ++i)
        for(int j = i + 1; j < n; ++j)
            if(ref_dist2(cloud[i], cloud[j]) <= t2) unite(i, j);

    std::vector<std::pair<int,int>> roots;
    for(int i = 0; i < n; ++i) roots.push_back({find(i), i});
    std::sort(roots.begin(), roots.end());

    std::vector<std::vector<int>> ret;
    for(int i = 0; i < n;) {
        int j = i;
        std::vector<int> comp;
        while(j < n && roots[j].first == roots[i].first) {
            comp.push_back(roots[j].second);
            ++j;
        }
        if((int)comp.size() >= minSize && (int)comp.size() <= maxSize)
            ret.push_back(comp);
        i = j;
    }
    std::sort(ret.begin(), ret.end(),
              [](const std::vector<int>& a, const std::vector<int>& b){ return a[0] < b[0]; });
    return ret;
}

static std::vector<std::vector<int>> normalize(std::vector<std::vector<int>> v) {
    for(auto& c : v) std::sort(c.begin(), c.end());
    std::sort(v.begin(), v.end(),
              [](const std::vector<int>& a, const std::vector<int>& b){ return a[0] < b[0]; });
    return v;
}

static std::vector<pointSingle> random_cloud(int n, double scale) {
    std::vector<pointSingle> c(n);
    for(auto& p : c) {
        p.x = (rand() % 1000) / 1000.0 * scale;
        p.y = (rand() % 1000) / 1000.0 * scale;
        p.z = (rand() % 1000) / 1000.0 * scale;
    }
    return c;
}

// cross-check KD radius search against brute force for `queries` random query points
static bool test_radius_search(const std::vector<pointSingle>& cloud,
                               double tol, int queries) {
    KD_Tree tree(cloud, tol);
    std::vector<int> got;
    for(int k = 0; k < queries; ++k) {
        int qi = rand() % (int)cloud.size();
        got.clear();
        tree.extract(tree.root_, cloud[qi], got);
        std::sort(got.begin(), got.end());
        std::vector<int> ref = brute_radius(cloud, cloud[qi], tol);
        if(got != ref) {
            std::cout << "  [FAIL] radius search @point " << qi
                      << ": kd=" << got.size() << " brute=" << ref.size() << "\n";
            return false;
        }
    }
    return true;
}

// cross-check clustering partition against brute-force union-find
static bool test_cluster(const std::vector<pointSingle>& cloud, double tol,
                         int minSize, int maxSize) {
    EuclideanCluster ec((float)tol, minSize, maxSize);
    std::vector<pointSingle> tmp = cloud;   // cluster() takes non-const ref
    std::vector<std::vector<int>> got = normalize(ec.cluster(tmp));
    std::vector<std::vector<int>> ref = brute_cluster(cloud, tol, minSize, maxSize);
    if(got != ref) {
        std::cout << "  [FAIL] cluster (tol=" << tol << ", min=" << minSize
                  << "): kd=" << got.size() << " clusters, brute=" << ref.size() << "\n";
        return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    std::vector<pointSingle> trps;
    for(int i = 0; i < 100; i++) {
        trps.push_back({rand()%10, rand()%10, rand()%10});
    }
    for(int i = 0; i < 200; i++) {
        trps.push_back({rand()%10 + 30, rand()%10 + 30, rand()%10 + 30});
    }
    for(int i = 0; i < 300; i++) {
        trps.push_back({rand()%10 + 60, rand()%10 + 60, rand()%10 + 60});
    }
    for(int i = 0; i < 500; i++) {
        trps.push_back({rand()%10 + 90, rand()%10 + 90, rand()%10 + 90});
    }
    for(int i = 0; i < 10; i++) {
        trps.push_back({rand()%10 + 200, rand()%10 + 200, rand()%10 + 200});
    }
    EuclideanCluster eclsk = EuclideanCluster(23, 2, 5000);
    std::vector<std::vector<int>> ret = eclsk.cluster(trps);
    std::cout << "clusters = " << ret.size() << "\n";   // 5
    for (auto& c : ret) std::cout << c.size() << " ";   // 100 200 300 500 10
    std::cout << "\n";

    // ---------------- appended verification cases ----------------
    int fail = 0;

    // [V1] radius search vs O(n) brute force on random clouds
    {
        bool ok = true;
        for(int t = 0; t < 5 && ok; ++t) {
            std::vector<pointSingle> c = random_cloud(300, 20.0);
            ok = test_radius_search(c, 3.0, 50);
        }
        std::cout << "[V1] radius search vs brute force: " << (ok ? "[OK]\n" : "[FAIL]\n");
        if(!ok) ++fail;
    }

    // [V2] cluster partition vs O(n^2) union-find on random clouds
    {
        bool ok = true;
        for(int t = 0; t < 10 && ok; ++t) {
            std::vector<pointSingle> c = random_cloud(400, 15.0);
            ok = test_cluster(c, 3.0, 1, 100000)
              && test_cluster(c, 5.0, 2, 100000)
              && test_cluster(c, 8.0, 3, 100000);
        }
        std::cout << "[V2] cluster partition vs brute force: " << (ok ? "[OK]\n" : "[FAIL]\n");
        if(!ok) ++fail;
    }

    // [V3] edge cases: empty / single / far-apart / chain / noise
    {
        bool ok = true;
        {
            std::vector<pointSingle> empty;
            EuclideanCluster ec(1.0f, 1, 100);
            if(!ec.cluster(empty).empty()) { std::cout << "  [FAIL] empty cloud\n"; ok = false; }
        }
        {
            std::vector<pointSingle> one = {{0,0,0}};
            EuclideanCluster ec(1.0f, 2, 100);
            if(ec.cluster(one).size() != 0) { std::cout << "  [FAIL] single point\n"; ok = false; }
        }
        {
            std::vector<pointSingle> two = {{0,0,0}, {100,0,0}};
            EuclideanCluster ec(1.0f, 1, 100);
            if(ec.cluster(two).size() != 2) { std::cout << "  [FAIL] far apart points\n"; ok = false; }
        }
        // chain: spacing 0.9 < tol 1.0 but total span 9.0 > tol -> must be ONE cluster
        {
            std::vector<pointSingle> chain;
            for(int i = 0; i < 11; ++i) chain.push_back({i * 0.9, 0, 0});
            EuclideanCluster ec(1.0f, 1, 100);
            std::vector<std::vector<int>> r = ec.cluster(chain);
            if(r.size() != 1 || r[0].size() != 11) {
                std::cout << "  [FAIL] chain growth: got " << r.size() << " clusters\n";
                ok = false;
            }
        }
        // 3 tight points + 1 isolated, minSize=2 -> 1 cluster of 3
        {
            std::vector<pointSingle> mix = {{0,0,0},{0.5,0,0},{0,0.5,0}, {50,50,50}};
            EuclideanCluster ec(1.0f, 2, 100);
            std::vector<std::vector<int>> r = ec.cluster(mix);
            if(r.size() != 1 || r[0].size() != 3) {
                std::cout << "  [FAIL] noise filter: got " << r.size() << " clusters\n";
                ok = false;
            }
        }
        std::cout << "[V3] edge cases (empty/single/far/chain/noise): "
                  << (ok ? "[OK]\n" : "[FAIL]\n");
        if(!ok) ++fail;
    }

    // [V4] MaxClusterSize cap demo (truncation splits the component)
    {
        std::vector<pointSingle> blob;
        for(int i = 0; i < 100; ++i)
            blob.push_back({(double)(rand()%5), (double)(rand()%5), (double)(rand()%5)});
        EuclideanCluster ec(20.0f, 1, 50);
        std::vector<std::vector<int>> r = ec.cluster(blob);
        std::cout << "[V4] cap demo: 100 connected points, MaxSize=50 -> "
                  << r.size() << " accepted clusters, sizes:";
        for(auto& c : r) std::cout << " " << c.size();
        std::cout << "\n     (truncation splits the component; large MaxSize = full-growth)\n";
    }

    // [V5] chain / path-graph cases (transitive region growing)
    {
        bool ok = true;
        auto expect = [&](const char* name, const std::vector<std::vector<int>>& r,
                          int nclusters, const std::vector<int>& sizes) {
            std::vector<int> got;
            for(auto& c : r) got.push_back((int)c.size());
            std::sort(got.begin(), got.end());
            std::vector<int> exp = sizes;
            std::sort(exp.begin(), exp.end());
            if((int)r.size() != nclusters || got != exp) {
                std::cout << "  [FAIL] " << name << ": got " << r.size()
                          << " clusters, sizes:";
                for(int s : got) std::cout << " " << s;
                std::cout << "\n";
                ok = false;
            }
        };

        // C1: spacing exactly == tol -> connected (<= is inclusive) -> 1 cluster
        {
            std::vector<pointSingle> ch;
            for(int i = 0; i < 10; ++i) ch.push_back({i * 1.0, 0, 0});
            EuclideanCluster ec(1.0f, 1, 100000);
            expect("C1 exact-boundary chain (d == tol)", ec.cluster(ch), 1, {10});
        }
        // C2: spacing just over tol -> no edges -> n singletons; minSize filters all
        {
            std::vector<pointSingle> ch;
            for(int i = 0; i < 10; ++i) ch.push_back({i * 1.01, 0, 0});
            EuclideanCluster ec(1.0f, 1, 100000);
            expect("C2 broken chain, minSize=1", ec.cluster(ch), 10,
                   std::vector<int>(10, 1));
            EuclideanCluster ec2(1.0f, 2, 100000);
            expect("C2 broken chain, minSize=2 (all filtered)", ec2.cluster(ch), 0, {});
        }
        // C3: 3D diagonal chain, per-step Euclidean distance 0.9 < tol -> 1 cluster
        {
            std::vector<pointSingle> ch;
            double s = 0.9 / std::sqrt(3.0);   // step per axis; 3*s^2 = 0.9^2
            for(int i = 0; i < 20; ++i) ch.push_back({i * s, i * s, i * s});
            EuclideanCluster ec(1.0f, 1, 100000);
            expect("C3 3D diagonal chain", ec.cluster(ch), 1, {20});
        }
        // C4: zigzag staircase (direction alternates x/y) -> still 1 cluster
        {
            std::vector<pointSingle> ch;
            for(int i = 0; i < 12; ++i) {
                double x = (i / 2) * 0.9 + (i % 2 ? 0.9 : 0.0);
                double y = (i / 2) * 0.9;
                ch.push_back({x, y, 0});
            }
            EuclideanCluster ec(1.0f, 1, 100000);
            expect("C4 zigzag staircase chain", ec.cluster(ch), 1, {12});
        }
        // C5: two parallel chains far apart -> 2 clusters
        {
            std::vector<pointSingle> ch;
            for(int i = 0; i < 10; ++i) ch.push_back({i * 0.9, 0, 0});
            for(int i = 0; i < 10; ++i) ch.push_back({i * 0.9, 100, 0});
            EuclideanCluster ec(1.0f, 1, 100000);
            expect("C5 two separated chains", ec.cluster(ch), 2, {10, 10});
        }
        // C6: one chain with a single big gap in the middle -> 2 clusters
        {
            std::vector<pointSingle> ch;
            for(int i = 0; i < 5; ++i) ch.push_back({i * 0.9, 0, 0});
            for(int i = 0; i < 6; ++i) ch.push_back({100 + i * 0.9, 0, 0});
            EuclideanCluster ec(1.0f, 1, 100000);
            expect("C6 chain with mid gap", ec.cluster(ch), 2, {5, 6});
        }
        // C7: closed loop (square ring) -> 1 cluster, and growth must TERMINATE
        //     (cycle in the radius graph; visited-marking prevents infinite loop)
        {
            std::vector<pointSingle> ch = {
                {0,    0,    0}, {0.45, 0,    0}, {0.9,  0,    0},
                {0.9,  0.45, 0}, {0.9,  0.9,  0}, {0.45, 0.9,  0},
                {0,    0.9,  0}, {0,    0.45, 0}
            };
            EuclideanCluster ec(1.0f, 1, 100000);
            expect("C7 closed loop (cycle, must terminate)", ec.cluster(ch), 1, {8});
        }
        // C8: long chain (500 points, spacing 0.5) -> 1 cluster of 500
        {
            std::vector<pointSingle> ch;
            for(int i = 0; i < 500; ++i) ch.push_back({i * 0.5, 0, 0});
            EuclideanCluster ec(1.0f, 1, 100000);
            expect("C8 long chain (500 pts)", ec.cluster(ch), 1, {500});
        }
        // C9: chain along z axis (exercises KD-Tree axis rotation x->y->z)
        {
            std::vector<pointSingle> ch;
            for(int i = 0; i < 15; ++i) ch.push_back({0, 0, i * 0.9});
            EuclideanCluster ec(1.0f, 1, 100000);
            expect("C9 z-axis chain", ec.cluster(ch), 1, {15});
        }
        // C10: chain segments with minSize filter: sizes 2 / 6 / 1, minSize=3
        //      -> only the 6-point segment survives
        {
            std::vector<pointSingle> ch;
            for(int i = 0; i < 2; ++i) ch.push_back({i * 0.5, 0, 0});
            for(int i = 0; i < 6; ++i) ch.push_back({100 + i * 0.5, 0, 0});
            ch.push_back({200, 0, 0});
            EuclideanCluster ec(1.0f, 3, 100000);
            expect("C10 segmented chain + minSize filter", ec.cluster(ch), 1, {6});
        }
        // C11: random-walk chain cross-checked against brute-force union-find
        {
            bool bok = true;
            for(int t = 0; t < 5 && bok; ++t) {
                std::vector<pointSingle> w;
                pointSingle p{0, 0, 0};
                for(int i = 0; i < 300; ++i) {
                    w.push_back(p);
                    double a = (rand() % 628) / 100.0;       // polar
                    double b = (rand() % 628) / 100.0;       // azimuth
                    p.x += std::sin(a) * std::cos(b) * 0.5;  // fixed step 0.5
                    p.y += std::sin(a) * std::sin(b) * 0.5;
                    p.z += std::cos(a) * 0.5;
                }
                bok = test_cluster(w, 1.0, 1, 100000)
                   && test_cluster(w, 1.0, 5, 100000);
            }
            if(!bok) { std::cout << "  [FAIL] C11 random-walk chain vs brute force\n"; ok = false; }
        }
        std::cout << "[V5] chain cases (boundary/broken/diagonal/zigzag/two/gap/"
                     "loop/long/z-axis/minSize/random-walk): "
                  << (ok ? "[OK]\n" : "[FAIL]\n");
        if(!ok) ++fail;
    }

    std::cout << "\n" << (fail == 0 ? "ALL VERIFICATIONS PASSED" : "SOME VERIFICATIONS FAILED") << "\n";
    return fail ? EXIT_FAILURE : EXIT_SUCCESS;
}