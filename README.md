# KD_Tree
kd-tree in c++, some code from doubao

# build
`g++ -O3 -Wall -Wextra KDTree_cluster.cpp -o kdt`

# `radius_search` Line-by-Line Breakdown: Why Is It Written This Way?

The core idea of this code is: **leverage the spatial partitioning of the KD-Tree to only search subtrees that the "search sphere" could possibly reach, and directly prune away entire regions that cannot contain any results.**

---

## 1. Geometric Background: Every Node Represents a "Half-Space"

When building a KD-Tree, each node splits the space into two halves along axis `node.axis` using a **splitting hyperplane**:

- Left subtree: all points whose coordinate on that axis is **≤** the node's coordinate
- Right subtree: all points whose coordinate on that axis is **>** the node's coordinate

```
            axis = x (vertical splitting line)
                |
   Left subtree |     Right subtree
   (far/near)   |     (near/far)
                |
          ● node.point
```

---

## 2. `diff`: The Perpendicular Distance from the Query Point to the Splitting Hyperplane

```python
diff = q[node.axis] - node.point[node.axis]
```

diff 就是查询点 q 到这个切分平面的 有符号距离 （只算切分轴方向上的差）：

- diff <= 0 → q 在平面 左侧 → 左子树是 near （q 所在的一侧）
- diff > 0 → q 在平面 右侧 → 右子树是 near

```
欧氏距离² = (Δaxis)² + (其他轴的差)² ≥ (Δaxis)² ≥ diff²
```
也就是说：far 子树中任何点到 q 的距离，都不可能小于 q 到切分平面的垂直距离 |diff| 。
  
`diff` is the **signed distance** from query point `q` to this splitting plane (counting only the difference along the splitting axis):

- `diff <= 0` → q is on the **left** of the plane → the left subtree is `near` (the side containing q)
- `diff > 0`  → q is on the **right** of the plane → the right subtree is `near`

```python
near = node.left if diff <= 0 else node.right
far  = node.right if diff <= 0 else node.left
```

---

## 3. The Key Pruning Test: `if diff * diff <= r2` — This Is the Heart of "Why It's Written This Way"

```python
search(near)                          # the near side must always be searched
if diff * diff <= r2:                 # the far side is searched only if this holds
    search(far)
```

**Why is `near` searched unconditionally?**
q itself lies inside the region governed by the near subtree, and q is inside the search sphere (it's the sphere's center). So the sphere **necessarily intersects** the near region — it may contain results and cannot be skipped.

**Why can the `far` side be conditionally skipped?**

Consider **any point p** in the far subtree: it and q lie on opposite sides of the splitting plane, so their coordinate difference along the splitting axis is at least `|diff|`:

```
Euclidean distance² = (Δaxis)² + (differences on other axes)² ≥ (Δaxis)² ≥ diff²
```

In other words:

> **No point in the far subtree can be closer to q than the perpendicular distance `|diff|` from q to the splitting plane.**

Therefore:

| Condition | Geometric meaning | Action |
|---|---|---|
| `diff² > r²` | The search sphere of radius r **cannot touch** the splitting plane at all; the entire far side lies outside the sphere | **Prune** — skip the whole far subtree |
| `diff² ≤ r²` | The search sphere **crosses** the splitting plane; part of the far side may lie inside the sphere | Must recursively search far |

2D illustration:

```
  Case diff > r:                 Case diff ≤ r:
       |  far                        |  far
   q   |     ● far point       q     |
  ( ●)-|----                  ( ●)--|----  the sphere crosses the plane!
   ↑   |   ↑                     \  |  /
  sphere|  distance > r,        sphere|  ● this far point may be inside
  can't |  cannot be in the        /  |
  reach |  sphere                     |
       |                             |
```

This is the essence of KD-Tree acceleration: **a single comparison eliminates an entire subtree (potentially thousands or tens of thousands of points)** instead of computing distances one by one.

---

## 4. Testing the Node Itself

```python
d2 = float(np.sum((q - node.point) ** 2))
if d2 <= r2:
    found.append((node.index, math.sqrt(d2)))
```

- Whenever a node is visited, check whether the point stored at that node lies inside the sphere (using the **full Euclidean distance** across all dimensions).
- The comparison uses `d2 <= r2` (squared distances) to avoid taking a square root for every point; `math.sqrt` is called only once when a point actually hits — a common performance optimization.
- Note: the pruning condition `diff² <= r²` is only a loose lower bound meaning "results **might** exist." Whether a point truly counts as a hit is decided here by the full-distance test.

---

## 5. Comparison with KNN Pruning

The same geometric idea applies; the pruning condition in the KNN version is:

```python
if len(heap) < k or diff * diff < -heap[0][0]:
    search(far)
```

The difference is the search radius:

- **Radius search**: the radius r is **fixed**, so the pruning threshold is always `r²`
- **KNN**: the radius is **dynamic** — the distance to the current k-th nearest neighbor (the heap top, `-heap[0][0]`). As the search proceeds and closer points are found, this "sphere" keeps shrinking, making pruning increasingly aggressive

---

## 6. A Small Detail: The `diff == 0` Boundary Case

`diff <= 0` assigns the case where q lies exactly on the splitting plane to the left subtree (near). This is fine — the splitting point itself is stored at the current node (and already tested separately), while for any other point on the plane, regardless of which subtree it falls into, `diff² = 0 ≤ r²` always holds, so the far side is guaranteed to be searched as well. **No points are missed.**

---

### One-Sentence Summary

> **`near` must always be searched (the sphere's center is on that side); the `far` side is first tested with the lower bound given by the "perpendicular distance to the splitting plane" — if even that minimum possible distance exceeds the radius, the entire far subtree cannot contain an answer and is pruned away; otherwise the sphere crosses the plane and we recurse into it.** This is how the O(n) brute-force search drops to an average of O(log n).
