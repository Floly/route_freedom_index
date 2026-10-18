# What the columns in `validation.ipynb` mean

Notation: s, t are the origin and destination nodes of a pair; d* is the shortest network path s→t;
euclid is the straight-line s–t distance; ε is the allowed detour (0.1); B = (1+ε)·d* is the path
budget; U is the total length of usable edges; A is the total length of edges whose both endpoints
lie inside the ellipse.

## Per-graph summary table (A5, first table)

| column | meaning |
|---|---|
| `RFI` | the final number for the graph: mean of `median_r` over bins with `used=True`. Higher means more freedom to vary the route. NaN if no bin reaches `MIN_PAIRS_PER_BIN` |
| `bins_used` | how many of the 3 bins enter the mean |
| `circuity` | median d*/euclid over all kept pairs. 1.0 = straight path, ~1.3 = typical grid |
| `out_of_range` | pairs whose euclid is < 300 m or > 3000 m. They belong to no bin and do not take part in the metric |
| `r_gt_1` | number of pairs with r > 1 + 1e-6. Must be 0; anything else is a bug or a snapping artefact |

## Per-bin table (A5, second table)

A bin is a stratum of pairs by Euclidean distance, `[lo, hi)` in meters. The last bin also includes 3000.

| column | meaning |
|---|---|
| `graph` | name of the test graph |
| `bin_m` | Euclidean s–t distance range in meters |
| `n_pairs` | pairs left in the bin after exclusions and subsampling (at most `MAX_PAIRS_PER_BIN` = 300) |
| `used` | `n_pairs ≥ MIN_PAIRS_PER_BIN` (50). If False, the bin is left out of RFI |
| `median_r` | median of r = U/A over the bin's pairs: the share of streets you can actually use among the streets within geometric reach under the budget |
| `iqr` | interquartile range of r (q75 − q25). Spread between pairs: a large value means the result depends on which pairs were sampled |
| `circuity` | median d*/euclid over the bin's pairs |
| `candidates` | all pairs i<j among the sampled origins whose euclid falls in the bin, before any exclusion |
| `unreachable` | candidates whose s and t are in different connected components (d* = ∞). Checked first |
| `too_circuitous` | reachable candidates with d* > 1.5·euclid (`CIRCUITY_CAP`). Excluded because the Dijkstra cutoff does not guarantee exact distances for them |

Bookkeeping: `candidates = n_pairs + unreachable + too_circuitous + boundary + (pairs removed by
subsampling to 300)`. The `boundary` column is not shown because on synthetic graphs it is always 0
(there is no download polygon until M2).

Why `n_pairs` is small for the tree: the tree path is almost always longer than 1.5·euclid, so
`too_circuitous` removes most candidates (only 6 pairs remain in the 1500–3000 bin).

## Pairs table (`result["pairs"]`)

| column | meaning |
|---|---|
| `s`, `t` | node indices in the graph |
| `row_s`, `row_t` | rows of these nodes in the Dijkstra matrix (indices into the origins list) |
| `bin` | bin number (0, 1, 2) |
| `euclid` | straight-line s–t distance, m |
| `dstar` | shortest network path, m |
| `U` | total length of usable edges: edges (u,v) with `min(d_s(u)+w+d_t(v), d_s(v)+w+d_t(u)) ≤ B` |
| `A` | total length of edges whose both endpoints lie inside the ellipse with foci s, t and focal sum B |
| `r` | U / A. Always ≤ 1 because U ⊆ A. Never clipped |

## `edge_freedom`

Per edge: the number of pairs for which the edge is usable, divided by the number of pairs for which
it lies inside the ellipse. NaN if no pair covers the edge. A high value means the edge almost always
serves detour routes; a low value (cul-de-sacs) means the edge is within reach but cannot be walked
within the budget.

## Reading them together

- High `RFI` with small `iqr`: route freedom is consistently high.
- Low `n_pairs` or `used=False` in a long bin: too few pairs, do not trust that bin (typical for trees and heavily fragmented graphs).
- Large `unreachable`: the graph is disconnected and part of the network does not take part in the metric.
- Removing edges raises `circuity` and lowers `RFI`: the two signals agree.
