# CaV-EMPS methodology

**CaV-EMPS** means *Cavity detection via Electrostatic Map Pocket Scoring*.
It proposes docking boxes from receptor geometry and receptor-derived AutoGrid
signals when a co-crystal ligand center is unavailable. The implementation is
`autogenerate_centers_tsv` in `docking/make_grids.py`, orchestrated by
`docking/dock_v02.py`. A withheld crystal ligand can be used afterward as
a benchmark reference; it is not a prediction input. The site score orders
search hypotheses and is not a ligand-binding affinity.

## Evidence and candidate families

The geometric path rasterizes receptor atoms onto the AutoGrid lattice.
An Euclidean distance transform (EDT) describes open volume and identifies
internal cavities and channels. Surface-cleft candidates cover another
receptor geometry. The map path uses AutoGrid carbon, electrostatic and
desolvation signals: favorable energies are clipped to physically useful
negative values, smoothed over a ligand-sized region and constrained to
receptor-proximal pocket shells. This avoids treating one isolated favorable
voxel as a complete binding pocket.

The default `receptor_search` policy assembles internal, surface and
consensus candidates into a compact portfolio. Å-scale non-maximum
suppression removes near-duplicates, and the retained centers and grid
dimensions define downstream docking boxes. The cavity-radius threshold can
be estimated from each receptor's EDT peaks instead of fixed globally.
Inter-site separation is clamped relative to box side length to keep boxes
distinct without becoming needlessly sparse.

## Search policies and controls

| Generated `SITE_POLICY` | Candidate selection |
| --- | --- |
| `receptor_search` | Compact complementary receptor-derived portfolio; default |
| `exhaustive_search` | Broader search across candidate families |
| `internal` | Internal cavities |
| `surface` | Surface clefts (`maps` is an internal alias) |
| `hybrid` | Consensus family |

`AUTOSITES` is a requested site budget; geometry and filtering can yield
fewer actual proposals. The `centers.tsv` metadata records policy and
parameters for reuse checks. Preserve it when comparing runs. Site IDs
`S1`, `S2`, etc. are output identifiers, not scores, priorities or
quality labels. [Configuration](../reference/configuration.md) lists the
adaptive EDT, surface-shell, map and separation settings. Edit non-CLI dials
only with a declared evaluation plan and preserve the actual generated config.

## Validation boundary

Map resolution, receptor extent, protonation and box settings affect
proposals. A site omitted from the proposal set cannot be recovered by local
docking. Compare geometry-only, map-based and combined settings under a fixed
site budget with held-out reference sites. Measure localization explicitly
(for example DCC or centroid distance), report successes and failures across
all generated sites, and avoid tuning thresholds on final evaluation targets.
See [Binding-site prediction](binding-site-prediction.md) and
[Benchmarking](../benchmarks/index.md).
