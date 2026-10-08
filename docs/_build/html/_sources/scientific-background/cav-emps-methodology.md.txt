# CaV-EMPS methodology

CaV-EMPS is Ultidock's receptor-derived site-proposal implementation. Its entry
point is `autogenerate_centers_tsv` in `docking/make_grids.py`. It combines geometric
cavity evidence and energetic map evidence to construct candidate docking boxes.
This page describes the implementation, not a claim of universal predictive accuracy.

## Evidence and candidate families

The geometric path rasterizes receptor occupancy and uses a Euclidean distance
transform to describe available space. Internal-cavity candidates and surface-cleft
candidates address different receptor geometries. The map-based path uses AutoGrid
fields, including carbon, electrostatic and desolvation information, with favorable
contributions, spatial smoothing and receptor-proximity constraints. Neighborhood
aggregation prevents the method from relying solely on an isolated favorable voxel.

Candidates are filtered, grouped and selected with spatial separation constraints.
The resulting centers and grid dimensions define the downstream search regions.
A proposal's internal score is a site-selection quantity, not a ligand binding energy.

## Search policies

| Generated `SITE_POLICY` | Candidate selection |
| --- | --- |
| `receptor_search` | Compact portfolio of complementary receptor-derived sites; default |
| `exhaustive_search` | Broader search across candidate families |
| `internal` | Internal-cavity family |
| `surface` | Surface-cleft family (`maps` is an internal alias) |
| `hybrid` | Consensus family |

The requested site count is a search budget; geometry and filtering affect the
actual proposals. Metadata in the centers file records the policy and parameters
used for reuse checks. Preserve it when comparing runs. For controlled changes to
non-CLI settings, see [Configuration](../reference/configuration.md).

## Validation boundary

Map resolution, receptor extent, protonation and box parameters can alter results.
A site omitted from the proposal set cannot be recovered by downstream local docking.
Evaluate geometry-only, map-based and combined settings under a fixed site budget
before attributing an improvement to a particular component. Use held-out reference
sites and report failures as well as successes. See [Binding-site prediction](binding-site-prediction.md).
