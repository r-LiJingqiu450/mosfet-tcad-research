# Si/SiO₂/Lateral-Contact Triple-Point Defect — Investigation and Correction

**Date:** 2026-10-30. **Branch:** `phase3/C2-2026-10-30`. Code and results at `fe44d68`; bit-identical reproduction is in `results/phase3/C2/repro_fe44d68/`.
**Scope:** construction of the R-I(a)-type C-1 deck at V_D = 0, F0 physics (R05 §0.2), DEVSIM 2.11.2.
**The identity numbers are previews.** The T1 error budget (E_eq, E_lin, E_freq; D-2/D-3 analysis; R05 mesh N_y,0 = 50; L ∈ {1, 3, 10} µm) has not been executed. **T1 is not claimed.**

## 1. Summary

- The original C-1 deck lost potential continuity at the two points where silicon, oxide and a lateral reservoir contact meet.
- Cause: a **DEVSIM topology rule in the built-in 2-D mesher**, combined with a **DEVSIM assembly rule for interfaces**. It is not geometry, not physics, and not a bug in the project code.
- Correction: an explicit mesh in which the oxide node *at* the triple point carries a one-node contact with the Si corner's Dirichlet value. Because the Si twin node is itself Dirichlet, this is exactly the discrete potential-continuity condition, not an extra constraint.
- Result: the continuity jump is 0, discrete Gauss closes to ≤ 8.6e-15, and no other oxide node is constrained. The T1 identity preview drops from −0.90 (original) to |r̂₀| ≤ 1.9e-9 on three refinement levels, which is the linear-solve noise floor.

## 2. Root cause (from DEVSIM source, tag v2.11.2rc5 = e3a25f08)

1. **Mesher ownership.** In `src/meshing/Mesh2d.cc` (≈ l. 302–385), every coordinate shared by more than one region goes to **at most one** contact or interface:
   - contacts are searched first, in declaration order, with `break` at the first match;
   - interfaces are searched only `if (!found)`.

   The Si/SiO₂ corner point belongs to the Si lateral contact, so the interface pair there is never created. The mesh log shows 19 of 21 pairs.
2. **The device layer allows sharing.** `src/Geometry/Device.cc::AddContact/AddInterface` only emits "shares a node" INFO warnings. The gmsh/array loader (`GmshLoader.cc`) builds contacts and interfaces independently.
3. **The assembly rule still drops the pair.** `src/Equation/InterfaceEquation.cc::RemoveContactNodesFromList` removes an interface node in a region when a contact in that region solves the same equation at that node. The `continuous` assembly (`NodeVolumeType1Assemble`) skips a pair unless both nodes are active.

   Even with shared ownership (explicit mesh), the pair is skipped because the Si corner is a ψ-Dirichlet contact node. The oxide corner then keeps an ordinary Laplace row, i.e. it floats.
4. **Declaration-order test.** Four orderings were tried: the first-declared contact wins, and the interface never gets the point. Evidence: `order.py` (C-1) and the mesh logs.

**Conclusion:** neither contact geometry nor declaration order can fix this. The oxide corner node needs its own equation.

## 3. Physically correct discrete condition

R-I(a) (R04 §4.2) prescribes:
- on the Si lateral boundary (full depth, *including* the interface point): ψ = ψ₁(x) and φ_n Dirichlet;
- on the oxide lateral boundary: homogeneous Neumann;
- at the interface: continuity of ψ and of ε∂ψ/∂n.

At the triple point, continuity forces ψ_ox(corner) = ψ_Si(corner) = ψ₁(0). In the discrete problem the Si corner node's Poisson row is already replaced by its Dirichlet value. The oxide corner node's correct equation is therefore ψ_ox,c = ψ_Si,c = PsiLat(0) (+ s_couple·bias when the lateral data follow V_D).

This is the continuity statement restricted to a node whose twin is Dirichlet. It does **not** constrain any other oxide node, so oxide lateral Neumann holds everywhere else. The displacement flux into that node becomes the charge of a terminal that is electrically the same as the Si lateral contact:
- `oxc_d` is attached to the drain circuit node `drain_bias`;
- `oxc_s` follows `source_bias`.

## 4. Candidate formulations compared

`phase3/c1/corner_study.py`; outputs in `results/phase3/C2/corner_study/` (v3; v1 and v2 are kept as superseded, see §7). All cases: V_G ∈ {0.8, 1.25, 2.0} V, V_D = 0, L = 1 µm.

| ID | Construction | Mesh | Constraint added | Corner ψ jump (V) | r̂₀ preview L1, 0.8 / 1.25 / 2.0 V | Gauss closure | Verdict |
|---|---|---|---|---|---|---|---|
| A0 | Original C-1 deck: lateral data from the 2-D Neumann run | built-in | none | 0.40 / 0.60 / 1.25 | −0.897 / −0.921 / −0.940 | ≤ 5.9e-16 | defect (regression evidence) |
| A0b | 1COL lateral data, no oxide treatment | built-in | none | 0.22 / 0.40 / 1.02 | +2.06e-2 / +9.6e-3 / +5.8e-3 | ≤ 2.2e-15 | defect |
| A1 | "Legacy pin": oxide node *adjacent* to the triple point (y = −0.243 nm) pinned to 1COL data | built-in | **artificial Dirichlet on an oxide lateral node** | 3.0e-3 / 5.4e-3 / 1.4e-2 | +5.24e-6 / +2.52e-6 / +1.50e-6 | **5.0e-2** (flux into the pinned node is not attributed to any terminal) | **rejected**: violates continuity and oxide Neumann |
| X0 | Explicit mesh, shared ownership, no treatment (control) | explicit L0 | none | 0.22 / 0.40 / 1.02 | +2.08e-2 / +9.7e-3 / +5.8e-3 | ≤ 1.7e-16 | defect reproduced: **not geometric** |
| X2 | Endpoint separation: Si lateral contacts start below the interface node | explicit L0–L2 | none, but the Si corner loses its ψ and φ_n Dirichlet data | 0 | L0 −4.9e-8 / −7.6e-8 / −1.6e-7; L1 −2.4e-8 / −3.6e-8 / −7.4e-8; L2 −1.2e-8 / −1.9e-8 / −3.9e-8 | ≤ 1.1e-15 | **rejected**: changes the R-I(a) BC at one node; error O(h²) (halves per √2 level) |
| **X1** | **Explicit mesh, one-node oxide contact ON the triple point = continuity tie** | **explicit L0–L2** | **none beyond continuity** | **0** | **L0 −7.3e-11 / −1.8e-10 / −3.8e-10; L1 +9.9e-11 / +3.0e-10 / +1.5e-9; L2 +5.8e-10 / −1.4e-9 / −1.9e-9** | **≤ 8.6e-15** | **adopted** |

Checks on X1:

| Item | Requirement | Result |
|---|---|---|
| 1. Mathematical consistency | Oxide corner row = continuity statement | Exported Newton matrix: the oxide-corner row has exactly 1 nonzero (its own column); 0 interface pairs dropped except at the two corners, where both twins are Dirichlet at the same value |
| 2. Potential continuity | ψ_ox = ψ_Si at every coincident pair | max jump **0.0 V** at all pairs, corners included, on all levels and V_G |
| 3. Displacement flux | Oxide-corner flux attributed to the lateral terminal | `oxc_s`/`oxc_d` charge −1.9076e-13 C/cm each at 1.25 V. This equals the gate flux terminating on the space charge of the two Dirichlet contact columns. `oxc_d` is on the drain circuit node, so the AC terminal current includes it |
| 4. Charge conservation | Discrete Gauss: Σ contact charges = Σ non-contact space charge | Closure ≤ 8.6e-15 on all levels. Contact rows identified from the exported permutation (dropped, or moved to a circuit row); interface rows are moved into oxide rows and are not contacts |
| 5. Mesh refinement | Residual must not depend on h | \|r̂₀\| ≤ 1.9e-9 on L0, L1, L2 (9 284 → 36 123 unknowns), with no trend in h. It is the same size as \|G_L1/G_L2 − 1\| (≤ 1.9e-9) and the decoupling changes (≤ 1.9e-9): **noise floor of the linear solves**, growing with system size |
| 6. Artificial Dirichlet | None | Oxide lateral rows with a single entry: only the two triple-point nodes (y = 0) on every level. A1's single-entry rows are at y = −0.243 nm, which is artificial |
| 7. R-I(a) assumptions preserved | Si lateral Dirichlet ψ and φ_n over the full depth incl. corner; oxide lateral Neumann; bottom ψ Dirichlet with zero electron flux; gate Dirichlet | Yes. Lateral ψ uniformity of the V_D = 0 solution 5.6e-17 to 1.1e-16 V |

**R-I(b).** In the explicit mesh, the Si lateral contacts do not solve PotentialEquation in R-I(b) mode, so the corner interface pairs are active and continuity holds by DEVSIM's own rule. This follows from the source reading in §2. **It was not exercised in this package.**

## 5. Not adopted, recorded for the reviewer

**Robin gate (oxide replaced by D·n = C_ox(V_G + ψ_b0 − ψ_s)).** This removes the oxide region and therefore the triple point. R04 §3.2 permits it, but it changes the 2-D reference by removing lateral oxide electrostatics. That alters what ε_red measures, so it is a **modeling decision** rather than a numerical fix. It was not implemented.

## 6. Status for C-3

The corrected construction is X1 (`deck.build_mesh_explicit(..., corner='tie')` + `deck.equilibrium_1col`). It is **ready for a separately authorized C-3 T1 study**, which was not started. C-3 must still:
- build the deck on the R05 meshes (N_y,0 = 50; L ∈ {1, 3, 10} µm);
- provide finite-V_D lateral data ψ₁(x; V_G, V_D) (available from the C-2 1COL deck with φ_n = V_D);
- compute every E component;
- run R-I(b) and verify its corner treatment.

## 7. Failed and superseded attempts (kept)

- `corner_study/v1_superseded`: valid r̂₀ and continuity data. The Gauss check there was ill-defined: it included the contact-node volume charge, giving a mismatch equal to dx/L.
- `corner_study/v2_superseded`: still mis-specified (Si interface rows counted as contact rows).
- `v1_superseded/logs/*_attempt1_FAILED.log`: X0/X2 first attempt, exit 1. The mesh builder declared physical names with no elements.
- The C-1 attempts 1–3 (`results/phase3/C1/c1a/`) are the original defect evidence.
