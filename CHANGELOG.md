# Changelog

## 1.2.0 - validation, sweeps, and publishing

- **Datasheet check:** transcribed all 57 points of the Copeland ZR34K3E-TFD sheet
  (`examples/copeland_zr34k3e_tfd_datasheet.csv`, every published column). The model
  reproduces published mass flow within 0.14 % (worst 0.40 %) and isentropic
  efficiency within 0.09 points. The previous bubble-point convention gave 4.6 %.
  The check is shown in the Test Validation tab with parity plots.
- **Full compressor map:** the example map grows from 27 to 57 points (condensing
  up to 150 °F). The map envelope is now the convex hull of the data, not a
  rectangle, so the fit cannot extrapolate into unpublished corners.
- **Design Sweep tab** (`sweep.py`): re-solves the full system across a range of
  one variable, plots COP/capacity and saturation temperatures, quantifies
  diminishing returns, and keeps infeasible points with their reason.
- **UI:** fixed button labels becoming unreadable (dark text on teal) after a click.
- **Publishing:** the README is rewritten for GitHub, with screenshots in
  `docs/images`. Added a GitHub Actions test workflow and a Streamlit Community
  Cloud-ready layout (the app runs without `pip install`).
- **Docs:** consolidated into `docs/`. Certification-readiness, implementation-status,
  verified-run, and model-fidelity notes were folded into `docs/VALIDATION.md` or
  dropped as out of date.

## 1.1.0 - robustness review

### Physics fixes
- **Plate secondary side:** laminar Nu changed from 4.36 (round tube) to 7.54
  (parallel plates). Laminar friction changed from 64/Re to 96/Re. Dittus-Boelter
  replaced by Gnielinski, with interpolation through the transition range.
  Default plate case: h 508 → 866 W/m²K, UA 754 → 1104 W/K.
- **Refrigerant-state check:** an outlet below saturation was labelled
  "superheated vapor" and gave a negative duty. It is now rejected. An outlet
  exactly at saturation is now handled. Blends are detected from glide, not a
  hard-coded list.
- **Saturation convention:** the condensing pressure now uses the dew point
  (AHRI 540 / EN 12900), with subcooling from the bubble point. This matters for R407C.
- **Zero superheat or subcooling** no longer crashes (CoolProp cannot fix a
  saturated state from P and T).
- **Compressor-map closure** now converts map capacity with the map's rating
  superheat and subcooling, then corrects to operating conditions.

### Data fix
- `examples/copeland_zr34k3e_tfd_map.csv` stored Btu/h in the `capacity_w`
  column, which implied COPs above Carnot. Converted to W. The original values
  are kept in `capacity_btu_h`.
- Map validation now rejects any row with COP above the Carnot limit.

### Solver robustness
- New `solver.py`: multi-start least squares plus a bracketed fallback.
  R290 on default inputs previously crashed with a temperature-cross error; it now solves.
- Infeasible cases raise `InfeasibleOperatingPoint` with the limiting exchanger and suggested fixes.
- Bounds respect the refrigerant's critical temperature and the map envelope.

### App
- Results persist when switching tabs, and a note appears if inputs changed after a run.
- Solver convergence warnings shown; Carnot COP and second-law efficiency reported.
- P-h diagram, heat-exchanger profiles, and state-point table added to the solver tab.
- Air uses CoolProp dry-air properties. Brine/custom fluids require an explicit specific heat (they previously used water silently).
- Bundled Copeland R407C map loads with one click; its rating conditions pre-fill the solver.
- A stale compressor fit is cleared when a different map is loaded.
- Superheat, subcooling and map rating conditions exposed; help text on inputs.
- "Project A" label removed, misleading active-card highlight removed, usage note matches tab order.
- Compatible with Streamlit versions before and after the `use_container_width` removal.

### Tests
- 25 → 53 tests, including hand-calculation references and a headless app smoke test (59 after 1.2.0).
