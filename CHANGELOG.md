# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/)
and this project follows semantic versioning principles.

---

## Version: feature/resolution-survival-optimization

### Date 2026-01-29

Wire resolution survival analysis into Dash UI

- Run resolution_survival_analysis on study selection
- Generate survival-based resolution outputs automatically
- Add Resolution → Survival page to Dash
- Display survival ranking table and score curve
- Add-only changes, no breaking modifications

### Date 2026-01-28


Add multi-resolution survival analysis with RSS v1

- Implement survival-based resolution selection (OS, PFS)
- Support robustness criteria: min_patients and min_events
- Compute Resolution Survival Score (RSS v1)
- Write survival rankings, plots, and best resolution per endpoint
- Add-only change, no impact on existing pipeline



## Version: pre-survival-refactor-stable

### Date: 2026-01-27

### Added

- Introduced structured survival handling with OS/PFS separation
- Added parent Survival Analysis page with tab-based navigation
- Added dedicated OS and PFS survival subpages
- Added automatic migration from legacy clinical survival config
- Added normalization and validation of survival event columns
- Added creation of OS_EVENT_string and PFS_EVENT_string for legacy plots
- Enabled suppress_callback_exceptions for dynamic multipage layout

### Changed
- Refactored routing for survival pages:
  - /survival
  - /survival/os
  - /survival/pfs
- Updated sidebar to support hierarchical survival navigation
- Updated study loading pipeline to support multiple survival endpoints
- Enforced numeric casting for survival time/event columns (lifelines compatibility)

### Fixed

- Fixed survival normalization producing object dtypes
- Fixed lifelines crashes due to invalid event formats
- Fixed broken cluster comparison callbacks
- Fixed duplicated/invalid Dash callback registrations
- Fixed routing conflicts between legacy and new survival pages
- Restored functionality of Cluster Comparison page after routing refactor

## Known Issues

- OS and PFS survival pages currently render empty (UI lifecycle issue under investigation)
- Survival subpage callbacks may not trigger on initial load
- Survival visualization refactor not yet complete

### Technical Notes

- Survival data preparation centralized in prepare_survival_columns
- Legacy survival fields migrated to config["survival"]
- Tab state synchronized with URL routing
- Global clinical dataframe updated at study selection time


## [Unreleased]

### Added
- Multi-resolution clustering analysis across resolutions from 1.0 to 0.1.
- Automatic identification of centroid genes for each cluster and resolution.
- Export of `resolution_centroids.tsv` containing all cluster centroid memberships.
- Gene-level centroid stability analysis across resolutions.
- Cluster-level centroid stability analysis across resolutions.
- Automatic selection of the best resolution based on stability criteria.
- Sankey plot (gene-level) showing patient flow across resolutions.
- Heatmap visualizations:
  - Gene-centroid vs resolution stability.
  - Cluster-centroid vs resolution stability.
- Model complexity analysis:
  - Resolution vs number of centroid genes curve.

### Dash UI
- New **Resolution Analysis** section in the dashboard.
- Added sub-pages:
  - **Overview**: best resolution summary and complexity curve.
  - **Gene Stability**: gene-level heatmap and stability table.
  - **Cluster Stability**: cluster-level heatmap.
  - **Sankey Flow**: interactive Sankey visualization of patient transitions.
- Embedded interactive Plotly Sankey visualization via iframe.
- Robust routing using URL-based callbacks.
- Clean separation between computation (backend) and visualization (frontend).

### Changed
- Improved internal pipeline organization for multi-resolution analyses.
- Standardized output generation under study output directories.

### Fixed
- Dash routing and callback triggering issues related to page navigation.
- Improved robustness of file loading and missing-output handling in the dashboard.

---

## [0.1.0] – Initial Release
- Initial implementation of Tortoise clustering and pathway analysis framework.
