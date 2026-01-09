# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/)
and this project follows semantic versioning principles.

---

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
