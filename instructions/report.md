# Report stream instructions

## Scope
- `workflow/report_step2.py` is the main script for this stream.
- The scope here is the interactive HTML report generated from step-2 outputs.
- Restrict attention to report-specific code paths unless the issue clearly crosses into `pipeline`.

## Stream boundaries
- Start with:
  - `workflow/report_step2.py`
  - `workflow/gather_annotations.py` when the report depends on annotation shaping
  - `step2.nf`
  - `step2.smk`
  - `data/species_tree.full.newick`
  - `genefam.csv` or `data/gene_families_searchinfo.csv`, depending on the entrypoint being exercised
  - related tests, demo helpers, and existing generated report outputs

## Regeneration
- After substantial report-side changes, always regenerate the report with `make_report.sh` before closing the task.
```bash
bash make_report.sh
```
- Keep `workflow/report_step2.py` compatible with the Python interpreter reached by `python` in `make_report.sh`; avoid type-hint syntax that breaks import-time evaluation on older interpreters unless annotations are deferred safely.
- If `make_report.sh` is temporarily unsuitable for the specific debugging run, note that explicitly and use the equivalent direct `workflow/report_step2.py` command instead:
```bash
python workflow/report_step2.py \
    --possvm_dir results/possvm \
    --possvm_prev_dir results/possvm_prev \
    --search_dir results/search \
    --cluster_dir results/clusters \
    --family_info data/gene_families_searchinfo.csv \
    --refsps Mmus \
    --refnames data/Mmus_gene_names.csv \
    --species_tree data/species_tree.full.newick \
    --align_dir results/align \
    --output report2.html
```
- If `results/possvm_prev/` is absent for the run you are testing, omit `--possvm_prev_dir`.

## Performance
- Keep the report self-contained as one HTML file, but treat file size, initial browser memory pressure, and on-demand decompression cost as first-class concerns.
- When `family_info` / `genefam.csv` limits the report to a subset of families, loaders should push that filter down before opening/parsing result files whenever the family can be inferred from the filename.
- Prefer lazy payloads over eager JSON when the data is only needed after a user action.
- Keep exact-domain and architecture payloads compact by interning repeated strings and storing them as gzip+base64 blobs.
- The current dominant report footprint is the per-HG tree detail payload, not the domain-architecture features.
- When reducing footprint, attack in this order unless the task clearly targets something else:
  - `treedata-*` payloads
  - `alndata-*` payloads
  - duplicated per-gene metadata inside per-HG payloads
  - only then smaller report-side catalogs such as exact-domain and architecture data

### Current performance contract
- `treedata-*` should be stored as gzip+base64 blobs in the HTML, not raw JSON text.
- `loadDetail()` and the alignment tree loader should inflate tree detail on demand and cache the parsed object so the same HG is not reparsed repeatedly.
- The Alignments tab `Alignment` toggle should hide the sequence/ruler pane and avoid drawing those canvases when hidden, so tree + compact domain-architecture browsing remains lighter.
- Lazy payload compression is intended to reduce initial page weight and initial memory pressure; it does not eliminate the live-memory cost of a tree once that HG has been opened and parsed.

### What to measure
- `python -m py_compile workflow/report_step2.py`
- full report regeneration with `bash make_report.sh`
- final HTML size:
```bash
stat -c '%y %s' report2.html
```

## report_step2.py architecture
- Single large file: Python data-loading layer plus one giant `HTML_TEMPLATE` raw string with HTML/CSS/JS inline.
- `main()` calls `build_report_context(args)` then `render_report_html(context)` which replaces `%%...%%` placeholders and writes one self-contained HTML file.
- Alignment data and tree detail are embedded as lazy `<script type="application/json">` payloads and decompressed in-browser with `pako.js`.

## Imperative source-awareness rule
- Treat each HG as potentially having two distinct annotation sources:
  - GeneRax-backed POSSVM outputs in `results/possvm/`
  - non-GeneRax / original IQ-TREE2 POSSVM outputs in `results/possvm_prev/`
- Do not assume orthogroup names, gene-to-OG assignments, reference orthologs, OG support, or reference support are identical across those two sources.
- Any report utility that displays, caches, filters, or navigates by OG names or reference-ortholog metadata must switch those annotations together with the active source.
- Flattening source-specific annotation metadata into one HG-global view is a bug unless the caller explicitly wants a merged summary.

## Python loaders and data shaping
- `get_species_prefix(gene_id)` strips pipe annotations and returns the species prefix.
- `load_family_info(csv)` and `load_family_details(csv)` load family metadata.
- `build_family_records(search_dir, family_info)` and `build_hg_records(cluster_dir, family_info)` assemble records for report tables.
- `load_possvm_trees(possvm_dir, source)` loads POSSVM outputs and per-gene metadata.
- `load_domain_hits(search_dir)` loads `*.domains.csv` domain calls.
- `load_reference_names(path, refsps)` loads reference names.
- `load_gene_lengths(cluster_dir)` loads protein lengths from FASTA.
- `parse_clade_groupings(species_tree_path)` derives clade color groupings from the species tree.
- `load_species_info(path)` returns `({prefix: full_name}, {prefix: group})` from `data/species_info.tsv`.
- `load_species_images(img_dir)` returns `{species_prefix: data_uri}` from `img/phylo/{prefix}.png`; it should silently return `{}` if the directory is absent.

## Species data and colours
- `SPECIES_INFO` is `{species_prefix: full_name}` and is shown in italic in species-facing UI.
- `SPECIES_GROUPS` is `{species_prefix: group_name}` from the optional third column of `species_info.tsv`.
- `SPECIES_IMAGES` is `{species_prefix: data_uri}` and powers inline species silhouettes in the report.
- `SPECIES_COLOR_PRESETS` and `applySpeciesPalette(name)` define the species-colour system used in the Species Tree tab.
- `refreshSpeciesColorViews()` is the shared redraw path after palette changes or single-species recolouring.
- `groupColors` / `recomputeGroupColors()` derive group colours from current species colours.
- `_cssToRgb(css)` and `_rgbToHex(r,g,b)` are the shared colour helpers.

## Species-group startup behavior
- Species groups must not be auto-collapsed on startup.
- `SPECIES_GROUPS` is metadata for colouring, badges, and optional alignment columns; it should not implicitly change the initial expanded/collapsed state of the Species Tree or Counts cladogram.

## Six tabs
- `families`: `drawFamilyTable()` renders the family table and category/family navigation.
- `architectures`: the family-centric domain-architecture browser lives here and should remain family-first rather than HG-first.
- `sptree`: `drawSpeciesTree()` renders the full species tree SVG with collapsible clades, prune-to-data, and heatmap split-group interactions.
- `heatmap`: `drawHeatmap()` plus `drawCladogram()` render the species-by-HG/OG matrix and cladogram.
- `trees`: `renderSidebar()`, `selectTree()`, and `renderTree()` render gene trees and their controls.
- `align`: `renderAlnSidebar()`, `selectAlignment()`, and `renderAlignment()` render the per-HG alignment viewer.

## Species Tree tab contract
- Full species tree SVG with collapsible clades (`spCollapsed`), prune-to-data toggle, annotation TSV / Newick download, and shift-click on internal nodes to add heatmap row groups.
- The Species Tree tab currently uses the rectangular renderer only.
- The layout is split between the tree on the left and a pinned species-details panel on the right.
- Clicking a species tip, label, or logo in the rectangular tree should pin that species into the right-side panel instead of opening a popup.
- The right-side species panel should show the species image, code, full name, group badge, summary counts, a species-colour picker, and annotation download actions.
- The annotation download actions in the species panel should include `All genes` plus one button per class, sorted by descending gene count.
- A shared `Logos` toggle controls whether the species-image column is shown, and toggling it should keep all label starts aligned.
- A `Names` toggle should switch tree labels between species prefixes and full names when `SPECIES_INFO` provides them; full names should render in italic without switching to bold.
- Tip labels show inline 24x18 species images from `SPECIES_IMAGES` when logos are enabled.
- Optional group pill badges should stay right-aligned, use the averaged group colour from `groupColors`, and reserve enough right-side space that they never overlap species labels.
- Collapsed clades still render each leaf tip aligned to the far-right tip column; the triangle base should stop slightly left of those tip circles.
- Species tab includes a palette selector plus Apply button, with Tableau as the default startup palette.
- Clicking an internal node should continue to open `showSpNodeActionPopup()` with `Collapse` and `Flip`.
- Root node keeps the flip-only dot behavior.

## Families tab contract
- Category rows stay expandable/collapsible, and each family row should further expand into a per-family HG browser.
- Classes, families, HGs, and OGs should default to largest-to-smallest ordering rather than alphabetical ordering.
- Expanded family details should list HGs showing:
  - HG id
  - number of sequences
  - number of species containing the HG
  - number of orthogroups
- each HG line should also show a trimmed grey comma-separated summary of reference sequence names when available
- If an HG has both GeneRax and original IQ-TREE annotations, the HG row should show both OG counts, labeled distinctly rather than flattening them into one value.
- The per-family HG list should be sortable by its visible columns.
- Expanding an HG should lazy-load and list that HG's orthogroups.
- The per-HG OG list should be sortable by its visible columns.
- Per-OG metadata extraction should be factored into a shared helper rather than embedded directly in the Families tab row renderer, because the same OG-level annotations will be extended later.
- That shared per-OG metadata should currently include:
  - MRCA node name derived from the species tree
  - dominant domain architecture for the OG, computed from its member genes and rendered with the same chip style used in the Domain Architectures tab
- Clicking an OG row in the Families tab should offer:
  - `Show counts`: open the Counts heatmap focused on that single OG within its HG
  - `Show gene tree`: open the Gene Trees tab with that OG focused
  - `Show alignment`: open the Alignments tab with other OGs collapsed
- OG popup actions are source-aware:
  - if the visible OG row comes from GeneRax annotations, `Show gene tree` and `Show alignment` must open the GeneRax source
  - if the visible OG row comes from original IQ-TREE annotations, those actions must open the IQ-TREE source
- the popup labels should make that source explicit rather than leaving the destination ambiguous
- The old explicit `Action` column in the OG list should not be shown; the OG row itself is the click target for the popup actions.
- Families search must stay lightweight while typing:
  - matching families may be filtered live, and OG hits may appear in the dropdown
  - typing a broad query must not auto-expand every matched family/HG subtree
  - only the explicitly selected OG result should expand its family/HG path and scroll to the matching row

## Alignment viewer contract
- The alignment viewer is canvas-based and uses lazy gzip+base64 payloads embedded in the HTML.
- Branch-length mode should be on by default.
- The `Alignment` toggle should hide the sequence/ruler pane while leaving the left tree/names/domain-architecture side visible.
- The Alignment tab species filter should support the same clade/species selection model as the Gene Trees tab:
  - free-text tags can resolve named species-tree clades as well as species substrings
  - a `Species tree` button should open the shared mini species-tree picker and add the selected clade/species filter tags
- When `SPECIES_GROUPS` metadata exists, the alignment controls should expose a `Species Group` toggle that adds a dedicated metadata column showing the group name for each sequence.
- Amino-acid residues in the sequence pane should render as plain monospace text, not bold, and slightly smaller than the previous styling.
- Column headers should be explicit, draggable, and have visible resize grips in the header row itself.
- The compact domain-architecture column should draw exact hits at real relative coordinates within the phylogeny span, not as stacked multiplicity chips.
- When most OGs are collapsed to single rows, the alignment gene-tree column must simplify collapsed subtrees instead of overplotting every original branch onto the same row.
- When a current-family lookup misses but the protein has exactly one domain track, the viewer should fall back to that single track.
- Exact-domain track matching must accept both bare family names like `Unc13` and full track ids like `neu.Unc13`; otherwise shared-domain genes can incorrectly collapse to `No domains`.
- Hover on the compact domain-architecture column should use the richer per-protein popup built from the lazy exact-domain catalog.
- Wheel scrolling over the left tree/names/header side should still move the alignment vertically.

## Alignment source-toggle rule
- When the alignment viewer toggles between GeneRax and IQ-TREE2 for a given HG, it must update all source-dependent gene metadata fields together:
  - gene tree / cladogram
  - OG assignments and OG-derived grouping / MRCA summaries
  - reference ortholog labels
  - OG support / reference support if displayed
- If a tooltip, popup, search helper, or grouping feature depends on OG names or reference-ortholog annotations, it must read from the currently active source-specific metadata.
- Direct navigation into the Alignment tab from Families/Counts OG actions must pass an explicit source preference; reusing whatever source the alignment viewer happened to have selected before is a bug.

## Gene tree viewer notes
- `renderTree()` lays out the tree as either a cladogram or a phylogram.
- Internal-node collapse uses `_children` and triangle/circle collapsed renderings.
- `annotateOGNodes()` marks MRCA nodes for OG labels.
- The shared tooltip layer is used by tree hovers and must stay above popups.
- Protein hover in gene trees should reuse the lazy exact-domain catalog rather than duplicating eager per-gene payloads.
- The Gene Trees tab should start with `Width = 0.5`.
- The default `Fit` / initial viewport behavior in the Gene Trees tab should left-align the fitted tree rather than centering it horizontally.
- In `Highlight OGs` mode, clade-highlight labels should scale with the same zoom-aware label sizing used by tree labels rather than staying near a fixed pixel size.
- The mini species-tree picker used for clade/species selection should be shared with the Alignment tab rather than implemented twice.

## Utility scripts
- `workflow/download_phylopic.py` is the utility for populating `img/phylo/` from PhyloPic: search by taxon name, resolve the chosen node's `primaryImage`, and download the `192x192` PNG thumbnail. Prefer this over ad hoc scraping when adding species silhouettes for the report.

## Alignment defaults and separators
- The Alignment tab should start with a wider gene-tree column than before; the current default is `240px`.
- OG separators in the Alignment tab should render as black horizontal lines across metadata and sequence panes, while the gene-tree subcolumn keeps grey separators.
