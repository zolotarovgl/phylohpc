# Claude instructions

## Project overview
PhyloHPC is a full phylogenomics pipeline that goes from multi-species proteomes to homology groups, gene trees, ortholog calls, species-level annotation tables, and a final interactive HTML report for hierarchical orthogroups (hOGs) across nested taxonomic clades.
The current step-4 source of truth is `step4_ancestry.smk` plus its helper scripts. If `README.md`, `step4.ancestry.nf`, and the code disagree, follow the `.smk` workflow and the Python scripts it calls.

## Manual summary
- The older `docs/manual.md` frames the intended repository contract as:
  - `step1.nf` for PFAM/domain search plus HG clustering
  - `step2.nf` for alignment, gene trees, POSSVM, and optional GeneRax
  - `workflow/gather_annotations.py` for per-species ortholog annotation tables
  - local and SLURM execution as first-class modes
  - `ids.txt` as the main handoff between step 1 and step 2
  - optional `resources.tsv` resource prediction for large step-2 runs
- Treat that manual as high-level intent, not exact implementation detail. Use it to preserve the overall pipeline purpose and user workflow, but prefer the current code when behavior differs.

## Stream selection
- Before doing substantial work, ask the user which stream to work on unless they already made it explicit.
- The stream names are:
  - `full_pipeline`
  - `hogs`
- If the user does not specify a stream but the request clearly maps to one of them, say which stream you are assuming and proceed.
- Once a stream is chosen, restrict attention to the relevant files and workflows unless the user explicitly asks for cross-stream work.

## Current workstreams
- `full_pipeline`
  - scope: improve the main pipeline itself
  - focus on `step1.nf`, `step2.nf`, `workflow/report_step2.py`, and supporting scripts directly used by those stages
  - preserve the end-to-end contract from proteomes to HGs, trees, POSSVM outputs, annotations, and report generation
  - prefer improvements that make local/SLURM runs more robust, outputs more consistent, and reports easier to interpret
  - for report work, see `## report_step2.py — architecture` section
- `hogs`
  - scope: test and stabilize step 4 / hOG reporting
  - treat `step4_ancestry.smk` as the active implementation
  - prioritize correctness of clade extraction, clade-level POSSVM outputs, cross-level OG linking, output layout, and `hog_hierarchy.html`
  - when changing step-4 logic, verify that downstream report-building still works on the expected file names and directory structure

## Pipeline layout
- Inputs:
  - `data/input.fasta` — concatenated proteomes
  - `species_list` — species prefixes included in the run
  - `genefam.csv` — family definitions for step 1 / step 2
  - `data/species_tree.full.newick` — named internal nodes; canonical full species tree
  - `data/species_tree.newick` — binary/pruned tree used in some downstream contexts
  - `data/Mmus_gene_names.csv` — reference gene-name mapping for POSSVM naming
- `step1.nf`:
  - `SEARCH` runs PFAM/HMMER domain search through `phylogeny/main.py hmmsearch`
  - `CLUSTER` runs MCL clustering through `phylogeny/main.py cluster`
  - Main outputs: `results/search/` and `results/clusters/`
- HG selection:
  - `workflow/select_hgs.py` selects the HGs that proceed downstream
  - `ids.txt` is the main step-2 input list
  - `ancestry_ids.txt` is the step-4/report subset when focusing on ancestry/hOG output
- `step2.nf`:
  - `ALN` creates trimmed MSAs in `results/align/`
  - `PHY` creates per-HG gene trees in `results/gene_trees/`
  - optional `GR_watcher` runs GeneRax and writes reconciled trees to `results/generax/`
  - `PVM` / `PVM_PREV` run POSSVM and write ortholog outputs to `results/possvm/` and `results/possvm_prev/`
  - `REPORT` builds `results/report_step2.html`
- `workflow/gather_annotations.py`:
  - combines search results, cluster membership, and POSSVM outputs into per-species annotation tables
- `step4_ancestry.smk`:
  - authoritative step-4 pipeline for the current hOG report
  - extracts clades, reruns POSSVM at each clade level, and builds the hierarchy report
  - main outputs live under `results/ancestry/`

## Desired functionality
- Preserve the full pipeline contract from `step1.nf` through `step4_ancestry.smk`; step 4 is downstream of HG selection and step-2 tree/orthology outputs.
- The step-4 goal is hierarchical orthogroup inference across nested clades, not just standalone per-clade OG tables.
- Keep the old-manual user workflow recognizable: install env, run step 1, filter HGs, run step 2, gather annotations, then run/test step 4 reporting.
- For each HG, step 4 should:
  - extract in-clade and out-of-clade species for every requested node
  - rerun POSSVM separately for each clade level
  - filter out-of-clade genes before cross-level comparisons
  - link broader-level OGs to narrower-level OGs by shared gene membership
  - summarize splits/retentions across consecutive levels in a single report
- The preferred final artifact is `results/ancestry/hog_hierarchy.html`, a self-contained HTML report for exploring hOG relationships across levels.
- When step-4 behavior is ambiguous, prefer `step4_ancestry.smk`, `workflow/build_hog_report.py`, `workflow/link_hog_levels.py`, and `workflow/visualize_hog_hierarchy.py` over older docs or alternate pipeline drafts.

## Key files
- `step1.nf` — step 1 Nextflow pipeline (PFAM/HMMER search → MCL clustering)
- `step2.nf` — step 2 Nextflow pipeline (alignment → gene tree inference → POSSVM, with optional GeneRax)
- `workflow/gather_annotations.py` — builds per-species annotation tables from search, cluster, and orthology outputs
- `step4_ancestry.smk` — authoritative step 4 Snakemake pipeline (extract clades → clade-specific POSSVM → build hOG report)
- `step4.ancestry.nf` — alternate/in-progress Nextflow version of step 4; not the main source of truth right now
- `workflow/visualize_hog_hierarchy.py` — **main visualization script**; contains `HTML_TEMPLATE` (all CSS/HTML/JS inline), `build_data()`, `parse_newick_tree()`
- `workflow/build_hog_report.py` — assembles the final HTML directly from ancestry outputs and imports `HTML_TEMPLATE`
- `workflow/link_hog_levels.py` — computes parent/child OG links across adjacent clade levels for one HG
- `workflow/extract_clade.py` — derives `{node}.in_species.txt`, `{node}.ignore_species.txt`, and `{node}.pruned.tree`
- `workflow/report_step2.py` — self-contained interactive HTML report for step-2 outputs; 4 tabs: Families, Species Tree, Counts (heatmap), Gene Trees; see `## report_step2.py — architecture` section for full detail
- `config_ancestry.yaml` — pipeline config (node names, paths, ref species)
- `ids.txt` / `ancestry_ids.txt` — HG lists used for step 2 and step 4
- `data/species_tree.full.newick` — full species tree with named internal nodes
- `results/ancestry/` — step-4 outputs; `{node}/possvm/{hg}.ortholog_groups.csv`, clade species lists, pruned trees, and `hog_hierarchy.html`
- `results/generax/` — GeneRax-reconciled gene trees (`{hg}.generax.tree`); fed into POSSVM

## Stream file boundaries
- For `full_pipeline`, start with:
  - `step1.nf`
  - `step2.nf`
  - `workflow/report_step2.py`
  - `workflow/gather_annotations.py`
  - `workflow/select_hgs.py`
  - `workflow/predict_resources.py`
  - related tests under `tests/` for those components
- For `hogs`, start with:
  - `step4_ancestry.smk`
  - `workflow/build_hog_report.py`
  - `workflow/visualize_hog_hierarchy.py`
  - `workflow/link_hog_levels.py`
  - `workflow/extract_clade.py`
  - `config_ancestry.yaml`
  - `ancestry_ids.txt`
  - related tests and `results/ancestry/` outputs
- Do not broaden into the other stream unless:
  - the user explicitly asks
  - the bug clearly crosses the stream boundary
  - a dependency in the other stream is the direct cause of the issue

## Regeneration command
Always run after editing `visualize_hog_hierarchy.py` or `build_hog_report.py`:
```
python workflow/build_hog_report.py \
    --ancestry_dir results/ancestry \
    --nodes "Metazoa,Bilateria,Euarchontoglires" \
    --ids ancestry_ids.txt \
    --output results/ancestry/hog_hierarchy.html
```

## hOG hierarchy HTML — architecture
- Single self-contained HTML file with embedded D3.js (v7, CDN)
- Data embedded as `%%DATA_JSON%%` → replaced at build time
- Main views:
  1. **Sankey/metro diagram** — hOG flow across taxonomic levels (main view)
  2. **OG detail popup** (`#popup-backdrop`) — species tree + presence/absence for a selected OG
  3. **Gene tree popup** (`#gtpopup-backdrop`) — full gene tree with OG membership boxes

## OG detail popup — species tree (renderPopupTree)
Implemented in `renderPopupTree()` inside `HTML_TEMPLATE`. Key behaviors:
- Shows **full species tree** for the OG's taxonomic level (`TREES[node.level]`)
- **Cladogram layout**: leaves aligned to right, internal nodes at their depth; uses `layout(n, yStart)` returning next yStart
- **Collapsible subtrees**: click internal node → collapses to triangle showing `[nPresent/nTotal]`; click triangle to expand. Collapse state in `popupState.collapsed` (Set of nodeKey strings). nodeKey = sorted leaf names joined by ","
- **Row height slider**: `#popup-rh-slider` controls `popupState.rowH`; `popupRhChanged()` re-renders
- **Presence/absence**: computed bottom-up from `popupState.presentSpecies` (Set of species names in the OG)
- **MRCA** (Dollo parsimony): deepest node where ≥2 children are present, or a single present leaf. Drawn as blue circle (`#2980b9`), always-visible label, shown in side panel
- **Secondary losses**: absent tips that are descendants of the MRCA. Colored red (`#e74c3c`). Internal loss-root nodes also drawn red
- **Hover**: per-tip mouseover shows species name (shared floating `<text>` element moved to end of SVG so it renders on top of circles)
- **Side panel**: species counts (with OG, without OG, secondary losses), MRCA clade name, loss node list

## Gene tree popup — renderGeneTree
- Parses POSSVM newick (`meta.gene_tree`) with `parseNewickJS()`
- Cladogram or branch-length layout
- **OG membership boxes**: colored rects behind tree (one column per level) + bar panel columns to right
- **Legend** (level color swatches): clicking a swatch **toggles** that level's boxes (both background highlight AND bar column) via `gtState.hiddenLevels` Set
- Column header abbreviations ("Met", "Bil", etc.) click to isolate that level
- Hovering any box shows OG name tooltip (`showTooltip()`, z-index 200 to appear above popup)
- Tax level header font-size: 10

## report_step2.py — architecture

Single ~6500-line file: Python data-loading layer + one giant `HTML_TEMPLATE` raw string (all HTML/CSS/JS inline, D3 v7 CDN). `main()` calls `build_report_context(args)` then `render_report_html(context)` which does `.replace("%%PLACEHOLDER%%", json.dumps(value))` substitutions and writes a self-contained HTML file.

### Python layer (lines 1–490, 6280–6456)

**Key loaders** (all return plain dicts/lists, no side effects):
- `get_species_prefix(gene_id)` — strips pipe annotations, returns prefix before first `_` or `.`
- `load_family_info(csv)` → `{family: cls}` from `genefam.csv`
- `load_family_details(csv)` → `{family: {pfam, category, cls}}`
- `build_family_records(search_dir, family_info)` → list of `{id, family, class, species_counts, total}` from `*.genes.list`
- `build_hg_records(cluster_dir, family_info)` → same shape from `*.fasta`
- `load_possvm_trees(possvm_dir, source)` → `(records, all_species, gene_meta)` from `*.ortholog_groups.newick` + `.csv`; each record has `{id, hg, family, tree_dict, ogs, og_names, species, ...}`
- `load_domain_hits(search_dir)` → `{gene_id: [{name,start,end}]}` from `*.domains.csv`
- `load_reference_names(path, refsps)` → `{gene_id: ref_name}`
- `load_gene_lengths(cluster_dir)` → `{gene_id: int}` from FASTA sequences
- `parse_clade_groupings(species_tree_path)` → list of `{name, groups:{sp: clade_name}}` for color groupings
- `load_species_info(path)` → `({prefix: full_name}, {prefix: group})` tuple from a 2–3 column TSV (`data/species_info.tsv`); 3rd column is the optional group name
- `load_species_images(img_dir)` → `{species_prefix: data_uri}` base64-encoded PNGs from `img/phylo/{prefix}.png`; silently returns `{}` if dir absent

**Assembly**:
- `build_report_context(args)` — runs all loaders, merges, returns one dict
- `render_report_html(context)` — injects 14 JSON blobs into `HTML_TEMPLATE` via `%%PLACEHOLDER%%`
- `_build_lazy_scripts(records, prev_records)` — per-HG `{tree, ogs, prev_tree?, prev_ogs?}` as `<script type="application/json" id="treedata-{id}">` tags (lazy-loaded on HG selection, not in initial JS constants)

### JS data constants (injected at build time)

```js
SPECIES_ORDER   // [str] — ordered leaf names from species tree
SP_TREE_DATA    // nested dict {name,dist,children?,leaf?}
FAMILY_DATA     // [{id,family,class,species_counts,total}]
HG_DATA         // [{id,family,hg,class,species_counts,total}]
TREE_INDEX      // [{id,hg,family,species,og_names,n_ogs,has_prev,source,...}] — no tree_dict
ALL_SPECIES     // [str] all species seen in gene trees
CLADE_DATA      // [{name,groups:{sp:clade}}] from parse_clade_groupings
NEWICK_RAW      // str — raw newick for download
FAMILY_INFO     // [{family,pfam,category,cls,n_hgs,n_trees,n_generax,total,n_species}]
HAVE_GENERAX    // bool
NO_TREE_GENES   // {hg_id: {species:[gene_ids]}} for HGs without trees
DOMAIN_DATA     // {gene_id: [{name,start,end}]}
GENE_META       // {gene_id: {length,og_support,ref_ortholog,ref_support}}
REFNAME_MAP     // {gene_id: ref_name}
SPECIES_INFO    // {species_prefix: full_name} — shown in italic in tip hover
SPECIES_GROUPS  // {species_prefix: group_name} — optional 3rd col of species_info.tsv; drives group pill badges
SPECIES_IMAGES  // {species_prefix: data_uri} — base64 PNG; shown inline (24×18px) next to tip label and full-size in hover
```

**Key species-tree runtime state / functions**:
- `speciesOrder` — `let` array derived from `SPECIES_ORDER`; mutable; used by heatmap and cladogram for row order
- `recomputeSpeciesOrder()` — re-walks `SP_TREE_DATA` leaves to rebuild `speciesOrder`; call after any flip so heatmap rows and cladogram stay in sync
- `showSpNodeActionPopup(ev, n)` — shows `#sptree-node-popup` at click; **Collapse** → `spCollapsed.add`; **Flip** → `spNodeById.get(n._spId).children.reverse()` + `recomputeSpeciesOrder()` + redraws all views; root node shows flip-only dot (no collapse option)
- `groupColors` — module-scope IIFE `{group_name: hex}`; computed once at startup as average RGB of member species colors using `_cssToRgb`/`_rgbToHex` helpers; referenced in `drawSpeciesTree()` for pill badges and in `leafColor()` for `"group"` color mode
- `_cssToRgb(css)` / `_rgbToHex(r,g,b)` — module-scope color helpers; handle both `#rrggbb` hex and `rgb(r,g,b)` strings (D3 interpolators return the latter)

### Four tabs

| Tab id | JS entry | Content |
|--------|----------|---------|
| `families` | `drawFamilyTable()` | Two-level accordion: **category rows** (aggregate stats — fam count, HGs, trees, species) expand to show per-family rows. State in `_famExpandedCats` Set; `famToggleCat(cat)` toggles. Auto-expands all categories when `_famFilter` is non-empty. Columns: category/family, PFAM links, class badge, species frac, Fam. count, HGs, trees, GeneRax, genes. `famGoToFamily()`/`famGoToClass()` navigate to Counts tab |
| `sptree` | `drawSpeciesTree()` | Full species tree SVG; collapsible clades (`spCollapsed` Set); prune-to-data toggle; shift+click internal node → adds to `hmSplitSets` for heatmap row groups; annotation TSV / newick download; tip labels show inline **24×18px species image** (base64 PNG from `SPECIES_IMAGES`) + full species name in italic on hover; optional 3rd-col **group pill badge** (right-aligned, colored = average RGB of member species via `groupColors`); click internal node → `showSpNodeActionPopup()` offers **Collapse** or **Flip** (reverses `children` in `SP_TREE_DATA` via `spNodeById`, calls `recomputeSpeciesOrder()` then redraws all views); root node has flip-only dot |
| `heatmap` | `drawHeatmap()` + `drawCladogram()` | Species × HG/OG matrix; drill-down state machine (`hmViewMode`: class→family→hg→og); z-score RdBu or absolute coloring; column drag-reorder (`hmColOrderOverride`); global search dropdown (navigate + custom OG selection via `hmCustomOGs`/`hmCustomGroups`); "Expand all to OGs" (`hmExpandToOGs`); col-sum bar chart; colorbar legend; **collapsible custom-selection bar** (`hmCustomBarExpanded`, toggle button collapses chip row); **Logos button** (`hmShowSpLogos`) toggles 12px species images in cladogram; **Group HG button** (`hmGroupByHG`) groups columns by parent HG with black separator lines and draggable group headers; always calls `drawCladogram()` last |
| `trees` | `renderSidebar()` + `selectTree()` + `renderTree()` | HG list sidebar (family groups, search, class filter); main tree SVG with pan/zoom; controls panel; `colorMode` includes `"group"` for tip coloring by species group |

### Gene tree viewer — key state & functions

**Global state**: `rootNode` (d3.hierarchy), `currentIndex`/`currentDetail` (loaded from lazy script tag), `useBranchLen`, `colorMode` (species/clade/og/**group**), `hlSet`/`hlQueries` (species highlights), `ogHlSet`/`ogHlQueries` (OG highlights), `collapsedFraction`, `cladeHighlights` Map (_uid→{color,label}), `nodeTriColors` Map

**Tree rendering pipeline**:
1. `selectTree(indexRec)` — fetches lazy JSON, sets `currentDetail`, calls `loadTree()`
2. `loadTree()` — builds `d3.hierarchy`, assigns `_uid` to every node, initializes collapse state
3. `renderTree(relayout)` — calls `annotateOGNodes()`, then `layoutTree()` (assigns `.x`/`.y`), then draws SVG
4. `layoutTree()` — d3 cluster layout (cladogram) or `assignBranchLenPos()` + scale (phylogram)
5. `elbowPath(s,t,mg,r)` — orthogonal elbow `M${sx},${sy}V${ty}H${tx}` for all links
6. `nodeX(d,mg)` — returns pixel x honouring branch-length mode

**Collapse mechanics**: left-click internal node → `showCollapseChoicePopup()` → triangle (`_isOgCol=true`, drawn as `col-tri` polygon) or circle; `_children` holds collapsed subtree. Right-click triangle → `showTriActionPopup()` (expand, focus, rename, color, compare).

**OG features**: `annotateOGNodes()` sets `._og_label` at MRCA of each OG's leaves. `toggleCollapseToOGs()` auto-collapses each OG clade to triangle. `toggleHighlightOGs()` draws colored rects behind OG subtrees.

**Species highlights**: `hlSet` (union Set of species), `hlGroupIndex` Map (sp→group i), per-group color via `hlTagColor(i)`. `focusHighlighted()` collapses non-highlighted branches to MRCA-triangle or circle per `focusCollapseAsTri`.

**POSSVM panel** (`#possvm-panel`): full in-browser orthogroup re-calling. `computePossvmAssignments(treeRoot, ingroupSps, sos)` implements species-overlap LPA: `collect()` gathers leaf sets bottom-up, `classify()` labels D/S events, then label propagation assigns OG ids. `pvmMidpointRoot()` — two-sweep diameter algorithm. `runHierPossvm()` — runs POSSVM for each selected clade node, builds `hogModel`, calls `renderHogMap()`.

**hOG map** (`renderHogMap(model)`): D3 SVG Sankey/metro. `model = {levels:[{name,ogs:[{id,og,genes,species,size}],species_count}], links:[{source,target,source_og,target_og,weight,source_level,target_level}]}`. `assignHogOrders()` sorts OG rows to minimize link crossings. Links drawn as cubic bezier paths colored by source OG. Hover highlights related OG chain via `relatedMap`.

### Heatmap detail

`drawHeatmap()` dispatches on `hmViewMode`:
- `class`: columns = unique classes (from `FAMILY_DATA`)
- `family`: columns = HGs in `hmActiveClass`
- `hg`: columns = HGs in `hmActiveFamily` (with per-species z-score or absolute cell coloring)
- `og`: columns = OGs in `hmActiveHG` (loaded from lazy tree data)
- custom: `getEffectiveCustomOGs()` union of `hmCustomOGs` + `hmCustomGroups[].ogs`

Cells are colored by `colorScale(z)` (RdBu) or `absoluteScale(count)`. Click cell → `selectTree()` + highlights that species. Shift+click column header (at HG level) → adds OG expansion. `drawCladogram()` draws species tree left panel synchronized with heatmap rows.

**Heatmap controls** (toolbar above heatmap):
- Label size / Rotation sliders — font size and column label rotation angle (default 90°)
- Colour dropdown — z-score (RdBu) or absolute count
- Sort cols dropdown — rank columns by a species' count
- PNG / SVG download buttons
- **Logos** button (`hmToggleSpLogos`) — toggles 12px species images in the left cladogram panel (`hmShowSpLogos` bool)
- **Group HG** button (`hmToggleGroupByHG`) — groups columns by parent entity: OGs→HG in `custom`/`og` modes, HGs→family in `hg` mode, families→class in `family` mode; highlights blue when active; resets `hmGroupOrderOverride` on toggle

**Group HG grouping internals**:
- `hmGroupByHG` (bool), `hmGroupOrderOverride` (array of group keys) — state vars
- `_colGrpKey(d)` — computed per drawHeatmap call; returns group key for each column record
- After column sort/filter, columns are re-sorted by group rank (stable, preserving `hmColOrderOverride` within group); `hmGrpMeta` array `[{key, startIdx, count}]` and `hmGrpBoundaries` Set track layout
- Group header band (20px, `HM_GROUP_H`) sits between column label bases and cell rows; contains a draggable `rect` per group labeled with the group key and column count
- Dragging a group header updates `hmGroupOrderOverride` and redraws; individual column drag (`hmColOrderOverride`) still works within groups
- Black vertical separator lines (`stroke #1a1a1a`, width 1.5) span from group header top through cells and bar chart

**Families tab accordion internals**:
- `_famExpandedCats` (Set of category keys), `famToggleCat(cat)` toggle function
- Category rows show aggregate stats: species union (not sum, to avoid double-counting), total HG/tree/gene counts across all families in the category
- Family rows are shown/hidden via `display:none` based on `_famExpandedCats`; auto-expand all when filter is non-empty

### Patterns to follow when modifying

- Tree node dict shape: `{name, dist, leaf?, gene_id?, species?, og?, ref?, children?, support?}`
- Collapsed node: `node._children` holds subtree; `node.children=null`; `node._isOgCol` = triangle vs circle
- `_uid` on every d3 node: assigned once in `loadTree()`, stable across re-renders, used as Map key
- OG color: `ogBaseColor(ogName, fallbackIndex)` checks `ogName2Color` then `palette[i%15]`
- Leaf color: `leafColor(sp)` respects `hlSet`, `colorMode` (`"species"` → `spColor`, `"group"` → `groupColor(sp)` via module-scope `groupColors` map, `"clade"`/`"og"` → `cladeSp2Color`)
- Label text: `leafLabelText(d)` = concatenation of gene_id · og · ref per toggle flags; `isLeafLabelVisible(d)` gates `hideNonHl`
- Tooltip: `showTip(event, htmlOrNode)` / `moveTip(event)` / `hideTip()` — shared `#tooltip` div
- Download: `downloadTreeSVG()` / `downloadTreePNG()` / `downloadHeatmapPNG()` / `downloadHeatmapSVG()`

## Coding conventions / preferences
- No permission prompts before editing — just proceed
- Always regenerate HTML after any script edit (see command above)
- JS style: concise, inline functions, D3 v7 patterns
- Follow `report_step2.py` patterns for tree layout and collapse (see `drawSpeciesTree`, `drawCladogram`)
- No species colors in the OG popup tree (only present=black, absent=grey, loss=red, MRCA=blue)
- Do not add docstrings or comments unless logic is non-obvious
- The `HTML_TEMPLATE` is a Python raw string in `visualize_hog_hierarchy.py`; JS/CSS is embedded there

## Pipeline notes
- `docs/manual.md` is useful mainly for intended user workflow, profile usage, and major outputs; it may be outdated in specifics
- `ids.txt` remains the main selected-HG handoff into step 2; `ancestry_ids.txt` is the focused subset for step 4/report work
- `resources.tsv` is an optional step-2 optimization input, not part of the biological contract
- `gather_annotations.py` is part of the intended pipeline functionality and should keep working against current POSSVM/search outputs
- Step 4 currently hangs off the Snakemake layout in `step4_ancestry.smk`, especially the file naming expected by `build_hog_report.py`
- Preferred step-4 output layout:
  - `results/ancestry/{node}/{node}.in_species.txt`
  - `results/ancestry/{node}/{node}.ignore_species.txt`
  - `results/ancestry/{node}/{node}.pruned.tree`
  - `results/ancestry/{node}/possvm/{hg}.ortholog_groups.csv`
  - `results/ancestry/hog_hierarchy.html`
- POSSVM is called with `-method lpa` (label propagation); support values do not affect OG inference
- GeneRax trees (`{hg}.generax.tree`) are fed into POSSVM; they carry posterior probability support values (0–1 scale) at internal nodes, not stripped
- Species tree node names must match exactly; "Filosoa" does not exist — use "Choanozoa" or other named nodes
