# TODOs:  


Report2: 
- [ ] species tree tab:  allow exporting species annotations by class in the species tree tab   
- [ ] species tree tab:  When hovering over a species label, display the number of annotated genes and some other metadata you find useful 
- [ ] species tree tab:  Make the triangles black outline and grey but also add an option to select the fill color  
- [ ] counts tab: allow full text search via a special bar in the heatmap pane so one can be immediately navigated to the relevant genes  
- [ ] counts tab: The column labels control part occupies too much vertical space - make it horizontal above the page instead  
- [ ] counts tab: The "prefix" selection button should be moved from the right part of the page to the left 
- [ ] counts tab: when clicking the counts heatmap - the whole HG with all the sequences from a particular species is shown. This behaviour shoud happen when clicked on a HG counts heatmap cell. On the OG level - only the sequences from the selected OG and species should be highlighted
- [ ] counts tab - when the clade of species is collapsed into the triangle, keep the node cirle at the tip of the triangle so it's easier to uncollapse it back   
- [ ] counts tab - when the clade of species is collapsed into the triangle, keep the node cirle at the tip of the triangle so it's easier to uncollapse it back   
- [ ] counts tab - when the clade of species is collapsed into the triangle, add the horizontal lines in the heatmap or (even better) a background shading behind the heatmap so the rows from the collapsed clade are visually separated from the rest. This should be possible to do for multiple clades at the same time. 
- [ ] counts tab - add an option to color not by the Z-score but by the absolute counts  
- [ ] counts tab - Move the explainer of the functions from the column labels control dropdown to the top of the page so it doesn't interphere with the heatmap or the cladogram   
- [ ] gene trees tab - selecting the OG to highlight works but it has no effect on which genes are highlighed. Fix this. 
- [ ] gene trees tab - add a button that toggles the node support labels on (useful if they correspond to bootstrap supports)   
- [ ] gene trees tab - when the terminal tip is clicked allow to copy the tip name to the clipboard 
- [ ] gene trees tab - when in the "focus" mode, the collapsed branches occupy too much vertical space - each tip occupies an equal space. Add an option to control the fraction of space occupied by the collapsed tips (from 1 to, say 0.5), so one can adjust this dynamically  
- [ ] gene trees tab - when in the "focus" mode, the collapsed branches occupy too much vertical space - each tip occupies an equal space. Add an option to control the fraction of space occupied by the collapsed tips (from 1 to, say 0.5), so one can adjust this dynamically  
- [ ] gene trees tab - add a small button with a tree symbol that would expand into a species tree with named nodes (only when hovering over them, that would allow to click on the node to be selected for highlight (in addition to be able to text-search the species / clades))
- [ ] gene trees tab - Color: tab - when the specific species / clade has been selected (or multiple), the value of that one should go to the "by species" to not create a confusion   
- [ ] gene trees tab - node collapsing does not work (triangle / circle )


Create the report   
```bash
python workflow/report_step2.py --species_tree data/species_tree.full.newick --family_info genefam.csv --possvm_dir results/possvm --possvm_prev_dir results/possvm_prev/ --output report2.html
```