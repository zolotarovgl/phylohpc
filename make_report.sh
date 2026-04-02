cat data/gene_families_searchinfo.csv | grep -w -E 'rbp|neu|chr|tfs' > tmp/genefam_test.csv
python workflow/report_step2.py --results_dir results/ --species_tree data/species_tree.full.newick --family_info tmp/genefam_test.csv  --output report2.html --species_info data/species_info.tsv  --refnames data/Mmus_gene_names.csv  --refsps Mmus --align_dir results/align/
