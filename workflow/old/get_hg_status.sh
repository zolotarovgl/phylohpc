mkdir -p tmp/
source workflow/functions.sh
read_config configs/config.txt
grep -c '>' results/clusters/*.fasta \
  | awk -F ':' '$2>=5{print $1}' \
  | xargs -n1 basename -s .fasta > tmp/all_ids

for FILEPREF in $(cat tmp/all_ids); do  ./workflow/check_status.sh $FILEPREF; done | grep -v '#' | awk '{print $1"\t"$2$3$4$5}'
