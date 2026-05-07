#!bin/bash

# Join all ORSUM results into a single summary file

python3 code/orsum/orsum.py \
  --gmt data/hsapiens.GO.BP.name.gmt \
  --fileAliases Up_acute Up_intermediate Down_intermediate Up_prolonged Up_acute_intermediate Up_intermediate_prolonged Down_intermediate_prolonged \
  --files "results/03_enrichement/control vs acute/orsum.go.bp.up.results.csv" \
           "results/03_enrichement/control vs intermediate/orsum.go.bp.up.results.csv" \
           "results/03_enrichement/control vs intermediate/orsum.go.bp.down.results.csv" \
            "results/03_enrichement/control vs prolonged/orsum.go.bp.up.results.csv" \
            "results/03_enrichement/acute vs intermediate/orsum.go.bp.up.results.csv" \
            "results/03_enrichement/intermediate vs prolonged/orsum.go.bp.up.results.csv" \
            "results/03_enrichement/intermediate vs prolonged/orsum.go.bp.down.results.csv" \
  --outputFolder results/03_enrichement/Summary_All_Contrasts \
  --minTermSize 5 
