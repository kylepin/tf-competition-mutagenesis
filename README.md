The best way to run this code is in a conda environment.
```
conda create --name tf_comp_analysis_env python==3.10
conda activate tf_comp_analysis_env
pip install -r requirements.txt
```

Then the several scripts can be run (in order) from the commandline via:
`bash get_genomic_site_matches.sh`, `python process_mutation_counts.py` and `python process_mutation_counts_w_sampling.py`.
`jupyter notebook` (also from the commandline) starts a jupyter notebook server, which can then be used to run the methylation_investigation.ipynb file.
