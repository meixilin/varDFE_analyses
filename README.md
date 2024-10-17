# varDFE_analyses

Last updated: 2024-10-17

Contact: Meixi Lin

Scripts for "The distribution of fitness effects varies phylogenetically across animals" by Lin et al. 2024

## File structure

```
.
├── README.md
├── step1_preparation
│   ├── customVCFfilter.py
│   ├── direct1DSFS.py
│   ├── easySFS.py
│   ├── step1_annotate_filter_vcf.sh
│   ├── step2_PCA_kinship_filter_samples.R
│   ├── step3_extract_synmis_vcf.sh
│   └── step4_gen_synmis_sfs.sh
├── step2_inference
│   ├── DFE1D_gridsearch_final.sh
│   ├── DFE1D_inferenceFIM_final.sh
│   ├── DFE1D_refspectra_final.sh
│   └── Demog1D_sizechangeFIM_final.sh
└── step3_plotting
    ├── FigSX_AICll_dfefunc.R
    ├── FigSX_DFE_props.R
    ├── FigSX_MIS_dfefunc.R
    ├── FigSX_PGLS.R
    ├── FigSX_PGLSllsurface.R
    ├── FigSX_SYN_MIS-SFS.R
    ├── FigSX_SYN_demog.R
    ├── FigSX_grid_scaled.R
    ├── FigSX_grid_unscaled.R
    ├── FigSX_lourenco_msig.R
    ├── FigSX_lourenco_params.R
    ├── Fig_GammaDFE.R
    ├── Fig_Gridsearch.R
    ├── Fig_Robustness.R
    ├── Fig_TraitCor.R
    ├── Fig_lourenco_DFE.R
    ├── Fig_theory.R
    ├── dfe_plot_util.R
    └── make_legend.R
```

