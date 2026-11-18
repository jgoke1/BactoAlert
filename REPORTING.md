# BACTOALERT reporting
Overlay this source on your existing pipeline checkout. Keep the local adapter FASTA, 16S BLAST database, reference data and samplesheets; they were not part of the source upload. Use Nextflow >=24.10.3 with your existing DSL2 settings and a local copy of jgoke1/bactoalert-report:0.1.0 built for the workstation architecture. A Docker Hub push was not verified.
    
Main entry point: main.nf. Add --run_id RUN_001 to the existing command. To add external Klebsiella context, supply --report_reference_dir /path/to/reference. All runs validate the bundled offline QualiBact snapshot; --qualibact_dir can point to an updated snapshot. A 16S-only run has no genome QC verdicts. Reporting is mandatory; --collect_results continues to control the separate raw collector only.

The report-data preparation step builds tables directly from the channels for this run, then BACTOALERT_REPORT renders one report into <outdir>/reports/<run_id>/bactoalert_report.html. Use a distinct run ID for each cohort. An absent general result is missing, not a negative detection. Mash and SKA matrices and Mash trees use the formatted species reports. The current HTML renderer does not display a separate reference-based phylogeny tab.

Changed: main.nf, nextflow.config, nextflow_schema.json, modules/messages.nf. Added: modules/bactoalert_report.nf, bin/prepare_report_input.py, assets/qualibact_20261030T235912Z/ and this guide. The startup and help banners use an ASCII alert beacon. GHRU Assembly remains credited for its upstream assembly workflows.
