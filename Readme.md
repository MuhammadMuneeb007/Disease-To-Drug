```markdown
# Drug Database Downloader 💊

Download drugs and compounds for any disease from major biomedical databases.

## Quick Start

```bash
# Install dependencies
pip install requests pandas

# Download drugs for any disease
python download_drugs.py "migraine"
python download_drugs.py "Alzheimer's Disease"
python download_drugs.py "Type 2 Diabetes"
```

## What It Does

Downloads drug information from:
- **ClinicalTrials.gov** - Clinical trials (all phases)
- **ChEMBL** - Chemical/drug database
- **Open Targets** - Disease-drug associations
- **PubChem** - Chemical compounds

## Example Output

```bash
$ python download_drugs.py migraine

🎯 COMPREHENSIVE DRUG DATABASE DOWNLOADER
Started: 2025-12-28 04:23:18

Resolving ontology IDs for 'migraine'...

[1/4] ClinicalTrials.gov... ✅ Found 1388 drugs
[2/4] ChEMBL... ⚠️ No EFO IDs (skipped)
[3/4] Open Targets... ✅ Found 25 drugs
[4/4] PubChem... ✅ Found 7 compounds

📊 SUMMARY
clinicaltrials      : 1388 drugs
chembl              :    0 drugs
opentargets         :   25 drugs
pubchem             :    7 drugs
TOTAL               : 1420 records
Unique drug names   :  854

✅ Files saved: 4
```

## Output Files

```
migraine_drugs_clinicaltrials_20251228_042403.csv    (1388 records)
migraine_drugs_opentargets_20251228_042403.csv       (25 records)
migraine_drugs_pubchem_20251228_042403.csv           (7 records)
migraine_drugs_ALL_20251228_042403.csv               (1420 total)
```

## Options

```bash
python download_drugs.py <disease_name> [--outdir PATH]

Examples:
  python download_drugs.py "migraine"
  python download_drugs.py "Breast Cancer" --outdir ./results
  python download_drugs.py "Parkinson's Disease" --outdir ./neuro_drugs
```

## Requirements

```
Python 3.7+
requests
pandas
```

## How It Works

1. Takes disease name (e.g., "migraine")
2. Resolves to ontology IDs (EFO/MONDO) automatically
3. Queries 4 databases in parallel
4. Saves individual CSVs + combined dataset

## Notes

- **ClinicalTrials.gov** always works (text search)
- **ChEMBL/Open Targets** need ontology IDs (may skip for rare diseases)
- **PubChem** searches drugs found in Open Targets
- All APIs are free and public (no authentication needed)

## License

MIT

## Author

Muhammad Muneeb ([@muhammadmuneeb007](https://github.com/muhammadmuneeb007))

## Data Sources

- [ClinicalTrials.gov](https://clinicaltrials.gov/) - U.S. National Library of Medicine
- [ChEMBL](https://www.ebi.ac.uk/chembl/) - EMBL-EBI
- [Open Targets](https://platform.opentargets.org/) - EMBL-EBI & Partners
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - NCBI
```

Simple, clean, and shows the exact output! 🚀
