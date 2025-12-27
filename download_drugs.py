#!/usr/bin/env python3
"""
Fixed: Comprehensive Migraine Drug Downloader
Downloads from: ClinicalTrials.gov, ChEMBL, Open Targets, PubChem
"""

import argparse
import os
import requests
import pandas as pd
import json
import time
from datetime import datetime

# ============================================================================
# PART 1: ClinicalTrials.gov (WORKING ✅)
# ============================================================================

def download_clinicaltrials(disease_name):
    """Download drugs for `disease_name` from ClinicalTrials.gov"""
    
    print("\n" + "="*70)
    print("📊 CLINICALTRIALS.GOV")
    print("="*70)
    
    BASE_URL = "https://clinicaltrials.gov/api/v2/studies"
    
    all_drugs = []
    page_token = None
    page_num = 1
    
    while True:
        print(f"📄 Page {page_num}...")
        
        params = {
            "query.cond": disease_name,
            "fields": "NCTId,BriefTitle,OverallStatus,Phase,InterventionName,InterventionType",
            "pageSize": 100
        }
        
        if page_token:
            params["pageToken"] = page_token
        
        try:
            response = requests.get(BASE_URL, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            studies = data.get("studies", [])
            
            for study in studies:
                protocol = study.get("protocolSection", {})
                identification = protocol.get("identificationModule", {})
                status = protocol.get("statusModule", {})
                design = protocol.get("designModule", {})
                arms = protocol.get("armsInterventionsModule", {})
                
                nct_id = identification.get("nctId", "")
                overall_status = status.get("overallStatus", "")
                phases = design.get("phases", [])
                phase = ", ".join(phases) if phases else "N/A"
                
                interventions = arms.get("interventions", [])
                
                for intervention in interventions:
                    drug_name = intervention.get("name", "")
                    intervention_type = intervention.get("type", "")
                    
                    if drug_name and intervention_type in ["DRUG", "BIOLOGICAL"]:
                        all_drugs.append({
                            "drug_name": drug_name,
                            "nct_id": nct_id,
                            "phase": phase,
                            "status": overall_status,
                            "source": "ClinicalTrials.gov"
                        })
            
            page_token = data.get("nextPageToken")
            if not page_token:
                break
                
            page_num += 1
            time.sleep(0.3)
            
        except Exception as e:
            print(f"❌ Error: {e}")
            break
    
    df = pd.DataFrame(all_drugs)
    if not df.empty:
        df = df[~df["drug_name"].str.contains("placebo", case=False, na=False)]
        print(f"✅ Found {len(df)} drugs")
    return df


# ============================================================================
# PART 2: ChEMBL (FIXED)
# ============================================================================

def download_chembl(efo_ids):
    """Download drugs from ChEMBL for EFO IDs in `efo_ids`"""
    
    print("\n" + "="*70)
    print("📊 ChEMBL")
    print("="*70)
    
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
    
    all_drugs = []
    
    # ChEMBL uses EFO IDs for diseases
    # Expect caller to provide `efo_ids` (list); if empty, function will return no results
    migraine_efo_ids = efo_ids

    if not migraine_efo_ids:
        print("⚠️ No EFO IDs provided for ChEMBL — skipping")
        return pd.DataFrame()

    print("🔍 Searching ChEMBL for provided EFO IDs...")

    for efo_id in migraine_efo_ids:
        print(f"\n  Querying EFO: {efo_id}")
        
        # Get drug indications for this disease
        drug_indication_url = f"{BASE_URL}/drug_indication.json"
        params = {
            "efo_id": efo_id,
            "limit": 1000
        }
        
        try:
            response = requests.get(drug_indication_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            indications = data.get("drug_indications", [])
            print(f"    Found {len(indications)} drug indications")
            
            for indication in indications:
                molecule_chembl_id = indication.get("molecule_chembl_id")
                max_phase = indication.get("max_phase_for_ind")
                
                # Get molecule details
                molecule_url = f"{BASE_URL}/molecule/{molecule_chembl_id}.json"
                
                try:
                    mol_response = requests.get(molecule_url, timeout=30)
                    mol_response.raise_for_status()
                    mol_data = mol_response.json()
                    
                    all_drugs.append({
                        "drug_name": mol_data.get("pref_name", molecule_chembl_id),
                        "chembl_id": molecule_chembl_id,
                        "max_phase": max_phase,
                        "efo_id": efo_id,
                        "molecule_type": mol_data.get("molecule_type"),
                        "molecular_formula": mol_data.get("molecule_properties", {}).get("full_molformula"),
                        "source": "ChEMBL"
                    })
                    
                    time.sleep(0.2)
                    
                except Exception as e:
                    print(f"      ⚠️ Error fetching {molecule_chembl_id}: {e}")
            
        except Exception as e:
            print(f"    ⚠️ Error: {e}")
        
        time.sleep(0.5)
    
    df = pd.DataFrame(all_drugs)
    if not df.empty:
        print(f"\n✅ Total ChEMBL drugs: {len(df)}")
    else:
        print("\n⚠️ No drugs found in ChEMBL")
    return df


# ============================================================================
# PART 3: Open Targets (FIXED)
# ============================================================================

def download_opentargets(disease_ids):
    """Download drugs from Open Targets for disease IDs in `disease_ids`"""
    
    print("\n" + "="*70)
    print("📊 OPEN TARGETS")
    print("="*70)
    
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    
    # Expect caller to provide `disease_ids` (list); if empty, function will return no results
    migraine_diseases = disease_ids

    if not migraine_diseases:
        print("⚠️ No disease IDs provided for Open Targets — skipping")
        return pd.DataFrame()
    
    all_drugs = []
    
    for disease_id in migraine_diseases:
        print(f"\n🔍 Querying {disease_id}...")
        
        # Fixed GraphQL query
        query = """
        query KnownDrugsQuery($diseaseId: String!) {
          disease(efoId: $diseaseId) {
            id
            name
            knownDrugs {
              count
              rows {
                drug {
                  id
                  name
                  drugType
                  maximumClinicalTrialPhase
                  hasBeenWithdrawn
                }
                phase
                status
                urls {
                  name
                  url
                }
              }
            }
          }
        }
        """
        
        variables = {"diseaseId": disease_id}
        
        try:
            response = requests.post(
                url,
                json={"query": query, "variables": variables},
                headers={"Content-Type": "application/json"},
                timeout=30
            )
            response.raise_for_status()
            data = response.json()
            
            # Check for errors
            if "errors" in data:
                print(f"  ⚠️ GraphQL errors: {data['errors']}")
                continue
            
            disease_data = data.get("data", {}).get("disease")
            
            if not disease_data:
                print(f"  ⚠️ No data for {disease_id}")
                continue
            
            disease_name = disease_data.get("name", disease_id)
            known_drugs = disease_data.get("knownDrugs", {})
            rows = known_drugs.get("rows", [])
            
            print(f"  ✓ Found {len(rows)} drugs for {disease_name}")
            
            for item in rows:
                drug_info = item.get("drug", {})
                
                all_drugs.append({
                    "drug_name": drug_info.get("name"),
                    "drug_id": drug_info.get("id"),
                    "drug_type": drug_info.get("drugType"),
                    "max_phase": drug_info.get("maximumClinicalTrialPhase"),
                    "phase": item.get("phase"),
                    "status": item.get("status"),
                    "withdrawn": drug_info.get("hasBeenWithdrawn", False),
                    "disease": disease_name,
                    "disease_id": disease_id,
                    "source": "Open Targets"
                })
            
        except Exception as e:
            print(f"  ❌ Error: {e}")
        
        time.sleep(0.5)
    
    df = pd.DataFrame(all_drugs)
    if not df.empty:
        print(f"\n✅ Total Open Targets drugs: {len(df)}")
    else:
        print("\n⚠️ No drugs found in Open Targets")
    return df


# ============================================================================
# PART 4: PubChem (Added as alternative)
# ============================================================================

def download_pubchem(search_terms):
    """Download compounds from PubChem for `search_terms` list"""
    
    print("\n" + "="*70)
    print("📊 PubChem")
    print("="*70)
    
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    # Expect caller to provide `search_terms` (list of strings); if empty, function will return no results
    if not search_terms:
        print("⚠️ No search terms provided for PubChem — skipping")
        return pd.DataFrame()

    all_compounds = []
    
    for term in search_terms:
        print(f"🔍 {term}...", end=" ")
        
        try:
            # Search for compound
            search_url = f"{BASE_URL}/compound/name/{term}/cids/JSON"
            response = requests.get(search_url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                cids = data.get("IdentifierList", {}).get("CID", [])
                
                if cids:
                    cid = cids[0]  # Take first result
                    
                    # Get compound details
                    detail_url = f"{BASE_URL}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
                    detail_response = requests.get(detail_url, timeout=10)
                    
                    if detail_response.status_code == 200:
                        compound_data = detail_response.json()
                        properties = compound_data["PropertyTable"]["Properties"][0]
                        
                        all_compounds.append({
                            "drug_name": term,
                            "pubchem_cid": properties.get("CID"),
                            "iupac_name": properties.get("IUPACName", "")[:200],
                            "molecular_formula": properties.get("MolecularFormula"),
                            "molecular_weight": properties.get("MolecularWeight"),
                            "smiles": properties.get("CanonicalSMILES"),
                            "source": "PubChem"
                        })
                        
                        print(f"✓ CID:{cid}")
                
            time.sleep(0.2)
            
        except Exception as e:
            print(f"✗ {e}")
    
    df = pd.DataFrame(all_compounds)
    if not df.empty:
        print(f"\n✅ Total PubChem compounds: {len(df)}")
    return df


def resolve_ontology_ids(name, ontology):
    """Try to resolve ontology IDs (EFO/MONDO) for a disease name using OLS.

    Returns a list like ['EFO_0003013'] or empty list on failure.
    """
    try:
        url = "https://www.ebi.ac.uk/ols/api/search"
        params = {"q": name, "ontology": ontology, "size": 10}
        r = requests.get(url, params=params, timeout=10)
        r.raise_for_status()
        data = r.json()

        docs = []
        # tolerate different response shapes
        if isinstance(data.get('response'), dict):
            docs = data['response'].get('docs', [])
        if not docs:
            docs = data.get('_embedded', {}).get('terms', []) if isinstance(data.get('_embedded'), dict) else []
        if not docs:
            docs = data.get('docs', []) or data.get('terms', []) or []

        ids = []
        for d in docs:
            if not isinstance(d, dict):
                continue
            val = d.get('obo_id') or d.get('short_form') or d.get('iri')
            if not val:
                # sometimes label contains the id embedded
                continue
            val = str(val).split('/')[-1]
            val = val.replace(':', '_')
            if val.upper().startswith(ontology.upper()):
                ids.append(val)
        # deduplicate while preserving order
        seen = set()
        out = []
        for x in ids:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out
    except Exception:
        return []


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    print("\n" + "🎯"*35)
    print("COMPREHENSIVE MIGRAINE DRUG DATABASE DOWNLOADER (FIXED)")
    print("🎯"*35)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    parser = argparse.ArgumentParser(description="Download drug lists from public resources for a disease")
    parser.add_argument('disease_name', help='Human-readable disease name (e.g. "migraine")')
    parser.add_argument('--outdir', default='.', help='Output directory for CSV files')

    args = parser.parse_args()

    disease = args.disease_name
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Resolve EFO / MONDO IDs automatically via OLS
    print(f"\nResolving ontology IDs for '{disease}'...")
    efo_ids = resolve_ontology_ids(disease, 'efo')
    mondo_ids = resolve_ontology_ids(disease, 'mondo')
    disease_ids = (efo_ids or []) + (mondo_ids or [])

    # Download from all sources
    results = {}

    print(f"\n[1/4] Downloading from ClinicalTrials.gov for '{disease}'...")
    results['clinicaltrials'] = download_clinicaltrials(disease)

    print("\n[2/4] Downloading from ChEMBL (using resolved EFO IDs)...")
    results['chembl'] = download_chembl(efo_ids)

    print("\n[3/4] Downloading from Open Targets (using resolved IDs)...")
    results['opentargets'] = download_opentargets(disease_ids)

    # Derive PubChem search terms from Open Targets known drugs (if available)
    pubchem_terms = []
    if not results['opentargets'].empty:
        if 'drug_name' in results['opentargets'].columns:
            pubchem_terms = results['opentargets']['drug_name'].dropna().astype(str).unique().tolist()
    # fallback: use disease name as a single search term
    if not pubchem_terms:
        pubchem_terms = [disease]

    print("\n[4/4] Downloading from PubChem (derived search terms)...")
    results['pubchem'] = download_pubchem(pubchem_terms)
    
    # Save individual files
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    saved_files = []
    
    print("\n" + "="*70)
    print("💾 SAVING FILES")
    print("="*70)
    
    # Make filenames disease-aware
    safe_name = ''.join(c if c.isalnum() else '_' for c in disease).lower()
    for source, df in results.items():
        if not df.empty:
            filename = os.path.join(outdir, f"{safe_name}_drugs_{source}_{timestamp}.csv")
            df.to_csv(filename, index=False)
            print(f"✅ {filename} ({len(df)} records)")
            saved_files.append(filename)
    
    # Create combined dataset
    all_dfs = [df for df in results.values() if not df.empty]
    
    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        combined_file = os.path.join(outdir, f"{safe_name}_drugs_ALL_{timestamp}.csv")
        combined_df.to_csv(combined_file, index=False)
        saved_files.append(combined_file)
        print(f"\n✅ COMBINED: {combined_file} ({len(combined_df)} total records)")
    
    # Summary
    print("\n" + "="*70)
    print("📊 SUMMARY")
    print("="*70)
    for source, df in results.items():
        print(f"{source:20s}: {len(df):4d} drugs")
    
    if all_dfs:
        print(f"{'TOTAL':20s}: {len(combined_df):4d} records")
        
        # Unique drugs
        if 'drug_name' in combined_df.columns:
            unique_drugs = combined_df['drug_name'].nunique()
            print(f"{'Unique drug names':20s}: {unique_drugs:4d}")
    
    print("\n✅ DOWNLOAD COMPLETE!")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\n📁 Files saved: {len(saved_files)}")


if __name__ == "__main__":
    try:
        import requests
        import pandas
    except ImportError:
        print("❌ Missing required packages!")
        print("Install with: pip install requests pandas")
        exit(1)
    
    main()
