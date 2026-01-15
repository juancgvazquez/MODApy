
import pandas as pd
import numpy as np
import pytest
from MODApy.variantsdb import VariantsDB

def test_remove_patient():
    # Create sample data
    data = {
        'CHROM': ['1', '1'],
        'POS': [100, 200],
        'REF': ['A', 'G'],
        'ALT': ['T', 'C'],
        'GENE_NAME': ['Gene1', 'Gene2'],
        'HGVS.C': ['c.1A>T', 'c.1G>C'],
        'HGVS.P': ['p.Asp1Val', 'p.Gly1Ala'],
        'Pat1': ['HET', 'HOM'],
        'Pat2': ['REF', 'HET'],
        'Pat3': ['HOM', 'REF']
    }
    df = pd.DataFrame(data)
    df.set_index(['CHROM', 'POS', 'REF', 'ALT', 'GENE_NAME', 'HGVS.C', 'HGVS.P'], inplace=True)
    
    # Initialize VariantsDB
    vdb = df.pipe(VariantsDB)
    
    # Calculate initial frequencies
    # Mocking that columns with "Pat" are patients
    # The real calcfreqs uses all columns except FREQ/ALLELE_FREQ as patients
    vdb = vdb.calcfreqs()
    
    assert 'Pat1' in vdb.columns
    assert 'Pat2' in vdb.columns
    assert 'Pat3' in vdb.columns
    assert 'FREQ' in vdb.columns
    assert 'ALLELE_FREQ' in vdb.columns
    
    # Helper to check freq
    # Pat1: 0/1 (HET), Pat2: 0/0 (Ref), Pat3: 1/1 (HOM) -> 3 patients. 
    # FREQ: 3/3 = 1.0 (all present/called) ? 
    # calcfreqs logic: self["FREQ"] = self[patients].notnull().sum(axis=1) / len(patients)
    # ALLELE_FREQ: ((HOM*2 + HET*1).sum()) / len(patients*2)
    # Row 1: Pat1(HET), Pat2(Ref), Pat3(HOM). Total alleles = 6. Alt alleles = 1 + 0 + 2 = 3. AF = 3/6 = 0.5
    
    # Remove Pat3
    vdb = vdb.remove_patient('Pat3')
    
    assert 'Pat3' not in vdb.columns
    assert 'Pat1' in vdb.columns
    assert 'Pat2' in vdb.columns
    
    # Check frequencies updated
    # Row 1: Pat1(HET), Pat2(Ref). Total alleles = 4. Alt alleles = 1 + 0 = 1. AF = 1/4 = 0.25
    # Let's check the first row value for ALLELE_FREQ
    
    # Accessing by iloc. First row is index 0.
    # We expect newly calculated AF
    
    expected_af = 0.25
    # The formatting in calcfreqs might leave it as float or something else?
    # It seems to return the dataframe.
    
    actual_af = vdb.iloc[0]['ALLELE_FREQ']
    np.testing.assert_almost_equal(actual_af, expected_af)

    # Test removing non-existent patient (should check logs or just not crash?)
    # The plan said "log error if doesn't", usually we verify it doesn't crash or raise.
    try:
        vdb = vdb.remove_patient('NonExistentPat')
    except SystemExit:
         # The original code tends to exit(1) on errors, let's see if we implement that or just log.
         # Ideally new code shouldn't exit(1), but match existing style.
         pass
