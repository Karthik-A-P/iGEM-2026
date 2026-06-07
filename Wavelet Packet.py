"""
Compare Hip1 (M. tuberculosis, P9WHR2) against a list of other proteases
using c-p-v encoding + WPT + cross-correlation similarity.
 
Outputs a ranked table sorted by C (max abs cross-correlation across 8 components).
"""
 
import os
import time
import requests
import numpy as np
import pywt
 
# ====================================================================
# CONFIG
# ====================================================================
HIP1_ACCESSION = "P9WHR2"
 
OTHER_PROTEASES = [
    "O04073", "O15393", "O35002", "O43464", "O53896",
    "P00740", "P00741", "P00742", "P00743", "P00750",
    "P04070", "P05156", "P08594", "P08709", "P0C0V0",
    "P24228", "P39099", "P42790", "P54493", "P83110",
    "P98119", "P9WHR3", "Q2HRB6", "Q56VR3", "Q8GB88",
    "Q8RR56", "Q92743", "Q9BYE2", "O06291", "O40922",
    "P42380", "P45161", "P88911", "Q83417", "C3SRW2",
    "O53945", "P42425", "Q9AQD1", "A0A0H2WX20", "Q5ZVV9",
    "Q8GB5",
]
 
CACHE_DIR = "uniprot_cache"
WAVELET = "bior3.3"
LEVEL = 3
PATHS = ['aaa', 'aad', 'ada', 'add', 'daa', 'dad', 'dda', 'ddd']
 
# ====================================================================
# c-p-v encoding (built once)
# ====================================================================
cpvh = {
    'A': (0.00, 8.1, 31.0, 1.8),
    'R': (0.65, 10.5, 124.0, -4.5),
    'N': (1.33, 11.6, 56.0, -3.5),
    'D': (1.38, 13.0, 54.0, -3.5),
    'C': (2.75, 5.5, 55.0, 2.5),
    'Q': (0.89, 10.5, 85.0, -3.5),
    'E': (0.92, 12.3, 83.0, -3.5),
    'G': (0.74, 9.0, 3.0, -0.4),
    'H': (0.58, 10.4, 96.0, -3.2),
    'I': (0.00, 5.2, 111.0, 4.5),
    'L': (0.00, 4.9, 111.0, 3.8),
    'K': (0.33, 11.3, 119.0, -3.9),
    'M': (0.00, 5.7, 105.0, 1.9),
    'F': (0.00, 5.2, 132.0, 2.8),
    'P': (0.39, 8.0, 32.5, -1.6),
    'S': (1.42, 9.2, 32.0, -0.8),
    'T': (0.71, 8.6, 61.0, -0.7),
    'W': (0.13, 5.4, 170.0, -0.9),
    'Y': (0.20, 6.2, 136.0, -1.3),
    'V': (0.00, 5.9, 84.0, 4.2),

}
_aas = list(cpvh.keys())
_arr = np.array([cpvh[a] for a in _aas])
_z = (_arr - _arr.mean(axis=0)) / _arr.std(axis=0)
A_VALUE = {aa: _z[i].sum() for i, aa in enumerate(_aas)}
 
 
# ====================================================================
# UniProt fetching with disk cache
# ====================================================================
def fetch_sequence(accession):
    """Fetch a sequence from UniProt, caching to disk. Returns None on failure."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache_path = os.path.join(CACHE_DIR, f"{accession}.fasta")
 
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            text = f.read()
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        try:
            response = requests.get(url, timeout=30)
            if response.status_code != 200 or not response.text.strip():
                return None
            text = response.text
            with open(cache_path, "w") as f:
                f.write(text)
            time.sleep(0.1)  # be polite
        except requests.RequestException:
            return None
 
    lines = text.strip().split("\n")
    if len(lines) < 2 or not lines[0].startswith(">"):
        return None
    return "".join(lines[1:])
 
 
# ====================================================================
# Pipeline
# ====================================================================
def protein_to_signal(seq):
    return np.array([A_VALUE[aa] for aa in seq if aa in A_VALUE])
 
 
def get_wpt_components(signal):
    wp = pywt.WaveletPacket(data=signal, wavelet=WAVELET, maxlevel=LEVEL, mode='symmetric')
    return {path: wp[path].data for path in PATHS}
 
 
def normalized_cross_correlation(s1, s2):
    """Max absolute value of normalized cross-correlation across all lags."""
    if len(s1) == 0 or len(s2) == 0:
        return 0.0
    raw = np.correlate(s1, s2, mode='full')
    norm = np.sqrt(np.sum(s1**2) * np.sum(s2**2))
    if norm == 0:
        return 0.0
    return float(np.max(np.abs(raw / norm)))
 
 
def similarity_vector(signal1, signal2):
    """Return dict {path: similarity} for all 8 components at level 3."""
    c1 = get_wpt_components(signal1)
    c2 = get_wpt_components(signal2)
    return {p: normalized_cross_correlation(c1[p], c2[p]) for p in PATHS}
 
 
# ====================================================================
# Main
# ====================================================================
def main():
    print(f"Fetching Hip1 ({HIP1_ACCESSION})...")
    hip1_seq = fetch_sequence(HIP1_ACCESSION)
    if hip1_seq is None:
        print(f"  ERROR: could not fetch Hip1. Aborting.")
        return
    print(f"  Length: {len(hip1_seq)} residues")
    hip1_signal = protein_to_signal(hip1_seq)
 
    print(f"\nFetching {len(OTHER_PROTEASES)} other proteases...")
    proteases = {}
    failed = []
    for acc in OTHER_PROTEASES:
        seq = fetch_sequence(acc)
        if seq is None:
            failed.append(acc)
            print(f"  {acc}: FAILED")
        else:
            proteases[acc] = seq
            print(f"  {acc}: {len(seq)} residues")
 
    if failed:
        print(f"\nWarning: {len(failed)} accession(s) failed: {', '.join(failed)}")
 
    print(f"\nComputing similarity vectors against Hip1...")
    results = []
    for acc, seq in proteases.items():
        signal = protein_to_signal(seq)
        sim_vec = similarity_vector(hip1_signal, signal)
        C = max(sim_vec.values())
        best_band = max(sim_vec, key=sim_vec.get)
        results.append({
            'accession': acc,
            'length': len(seq),
            'C': C,
            'best_band': best_band.upper(),
            'vector': sim_vec,
        })
 
    # Sort by C descending
    results.sort(key=lambda r: r['C'], reverse=True)
 
    # ===== Ranked table =====
    print(f"\n{'='*70}")
    print(f"Ranked similarity to Hip1 ({HIP1_ACCESSION}, length {len(hip1_seq)})")
    print(f"{'='*70}")
    print(f"{'Rank':<6}{'Accession':<14}{'Length':>8}{'C':>10}{'Best band':>12}")
    print("-" * 50)
    for i, r in enumerate(results, 1):
        print(f"{i:<6}{r['accession']:<14}{r['length']:>8}{r['C']:>10.4f}{r['best_band']:>12}")
 
    # ===== Full similarity vectors =====
    print(f"\n{'='*100}")
    print("Full similarity vectors (all 8 WPT components)")
    print(f"{'='*100}")
    header = f"{'Accession':<14}" + "".join(f"{p.upper():>9}" for p in PATHS) + f"{'C':>9}"
    print(header)
    print("-" * len(header))
    for r in results:
        row = f"{r['accession']:<14}"
        for p in PATHS:
            row += f"{r['vector'][p]:>9.4f}"
        row += f"{r['C']:>9.4f}"
        print(row)
 
    # ===== Save CSV for downstream use =====
    csv_path = "hip1_similarity.csv"
    with open(csv_path, "w") as f:
        f.write("rank,accession,length,C,best_band," + ",".join(p.upper() for p in PATHS) + "\n")
        for i, r in enumerate(results, 1):
            row = [str(i), r['accession'], str(r['length']), f"{r['C']:.6f}", r['best_band']]
            row.extend(f"{r['vector'][p]:.6f}" for p in PATHS)
            f.write(",".join(row) + "\n")
    print(f"\nFull results saved to {csv_path}")
 
 
if __name__ == "__main__":
    main()