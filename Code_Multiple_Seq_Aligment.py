from collections import OrderedDict
import sys
import csv

# --- Input file handling ---
if len(sys.argv) > 1:
    alignment_file = sys.argv[1]  # take filename from command line
else:
    alignment_file = "three_seq.txt"  # default if none provided

default_start = 1  # starting residue number

# --- Read alignment and concatenate blocks ---
sequences = OrderedDict()
try:
    with open(alignment_file, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("CLUSTAL") or line.startswith("*") or set(line) <= set("*. :"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            seq_id, seq_chunk = parts[0], parts[1]
            pdb_id = seq_id.split("_")[0]  # ✅ keep only PDB ID
            if pdb_id not in sequences:
                sequences[pdb_id] = []
            sequences[pdb_id].append(seq_chunk)
except FileNotFoundError:
    print(f"ERROR: alignment file '{alignment_file}' not found.")
    sys.exit(1)

# --- Join sequence blocks into full alignment ---
for sid in sequences:
    sequences[sid] = "".join(sequences[sid])

# --- Ensure equal length (pad with gaps if needed) ---
max_len = max(len(s) for s in sequences.values())
for sid in sequences:
    sequences[sid] = sequences[sid].ljust(max_len, "-")

# --- Initialize residue counters ---
counters = {sid: default_start for sid in sequences}

# --- Build correspondence table ---
seq_ids = list(sequences.keys())
table = []
for i in range(max_len):
    row = []
    for sid in seq_ids:
        aa = sequences[sid][i]
        if aa == "-":
            row.append("-")
        else:
            row.append(f"{counters[sid]}_{aa}")
            counters[sid] += 1
    table.append(row)

# --- Print table to console ---
print("\t".join(seq_ids))
for r in table:
    print("\t".join(r))

# --- Save table to CSV (Excel-ready) ---
csv_filename = "alignment_table.csv"
with open(csv_filename, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(seq_ids)
    writer.writerows(table)

