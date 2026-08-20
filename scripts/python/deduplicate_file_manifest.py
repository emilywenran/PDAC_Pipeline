import os

manifest_file = "/Users/emily/thesis/PDAC/01_metadata/file_manifest.tsv"

with open(manifest_file, "r") as f:
    lines = f.read().splitlines()

header = lines[0]
unique_lines = [header]
seen_ids = set()

duplicate_count = 0
for line in lines[1:]:
    if not line.strip():
        continue
    parts = line.split("\t")
    file_id = parts[0]
    if file_id in seen_ids:
        duplicate_count += 1
        continue
    seen_ids.add(file_id)
    unique_lines.append(line)

with open(manifest_file, "w") as f:
    f.write("\n".join(unique_lines) + "\n")

print(f"Deduplicated file manifest: removed {duplicate_count} duplicate lines.")
