nextflow.enable.dsl = 2

def loadEntrezEmail() {
  def envFile = new File("${projectDir}/.env")
  if (envFile.exists()) {
    def match = envFile.readLines().find { it.startsWith('ENTREZ_EMAIL=') }
    if (match) {
      return match.split('=', 2)[1].trim()
    }
  }
  return System.getenv('ENTREZ_EMAIL') ?: 'you@example.com'
}

params {
  email = loadEntrezEmail()
  seasons = ['2021_22', '2022_23', '2023_24', '2024_25', '2025_26']
  n_per_season = 10
  seed = 42
  min_length = 566
  // For future seasons, store the season -> HA structure mapping in a config file or command-line override.
  season_structure_map = [
    '2025_26': '7KNA'
  ]
  pdb_id = '7KNA'
  chain = 'A'
  expected_protein = 'Hemagglutinin'
  expected_organism = 'Influenza A virus'
  reference = 'data/ref/ref_A_Missouri_11_2025.faa'
}

workflow {
  fetch_ch = Channel.fromList(params.seasons)
  fetched = fetch(fetch_ch)

  selected = select(fetched.collect())
  selected_fasta = selected.map { it[0] }
  selected_meta = selected.map { it[1] }

  ref_path = file("${projectDir}/${params.reference}")
  aln = align(selected_fasta, Channel.fromPath(ref_path))
  variability_out = variability(aln, selected_meta)
  numbering_out = numbering(Channel.fromPath(ref_path))
  pymol_out = pymol(variability_out.map { it[1] }, numbering_out.map { it[0] })
  validate(variability_out.map { it[0] }, variability_out.map { it[1] }, numbering_out.map { it[0] }, numbering_out.map { it[1] }, pymol_out)
}

process fetch {
  tag { season }
  maxForks 1
  publishDir 'data/raw', mode: 'copy'

  input:
  val season

  output:
  tuple path("${season}.faa"), path("${season}_meta.csv")

  script:
  """
  python3 ${projectDir}/scripts/01_fetch.py \
    --season ${season} \
    --email "${params.email}" \
    --output-dir . \
    --meta-output ${season}_meta.csv
  """
}

process select {
  publishDir 'data/selected', mode: 'copy'

  input:
  path raw_files

  output:
  tuple path('selected_*.faa'), path('selected_meta.csv')

  script:
  """
  python3 ${projectDir}/scripts/02_select.py \
    --n-per-season ${params.n_per_season} \
    --seed ${params.seed} \
    --min-length ${params.min_length} \
    --input-dir . \
    --output-dir .
  """
}

process align {
  publishDir 'results', mode: 'copy'

  input:
  path selected_fa
  path reference_fa

  output:
  path 'aln.faa'

  script:
  """
  cat ${reference_fa} ${selected_fa} > all.faa
  mafft --auto all.faa > aln.faa
  """
}

process variability {
  publishDir 'results', mode: 'copy'

  input:
  path aln_fa
  path meta_csv

  output:
  tuple path('variability.csv'), path('drift.csv')

  script:
  """
  python3 ${projectDir}/scripts/04_variability.py --alignment ${aln_fa} --meta ${meta_csv} --output-dir .
  """
}

process pymol {
  publishDir 'results', mode: 'copy'

  input:
  path drift_csv
  path numbering_csv

  output:
  path 'drift.pml'

  script:
  """
  python3 ${projectDir}/scripts/05_pymol.py \
    --drift-csv ${drift_csv} \
    --numbering-csv ${numbering_csv} \
    --pdb-id "${params.pdb_id}" \
    --chain "${params.chain}" \
    --output-dir .
  """
}

process validate {
  publishDir 'results', mode: 'copy'

  input:
  path variability_csv
  path drift_csv
  path numbering_map_csv
  path numbering_check_txt
  path pymol_file

  output:
  path 'validation_summary.txt'

  script:
  """
  python3 - <<'PY'
from pathlib import Path
import csv

required = [
    Path("${variability_csv}"),
    Path("${drift_csv}"),
    Path("${numbering_map_csv}"),
    Path("${numbering_check_txt}"),
    Path("${pymol_file}"),
]

missing_or_empty = []
for path in required:
    if not path.exists():
        missing_or_empty.append(f"missing: {path}")
    elif path.stat().st_size == 0:
        missing_or_empty.append(f"empty: {path}")

if missing_or_empty:
    for msg in missing_or_empty:
        print(msg)
    raise SystemExit(1)

for csv_path in [Path("${variability_csv}"), Path("${drift_csv}")]:
    with csv_path.open() as fh:
        rows = list(csv.reader(fh))
    if len(rows) <= 1:
        raise SystemExit(f"{csv_path.name} has no data rows")

numbering_map = Path("${numbering_map_csv}")
with numbering_map.open() as fh:
    reader = csv.DictReader(fh)
    if reader.fieldnames is None:
        raise SystemExit("numbering_map.csv is missing a header")
    required_cols = {"ref_index", "mature_index", "pdb_resi"}
    missing_cols = required_cols - set(reader.fieldnames)
    if missing_cols:
        raise SystemExit(f"numbering_map.csv missing required columns: {sorted(missing_cols)}")

check_path = Path("${numbering_check_txt}")
with check_path.open() as fh:
    check_text = fh.read().strip()
if not check_text:
    raise SystemExit("numbering_check.txt is empty")

summary = "validation_summary.txt"
with Path(summary).open("w") as out:
    out.write("VALIDATION PASSED" + chr(10))
    out.write("variability_csv=" + Path("${variability_csv}").name + chr(10))
    out.write("drift_csv=" + Path("${drift_csv}").name + chr(10))
    out.write("numbering_map_csv=" + Path("${numbering_map_csv}").name + chr(10))
    out.write("pymol_file=" + Path("${pymol_file}").name + chr(10))

print("VALIDATION PASSED")
PY
  """
}

process numbering {
  publishDir 'results', mode: 'copy'

  input:
  path ref_fa

  output:
  tuple path('numbering_map.csv'), path('numbering_check.txt')

  script:
  """
  python3 ${projectDir}/scripts/03_numbering.py \
    --reference ${ref_fa} \
    --pdb-id "${params.pdb_id}" \
    --chain "${params.chain}" \
    --expected-protein "${params.expected_protein}" \
    --expected-organism "${params.expected_organism}" \
    --output-dir .
  """
}
