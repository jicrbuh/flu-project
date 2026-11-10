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
  term = 'Influenza A H1N1 hemagglutinin'
  seasons = ['2025_26']
  n_per_season = 10
  seed = 42
  min_length = 200
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

  aln = align(selected_fasta, Channel.fromPath(params.reference))
  variability_out = variability(aln, selected_meta)
  numbering_out = numbering(Channel.fromPath(params.reference))
  pymol(variability_out.map { it[1] }, numbering_out.map { it[0] })
}

process fetch {
  tag { season }
  publishDir 'data/raw', mode: 'copy'

  input:
  val season

  output:
  tuple path("${season}.faa"), path("${season}_meta.csv")

  script:
  """
  python3 ${projectDir}/scripts/01_fetch.py \
    --season ${season} \
    --term "${params.term}" \
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
  python3 ${projectDir}/scripts/05_pymol.py --drift-csv ${drift_csv} --numbering-csv ${numbering_csv} --output-dir .
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
