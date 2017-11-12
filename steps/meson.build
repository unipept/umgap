steps = [
  'preprocess.sh',
  'gene_prediction.sh',
  'prot2pept.sh',
  'pept2lca.sh',
  'taxa2agg.sh',
  'agg2freq.sh',
]

foreach step: steps
  configure_file(
    input: step + '.in',
    output: step,
    configuration: conf,
  )
endforeach
