project('UMGAP', [ 'rust' ],
  version: '0.1.0',
)

# Makes all bash-script verbose
debug = false

#─────────────────────────────────────────────────────────────────────────────┐
# DEPENDENCIES                                                                │
#─────────────────────────────────────────────────────────────────────────────┘
# First of all, find a recent enough version of Bash
bash = find_program('bash')
bash_required_version = '4.0.0'
bash_out = run_command(bash, '-c', 'IFS=.;printf "${BASH_VERSINFO[*]: 0:3}"').stdout()
if bash_out.version_compare('< ' + bash_required_version)
  error('Found Bash version @0@, but at least version @1@ is required.'
      .format(bash_out, bash_required_version)
  )
endif

# Check whether we have GNU getopt
getopt = find_program('getopt')
if run_command(getopt, '-T').returncode() != 4
  error('The GNU version of getopt is required.')
endif

# Coreutils
sed = find_program('gsed', 'sed')
awk = find_program('gawk', 'awk')
grep = find_program('ggrep', 'grep')
sort = find_program('gsort', 'sort')
# Check for --parallel support in `sort`
sort_check = run_command(sort, '--parallel=2', 'meson.build')
if sort_check.returncode() == 0
  sort_command = '"@0@" --parallel=8'.format(sort.path())
else
  message('''sort does not support "--parallel".
    This will condiderably slow down the pipeline.
    Please consider installing a recent version of GNU coreutils.'''
  )
  sort_command = sort.path()
endif

# Other
zcat = find_program('zcat')
gunzip = find_program('gunzip')

#─────────────────────────────────────────────────────────────────────────────┐
# UNIPEPT (UMPAP)                                                             │
#─────────────────────────────────────────────────────────────────────────────┘
ruby = find_program('ruby')
unipept = find_program('unipept', required: false)
prot2pept = find_program('prot2pept', required: unipept.found())
peptfilter = find_program('peptfilter', required: unipept.found())
if not unipept.found()
  gem = find_program('gem', required: false)
  if gem.found()
    error('Please install unipept using `gem install unipept`')
  else
    error('Please install the Ruby "gem" program')
  endif
endif


#─────────────────────────────────────────────────────────────────────────────┐
# UMGAP RUST TOOLS                                                            │
#─────────────────────────────────────────────────────────────────────────────┘
cargo = find_program('cargo')
umgap_bin = custom_target('umgap',
  input: join_paths('src', 'main.rs'),
  output: 'umgap',
  command: [ cargo,
    'install',
    '--force',
    '--root=@0@'.format(meson.build_root()),
    '--path=@0@'.format(meson.source_root()),
  ],
  build_by_default: true,
  install: true,
  install_dir: join_paths(get_option('prefix'), get_option('bindir')),
)
umgap_path = join_paths(meson.build_root(), 'bin', 'umgap')

#─────────────────────────────────────────────────────────────────────────────┐
# GENE PREDICTORS                                                             │
#─────────────────────────────────────────────────────────────────────────────┘
FGS = find_program('FGS', required: false)
FGSpp = find_program('FGS++', required: false)
if FGSpp.found()
  FGSpp_path = FGSpp.path()
else
  FGSpp_subproj = subproject('FGSpp')
  FGSpp_bin = FGSpp_subproj.get_variable('FGSpp')
  FGSpp_path = FGSpp_bin.full_path()
endif
maven = find_program('mvn')


#─────────────────────────────────────────────────────────────────────────────┐
# UMGAP PIPELINE                                                              │
#─────────────────────────────────────────────────────────────────────────────┘
# Configuration of the shell scripts
conf = configuration_data()

conf.set('CONFIG_DIR', join_paths(meson.source_root(), 'config'))

if debug
  conf.set('HASHBANG', '#!@0@ -x'.format(bash.path()))
else
  conf.set('HASHBANG', '#!@0@'.format(bash.path()))
endif
conf.set_quoted('SED', '@0@'.format(sed.path()))
conf.set_quoted('AWK', '@0@'.format(awk.path()))
conf.set_quoted('GREP', '@0@'.format(grep.path()))
conf.set('SORT', sort_command)
conf.set_quoted('UNIPEPT_GEM', '@0@'.format(unipept.path()))
conf.set_quoted('PROT2PEPT', '@0@'.format(prot2pept.path()))
conf.set_quoted('PEPTFILTER', '@0@'.format(peptfilter.path()))
conf.set_quoted('UMGAP', '@0@'.format(umgap_path))
conf.set_quoted('GUNZIP', '@0@'.format(gunzip.path()))
conf.set_quoted('ZCAT', '@0@'.format(zcat.path()))
conf.set_quoted('FGS', '@0@'.format(FGS.path()))
conf.set_quoted('FGSPP', '@0@'.format(FGSpp_path))
conf.set('FGSPP_INSTALLED', '@0@'.format(FGSpp.found()))
conf.set_quoted('MAVEN', '@0@'.format(maven.path()))

umgap_pipeline = configure_file(
  input: 'pipeline.sh.in',
  output: 'pipeline.sh',
  configuration: conf,
)


#─────────────────────────────────────────────────────────────────────────────┐
# SUBDIRECTORIES                                                              │
#─────────────────────────────────────────────────────────────────────────────┘
subdir('steps')
# if build_prototypes
#   subdir('prototypes')
# endif
