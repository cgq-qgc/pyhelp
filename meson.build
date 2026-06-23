project('pyhelp', 'fortran', 'c',
  version: '0.5.0',
  default_options: ['warning_level=1']
  )

add_project_arguments('-Wno-format', language: 'c')
add_project_arguments('-Wno-conversion', language: 'fortran')

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true
  ).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

# Generate the C wrapper using f2py
help3o_c = custom_target('help3o_c',
  input : 'pyhelp/HELP3O.FOR',
  output : ['HELP3Omodule.c', 'HELP3O-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'HELP3O',
             '--lower', '--build-dir', '@OUTDIR@',
             'only:', 'run_simulation', ':'
             ]
)

# Windows uses static libs directly
# For other OS this is done in the CI
link_args = []
if host_machine.system() == 'windows'
  link_args = ['-static-libgfortran', '-static-libgcc']
endif

# Build extension using the generated C wrapper AND the original Fortran file
py.extension_module(
  'HELP3O',
  ['pyhelp/HELP3O.FOR', help3o_c, incdir_f2py + '/fortranobject.c'],
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
  subdir: 'pyhelp',
  link_args: link_args
)

# Install all Python files.
install_subdir('pyhelp',
  install_dir: py.get_install_dir(),
  exclude_files: ['*.FOR', '__pycache__', '*.pyc']
  )
