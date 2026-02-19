project('pyhelp', 'fortran',
  version: '0.5.0',
  default_options: ['warning_level=1']
  )

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true
  ).stdout().strip()

inc_np = include_directories(incdir_numpy)

# Build extension
py.extension_module(
  'HELP3O',
  'pyhelp/HELP3O.FOR',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
  subdir: 'pyhelp',
  fortran_args: ['-fdefault-real-8'],
  link_args: ['-static', '-static-libgfortran', '-static-libgcc']
  )

# Install all Python files
install_subdir('pyhelp',
  install_dir: py.get_install_dir(),
  exclude_files: ['*.FOR', '__pycache__', '*.pyc']
  )
