project('fmm3dpy', 'c', 'fortran',
  version: '2.1.0',
  license: 'Apache-2.0',
  default_options: [
    'optimization=3',
    'fortran_std=legacy'
  ]
)

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy)
inc_f2py = include_directories(incdir_f2py)

f2py_lib = static_library('_f2py',
  incdir_f2py / 'fortranobject.c',
  dependencies: py_dep,
  include_directories: [inc_np, inc_f2py],
)
f2py_dep = declare_dependency(
  include_directories: [inc_np, inc_f2py],
  link_with: f2py_lib,
  dependencies: py_dep
)
f2py_prog = find_program('f2py', 'f2py3', required: true)

fmm_src = [
  'src/Helmholtz/h3dcommon.f',
  'src/Helmholtz/hnumfour.f',
  'src/Helmholtz/hfmm3dwrap.f',
  'src/Helmholtz/hwts3e.f',
  'src/Helmholtz/hfmm3d_ndiv.f',
  'src/Helmholtz/h3dtrans.f',
  'src/Helmholtz/hpwrouts.f',
  'src/Helmholtz/helmrouts3d.f',
  'src/Helmholtz/h3dterms.f',
  'src/Helmholtz/projections.f',
  'src/Helmholtz/hfmm3dwrap_vec.f',
  'src/Helmholtz/hfmm3d.f',
  'src/Helmholtz/hndiv.f',
  'src/Helmholtz/helmkernels.f',
  'src/Helmholtz/hnumphys.f',
  'src/Laplace/lfmm3d_ndiv.f',
  'src/Laplace/l3dterms.f',
  'src/Laplace/l3dtrans.f',
  'src/Laplace/lpwrouts.f',
  'src/Laplace/lapkernels.f',
  'src/Laplace/lwtsexp_sep2.f',
  'src/Laplace/lfmm3d.f',
  'src/Laplace/lwtsexp_sep1.f',
  'src/Laplace/laprouts3d.f',
  'src/Laplace/lfmm3dwrap.f',
  'src/Laplace/lfmm3dwrap_vec.f',
  'src/Laplace/lndiv.f',
  'src/Stokes/stokkernels.f',
  'src/Stokes/stfmm3d.f',
  'src/Common/besseljs3d.f',
  'src/Common/prini.f',
  'src/Common/pts_tree3d.f',
  'src/Common/dfft.f',
  'src/Common/cdjseval3d.f',
  'src/Common/legeexps.f',
  'src/Common/yrecursion.f',
  'src/Common/rotproj.f',
  'src/Common/tree_routs3d.f',
  'src/Common/dlaran.f',
  'src/Common/cumsum.f',
  'src/Common/rotviarecur.f',
  'src/Common/hkrand.f',
  'src/Common/rotgen.f',
  'src/Common/fmmcommon.f',
  'src/Maxwell/emfmm3d.f90'
]

hfmm_src = [
  'src/Helmholtz/hfmm3dwrap.f',
  'src/Helmholtz/hfmm3dwrap_vec.f',
  'src/Helmholtz/helmkernels.f'
]

lfmm_src = [
  'src/Laplace/lfmm3dwrap.f',
  'src/Laplace/lfmm3dwrap_vec.f',
  'src/Laplace/lapkernels.f'
]

emfmm_src = [
  'src/Helmholtz/hfmm3dwrap.f',
  'src/Helmholtz/hfmm3dwrap_vec.f',
  'src/Helmholtz/helmkernels.f',
  'src/Maxwell/emfmm3d.f90'
]

stfmm_src = [
  'src/Laplace/lfmm3dwrap.f',
  'src/Laplace/lfmm3dwrap_vec.f',
  'src/Laplace/lapkernels.f',
  'src/Stokes/stfmm3d.f',
  'src/Stokes/stokkernels.f'
]

fmm3d_lib = static_library('fmm3d',
  fmm_src,
  fortran_args: ['-w'],
  gnu_symbol_visibility: 'hidden'
)

# These source lists are passed to f2py only for signature extraction.
# They are NOT compiled separately -- all Fortran compilation happens
# through fmm_src -> fmm3d_lib above.
modules_to_build = {
  'hfmm': hfmm_src,
  'lfmm': lfmm_src,
  'emfmm': emfmm_src,
  'stfmm': stfmm_src
}

foreach name, sources : modules_to_build
  mod_name = name + '3d_fortran'

  # Generate C and Fortran wrappers with f2py
  gen_pyf = custom_target(
    '@0@_gen'.format(name),
    output: ['@0@module.c'.format(mod_name), '@0@-f2pywrappers.f'.format(mod_name)],
    input: sources,
    command: [f2py_prog, '-m', mod_name, '@INPUT@', '--lower', '--quiet']
  )

  # Compile the Python extension module
  py.extension_module(
    mod_name,
    gen_pyf,
    dependencies: [f2py_dep],
    link_with: [fmm3d_lib],
    c_args: ['-Wno-format'],
    install: true,
    subdir: 'fmm3dpy'
  )
endforeach

py.install_sources(
  ['python/fmm3dpy/__init__.py', 'python/fmm3dpy/fmm3d.py'],
  subdir: 'fmm3dpy'
)
