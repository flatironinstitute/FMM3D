project('fmm3dpy', 'c', 'fortran',
  version: '2.1.0',
  license: 'Apache-2.0',
  default_options: ['warning_level=1']
)

py = import('python').find_installation(pure: false)
py_dep = py.dependency()
np_dep = dependency('numpy', required: true)

incdir_numpy = run_command(py, ['-c', 'import sys; sys.path = [p for p in sys.path if p != "."]; import numpy; print(numpy.get_include())'], check: true).stdout().strip()
inc_np = include_directories(incdir_numpy)
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)

fortranobject_c = incdir_f2py / 'fortranobject.c'
fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'],
  dependencies: py_dep,
  include_directories: [inc_np, inc_f2py],
  gnu_symbol_visibility: 'hidden',
)
fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [inc_np, inc_f2py],
)

fmm_src = [
  'src/Helmholtz/hfmm3dwrap_legacy.f',
  'src/Helmholtz/h3dcommon.f',
  'src/Helmholtz/hnumfour.f',
  'src/Helmholtz/hfmm3dwrap.f',
  'src/Helmholtz/hfmm3d_memest.f',
  'src/Helmholtz/hwts3e.f',
  'src/Helmholtz/hfmm3d_ndiv.f',
  'src/Helmholtz/h3dtrans.f',
  'src/Helmholtz/hpwrouts.f',
  'src/Helmholtz/helmrouts3d.f',
  'src/Helmholtz/h3dterms.f',
  'src/Helmholtz/projections.f',
  'src/Helmholtz/helmkernels_hess.f',
  'src/Helmholtz/hfmm3dwrap_vec.f',
  'src/Helmholtz/hfmm3d.f',
  'src/Helmholtz/hndiv.f',
  'src/Helmholtz/helmkernels_dr.f',
  'src/Helmholtz/helmkernels.f',
  'src/Helmholtz/hfmm3d_mps.f90',
  'src/Helmholtz/hnumphys.f',
  'src/Laplace/lfmm3d_ndiv.f',
  'src/Laplace/lfmm3dwrap_legacy.f',
  'src/Laplace/l3dterms.f',
  'src/Laplace/l3dtrans.f',
  'src/Laplace/lpwrouts.f',
  'src/Laplace/lapkernels_dr.f',
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
  fortran_args: ['-fPIC', '-O3', '-funroll-loops', '-std=legacy', '-w'],
  gnu_symbol_visibility: 'hidden'
)

f2py_prog = find_program('f2py', 'f2py3', required: true)

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
    name + '_genpyf',
    output: [mod_name + 'module.c', mod_name + '-f2pywrappers.f'],
    input: sources,
    command: [f2py_prog, '-m', mod_name, '@INPUT@', '--lower']
  )

  # Compile the Python extension module
  py.extension_module(
    mod_name,
    gen_pyf,
    link_with: [fmm3d_lib],
    dependencies: [np_dep, fortranobject_dep],
    link_language: 'fortran',
    install: true,
    subdir: 'fmm3dpy'
  )
endforeach

py.install_sources([
  'python/fmm3dpy/__init__.py',
  'python/fmm3dpy/fmm3d.py',
], subdir: 'fmm3dpy'
)
