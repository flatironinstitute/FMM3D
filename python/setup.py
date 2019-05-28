import setuptools
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

pkg_name = "fmm3dpy"

## TODO: this should be automatically populated using "read directory, or whatever"
## TODO: fix problem with relative location for executable
list_helm=['h3dcommon.f','helmrouts3d.f','hfmm3dwrap_vec.f','quadread.f',
'h3dterms.f','hfmm3d.f','hpwrouts.f','h3dtrans.f','hfmm3dwrap.f','hwts3.f','projections.f',
'numphysfour.f']
list_lap=['l3dterms.f','l3dtrans.f','laprouts3d.f','lfmm3d.f','lpwrouts.f','lwtsexp_sep1.f',
'lwtsexp_sep2.f','lfmm3dwrap.f','lfmm3dwrap_vec.f']
list_common=['besseljs3d.f','tree_lr_3d.f','dlaran.f','hkrand.f','prini.f','rotgen.f','rotviarecur.f',
'cdjseval3d.f','dfft.f','fmmcommon.f','legeexps.f','rotproj.f','yrecursion.f']

c_opts = ['_c','_d','_cd']
c_opts2 = ['c','d','cd']
st_opts = ['_s','_t','_st']
p_opts = ['_p','_g']
p_opts2 = ['p','g']

list_int_helm = []
list_int_helm_vec = []
list_int_helm_dir = []

list_int_lap = []
list_int_lap_vec = []
list_int_lap_dir = []

for st in st_opts:
    for cd in c_opts:
        for pg in p_opts:
            list_int_helm.append('hfmm3d'+st+cd+pg)
            list_int_helm_vec.append('hfmm3d'+st+cd+pg+'_vec')
            list_int_lap.append('lfmm3d'+st+cd+pg)
            list_int_lap_vec.append('lfmm3d'+st+cd+pg+'_vec')

for cd in c_opts2:
    for pg in p_opts2:
        list_int_helm_dir.append('h3ddirect'+cd+pg)
        list_int_lap_dir.append('l3ddirect'+cd+pg)

ext_helm = Extension(
    name='hfmm3d_fortran',
    sources=['./src/Helmholtz/'+item for item in list_helm]+['./src/Common/'+item for item in list_common],
    f2py_options=['only:']+list_int_helm+list_int_helm_vec+list_int_helm_dir+[':'],
    extra_f77_compile_args=['-O3', '-W','-fopenmp'],
    extra_link_args=['-O3','-lgomp']
)

ext_lap = Extension(
    name='lfmm3d_fortran',
    sources=['./src/Laplace/'+item for item in list_lap]+['./src/Common/'+item for item in list_common],
    f2py_options=['only:']+list_int_lap+list_int_lap_vec+list_int_lap_dir+[':'],
    extra_f77_compile_args=['-O3', '-W','-fopenmp'],
    extra_link_args=['-O3','-lgomp']
)

## TODO: fill in the info below
setup(
    name=pkg_name,
    version="0.1.0",
    author="Leslie Greengard, Zydrunas Gimbutas, Jeremy Magland, and Manas Rachh",
    author_email="mrachh@flatironinstitute.org",
    description="This pacakge contains basic routines for Laplace and Helmholtz fast multipole methods in three dimensions",
    url="",
    packages=setuptools.find_packages(),
    install_requires=[
        "pytest"
    ],
    ext_modules=[ext_helm,ext_lap],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )    
)
