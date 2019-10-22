import setuptools
import string
import os
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

pkg_name = "fmm3dpy"

## TODO: this should be automatically populated using "read directory, or whatever"
## TODO: fix problem with relative location for executable

list_helm=['hfmm3dwrap.f','hfmm3dwrap_vec.f']
list_lap=['lfmm3dwrap.f','lfmm3dwrap_vec.f']
list_common=[]

FAST_KER = os.getenv('FAST_KER')
FLIBS = os.getenv('FLIBS')
FFLAGS = os.getenv('FFLAGS')

if(FAST_KER=='ON'):
    list_helm.append('helmkernels_fast.f')
    list_lap.append('lapkernels_fast.f')
else:
    list_helm.append('helmkernels.f')
    list_lap.append('lapkernels.f')


FLIBS = FLIBS.rstrip().split(' ')
FLIBS = list(filter(None, FLIBS))
FLIBS.append('../lib-static/libfmm3d.a')

FFLAGS = FFLAGS.rstrip().split(' ')
FFLAGS = list(filter(None,FFLAGS))


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
    sources=['../src/Helmholtz/'+item for item in list_helm]+['../src/Common/'+item for item in list_common],
    f2py_options=['only:']+list_int_helm+list_int_helm_vec+list_int_helm_dir+[':'],
    extra_f77_compile_args=FFLAGS,
    extra_link_args=FLIBS
)

ext_lap = Extension(
    name='lfmm3d_fortran',
    sources=['../src/Laplace/'+item for item in list_lap]+['../src/Common/'+item for item in list_common],
    f2py_options=['only:']+list_int_lap+list_int_lap_vec+list_int_lap_dir+[':'],
    extra_f77_compile_args=FFLAGS,
    extra_link_args=FLIBS
)

## TODO: fill in the info below
setup(
    name=pkg_name,
    version="0.1.0",
    author="Leslie Greengard, Zydrunas Gimbutas, Libin Lu, Jeremy Magland, and Manas Rachh",
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
