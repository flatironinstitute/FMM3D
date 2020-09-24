import setuptools
import string
import os
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from sys import platform

pkg_name = "fmm3dpy"

## TODO: this should be automatically populated using "read directory, or whatever"
## TODO: fix problem with relative location for executable

list_helm=['hfmm3dwrap.f','hfmm3dwrap_vec.f','helmkernels.f']
list_lap=['lfmm3dwrap.f','lfmm3dwrap_vec.f','lapkernels.f']
list_common=[]

FLIBS = os.getenv('FMM_FLIBS')
FLIBS = FLIBS.rstrip().split(' ')
FLIBS = list(filter(None, FLIBS))
FLIBS.append('../lib-static/libfmm3d.a')

c_opts = ['_c','_d','_cd']
c_opts2 = ['c','d','cd']
st_opts = ['_s','_t','_st']
p_optsh = ['_p','_g']
p_optsh2 = ['p','g']

p_optsl = ['_p','_g','_h']
p_optsl2 = ['p','g','h']

list_int_helm = []
list_int_helm_vec = []
list_int_helm_dir = []

list_int_lap = []
list_int_lap_vec = []
list_int_lap_dir = []

for st in st_opts:
    for cd in c_opts:
        for pg in p_optsh:
            list_int_helm.append('hfmm3d'+st+cd+pg)
            list_int_helm_vec.append('hfmm3d'+st+cd+pg+'_vec')
        for pg in p_optsl:
            list_int_lap.append('lfmm3d'+st+cd+pg)
            list_int_lap_vec.append('lfmm3d'+st+cd+pg+'_vec')

for cd in c_opts2:
    for pg in p_optsh2:
        list_int_helm_dir.append('h3ddirect'+cd+pg)
    for pg in p_optsl2:
        list_int_lap_dir.append('l3ddirect'+cd+pg)

ext_helm = Extension(
    name='fmm3dpy.hfmm3d_fortran',
    sources=['../src/Helmholtz/'+item for item in list_helm]+['../src/Common/'+item for item in list_common],
    f2py_options=['only:']+list_int_helm+list_int_helm_vec+list_int_helm_dir+[':'],
    extra_link_args=FLIBS
)

ext_lap = Extension(
    name='fmm3dpy.lfmm3d_fortran',
    sources=['../src/Laplace/'+item for item in list_lap]+['../src/Common/'+item for item in list_common],
    f2py_options=['only:']+list_int_lap+list_int_lap_vec+list_int_lap_dir+[':'],
    extra_link_args=FLIBS
)

## TODO: fill in the info below
setup(
    name=pkg_name,
    python_requires='>=3.0.0',
    version="0.0.5",
    author="Zydrunas Gimbutas, Leslie Greengard, Libin Lu, Jeremy Magland, Dhairya Malhotra, Michael O'Neil, Manas Rachh, and Vladimir Rokhlin",
    author_email="mrachh@flatironinstitute.org",
    description="This pacakge contains basic routines for Laplace and Helmholtz fast multipole methods in three dimensions",
    long_description=open('../README.md').read(),
    long_description_content_type='text/markdown',
    url="",
    packages=['fmm3dpy'],
    install_requires=[
        "numpy",
        "pytest"
    ],
    ext_modules=[ext_helm,ext_lap],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ]
)
