import setuptools
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

pkg_name = "fmm3dpy"

## TODO: this should be automatically populated using "read directory, or whatever"
## TODO: fix problem with relative location for executable
list1=['h3dcommon.f','helmrouts3d.f','hfmm3dpartwrap_vec.f','quadread.f',
'h3dterms.f','hfmm3dpart.f','hpwrouts.f','h3dtrans.f','hfmm3dpartwrap.f','hwts3.f','projections.f',
'numphysfour.f']
list2=['besseljs3d.f','d3hplratree.f','dlaran.f','hkrand.f','prini.f','rotgen.f','rotviarecur.f',
'cdjseval3d.f','dfft.f','fmmcommon.f','legeexps.f','rotgen2.f','rotproj.f','yrecursion.f']

ext1 = Extension(
    name='fmm3d_fortran',
    sources=['src/Helmholtz/'+item for item in list1]+['src/Common/'+item for item in list2],
    f2py_options=['only:','hfmm3dpartwrap.f',':']
)

## TODO: fill in the info below
setup(
    name=pkg_name,
    version="0.1.0",
    author="Manas Rachh",
    author_email="",
    description="Something.",
    url="",
    packages=setuptools.find_packages(),
    install_requires=[
    ],
    ext_modules=[ext1],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )
)
