FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y build-essential

RUN apt-get update && apt-get install -y gfortran
#RUN ln -s /usr/bin/gfortran-4.8 /usr/bin/gfortran

#RUN apt-get update && apt-get install -y nano git
#RUN conda install jupyterlab

RUN conda install numpy

ADD src /source/src
ADD examples /source/examples
#RUN mkdir /source/build
#WORKDIR /source/examples
#RUN ulimit -s unlimited && make -f lfmm3dpart.make -j

WORKDIR /source/src
#RUN rm Common/prinm.f
RUN f2py -c -m lfmm3dpart Helmholtz/*.f Common/*.f

WORKDIR /source
