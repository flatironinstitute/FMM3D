FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y build-essential

#RUN apt-get update && apt-get install -y gfortran
#RUN ln -s /usr/bin/gfortran-4.8 /usr/bin/gfortran

#RUN apt-get update && apt-get install -y nano git
RUN conda install jupyterlab

ADD src /source/src
ADD examples /source/examples
#RUN mkdir /source/build
#WORKDIR /source/examples
#RUN ulimit -s unlimited && make -f lfmm3dpart.make -j

RUN conda install numpy

WORKDIR /source/src
#RUN rm Common/prinm.f
RUN echo "test" && f2py -c -m lfmm3dpart Laplace/*.f Common/*.f

WORKDIR /source
