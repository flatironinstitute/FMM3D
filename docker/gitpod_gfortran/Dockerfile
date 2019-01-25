FROM gitpod/workspace-full

USER root

RUN apt-get update && apt-get install -y build-essential
RUN apt-get update && apt-get install -y gfortran

USER gitpod

RUN pip install numpy
