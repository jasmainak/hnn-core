FROM jupyter/minimal-notebook:x86_64-2022-10-09

MAINTAINER Jupyter Help <jupyter-help@brown.edu>

USER root

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get -yq dist-upgrade \
    && apt-get install -yq --no-install-recommends \
    openssh-client \
    vim \ 
    curl \
    gcc \
    libopenmpi-dev \
    && apt-get clean

USER $NB_UID

RUN pip install --upgrade pip

USER $NB_UID

RUN pip install --upgrade pip

USER root

RUN apt-get install -yq --no-install-recommends \
    xvfb \
    x11-utils \
    libx11-dev \
    qtbase5-dev \
    qt5-qmake \
    && apt-get clean

ENV DISPLAY=:99

USER $NB_UID

RUN pip install vtk && \
    pip install numpy && \
    pip install scipy && \
    pip install pyqt5 && \
    pip install xvfbwrapper && \
    pip install pyvista && \
    pip install pyvistaqt && \
    pip install nibabel && \
    pip install ipyevents && \
    pip install darkdetect

USER ${NB_UID}

RUN pip install mpi4py

RUN git clone https://github.com/neuronsimulator/nrn && \
	cd nrn && \
	./build.sh && \
	./configure --with-nrnpython=python3 --without-iv --prefix=$HOME && \
	make -j4 && \
	make install -j4 && \
	cd src/nrnpython/ && \
	python3 setup.py install

ENV PATH=${PATH}:${HOME}/x86_64/bin

RUN pip install hnn_core

# Add an x-server to the entrypoint. This is needed by Mayavi
ENTRYPOINT ["tini", "-g", "--", "xvfb-run"]
