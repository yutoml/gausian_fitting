FROM quay.io/jupyter/base-notebook:latest
USER root
RUN apt-get update --yes && apt-get install --yes git build-essential &&\
    python -m pip install --upgrade pip

USER ${NB_UID}
# "ipywidgets" is not necessary for stemtool or MaximaFinder but need to show the result of matplotlib on notebook.
RUN git clone --recursive https://github.com/yutoml/gausian_fitting.git &&\
    cd gausian_fitting && pip install Cython numpy ipywidgets && pip install stemtool/ MaximaFinder/ 
WORKDIR /home/jovyan/gausian_fitting