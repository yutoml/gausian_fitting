FROM quay.io/jupyter/base-notebook:latest
USER root
RUN apt-get update --yes && apt-get install --yes git build-essential &&\
    python -m pip install --upgrade pip

USER ${NB_UID}
RUN git clone https://github.com/yutoml/MaximaFinder.git &&\
    git clone https://github.com/yutoml/stemtool.git 
# "ipywidgets" is not necessary for stemtool or MaximaFinder but need to show the result of matplotlib on notebook. 
RUN pip install Cython numpy ipywidgets && pip install stemtool/ MaximaFinder/ 
COPY --chown=${NB_UID} peak_fitting.ipynb ${HOME}/
COPY --chown=${NB_UID} test_3.tif ${HOME}/