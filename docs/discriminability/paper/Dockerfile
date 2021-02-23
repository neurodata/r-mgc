FROM rocker/tidyverse:latest

RUN R -e "install.packages(c('MASS', 'mvtnorm', 'parallel', 'ICC', 'abind', 'emdbook', 'gtools', 'kernelPSI', 'energy', 'reticulate', 'mltools', 'data.table', 'kernlab'))"
RUN R -e "require(devtools); install_github('muschellij2/I2C2')"
RUN R -e "require(devtools); install_github('neurodata/r-mgc'); install_github('neurodata/lol')"
RUN R -e "install.packages(c('igraph', 'reshape2', 'stringr', 'FNN', 'Metrics', 'randomForest', 'rARPACK'))"
RUN apt-get update && apt-get install -y \
    curl
CMD /bin/bash
RUN sudo apt-get install -y python3.8-venv
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python3.8 get-pip.py
RUN pip install -U pip
RUN pip install virtualenv
RUN mkdir $HOME/.virtualenvs
RUN python3.8 -m venv $HOME/.virtualenvs/hyppo
RUN git clone https://github.com/neurodata/hyppo.git /hyppo
RUN cd /hyppo; $HOME/.virtualenvs/hyppo/bin/python setup.py install
RUN sudo apt-get install -y vim libglpk-dev
