ARG BASE_CONTAINER=jupyter/scipy-notebook
ARG NB_USER=jovyan
ARG NB_UID=1000

FROM $BASE_CONTAINER
LABEL maintainer="Kevin Sawade <kevin.sawade@uni-konstanz.de>"

USER root

COPY binder/apt.txt /tmp/apt.txt
RUN apt-get update
RUN xargs -a /tmp/apt.txt sudo apt-get install -y


RUN pip install --no-cache-dir notebook
RUN pip install --no-cache-dir jupyterhub
COPY requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt
RUN python -m ipykernel install --sys-prefix
RUN jupyter contrib nbextension install --sys-prefix
RUN jupyter nbextension enable --py --sys-prefix widgetsnbextension
RUN jupyter nbextension enable --py --sys-prefix nglview
RUN jupyter nbextension enable scratchpad/main
RUN jupyter nbextension enable rubberband/main
RUN jupyter nbextension enable exercise2/main
RUN jupyter nbextension enable toggle_all_line_numbers/main
RUN jupyter nbextension enable execute_time/ExecuteTime
RUN jupyter nbextension enable latex_envs/latex_envs

ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

# RUN adduser --disabled-password \
#     --gecos "Default user" \
#     --uid ${NB_UID} \
#     ${NB_USER}

COPY . ${HOME}
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

ENTRYPOINT jupyter-notebook
