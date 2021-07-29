FROM ubuntu:20.04

COPY repo/ /repo
RUN useradd -ms /bin/bash jovyan
RUN apt-get update
RUN apt-get install -y curl git vim
USER jovyan
WORKDIR /home/jovyan
RUN mkdir -p git && \
    cd git && \
    cp /repo/master.txt ingredients.txt && \
    git config --global user.email "you@example.com" && \
    git config --global user.name "jovyan" && \
    git config --global alias.graph "log --all --graph --decorate --oneline" && \
    git init && \
    git add ingredients.txt && \
    git commit -m "Initial commit" && \
    git checkout -b like-cilantro && \
    cp /repo/like_cilantro.txt ingredients.txt && \
    git add ingredients.txt && \
    git commit -m "please more cilantro" && \
    git checkout master && \
    git checkout -b dislike-cilantro && \
    cp /repo/dislike_cilantro.txt ingredients.txt && \
    git add ingredients.txt && \
    git commit -m "reduce cilantro to 0.5" && \
    git checkout master
