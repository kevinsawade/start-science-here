FROM ubuntu:20.04

COPY repo/ /repo
RUN useradd -ms /bin/bash jovyan
RUN apt-get update
RUN apt-get install -y curl git vim
RUN apt-get install -y python-setuptools python-dev
USER jovyan
WORKDIR /home/jovyan
RUN mkdir -p git && \
    cd git && \
    cp /repo/master.txt ingredients.txt && \
    git config --global user.email "you@example.com" && \
    git config --global user.name "jovyan" && \
    git config --global alias.graph "log --all --graph --decorate --oneline" && \
    git config --global merge.tool vimdiff && \
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
    git checkout master && \
    cp /repo/result.txt . && \
    cp /repo/check_merge.py .
RUN cd git && \
    git merge like-cilantro && \
    git merge dislike-cilantro || true
WORKDIR /home/jovyan/git
ENTRYPOINT read -p "Make the resulting ingredients.txt have MORE cilantro. Continue with enter": some_var ; git mergetool ; python check_merge.py
