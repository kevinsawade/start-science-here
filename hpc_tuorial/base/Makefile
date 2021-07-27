SLURM_VERSION = 19.05.1
IMAGE = kevinsawade/slurmbase

DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

.PHONY: all build clean test

all: build

build:
	docker build -t $(IMAGE):$(SLURM_VERSION) .

clean:
	@[ -z $(docker images -q $(IMAGE):$(SLURM_VERSION)) ] || docker rmi $(IMAGE):$(SLURM_VERSION)
