MC = .mc
BASE_BIN = .mc/bin
PY2 = .mc/envs/py2
export BIN_PATH := $(abspath ${BASE_BIN})
export PATH := ${BIN_PATH}:${PATH}

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CCFLAGS += -D LINUX
	MC_LINK := https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
endif
ifeq ($(UNAME_S),Darwin)
	CCFLAGS += -D OSX
	MC_LINK := https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
endif
UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),x86_64)
	CCFLAGS += -D AMD64
endif
ifneq ($(filter %86,$(UNAME_P)),)
	CCFLAGS += -D IA32
endif
ifneq ($(filter arm%,$(UNAME_P)),)
	CCFLAGS += -D ARM
endif



all: ${BASE_BIN}/python

clean: 
	rm -rf ${MC} mc.sh

.PHONY: all clean
.SECONDARY:

${BASE_BIN}/python:
	mkdir -p raw_data/reads raw_data/reference_sequence analysis_results 
	wget -O - ${MC_LINK} > mc.sh
	bash mc.sh -bf -p ${MC}
	.mc/bin/conda config --system --add channels r --add channels bioconda --add channels conda-forge
	.mc/bin/conda config --system --set always_yes True
	conda install conda-build
	cd scripts && conda build gubbins
	conda install snippy raxml bcbiogff
	conda install --use-local gubbins
	chmod 755 run.sh scripts/*
	rm -fr mc.sh
	cp -a scripts/vcffirstheader ${BASE_BIN}
	sed -i 's~../vcflib/scripts/vcffirstheader~vcffirstheader~g' ${BASE_BIN}/freebayes-parallel
	sed -i 's~../vcflib/bin/vcfstreamsort~vcfstreamsort~g' ${BASE_BIN}/freebayes-parallel

${PY2}/bin/python: ${BASE_BIN}/python
	conda create -n py2 python=2

