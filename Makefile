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

all: ${BASE_BIN}/python ${MC}/.condarc ${BASE_BIN}/conda-build ${MC}/conda-bld/linux-64/gubbins-2.2.0-0.tar.bz2 ${BASE_BIN}/snippy ${BASE_BIN}/run_gubbins.py

mc_only: ${BASE_BIN}/python ${MC}/.condarc

after_mc_only: ${BASE_BIN}/snippy ${BASE_BIN}/run_gubbins.py

clean: 
	rm -rf ${MC} mc.sh

.PHONY: all clean mc_only after_mc_only
.SECONDARY:

${BASE_BIN}/python:
	mkdir -p raw_data/reads raw_data/reference_sequence analysis_results logs .temp
	wget -O - ${MC_LINK} > mc.sh
	bash mc.sh -bf -p ${MC}
	conda install -y python=3.5
	chmod 755 run.sh run-hpc.sh scripts/*
	rm -fr mc.sh

${MC}/.condarc: ${BASE_BIN}/python
	.mc/bin/conda config --system --add channels r --add channels bioconda --add channels conda-forge
	.mc/bin/conda config --system --set always_yes True

${BASE_BIN}/conda-build: ${BASE_BIN}/python ${MC}/.condarc ${BASE_BIN}/snippy
	conda install conda-build

${MC}/conda-bld/linux-64/gubbins-2.2.0-0.tar.bz2: ${BASE_BIN}/python ${MC}/.condarc ${BASE_BIN}/conda-build ${BASE_BIN}/snippy
	cd scripts && conda build gubbins

${BASE_BIN}/snippy: ${BASE_BIN}/python ${MC}/.condarc
	conda install biopython=1.65 raxml bcbiogff pandas sra-tools
	conda install snippy
	cp -a scripts/vcffirstheader ${BASE_BIN}
	sed -i 's~../vcflib/scripts/vcffirstheader~vcffirstheader~g' ${BASE_BIN}/freebayes-parallel
	sed -i 's~../vcflib/bin/vcfstreamsort~vcfstreamsort~g' ${BASE_BIN}/freebayes-parallel

${BASE_BIN}/run_gubbins.py: ${BASE_BIN}/python ${MC}/.condarc ${BASE_BIN}/conda-build ${MC}/conda-bld/linux-64/gubbins-2.2.0-0.tar.bz2 ${BASE_BIN}/snippy
	conda install --use-local gubbins

${PY2}/bin/python: ${BASE_BIN}/python
	conda create -n py2 python=2

