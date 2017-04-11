MC = .mc
BASE_BIN = .mc/bin
PY2 = .mc/envs/py2
export MINICONDA := $(abspath ${MC})
export BIN_PATH := $(abspath ${BASE_BIN})
scripts = scripts
export SCRIPTS := $(abspath ${scripts})
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

all: ${BASE_BIN}/python ${BASE_BIN}/snippy ${BASE_BIN}/run_gubbins.py

mc_only: ${BASE_BIN}/python

after_mc_only: ${BASE_BIN}/snippy ${BASE_BIN}/run_gubbins.py

clean: 
	rm -rf ${MC} mc.sh

.PHONY: all clean mc_only after_mc_only
.SECONDARY:

${BASE_BIN}/python:
	mkdir -p raw_data/reads raw_data/reference_sequence analysis_results logs .temp
	wget -O - ${MC_LINK} > mc.sh
	bash mc.sh -bf -p ${MC}
	.mc/bin/conda config --system --add channels r --add channels bioconda --add channels conda-forge
	.mc/bin/conda config --system --set always_yes True
	conda install -y python=3.5 fasttree raxml biopython bcbiogff pandas sra-tools
	conda create -yn gubbins-env python=3.5 autoconf automake libtool pkg-config nose reportlab dendropy certifi pillow libgcc zlib fasttree raxml biopython=1.65 
	chmod 755 run.sh run-hpc.sh ${SCRIPTS}/*
	rm -fr mc.sh

${BASE_BIN}/snippy: ${BASE_BIN}/python
	conda install snippy
	cp -a ${SCRIPTS}/vcffirstheader ${BASE_BIN}
	sed -i 's~../vcflib/scripts/vcffirstheader~vcffirstheader~g' ${BASE_BIN}/freebayes-parallel
	sed -i 's~../vcflib/bin/vcfstreamsort~vcfstreamsort~g' ${BASE_BIN}/freebayes-parallel

${BASE_BIN}/run_gubbins.py: ${BASE_BIN}/python ${BASE_BIN}/snippy
	git clone https://github.com/sanger-pathogens/gubbins.git
	source activate gubbins-env && cd gubbins && autoreconf -i && ./configure --exec-prefix=${MINICONDA} && make && make install && cd python && python setup.py install
	mv gubbins/python/scripts/gubbins_drawer.py .mc/bin
	rm -rf gubbins

${PY2}/bin/python: ${BASE_BIN}/python
	conda create -n py2 python=2

