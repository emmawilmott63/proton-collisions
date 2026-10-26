PYTHIA = $(HOME)/Desktop/projects/pythia8318
FASTJET = $(HOME)/Desktop/projects/fastjet-install

CXX = c++
CXXFLAGS = -O2 -std=c++17 $(shell root-config --cflags) -I$(PYTHIA)/include -I$(FASTJET)/include
LDFLAGS = $(shell root-config --libs) -L$(PYTHIA)/lib -lpythia8 -L$(FASTJET)/lib -lfastjet
RPATH = -Wl,-rpath,$(PYTHIA)/lib -Wl,-rpath,$(FASTJET)/lib

# Generation programs (share generateCommon.h)
GEN = jets-charm jets-inc

# Analysis programs (share analysisUtils.h)
ANA = sanity mass buildtree features

# TMVA programs (share analysisUtils.h and tmvaCommon.h); need -lTMVA, which
# root-config --libs does not include
TMVA_PROGS = trainTMVA applyTMVA

all: $(GEN) $(ANA) $(TMVA_PROGS) reconstruct-d0

$(GEN): %: macros/%.cc macros/generateCommon.h macros/analysisUtils.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

$(ANA): %: macros/%.cc macros/analysisUtils.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

$(TMVA_PROGS): %: macros/%.cc macros/analysisUtils.h macros/tmvaCommon.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) -lTMVA $(RPATH)

reconstruct-d0: macros/reconstructD0.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

clean:
	rm -f $(GEN) $(ANA) $(TMVA_PROGS) test reconstruct-d0

# ---------------------------------------------------------------------------
# Analysis pipeline:  generate -> candidates -> train -> apply
#
#   make study                      whole smearing study; only redoes what is
#                                   out of date, so it can be interrupted and
#                                   resumed
#   make study LEVELS="0 0.02 0.05" other smearing levels (relative pT)
#   make study METHOD=BDTG          other TMVA method
#   make study TRAIN_OPTS="--depth=5 --ntrees=600 --min-node=0.5"
#   make reapply                    rerun only applyTMVA (keeps the trained BDTs)
#   make clean-study                delete this study's candidates, trained
#                                   BDTs and plots (asks for nothing: be sure)
#   make generate EVENTS=1000000    create $(INPUT) if it does not exist
#
# Rebuilding a program does NOT redo finished results (the programs are
# order-only prerequisites). After changing the code, run `make clean-study`
# first if you want the results regenerated.
# ---------------------------------------------------------------------------
INPUT      ?= data/d0_inc_5_100.root
LEVELS     ?= 0 0.01 0.03 0.05 0.10
METHOD     ?= BDT
EVENTS     ?= 1000000
JOBS       ?= 2
BUILD_OPTS ?=
TRAIN_OPTS ?=
APPLY_OPTS ?=

# 0.03 -> smear0p03
tagof   = smear$(subst .,p,$(1))
TAGS    = $(foreach s,$(LEVELS),$(call tagof,$(s)))
RESULTS = $(foreach t,$(TAGS),plots/tmva_$(METHOD)_$(t)_significance.png)

# Particle file: only the default name has a rule (hours of generation, so
# a rebuilt jets-inc is order-only too). Other INPUT files must exist.
data/d0_inc_5_100.root: | jets-inc
	./scripts/run_generation.sh inc $(EVENTS) $(JOBS)

generate: | jets-inc
	./scripts/run_generation.sh inc $(EVENTS) $(JOBS)

# $(1) = smearing level, $(2) = tag
define LEVEL_RULES
data/candidates_$(2).root: $(INPUT) | buildtree
	@mkdir -p data logs
	./buildtree $(INPUT) --smear-pt=$(1) --tag=$(2) $(BUILD_OPTS) > logs/$(2)_build.log 2>&1

ds_$(2)/weights/TMVAClassification_$(METHOD).weights.xml: data/candidates_$(2).root | trainTMVA
	@mkdir -p logs
	./trainTMVA data/candidates_$(2).root --dataset=ds_$(2) --method=$(METHOD) $(TRAIN_OPTS) > logs/$(2)_train.log 2>&1

plots/tmva_$(METHOD)_$(2)_significance.png: ds_$(2)/weights/TMVAClassification_$(METHOD).weights.xml | applyTMVA
	@mkdir -p plots logs
	./applyTMVA data/candidates_$(2).root --dataset=ds_$(2) --method=$(METHOD) $(APPLY_OPTS) > logs/$(2)_apply.log 2>&1
endef
$(foreach s,$(LEVELS),$(eval $(call LEVEL_RULES,$(s),$(call tagof,$(s)))))

study: $(RESULTS)
	@echo "Summary (data/tmva_summary.csv):"
	@column -s, -t data/tmva_summary.csv

# Redo only the apply step (new cut scan / summary) without retraining
reapply:
	rm -f $(RESULTS)
	$(MAKE) study

clean-study:
	rm -f $(foreach t,$(TAGS),data/candidates_$(t).root data/tmva_training_$(t).root data/bdt_response_$(t).root plots/tmva_$(METHOD)_$(t)_*.png)
	rm -rf $(foreach t,$(TAGS),ds_$(t))

help:
	@sed -n '/^# Analysis pipeline/,/^# -----/p' Makefile

.PHONY: all clean study reapply clean-study generate help
