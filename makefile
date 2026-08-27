.ONESHELL:
SHELL = /bin/bash

CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
# Variables
BOOK_DIR=bibook
BUILD_DIR=$(BOOK_DIR)/exports
TEX_FILE=$(BUILD_DIR)/bibook.tex
PDF_FILE=$(BUILD_DIR)/bibook.pdf


all: build-book $(PDF_FILE)

imgs:
	@$(MAKE) -C bibook/msa/img

# Live preview. --execute is NOT the default for `jupyter book start`, and
# without it every alignviz figure renders as an empty box.
serve: imgs
	$(CONDA_ACTIVATE) jb; cd $(BOOK_DIR) && jupyter book start --execute

# Mirror the executable .md chapters as .ipynb, so GitHub renders them and
# Colab can open them. The .md stays the source of truth; the .ipynb is
# generated - never edit it by hand.
NOTEBOOK_PAGES = bibook/pairwise/nw_code bibook/pairwise/needleman bibook/pairwise/waterman bibook/pairwise/semi

notebooks:
	$(CONDA_ACTIVATE) jb
	for p in $(NOTEBOOK_PAGES); do \
	    jupytext --to ipynb --from myst "$$p.md" -o "$$p.ipynb" ; \
	done

# Fail if a generated .ipynb is missing or older than its .md, so a stale
# mirror cannot be committed unnoticed.
check-notebooks:
	@$(CONDA_ACTIVATE) jb
	@stale=0 ; \
	for p in $(NOTEBOOK_PAGES); do \
	    if [ ! -f "$$p.ipynb" ] || [ "$$p.md" -nt "$$p.ipynb" ]; then \
	        echo "stale or missing: $$p.ipynb" ; stale=1 ; \
	    fi ; \
	done ; \
	if [ $$stale -ne 0 ]; then echo "run: make notebooks" >&2 ; exit 1 ; fi ; \
	echo "notebook mirrors up to date"

build-book: imgs notebooks
	$(CONDA_ACTIVATE) jb; cd $(BOOK_DIR) && jupyter book build --html --execute

# Step 1: Build the book with Jupyter Book using the LaTeX builder
$(BUILD_DIR):
	cd $(BOOK_DIR) && jupyter book build --tex --execute

# Step 2: Replace all occurrences of "align*" with "aligned" in bibook.tex
$(TEX_FILE): $(BUILD_DIR)
	sed -i 's/align\*/aligned/g' $(TEX_FILE)

# Step 3: Compile the LaTeX file to PDF using latexmk
$(PDF_FILE): $(TEX_FILE)
	cd $(BUILD_DIR) && latexmk -f -pdf -dvi- bibook.tex

# Clean up build files
clean:
	cd $(BOOK_DIR) && jupyter book clean

# Clean up everything including the PDF
clean-all:
	cd $(BOOK_DIR) && jupyter book clean && rm -rf exports/

.PHONY: all serve notebooks check-notebooks clean clean-all


