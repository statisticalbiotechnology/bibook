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

build-book: imgs
	$(CONDA_ACTIVATE) jb; cd $(BOOK_DIR) && jupyter book build

# Step 1: Build the book with Jupyter Book using the LaTeX builder
$(BUILD_DIR):
	cd $(BOOK_DIR) && jupyter book build --tex

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

.PHONY: all clean clean-all


