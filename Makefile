.PHONY: setup setup-pip run clean

# Create the conda environment
setup:
	conda env create -f environment.yml

# Or install with pip into the current (virtual) environment
setup-pip:
	pip install -r requirements.txt

# Run the pipeline with the default config (example dataset out of the box)
run:
	python scrna_pipeline.py --config config/config.yaml

# Remove generated outputs
clean:
	rm -rf results/*.h5ad results/*.csv figures/*.png data/* figures
	touch results/.gitkeep figures/.gitkeep
