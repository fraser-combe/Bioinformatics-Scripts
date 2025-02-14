#!/bin/bash

# Change to the abricate_db directory
cd /abricate_db

# Download the vfdb database
abricate-get_db --db vfdb

# Set up the database
abricate --setupdb

# Verify the setup
abricate --list

# Create a symbolic link to the abricate_db directory for the environment
ln -s /abricate_db /opt/conda/envs/abricate
