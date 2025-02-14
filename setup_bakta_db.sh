#!/bin/bash

BAKTA_DB_DIR="$HOME/bakta_db"

mkdir -p "$BAKTA_DB_DIR/full"
mkdir -p "$BAKTA_DB_DIR/light"

if [ ! -d "$BAKTA_DB_DIR/$1/db" ]; then
  echo "Downloading Bakta $1 database..."
  bakta_db download --output "$BAKTA_DB_DIR/$1/db" --type "$1"
else
  echo "Bakta $1 database already exists, skipping download."
fi
