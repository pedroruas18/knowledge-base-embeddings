#!/usr/bin/env bash

kb=$1 # Target knowledge base. Ex: medic
kb_format=$2 # File format of the target knowledge base: tsv, obo, txt

echo '» preparing input to node2vec...' 
python input.py $kb $kb_format

echo '» Generating candidate embeddings with node2vec...'

mkdir -p data/embeddings/$kb

python node2vec/src/main.py \
    --input node2vec/graph/$kb.edgelist \
    --output ./data/embeddings/$kb/$kb.emb \
    --dimensions 200 \
    --directed
echo '» Done! knowledge base embeddings in dir data/embeddings/'$kb