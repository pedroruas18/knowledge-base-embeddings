# Generate node2vec embeddings for knowledge bases

Install required packages:

```
apt update && apt upgrade && xargs apt-get -y install < apt.txt && pip install -r requirements.txt 
```

Download the necessary data:

```
./get_data.sh
```

node2vec implementation adapted from https://github.com/aditya-grover/node2vec available at: https://github.com/pedroruas18/node2vec

Generate input for node2vec algorithm and run node2vec:

```
./gen_embeds.sh $kb $format
```

Change node2vec input arguments ('node2vec/src/main.py') in 'gen_embeds.sh', line 10.