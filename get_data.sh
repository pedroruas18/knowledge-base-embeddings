
mkdir -p data/


###############################################################################
#                        Get node2vec repository                              #
###############################################################################
git clone https://github.com/pedroruas18/node2vec.git


###############################################################################
#                           DOWNLOAD KB FILES                                 #
###############################################################################
mkdir -p data/kbs/
cd data/kbs/

# MEDIC 
mkdir -p medic/
cd medic/
wget ctdbase.org/reports/CTD_diseases.tsv.gz 
gzip -d CTD_diseases.tsv.gz 
cd ../

cd ../../