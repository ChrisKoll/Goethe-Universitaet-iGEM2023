#!\bin\sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
conda env create -f promotech_env.yml
conda activate promotech_env

mkdir model
wget http://www.cs.mun.ca/~lourdes/public/PromoTech_models/RF-HOT.zip
mv RF-HOT.zip /model/RF-HOT.zip
unzip model/RF-HOT.zip


wget http://www.cs.mun.ca/~lourdes/public/PromoTech_models/RF-TETRA.zip
mv RF-HOT.zip /model/RF-TETRA.zip
unzip model/RF-TETRA.zip
