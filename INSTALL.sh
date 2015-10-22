#!/usr/bin/env bash

echo "Installing git ..."

sudo apt-get install git

echo "Fetching TEES ..."

git clone https://github.com/jbjorne/TEES.git
sudo apt-get install ruby
sudo apt-get install make
sudo apt-get install g++-multilib
sudo apt-get install flex
sudo apt-get install libboost-all-dev
sudo apt-get install python-numpy
sudo apt-get install python-scipy
sudo apt-get install python-sklearn
sudo apt-get install python-nltk

cd TEES
python configure.py
sudo python setup.py install 
