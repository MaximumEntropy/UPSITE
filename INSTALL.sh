#!/usr/bin/env bash

echo "Installing git ..."

sudo apt-get install git

sudo apt-get -y install ruby
sudo apt-get -y install make
sudo apt-get -y install g++-multilib
sudo apt-get -y install flex
sudo apt-get -y install libboost-all-dev
sudo apt-get -y install python-numpy
sudo apt-get -y install python-scipy
sudo apt-get -y install python-sklearn
sudo apt-get -y install python-nltk
sudo apt-get -y install default-jre
sudo apt-get -y install default-jdk

cd TEES
python configure.py
sudo python setup.py install 
