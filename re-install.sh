#!bin/bash

#sudo ~/.anaconda3/bin/python3 setup.py install >> output
#sudo python3 setup.py install
python3 setup.py install --user && sh ./addon/install_addon.sh
