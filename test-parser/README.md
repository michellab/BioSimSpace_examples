
Collection of inputs to test that the Sire amber reader/writer gives consistent energies. 

bash driver.sh

Currently fails for
# Fail
~/sire.app/bin/python test-sireparser.py  inputs/1a.parm7 inputs/1a_equilibrated.rst7
# Fail
~/sire.app/bin/python test-sireparser.py  inputs/5.parm7 inputs/5_equilibrated.rst7
# Fail
~/sire.app/bin/python test-sireparser.py  inputs/THROM-1a.parm7 inputs/THROM-1a_equilibrated.rst7
# Fail
~/sire.app/bin/python test-sireparser.py  inputs/THROM-5.parm7 inputs/THROM-5_equilibrated.rst7

