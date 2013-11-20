'''
Created on 9 Jan 2013

@author: musselle
'''

# Relative path additions using import _addpaths.py in each script
import os, sys

this_dir = os.path.dirname(__file__)
util_dir = os.path.abspath(os.path.join(this_dir, '../lib'))
testdata_dir = os.path.abspath(os.path.join(this_dir, '../testdata'))

sys.path.insert(0, this_dir)
sys.path.insert(0, util_dir)
sys.path.insert(0, testdata_dir)