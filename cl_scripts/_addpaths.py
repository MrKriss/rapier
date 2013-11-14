'''
Created on 9 Jan 2013

@author: musselle
'''

# Relative path additions 
""" In a file called _mypath.py"""

import os, sys

#thisdir = os.path.dirname(__file__)
#
#utildir = os.path.join(thisdir, '../utils')
#runscriptdir = os.path.join(thisdir, '../runscripts')

# Work out where data is stored
cpu_name = os.uname()[1] 

if cpu_name == 'yildun':
    data_prefix = '/space/musselle/data'
    work_prefix = '/home/pgrad/musselle/ubuntu/workspace/popGen/'
        
elif cpu_name == 'luca':
    data_prefix = '/home/musselle/san/data'
    work_prefix = '/home/musselle/popGen/'
       
elif cpu_name == 'gg-pc6':
    data_prefix = '/home/musselle/data'
    work_prefix = '/home/musselle/popGen/'

# elif cpu_name == "Musselles-MacBook":
else:
    data_prefix = '/Users/chris/data'
    work_prefix = '/Users/chris/Dropbox/work/popGen/code/'

sys.path.append(work_prefix)
sys.path.append(os.path.join(work_prefix, 'utils'))
sys.path.append(os.path.join(work_prefix, 'plotscripts'))
sys.path.append(os.path.join(work_prefix, 'runscripts'))
sys.path.append(os.path.join(work_prefix, 'clf'))
    
del work_prefix
del data_prefix

""" then put import _mypath at the top of each script """