"""
Download and extract the data archive from the web.
"""

import tarfile
import os
from urllib.request import urlretrieve

data_archive_path = 'data.tar.gz'
data_archive_url = 'http://staff.washington.edu/bmmorris/docs/data.tar.gz'
p = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# Write local copy of the data archive:
print('Downloading (150 MB) and extracting data to "{0}"...'.format(p))
urlretrieve(data_archive_url, data_archive_path)

# Extract all files from the gzipped tar archive
with tarfile.open(data_archive_path, 'r|gz') as f:
    def is_within_directory(directory, target):
        
        abs_directory = os.path.abspath(directory)
        abs_target = os.path.abspath(target)
    
        prefix = os.path.commonprefix([abs_directory, abs_target])
        
        return prefix == abs_directory
    
    def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
    
        for member in tar.getmembers():
            member_path = os.path.join(path, member.name)
            if not is_within_directory(path, member_path):
                raise Exception("Attempted Path Traversal in Tar File")
    
        tar.extractall(path, members, numeric_owner=numeric_owner) 
        
    
    safe_extract(f)

# Remind the user to set an environment variable:
print('Now set an environment variable "MRSPOC_DATA_DIR" to the path to the '
      'data, directory {0}'.format(p))
