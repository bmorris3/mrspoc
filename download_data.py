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
    f.extractall()

# Remind the user to set an environment variable:
print('Now set an environment variable "MRSPOC_DATA_DIR" to the path to the '
      'data, directory {0}'.format(p))
