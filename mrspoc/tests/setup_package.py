# import os

# If this package has tests data in the tests/data directory, add them to
# the paths here, see commented example
paths = [] #[os.path.join('data', '*fits')]

def get_package_data():
    return {'mrspoc.tests': paths}
