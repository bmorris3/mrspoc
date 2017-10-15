from setuptools import setup


setup(name='mrspoc',
      version='0.1',
      description='M-dwarf Rotational Stellar Photocenter Offset Calculator',
      install_requires=['numpy', 'astropy', 'matplotlib', 'scipy'],
      author='Brett Morris',
      author_email='bmmorris@uw.edu',
      license='MIT',
      url='https://github.com/bmorris3/mrspoc',
      zip_safe=False,
      use_2to3=False,
)