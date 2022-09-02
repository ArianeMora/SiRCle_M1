from setuptools import setup
import os
import re


def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'sircle/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='sircle',
      version=read_version(),
      description='',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Ariane Mora',
      author_email='ariane.n.mora@gmail.com',
      url='https://github.com/ArianeMora/SiRCLe_M1',
      license='GPL3',
      project_urls={
          "Bug Tracker": "https://github.com/ArianeMora/SiRCLe_M1/issues",
          "Documentation": "https://github.com/ArianeMora/SiRCLe_M1",
          "Source Code": "https://github.com/ArianeMora/SiRCLe_M1",
      },
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords=['epigenetics', 'bioinformatics'],
      packages=['sircle'],
      entry_points={
          'console_scripts': [
              'sircle = sircle.__main__:main'
          ]
      },
      install_requires=['pandas', 'numpy', 'scibiomart', 'sciutil>=1.0.3', 'scimotf', 'scie2g'],
      python_requires='>=3.6',
      data_files=[("", ["LICENSE"])]
      )