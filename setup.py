from setuptools import setup
#from distutils.core import setup
from pathlib import Path

with open("README.rst", "r", encoding="utf-8") as f:
    long_description = f.read()
with open("influx_si/influx_version.txt", "r") as f:
    version = f.read().rstrip()
    
testdir = Path(__file__).parent / 'influx_si' / 'test'
tfiles = [str(p.relative_to(testdir)) for p in testdir.rglob('*')]
rdir = Path(__file__).parent / 'influx_si' / 'R'
rfiles = [str(p.relative_to(rdir)) for p in rdir.rglob('*')]

setup(
   name='influx_si',
   version=version,
   description='Metabolic flux and concentration estimation based on stable isotope labeling',
   keywords='metabolic flux analysis, least squares, isotope labeling, systems biology',
   license="GNU General Public License v2 or later (GPLv2+)",
   long_description=long_description,
   author='Serguei Sokol',
   author_email='sokol@insa-toulouse.fr',
   url="https://github.com/sgsokol/influx/",
   packages=['influx_si'],
   package_data={'influx_si': ['licence_en.txt', 'influx_version.txt', 'R/*.R', 'test/*/*', 'test/ok/*/*/*/*']},
   include_package_data=True,
   python_requires='>3.5',
   install_requires=['scipy', 'python-libsbml', 'pandas', 'kvh', 'packaging', 'asteval'],
   scripts=[
      'influx_si/bin/ff2ftbl.py',
      'influx_si/bin/ftbl2code.py',
      'influx_si/bin/ftbl2cumoAb.py',
      'influx_si/bin/ftbl2kvh.py',
      'influx_si/bin/ftbl2netan.py',
      'influx_si/bin/ftbl2optR.py',
      'influx_si/bin/ftbl2xgmml.py',
      'influx_si/bin/res2ftbl_meas.py',
      'influx_si/influx_s.py',
      'influx_si/influx_i.py',
      'influx_si/txt2ftbl.py',
      'influx_si/ftbl2mtf.py',
      'influx_si/ftbl2metxml.py',
      'influx_si/ftbl2labcin.py',
   ],
   entry_points={
      "console_scripts": [
         "ff2ftbl = influx_si.cli:cli",
         "ftbl2code = influx_si.cli:cli",
         "ftbl2cumoAb = influx_si.cli:cli",
         "ftbl2kvh = influx_si.cli:cli",
         "ftbl2netan = influx_si.cli:cli",
         "ftbl2optR = influx_si.cli:cli",
         "ftbl2xgmml = influx_si.cli:cli",
         "res2ftbl_meas = influx_si.cli:cli",
         "influx_s = influx_si.influx_s:main",
         "influx_i = influx_si.influx_i:main",
         "txt2ftbl = influx_si.txt2ftbl:main",
         "ftbl2mtf = influx_si.ftbl2mtf:main",
         "ftbl2metxml = influx_si.ftbl2metxml:main",
         "ftbl2labcin = influx_si.ftbl2labcin:main",
      ],
   },
   classifiers=[
      'Environment :: Console',
      'Intended Audience :: End Users/Desktop',
      'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
      'Operating System :: MacOS :: MacOS X',
      'Operating System :: Microsoft :: Windows',
      'Operating System :: POSIX',
      'Programming Language :: Python :: 3 :: Only',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
   ],
   project_urls={
      'Documentation': 'https://influx-si.readthedocs.io/',
      'Source': 'https://github.com/sgsokol/influx',
      'Tracker': 'https://github.com/sgsokol/influx/issues',
   },
)
