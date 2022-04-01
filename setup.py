from setuptools import setup
#from distutils.core import setup

with open("README.rst", "r", encoding="utf-8") as f:
    long_description = f.read()
with open("influx_si/influx_version.txt", "r") as f:
    version = f.read().rstrip()

setup(
   name='influx_si',
   version=version,
   description='Metabolic flux and concentration estimation based on stable isotope labeling',
   keywords='metabolic flux analysis, least squares, isotope labeling, systems biology',
   license="GNU General Public License v2 or later (GPLv2+)",
   long_description=long_description,
   author='Serguei Sokol',
   author_email='sokol@insa-toulouse.fr',
   url="https://metasys.insa-toulouse.fr/software/influx/",
   packages=['influx_si'],
   package_data={'influx_si': ['licence_en.txt', 'influx_version.txt', 'R/*.R', 'test/e_coli.ftbl',  'test/e_coli_growth.ftbl', 'test/e_coli_i.ftbl', 'test/e_coli_msen.txt', 'test/ok/*', 'test/prl_exp/*', 'doc/html/*', 'doc/*.pdf', 'doc/*.html']},
   install_requires=['scipy', 'python-libsbml', 'pandas'],
   scripts=[
      'influx_si/bin/ff2ftbl.py',
      'influx_si/bin/ftbl2code.py',
      'influx_si/bin/ftbl2cumoAb.py',
      'influx_si/bin/ftbl2kvh.py',
      'influx_si/bin/ftbl2netan.py',
      'influx_si/bin/ftbl2optR.py',
      'influx_si/bin/ftbl2xgmml.py',
      'influx_si/bin/influx_s.py',
      'influx_si/bin/influx_i.py',
      'influx_si/bin/res2ftbl_meas.py',
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
         "influx_s = influx_si.cli:cli",
         "influx_i = influx_si.cli:cli",
         "res2ftbl_meas = influx_si.cli:cli",
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
      'Documentation': 'https://metasys.insa-toulouse.fr/software/influx/doc/',
      'Source': 'https://github.com/sgsokol/influx',
      'Tracker': 'https://github.com/sgsokol/influx/issues',
   },
)
