#from setuptools import setup
from distutils.core import setup

with open("README.rst", 'r') as f:
    long_description = f.read()
with open("influx_version.txt", 'r') as f:
    version = f.read().rstrip()

setup(
   name='influx_si',
   version=version,
   description='Metabolic flux and concetration estimation based on stable isotope labeling',
   keywords='metabolic flux analysis, least squares, isotope labeling, systems biology',
   license="GPL2+",
   long_description=long_description,
   author='Serguei Sokol',
   author_email='sokol@insa-toulouse.fr',
   url="https://metasys.insa-toulouse.fr/software/influx/",
   packages=['influx_si'],
   package_data={'influx_si': ['licence_en.txt', 'influx_version.txt', 'R/*.R', 'test/e_coli.ftbl', 'test/e_coli_i.ftbl', 'test/e_coli_msen.txt', 'test/ok/*', 'doc/html/*', 'doc/*.pdf', 'doc/*.html']},
   requires=['scipy', 'kvh'],
   scripts=[
      'bin/ff2ftbl.py',
      'bin/res2ftbl_meas.py',
      'bin/ftbl2code.py',
      'bin/ftbl2cumoAb.py',
      'bin/ftbl2kvh.py',
      'bin/ftbl2metxml.py',
      'bin/ftbl2netan.py',
      'bin/ftbl2optR.py',
      'bin/ftbl2xgmml.py',
      'bin/influx_s.py',
      'bin/influx_i.py'
   ],
   classifiers=[
      'Environment :: Console',
      'Intended Audience :: End Users/Desktop',
      'License :: OSI Approved :: GPL2+',
      'Operating System :: MacOS :: MacOS X',
      'Operating System :: Microsoft :: Windows',
      'Operating System :: POSIX',
      'Programming Language :: Python :: Pyhton3',
      'Topic :: Science/Biology/Metabolism/Fluxes',
   ],
)
