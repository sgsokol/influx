#from setuptools import setup
from distutils.core import setup

with open("README.rst", 'r') as f:
    long_description = f.read()
with open("influx_version", 'r') as f:
    version = f.read()

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
   package_data={'influx_si': ['licence_en.txt', 'influx_version.txt', '../*.R', '../test/е_coli.ftbl', '../test/е_coli_i.ftbl', '../test/e_coli_msen.txt', '../test/ok/*', '../doc/html/*', '../doc/*.pdf']},
   requires=['scipy', 'kvh'],
   scripts=[
      'ftbl2code.py',
      'ftbl2cumoAb.py',
      'ftbl2kvh.py',
      'ftbl2metxml.py',
      'ftbl2netan.py',
      'ftbl2optR.py',
      'ftbl2xgmml.py',
   ],
   classifiers=[
      'Environment :: Console',
      'Intended Audience :: End Users/Desktop',
      'License :: OSI Approved :: GPL2+',
      'Operating System :: MacOS :: MacOS X',
      'Operating System :: Microsoft :: Windows',
      'Operating System :: POSIX',
      'Programming Language :: Python :: Pyhton3',
      'Topic :: Science/Biology/Metabolism',
   ],
)
