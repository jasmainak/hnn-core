#! /usr/bin/env python
import os.path as op
import os
import subprocess

from setuptools import setup, find_packages

from distutils.command.build_py import build_py
from distutils.cmd import Command

descr = """Experimental code for simulating evoked using Neuron"""

DISTNAME = 'hnn-core'
DESCRIPTION = descr
MAINTAINER = 'Mainak Jas'
MAINTAINER_EMAIL = 'mainakjas@gmail.com'
URL = ''
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'http://github.com/jonescompneurolab/hnn-core'
VERSION = '0.1.dev0'

# test install with:
# $ python setup.py clean --all install
#
# to make sure there are no residual mod files
#
# also see following link to understand why build_py must be overriden:
# https://stackoverflow.com/questions/51243633/python-setuptools-setup-py-install-does-not-automatically-call-build
class BuildMod(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        print("=> Building mod files ...")
        mod_path = op.join(op.dirname(__file__), 'hnn_core', 'mod')
        process = subprocess.Popen(['nrnivmodl'], cwd=mod_path,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        outs, errs = process.communicate()
        print(outs)


class build_py_mod(build_py):
    def run(self):
        build_py.run(self)
        self.run_command("build_mod")


if __name__ == "__main__":
    setup(name=DISTNAME,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          long_description=open('README.rst').read(),
          classifiers=[
              'Intended Audience :: Science/Research',
              'Intended Audience :: Developers',
              'License :: OSI Approved',
              'Programming Language :: Python',
              'Topic :: Software Development',
              'Topic :: Scientific/Engineering',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS',
          ],
          platforms='any',
          packages=find_packages(),
          package_data={'hnn_core':
            ['param/*.json',
             'mod/*',
             'mod/x86_64/*',
             'mod/x86_64/.libs/*']},
          include_package_data=True,
          cmdclass={'build_py': build_py_mod, 'build_mod': BuildMod}
          )
