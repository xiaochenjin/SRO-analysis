#https://cython.readthedocs.io/en/latest/src/quickstart/build.html
import distutils.core
import Cython.Build
distutils.core.setup(
    ext_modules = Cython.Build.cythonize("computeSRO_dist.pyx"))
