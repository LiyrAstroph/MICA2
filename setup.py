import os
from setuptools import setup
from setuptools.extension import Extension
from glob import glob
import pkgconfig
import numpy 

# if CC is not set, use the default value
if not os.environ.get("CC"):
  os.environ["CC"] = "mpicc"

# check Intel OneAPI MKL
if os.environ.get("MKLROOT"):
  mklroot = os.environ.get("MKLROOT")
  FlagIntelMKL=True
else:
  FlagIntelMKL=False

try:
  from Cython.Build import cythonize
except ImportError:
  raise RuntimeError('Cython not found.')

def configure_mpi():
  """
  get configurations of mpi
  """
  if pkgconfig.exists('mpich'):
    mpiconf = pkgconfig.parse('mpich')
  elif pkgconfig.exists('ompi'):
    mpiconf = pkgconfig.parse('ompi')
  else:
    raise SystemError("Not found MPICH or OpenMPI installed.")

  return mpiconf

mpiconf = configure_mpi()

basedir = os.path.dirname(os.path.abspath(__file__))
homedir = os.environ['HOME']
include_dirs = [basedir, os.path.join(basedir, "src"), numpy.get_include(),] + mpiconf['include_dirs'] \
              +[os.path.join(basedir, "cdnest"),]
library_dirs = [basedir] + mpiconf['library_dirs']

if os.name == 'nt':  # Windows, assumming MSVC compiler
  libraries = ['dnest']
  compiler_args = ['/Ox', '/fp:fast']
  link_args = []
elif os.name == 'posix':  # UNIX, assumming GCC compiler
  if FlagIntelMKL:
    include_dirs += [os.path.join(mklroot, "include"),]
    library_dirs += [os.path.join(mklroot, "lib"),]
    libraries = ['m', 'c', 'gsl', 'mkl_intel_ilp64', 'mkl_core', 'mkl_gnu_thread', 'gomp', 'pthread', ] 
    compiler_args = ['-DIntelMKL',]
  else:
    libraries = ['m', 'c', 'gsl', 'gslcblas', 'lapack', 'lapacke']
    compiler_args = []
  
  libraries += mpiconf['libraries']
  compiler_args += ['-O3', '-ffast-math', '-fcommon'] 
  link_args = []

src = [os.path.join(basedir, "python", "pymica", "pymica.pyx")] + glob(os.path.join(basedir, "src", "*.c")) \
     + glob(os.path.join(basedir, "cdnest", "*.c"))
headerfiles = [os.path.join(basedir, "python", "pymica", "pymica.pxd")] + glob(os.path.join(basedir, "src", "*.h")) \
     + glob(os.path.join(basedir, "cdnest", "*.h"))

extensions = cythonize([
    Extension("pymica.pymica", 
	  sources=src,
    depends=headerfiles,
	  extra_compile_args=compiler_args,
    extra_link_args=link_args,
    include_dirs=include_dirs,
    libraries=libraries,
    library_dirs=library_dirs
    ),
  ], annotate=False)

setup(
      name="pymica",
      version="2.1.1",
      author="Yan-Rong Li",
      author_email="liyanrong@mail.ihep.ac.cn",
      packages=["pymica", "pymica.utility"],
      package_dir={"pymica":"python/pymica", "pymica.utility":''},
      ext_modules = extensions,
)