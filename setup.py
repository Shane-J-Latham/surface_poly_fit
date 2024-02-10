from skbuild import setup
from setuptools import find_packages
import versioneer

setup(
    name="surface_poly_fit",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
    packages=find_packages(),
    package_data={'surface_poly_fit': ["*.txt"]},
    cmake_with_sdist=True
)
