import sys
from setuptools import setup, find_packages

if sys.version_info < (3, 8):
    sys.exit("Sorry, Python < 3.8 is not supported")

setup(
    name="nanoreactor_processing",
    packages=["nanoreactor_processing"],
    version="1.1.1",
    license="MIT",
    description="Automated evaluation of computational nanoreactor simulations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Alexandra Stan",
    author_email="alexandra.stan@cup.uni-muenchen.de",
    url="https://github.com/ochsenfeld-lab/nanoreactor_processing",
    download_url="https://github.com/ochsenfeld-lab/nanoreactor_processing/archive/refs/heads/main.zip",
    keywords=[
        "computational chemistry",
        "automated",
        "molecular dynamics",
        "chemical reactions",
        "ab initio nanoreactor",
	"rdkit"
    ],
    install_requires=[
        "numpy>=1.23.1",
        "scipy>=1.8.1",
	"pandas>=1.4.3",
        "networkx>=2.8.4",
        "matplotlib>=3.5.2",
        "seaborn>=0.11.2",
        "pyvis>=0.2.1",
        "rdkit-pypi==2021.3.4",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",  # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        "Environment :: Console",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
    ],
    zip_safe=False,
)
