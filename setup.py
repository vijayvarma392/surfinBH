'''
Standard setup.py to upload the code on pypi.
Make sure the fit data for all fits is in surfinBH/data before doing this.
Do:
    python setup.py sdist bdist_wheel
    twine upload dist/*
'''
import setuptools

with open("README.md", "rb") as fh:
    long_description = fh.read().decode("UTF-8")

# Extract code version from __init__.py
def get_version():
    with open('surfinBH/surfinBH.py') as f:
        for line in f.readlines():
            if "__version__" in line:
                return line.split('"')[1]

setuptools.setup(
    name="surfinBH",
    version=get_version(),
    author="Vijay Varma",
    author_email="vvarma@caltech.edu",
    description="Surrogate Final BH properties.",
    keywords='black-holes gravitational-waves Gaussian-process-regression',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vijayvarma392/surfinBH",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'h5py',
        'scikit-learn>=0.19.1,<0.23',
        'lalsuite>=6.70',
        'gwsurrogate',
    ],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Operating System :: OS Independent",
    ],
)
