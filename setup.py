import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

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
