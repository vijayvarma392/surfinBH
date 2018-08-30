import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="surfinBH",
    version="0.0.8.dev1",
    author="Vijay Varma",
    author_email="vvarma@caltech.edu",
    description="Surrogate Final BH properties.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vijayvarma392/surfinBH",
    packages=setuptools.find_packages(),
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Operating System :: OS Independent",
    ],
)
