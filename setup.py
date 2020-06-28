import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="feats", 
    version="1.0.1",
    author="Edwin Vans",
    author_email="vans.edw@gmail.com",
    description="A Python tool for performing downstream analysis on Single Cell RNA-seq datasets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/edwinv87/feats",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)