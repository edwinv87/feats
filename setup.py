import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="feats", # Replace with your own username
    version="1.0.0",
    author="Edwin Vans",
    author_email="vans.edw@gmail.com",
    description="A python tool for performing downstream analysis on Single Cell RNA-seq datasets",
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