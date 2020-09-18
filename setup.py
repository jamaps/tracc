import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tracc",
    version="0.0.2",
    author="Jeff Allen",
    author_email="jeff.allen@utoronto.ca",
    description="Transport accessibility measures in Python",
    long_description="Transport accessibility measures in Python",
    long_description_content_type="text/markdown",
    url="https://github.com/jamaps/tracc",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'pandas>=1.0',
        'numpy>=1.18.5',
        'geopandas>=0.7.0',
        "libpysal>=4.3.0"
    ]
)

# https://packaging.python.org/tutorials/packaging-projects/
