from setuptools import setup


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="paraphase",
    version="1.1.1",
    description="paraphase: HiFi-based SMN1/SMN2 variant caller",
    long_description=readme(),
    url="https://github.com/PacificBiosciences/paraphase",
    author="Xiao Chen",
    author_email="xchen@pacificbiosciences.com",
    license="BSD-3-Clause-Clear",
    packages=["paraphase"],
    package_data={"paraphase": ["data/*", "data/smn1/*"]},
    install_requires=["pysam", "numpy", "scipy", "networkx", "matplotlib"],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    entry_points={"console_scripts": ["paraphase=paraphase.paraphase:main"]},
    include_package_data=True,
    long_description_content_type="text/markdown",
)
