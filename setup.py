from setuptools import setup
import paraphase


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="paraphase",
    version=paraphase.__version__,
    description="paraphase: HiFi-based caller for highly homologous genes",
    long_description=readme(),
    url="https://github.com/PacificBiosciences/paraphase",
    author="Xiao Chen",
    author_email="xchen@pacificbiosciences.com",
    license="BSD-3-Clause-Clear",
    packages=["paraphase", "paraphase.genes"],
    package_data={"paraphase": ["data/*", "data/**/*"]},
    install_requires=["pysam", "numpy", "scipy", "networkx", "matplotlib", "pyyaml"],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    entry_points={"console_scripts": ["paraphase=paraphase.__main__:main"]},
    long_description_content_type="text/markdown",
)
