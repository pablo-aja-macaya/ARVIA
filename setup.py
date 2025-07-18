from setuptools import setup

setup(
    name="ARVIA",
    version="0.1.0",
    description="Antibiotic Resistance Variant Identifier for Pseudomonas aeruginosa",
    url='https://github.com/Pablo-Aja-Macaya/ARVIA',
    author="Pablo Aja",
    author_email="test@gmail.com",
    license="",
    packages=["arvia"],
    entry_points={
        "console_scripts": ["arvia=arvia.arvia:main"],
    },
    install_requires=[
        # "pandas",
        # "numpy",
        # "pytest",
        # "pytest-cov",
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        # "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.8",
    ],
)
