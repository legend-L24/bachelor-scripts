from distutils.core import setup

setup(
    name='kqeq',
    version='0.1',
    packages=['kqeq'],
    url='https://jmargraf.gitlab.io/kqeq/',
    license='The software is licensed under the MIT License',
    author='Johannes Margraf, ',
    author_email='margraf@fhi.mpg.de',
    description='A Python implementation of the kernel Charge Equilibration method, allows training and using physics-based machine-learning models for predicting charge distributions in molecules and materials.',
    #install_requires=['numpy<1.20.0','scipy','ase'],#I delete dscribe
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
