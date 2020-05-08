from setuptools import setup

setup(
    name='ARTDeco',
    version='0.4',
    packages=['ARTDeco'],
    url='https://github.com/sjroth/ARTDeco',
    license='MIT',
    author='Samuel J. Roth, M.A.',
    author_email='sjroth@eng.ucsd.edu',
    description='Automatic Readthrough DEteCtiOn (ARTDeco)-a pipeline for analyzing and quantifying transcriptional readthrough',
    entry_points={
        'console_scripts': [
            'ARTDeco=ARTDeco.main:main'
        ]
    },
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=[
        'bx-python>=0.8.2',
        'networkx>=2.2',
        'numpy>=1.16.2',
        'pandas>=0.24.2',
        'rpy2>=2.9.4',
        'RSeQC>=3.0.0',
    ],
    include_package_data = True
)
