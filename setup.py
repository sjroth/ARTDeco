from setuptools import setup

setup(
    name='ARTDeco',
    version='0.2',
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
        'Topic :: Software Development :: Libraries'
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=[
        'pandas',
        'rpy2',
        'numpy',
        'bx-python',
        'RSeQC',
    ],
)
