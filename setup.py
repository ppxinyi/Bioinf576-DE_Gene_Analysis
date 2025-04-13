from setuptools import setup, find_packages

setup(
    name='dge',
    version='0.1.0',
    author='Xinyi Deng',
    description='A toolkit for differential gene expression analysis',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'matplotlib',
        'seaborn',
        'rpy2',
        'scikit-learn',
        'numpy'
    ],
    entry_points={
        'console_scripts': [
            'dge=dge.cli:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.8',
)
