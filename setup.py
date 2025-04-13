from setuptools import setup, find_packages

setup(
    name="DGE",
    version="0.1.0",
    author="Xinyi Deng",
    author_email="your_email@example.com",
    description="A tool for RNA-seq differential gene expression analysis",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "seaborn",
        "matplotlib",
        "scipy",
        "scikit-learn",
    ],
    entry_points={
        "console_scripts": [
            "dge-run = DGE.main:main"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.7",
)
