#!/usr/bin/env python
import setuptools

setuptools.setup(
    name="stacpro",
    version="0.1.1",
    author="Dantong Wang, Mengyang XU",
    author_email="wangdantong@genomics.cn, xumengyang@genomics.cn",
    url="https://github.com/BGI-Qingdao/stacpro",
    #long_description=Path('README.md').read_text('utf-8'),
    python_requires=">=3.7,<3.11",
    packages=setuptools.find_packages(),
    install_requires=[
        "pandas",
        "matplotlib",
        "numpy",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    license="GPL-3.0+",
    description="Structure based alignment and clustering of protein.",
    platforms='any',
        ],
    },
)