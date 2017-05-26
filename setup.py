from setuptools import setup

setup(name='gsea',
    version='0.0.1',
    author = "Tristan R. Brown",
    author_email = "brown.tristan.r@gmail.com",
    description = ("An implementation of the Gene Set Enrichment Analysis "
                        "algorithm."),
    url = 'https://github.com/tristanbrown/gsea',
    license = "MIT",
    packages = ['gsea', 'test'],
    install_requires = ['numpy'],
    entry_points = {
        'console_scripts': [
            'gsea = my_project.__main__:main'
        ]
    },
    )