from setuptools import find_packages, setup

setup(
    author="Jie Li",
    description="Bifunctional Canonization",
    name='CanonizedRMSD',
    packages=find_packages(
        include=['CanonizedRMSD', 'CanonizedRMSD.*']),
    include_package_data=True,
    entry_points={'console_scripts': [
        'canonize = CanonizedRMSD.scripts.canonize:main',
    ]},
    version='1.0.0',
)
