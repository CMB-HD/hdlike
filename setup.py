from setuptools import setup

setup(
    name="hdlike",
    version="1.0",
    description="Likelihood for CMB-HD",
    url="https://github.com/CMB-HD/hdlike",
    author="CMB-HD Collaboration",
    python_requires=">=3",
    install_requires=["numpy"],
    packages=["hdlike"],
    include_package_data=True,
)
