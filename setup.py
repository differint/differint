from setuptools import setup

def readme():
      with open(README.rst) as f:
            return f.read()

setup(name = 'differint',
      version = '0.1',
      description = 'Collection of algorithms for numerically calculating fractional derivatives.',
      long_description = readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.3',
        'Intended Audience :: Scientists',
        'Topic :: Numerical Analysis :: Fractional Calculus',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
      ],
      url = 'http://github.com/snimpids/differint',
      author = 'Matthew Adams',
      author_email = 'Matthew.Adams@ucalgary.ca',
      license = 'MIT',
      packages = ['differint'],
      zip_safe = False,
      include_package_data = True,
      install_requires = ['numpy'])
