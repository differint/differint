from setuptools import setup

def readme():
      with open('README.md') as f:
            return f.read()

setup(name = 'differint',
      version = '1.0.0',
      description = 'Collection of algorithms for numerically calculating fractional derivatives.',
      long_description = readme(),
      long_description_content_type="text/markdown",
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
      ],
      url = 'http://github.com/differint/differint',
      author = 'Matthew Adams',
      author_email = 'Matthew.Adams@ucalgary.ca',
      license = 'MIT',
      packages = ['differint'],
      zip_safe = False,
      include_package_data = True,
      install_requires = ['numpy'],
      test_suite = 'nose.collector',
      tests_require = ['nose']
     )
