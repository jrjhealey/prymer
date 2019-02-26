from setuptools import setup

setup(name='prymer',
      version='1.0',
      description='Return primers for given input sequences',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.x',
        'Topic:: Scientific / Engineering:: Bio - Informatics'],
      url='https://github.com/jrjhealey/PrYmer',
      author='Joe R. J. Healey',
      author_email='jrj.healey@gmail.com',
      license='MIT',
      packages=['prymer'],
      install_requires=['biopython'],
      zip_safe=False)
