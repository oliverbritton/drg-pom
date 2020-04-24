try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup
	
config = {
	'description': 'Population of models simulations for dorsal root ganglion neurons',
	'author': 'Oliver Britton',
	'url': 'https://github.com/oliverbritton/drg-pom'`,
	'download_url': 'https://github.com/oliverbritton/drg-pom/archive/master.zip',
	'author_email': 'oliverjbritton@gmail.com',
	'version': '0.1',
	'install_requires': [
            'nose2',
            'numpy<1.16.3',
            'pandas',
            'pyDOE',
            'matplotlib',
            'seaborn',
            ],
	'packages': ['drgpom'],
	'scripts': [],
	'name': 'drgpom'
        
}

setup(**config)
