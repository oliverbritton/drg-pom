try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup
	
config = {
	'description': 'ProjectName',
	'author': 'Oliver Britton',
	'url': 'URL to get it at.',
	'download_url': 'Where to download it.',
	'author_email': 'oliverjbritton@gmail.com',
	'version': '0.1',
	'install_requires': ['nose2'],
	'packages': ['NAME'],
	'scripts': [],
	'name': 'projectname'
}

setup(**config)