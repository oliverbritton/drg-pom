import os

try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup
	
config = {
	'description': 'Population of models simulations for dorsal root ganglion neurons',
	'author': 'Oliver Britton',
	'url': 'https://github.com/oliverbritton/drg-pom',
	'download_url': 'https://github.com/oliverbritton/drg-pom/archive/master.zip',
	'author_email': 'oliverjbritton@gmail.com',
	'version': '1.0.4',
	'install_requires': [
            'nose2',
            'numpy<1.16.3',
            'pandas',
            'pyDOE',
            'matplotlib',
            'seaborn',
            ],
	'packages': ['drgpom', 'drgpom.methods', 'drgpom.methods.biomarkers',
            'drgpom.examples', 'drgpom.examples.data', 'drgpom.models'],
	'scripts': [],
	'name': 'drgpom',
        'include_package_data': True,
        'package_data': {'drgpom':[
                                   os.path.join('examples','data', '*.*'),
                                   os.path.join('examples','results','*.*'),
                                   os.path.join('models','*.mod'),
                                   os.path.join('examples','*.ipynb'),
                                   ]
                        }
        
}

setup(**config)
