 # OmegaLoMo

The aim of this project is to find new interactions between proteins in a proteome, based on sequence homology.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

This pipeline is currently using [Jupyter](http://jupyter.org/) to work.

These packages are necessary to make OmegaLoMo work fine

```python
import sys
import os
import subprocess
import fnmatch
import networkx as nx
import matplotlib.pyplot as plt
import copy
import collections
import numpy as np
import operator
import networkx as nx
```


### Installing

A step by step series of examples that tell you have to get a development env running

First of all, you have to download the full package.

The first step is to clone the git repository. In order to do that, open a terminal window and write this

command line :

```
git clone https://github.com/glaunay/omegaLoMo.git 
```

Import the package to your .py project

```
import PACKAGE_PATH.core as core
```






## Running the tests

After downloading the pipeline and the example files, here are some tests that can help you to familiarise with all data.

Examples files contains an index file, a few blast files and a serialized topology.

A Jupyter notebook is also join to help you at each step.



### Break down into end to end tests


Create a new cell and initialize a new HomegaSet. Provide a serialize JSON File or file containing BLAST files.
Provide also the R6 Index file that contain all R6 IDs.

Using serialize JSON File:
```
omegaSet = core.HomegaSet(bean= 'JSON_FILE_PATH.json' , queryIdList= INDEX_R6_FILE_PATH)
```

Using a directory containing BLAST Files:
```
omegaSet = core.HomegaSet(path= 'DIR_PATH' , queryIdList= INDEX_R6_FILE_PATH)
```

Then, create a OmegaMatrix
```
omegaMatrix = ca.OmegaMatrix(topo = newDic, omegaSet = omegaSet)
omegaMatrix.reduceAndVectorInject()



queryTopo = omegaMatrix.project()
```

At this point, you should be able to observe some graphs corresponding to the first neighbors of each proteins of interest.

These are the predicted interactions proteins-proteins in the genome of your organism.


## Deployment

If you want to use this pipeline with your own data set, make sure to .........

## Built With

* [Jupyter](http://jupyter.org/) - An open-source web application


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Loïc Barlet** - *Initial work* - [LoBarlet](https://github.com/LoBarlet)
* **Mohamed Bachir-Cherif** - *Initial work* - [mbachircherif](https://github.com/mbachircherif)
* **Guillaume Launay** - *Initial work* - [glaunay](https://github.com/glaunay)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

