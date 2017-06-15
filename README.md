<p align="center"><img src="./pictures/logomegalomo.png" alt="Drawing" style="width: 300px;"/></p>

Finding new interactions between proteins in a proteome, based on sequence homology.

<img src="./pictures/FullPipeline.png" alt="Drawing" style="width: 1000px;"/>

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites


This pipeline is currently using [Jupyter](http://jupyter.org/) to work.

It also use [Python 2.7](https://www.python.org/download/releases/2.7/) and more particularly these packages :


```python
import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import copy
import json
```


### Installing

First of all, you have to download the full package.

The first step is to clone the git repository. In order to do that, open a terminal console and write this

command line :

```
git clone https://github.com/glaunay/omegaLoMo.git 
```

Import the package to your .py project

```
import PACKAGE_PATH.core as core
import PACKAGE_PATH.createTopo as cT
import PACKAGE_PATH.graph as graph
```






## Running the tests

After downloading the package and the example files, here are some tests that can help you to familiarise with all data and inputs.

Example files contains :

>An index file

>>List of the full proteome of our organism

>A small database 

>>Known interactions in numerous organisms (from litterature, ...)

>A few blast files

>A serialized topology.

>A list of target proteins

>>Proteins of interest 

A Jupyter notebook is also join in the git to help you at each step.



### Break down into end to end tests

To help you run the pipeline, please refere to the jupyter notebook provide.



At the end, you should be able to observe some graphs corresponding to the first neighbors of each proteins of interest.

These are the predicted interactions proteins-proteins in the genome of your organism.


## Deployment

If you want to use this pipeline with your own data set, make sure to check all input format


## Built With

* [Jupyter](http://jupyter.org/) - An open-source web application



## Authors

* **Loïc Barlet** - *Initial work* - [LoBarlet](https://github.com/LoBarlet)
* **Mohamed Bachir-Cherif** - *Initial work* - [mbachircherif](https://github.com/mbachircherif)
* **Guillaume Launay** - *Initial work* - [glaunay](https://github.com/glaunay)



