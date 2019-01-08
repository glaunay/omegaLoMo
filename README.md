# Building Protein Interaction Network

Collection of notebooks

## External documentation
The [pyproteinsExt documentation](https://github.com/glaunay/pyproteinsExt) details the API of the  _mitabTopology_ and _database.uniprotFastaFS_ used in these notebooks.

* omega package

## Target Proteome

 * Specie: Pneumococcus Pneumoniae strain ATCC BAA-255 / R6
 * Size : 2033 proteins in fasta sequences

### Enriching target proteome sequence set

See [divisomeFactory repository](https://github.com/glaunay/divisomeFactory) for details.
This consists of the following steps

 1. Collecting Uniclust
 2. Setting Uniclust Blast Database
 3. Querying Uniclust for each target proteome sequence
 4. Merging the intial proteome fasta set with its uniclust homologs
 5. Building the corresponding enriched target proteome blast database

## Building The experimental protein-protein interaction set

### Provider database and version

* Intact
* MatrixDB
* Mint
* Biogrid

See [mitab culling notebook](https://github.com/glaunay/omegaLoMo/blob/master/interaction%20dataset%20building%20and%20uniprot%20FS%20API.ipynb) for details.

### Discard experiments not supported by valid uniprot identifiers

See [mitab culling notebook](https://github.com/glaunay/omegaLoMo/blob/master/interaction%20dataset%20building%20and%20uniprot%20FS%20API.ipynb) for details.


## Linking Target Proteome with experimental interaction

This produces the following file use in this analysis

1. **Molecular interaction File**Molecular interaction File. It respects the mitab format, all pairs of interactors do not necessarly have corresponding homologs in target proteome.
2. **homologyTree.json**. Two level dictionary binding experimental interactors to target proteome sequences based on homology.

Please see [divisomeFactory repository](https://github.com/glaunay/divisomeFactory) for details on procedure and output files format.

### Culling the interaction dataset

Using the provided notebook, [reference_MI_toplogy](https://github.com/glaunay/omegaLoMo/blob/master/reference_MI_toplogy.ipynb) discard interactions without homologs in target proteome

#### inputs

* _mitabFile_ (**Molecular interaction File**)
* _homologyFile_ (**homologyTree.json**)

#### outputs
A file pickle dump of an instance of a _mitabTopology_ Object which stores all experimental interactions with both interactors havong homologs in the target proteome



### Modeling interaction in target proteome
The complete network as an adjacency matrix along with edge information is always kept.
An edge element has the folowing strucure:


#### Implementation

Edge is a list of 2-tuples
2-tuples can be invalidated if they violate trim conditions


networkx like graph can be obtained by filtering the complete graph based on egde, node properties

#### Edge triming

An egdge with no valid 2-tuples is not Visible
A node not connected to a visible edge is hidden

#### Seed pruning

A **previously trimmed** network can be pruned using a set of seed nodes.
A path must exist between any node of the network and at least one seed node. Otherwise, the nodes is hidden.



### Enriching target proteome ppi with uniprot information

### Visualizing the network

#### The network litteral object

Top three keys :
* nodes
* links
* registry

```json
{
    "nodes": [
        {
            "id": "Q8DNW6",
            "group": 0,
            "val": 4
        },
        {
            "id": "Q8DQD4",
            "group": 0,
            "val": 25
        },
...
```

```json
"links": [
        {
            "source": "Q8DNW6",
            "target": "Q8DQD4",
            "data": {
                "lowQueryParam": [
                    [
                        "P75831",
                        "236",
                        "7",
                        "216",
                        "648",
                        "5",
                        "218",
                        "138",
                        "89",
                        "2.52549e-64"
                    ],
                    [
                        "P75831",
                        "236",
                        "7",
                        "216",
                        "648",
                        "5",
                        "218",
                        "138",
                        "89",
                        "2.52549e-64"
                    ]
                ],
                "highQueryParam": [
                    [
                        "P75831",
                        "244",
                        "1",
                        "232",
                        "648",
                        "1",
                        "235",
                        "141",
                        "95",
                        "4.18918e-74"
                    ],
                    [
                        "P75831",
                        "244",
                        "1",
                        "232",
                        "648",
                        "1",
                        "235",
                        "141",
                        "95",
                        "4.18918e-74"
                    ]
                ]
            }
        },
        ....
```

```json
"registry": {
        "P0A2U9": {
            "id": "P0A2U9",
            "name": "AMIE_STRR6",
            "fullName": "Oligopeptide transport ATP-binding protein AmiE",
            "geneName": "amiE",
            "GO": [
                {
                    "id": "GO:0005886",
                    "term": "C:plasma membrane",
                    "evidence": "ECO:0000501"
                },
                {
                    "id": "GO:0005524",
                    "term": "F:ATP binding",
                    "evidence": "ECO:0000501"
                },
                {
                    "id": "GO:0016887",
                    "term": "F:ATPase activity",
                    "evidence": "ECO:0000501"
                },
                {
                    "id": "GO:0015031",
                    "term": "P:protein transport",
                    "evidence": "ECO:0000501"
                }
            ]
        },
        "P33916": {
            "id": "P33916",
            "name": "YEJF_ECOLI",
            "fullName": "Uncharacterized ABC transporter ATP-binding protein YejF",
            "geneName": "yejF",
            "GO": [
                {
                    "id": "GO:0005886",
                    "term": "C:plasma membrane",
                    "evidence": "ECO:0000314"
                },
                {
                    "id": "GO:0005524",
                    "term": "F:ATP binding",
                    "evidence": "ECO:0000255"
                },
                {
                    "id": "GO:0042626",
                    "term": "F:ATPase activity, coupled to transmembrane movement of substances",
                    "evidence": "ECO:0000250"
                },
                {
                    "id": "GO:0042884",
                    "term": "P:microcin transport",
                    "evidence": "ECO:0000315"
                },
                {
                    "id": "GO:0035672",
                    "term": "P:oligopeptide transmembrane transport",
                    "evidence": "ECO:0000314"
                }
            ]
        }
    ...
```

#### The mitab litteral object

