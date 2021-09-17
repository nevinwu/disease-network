# disease-network
Some scripts I used during my master's final project:

- **GSEXXXX.r**: this files' names refer to GEO's entries that are being processed. Each file belongs to one of this datasets.
- **join_datasets.r**: this is the file used to integrate and join information and data among differeten datasets. Furthermore, it obtains DEGs that will be nodes in our network.
- **create_netwokr.r**: uses previously created interactome and previously obtained DEGs to map these DEGs over the interactome and so obtain out disease network.
- **node_analysis.r**: analysis results of Cytoscape to prioritaze certain amount of nodes which we will study as drug targets.
- **promiscuity.r**: studies our drug's promiscuity among our disease network.
- **network.cys**: Cytoscape file that contains the created disease network.
- **funtional_analysis**: directory containing files related to funtional analysis performed with BiNGO.
- **functions**: directory containing R functions used during this project.
- **interactome.rda**: interactome used to create our disease network. Obtained from Luck, K., Kim, DK., Lambourne, L. et al. A reference map of the human binary protein interactome. Nature 580, 402–408 (2020). Link: http://www.interactome-atlas.org/
