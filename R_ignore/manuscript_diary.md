# Manuscript of Interest

Diary of manuscripts for Tri-CePPI project and `fomes` model.  
**_TA_**: are take away points for consideration or useful insights.




## Agent-Based Stochastic Models
### [Shchur _et al._ 2022/VGsim](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010409)

**_TA_**:  Fast/efficient, somewhat flexible, compartment based simulation software: accounts for host population structure, pathogen evolution, host immunity, contact matrices.

##### Model Framework      
- Runs forward sims through a Gillespie-style algorithm in a SIS Compartment model framework (can vary S and I compartments to model waning immunity, vaccination, etc.).    
- code base is C++ w/ Python interface/wrapper

##### Coding Decisions
- Allows for host population (contact) structure.    
- Adaptive molecular evolution but only positive fitness effects.    
- Immunity is modeled as a Markovian (no accruing of immunity)
- _Within_ population transmission depends on a contact density parameter --> NPIs are modeled through this contact density parameter
- _Between_ population transmission is under a migration model
	- From/To assymetric probability matrix
		- Deal with extinction of demes by always having individuals return home/have short trips  
	- Cumulative upper bounds of migration --> assumes within >> between and authors note that it is suboptimal if a freely mixing population
- Populations models as _distrete demes_ (within transmission above) with individuals traveling between demes by migration above



### [Moshiri _et al._ 2018/FAVITES](https://academic.oup.com/bioinformatics/article/35/11/1852/5161084?login=false)

**_TA_**: Simulates the full end-to-end epidemic dataset (social contact network, transmission history, incomplete sampling, viral phylogeny, error-free sequences and real-world sequencing imperfections) - via a generative model = computationally expensive.     

- _NB_ his dissertation on FAVITES has additional [details](https://escholarship.org/uc/item/62s7q92d)


##### Model Framework    
- Agent based stochastic simulation
	- truly modular framework allow for calling wide spectrum
-  Code base in Python with API for different modules.   
- File formats somewhat unique to program (multiple of them)

##### Coding Decisions    
- Pays particular attention to contact network: uses `NetworkX` for various implementations
- [General workflow](https://github.com/niemasd/FAVITES/wiki/General-Workflow) or [Fig1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6931354/)
	- Initialization via "seeds" or individuals infected at t=0
		- For _seed_ step genetics, uses HMM to carry around realistic viral sequences (throughout? - comp expensive...)  
- Transmission seems to be on a scheduler versus loop? Says various models either assume exponential waiting times versus POMPs moving between states (Need to dig into modules to see how these differ)  
- Final set of steps involve making realistic sampling frameworks and phylogenetic trees (_i.e._ branching processes)


**Considerations**

- Large number of dependencies, some of which are author's packages (eg. niemasd/Dual-Birth-Simulator, niemasd/GEMF in python vs C, TreeSwift)
		- _N.B._ many related to simulating sequences/seq error rates   

- Notes the importance of epistasis for viral evolution/dynamics?  

- Transmission network file format: Self-edges (i.e., same node in columns 1 and 2) denote removal of infection, either via recovery or death

  

### [Lequime _et al._ 2020/nosoi](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13422)

**_TA_**:  

##### Model Framework     
##### Coding Decisions  

### [Campbell _et al._ 2019/outbreaker2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006930)

**_TA_**:  

##### Model Framework     
##### Coding Decisions  



## ID Dynamic Modeling Insights

### General
#### [Blenkinsop _et al._ 2022](https://elifesciences.org/articles/76487)

**_TA_**:: Clinical data of HIV VL and Tcell infxns to resolve the branch length/ancestral depth in the Coal Tree

#### [Lee _et al._ 2020](https://www.science.org/doi/10.1126/science.abd8755)

**_TA_**: Persepctive on various "engines" or factors driving SARS-Cov-2 transmission. Useful to think about predicitive variables. Type of transmission routes will differ depending on spatial scale and environmental conditions (e.g. household contacts, droplet vs fomite doesn't matter; aerosolized matters at community level; influenza in tropical regions more by fomites vs in low-humid regions, more by aerosolized)
Idea of overdispersion key at low incidence intersections. Takes very few inter-regional connections to produce "small world" network dynamics

### Phylodynamics
#### [Rasmussen _et al._ 2011/Nonlinear SIR Inference with Time Series](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3161897/)

**_TA_**: First example of a stochastic, mechanistic population dynamic models (with the bonus of being integrated with time series data). Framework uses state-space models (SSMs) for biological process/SIR mechanics.


### Hidden Geometries
#### [Brockmann & Helbing 2013](https://www.science.org/doi/10.1126/science.1245200)

**_TA_**: Canonical flu spread via airlines/networks manuscript. Hidden geometries and flux model


### Superspreading
#### [Loyd-Smith _et al._ 2005](https://www.nature.com/articles/nature04153)

**_TA_**:

### Dynamic Networks

#### [Volz & Meyer 2007](https://pubmed.ncbi.nlm.nih.gov/17878137/)

**_TA_**: Neighbor exchange model for including dynamicism into contact networks. 

#### [Volz 2008](https://link.springer.com/article/10.1007/s00285-007-0116-4)

**_TA_**: ... 





## Genomics
### Software
#### [Turakhia _et al._ 2021/UShER](https://www.nature.com/articles/s41588-021-00862-7)
**_TA_**: USHER algorithmm/framework is similar to a segregating sites MSA but uses precomputed data object storing the inferred histories of mutation events/seg sites on the tree (almost like a look up time table). Very fast and scalable to big data  

#### [Maio_et al._ 2016/SCOTTI](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5040440/)

**_TA_**: Uses the structured coalescence to overcome issues of within host evolution and non-sampled hosts. Models each host as it's own "deme" and transmission events are migration events between demes. Also brings up several transmission complexities (see figure 1)

### Methods/Methodological Insights
