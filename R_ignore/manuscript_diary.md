# Manuscript of Interest

Diary of manuscripts for Tri-CePPI project and `fomes` model.  
**_TA_**: are take away points for consideration or useful insights.



[TOC]




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



### [Cárdenas _et al._ 2023/Opqua](https://www.biorxiv.org/content/10.1101/2021.12.16.473045v2.full.pdf)

**_TA_**:  Fitness landscapes

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

**_TA_**: Neighbor exchange model for including dynamicism into contact networks. NE essentially leverages rewiring probability to toggle between a static and closer to dynamic network dynamics. 

#### [Volz 2008](https://link.springer.com/article/10.1007/s00285-007-0116-4)

**_TA_**: ... 

#### [Vernon & Keeling 2008](https://royalsocietypublishing.org/doi/10.1098/rspb.2008.1009)

**_TA_**: Dyanmic networks needed to capture nuances of disease spread. Final epidemic sizes have bimodal shape with dynamic network versus static 2/2 fact that sometimes a node this well connected to a "clique" is accessed during an infectious opportunity and sometimes it is not (infectious timing needs to be right to promote propagation). Interestingly, theoretically year-long weighted static networks do not follow the theoretical derivations given the fact that connections are not made or broken randomly but follow an underlying process related to seasonality, day of week, etc. that constitute complex behaviors: "positive and negative correlations [of connectivity]  at a range of temporal lags". Overall, caution against using approximations of the full dynamic network (as shown above, the intuitive weighted static network fails)

#### [Bansal et al. 2010](https://doi.org/10.1080/17513758.2010.503376)

**_TA_**: Review article. Interesting points: 

- acute infections w/in a small population appropriately captured by a static network --> acute means short lived (no time for dynamicism to take place).
- Prev work has use time-integrated networks, where "models consider all contacts that may have occurred within a relevant period of time (e.g. an average infectious period) and aggregate these contacts into a single static contact network." --> however, different dyanmicism can result in same static shapes (different processes result in same lower-dimensional representation)
- Justin's point about ground-hog day reiterated here: there are a number of "chance" encounters that happen in the majority of people's lives 
- Three types of extrinsic classes relevant to ID Epi that cause dynamicism: [1] extrinsic changes (demography: birth-death; rewiring; migration; socioeconomic changes); [2] pathogen mediated changes (infection -> immunity; change in behavior 2/2 symptoms = infected node changing behavior not susceptible, eg rabies makes more likely to be aggressive); [3] public health mediated changes (vaccination; treatment (antivirals, abx); avoidance behavior; closures and lockdowns)
- "The work of Colizza et al. [17] and Viboud et al. [61] has shown that global connectivity patterns drive the diffusion of infectious diseases on large scales. Most of these models, however, do not consider dynamic changes in individual-level network connectivity."
  - scales of connectedness 
- Social avoidance = edge rewiring not a random -->  "Using a pair-approximation framework, Gross et al. [27] show that such adaptive rewiring leads to modular network structure (consisting of highly intraconnected but only loosely interconnected sub-networks), a wide degree distribution and degree assortavity (a positive correlation in degree between pairs of linked nodes). This reveals two means by which avoidance behaviour can counterintuitively exacerbate epidemics: degree correlations can decrease the effectiveness of targeted vaccination; and high intraconnectivity of the susceptible cluster enables the persistence of epidemics even below the epidemic threshold"
- [Ferrari et. Al. 2006](https://royalsocietypublishing.org/doi/10.1098/rspb.2006.3636) "demonstrate that infection-acquired immunity due to a prior epidemic results in greater population protection in a heterogeneously structured population, while random vaccination leads to smaller epidemics in highly structured small-world populations.



#### [Read et al. 2008](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2607433/)

**_TA_**: Determined casual contact (conversations) and close contact (physical contact). Found that casual contacts were more transient/irregular and situation/location dependent (eg workspace) versus close contact occur more often/stable and are often situation/location dependent (eg home or social). Overall suggests that there is modularity in dynamic connections = "core groups". Although, arguably this modularity is in part due to the collection technique of diaries and focusing on the participant's connections (ego) and not all connections ego+alter. _Suggest that daily contacts do not have a long tailed distribution and therefore are not appropriately captured by the power law (as has been claimed in sexual contact networks)_ 





## Genomics
### Software
#### [Turakhia _et al._ 2021/UShER](https://www.nature.com/articles/s41588-021-00862-7)
**_TA_**: USHER algorithmm/framework is similar to a segregating sites MSA but uses precomputed data object storing the inferred histories of mutation events/seg sites on the tree (almost like a look up time table). Very fast and scalable to big data  

#### [Maio_et al._ 2016/SCOTTI](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5040440/)

**_TA_**: Uses the structured coalescence to overcome issues of within host evolution and non-sampled hosts. Models each host as it's own "deme" and transmission events are migration events between demes. Also brings up several transmission complexities (see figure 1)

### Methods/Methodological Insights
