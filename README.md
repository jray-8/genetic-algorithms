# __Genetic Algorithms (GA)__

A genetic algorithm is a type of search and optimization technique inspired by the process of natural selection from biological evolution. GAs are used to find near-optimal solutions to problems where brute-force or direct methods are infeasible. 

They operate on a population of potential solutions, iteratively improving them through means of selection, crossover and mutation, where better solutions pass on their traits to the next generation, while weaker ones are removed.

GAs are powerful tools designed to find solutions without exact knowledge of the problem. As long as you have a way of symbolically representing solutions to a problem, and understand what makes one solution better than another, the population will converge towards better approximations of the answer.

Consequently, GAs are well-suiting for:
- Opitmization problems
- Complex search-spaces
- Constraint satisfaction
- Combinatorial problems
- Problems that have no closed-form solutions

Every problem we solve using a genetic algorithm will require a similar setup before we can begin writing the evolution lifecycle.


## __Key Components:__

1. __Representation:__
A way to encode solutions to the problem.
	- binary strings
	- integers
	- real numbers
	- arrays of symbols or objects
	- permutations (lists of indices)

1. __Fitness Function:__
The function which evaluates how good a given solution is.

1. __Crossover Function:__
The method to combine two "parent" solutions to create "offspring".

1. __Mutation Function:__
Introduce random changes in the offspring to promote diversity.
	- injects new genetic information into the population
	- search different parts of the solution space

1. __Population Dynamics:__  

	__Population Size ($N$):__
	How many solutions will be tracked at once?  
	- the maximum capacity of the population
	- this will dictate the amount of memory used by the algorithm

	__Birth Rate ($n$):__
	The number of new solutions added per generation.

	__Number of Parents ($p$):__
	How many solutions are combined to produce offspring each generation.
	- this adds diversity to the solutions,
	- but may slow down convergence to what's optimal

1. __Initial Population:__
How the starting population is generated.
	- randomly
	- biased
	- selected archetypes

1. __Selection Function:__
How individuals are chosen for reproduction.
	- rank selection (choose the strongest)
	- tournament selection (best from a random subset)
	- roulette wheel (random with biased towards greater fitness)
	- random selection
	- mate strongest with randomly selected (for diversity)

1. __Replacement Strategy:__
Decide which individuals stay in the population for the next generation.
	- generational replacement (entire population replaced)
	- steady-state replacement (choose only a small number of individuals each generation)
	- elitist selection (keep the strongest)
	- dynamic replacement (start off by only replacing a few, and gradually replace more)
	- favor the removal of similar individuals (crowding)
	- maintain multiple subpopulations (niching)

1. __Termination Condition:__
When the algorithm should stop.
	- maximum number of generations
	- strongest member reached a satisfactory fitness level
	- _golden child_ is found (ideal solution)
	- convergence has plateaued
	- time elapsed


## __Evolution Lifecycle__

Every GA will follow this basic routine.

1. __Initialize Population__
	- Generate the first generation of candidate solutions

1. __Evaluate Starting Population__
	- Calculate the fitness of every member in the population and rank them

1. __Selection__
	- Select individuals from the population, based on their fitness, for reproduction

1. __Crossover (Recombination)__
	- Perform crossover on the selected parents to create offspring, combining their genetic material

1. __Mutation__
	- Introduce random changes to the offspring's genes to add new genetic information to the population

1. __Evaluate Fitness__
	- Assess the health of the offspring and rank the new population

1. __Replacement__
	- Replace members of the old population with the new children

1. __Check Termination Conditions__
	- Stop the evolution loop, if a condition is met

1. __Repeat__  
	- If no termination condition is met, go back to __step 3__ and continue the algorithm


## __Glossary__

> Here are some common terms used in GAs that are derived from the biology of the natural world.

__Genetics:__
The scientific study of genes, heredity, and variation in living organisms—forming the basis for how features of a species evolve over time.

__Heredity:__
The biological processes by which physical or mental characteristics are passed on from parents to their offspring.

__Species:__
A group of organisms that can interbreed and produce fertile offspring under natural conditions. Members of a species are genetically similar to each other, and are not typically capable of reproducing with members of another species.

__DNA:__
Deoxyribonucleic acid—the molecular structure that encodes genetic instructions for building and maintaining living organisms, often referred to as the _blueprint of life_.

> Think of DNA as the low-level building material for other genetic structures, like bricks to build houses.

__Gene:__
A functional sequence of DNA that contains the code to produce a specific protein or determine a particular trait.

> A gene can then be imagined as a singular house made of DNA.

__Chromosomes:__
A structured package of genes and non-coding regions, comprised of tightly-coiled strands of DNA. Chromosomes organize the genes that exist within an organism in a specific order. Each species typically has a fixed number of unique chromosomes. These chromosomes form a _set_, and most cells of an organism (except for _gametes_ like sperm or egg cells) will contain an exact copy of this set of chromosomes.

> If a gene is a house made of DNA (bricks), then a chromosome is like a street of houses.

__Gamete:__
A reproductive cell of a plant or animal.

__Ploidy:__
An organism may have more than one complete sets of chromosomes in each cell. For instance, humans have 23 chromosome pairs. This means humans have 23 unique chromosomes, but have 2 of each kind—one from the mother and the other from the father. Ploidy refers to the number of complete sets of chromosomes in a cell, or the cells of an organism.

__Diploid (2n):__
A diploid cell, organism or species is one that contains two complete sets of chromosomes (one from each parent). Humans are diploid.

__Haploid (n):__
Describes a cell that contains only a single set of chromosomes. Human gametes are haploid.

__Locus:__ (Plural: Loci.)
The fixed position on a chromosome where a particular gene (or other DNA sequence) is located.

> Like a home address for a gene.

__Genome:__
The complete collection of genetic material within an orgianism, including all of its chromosomes (which holds all its DNA). The _genetic library_ of an individual.

> If a chromosome were a street of houses (genes), then a genome would be a self-functioning city of streets and all the infrastructure necessary to keep the city running smoothly.

---
<br>

__Allele:__
A specific variant of a gene, contributing to the diversity of genetic traits in a population. Alleles of a particular gene and are found at the same locus.

> __Ex.__ &nbsp; 
A gene that determines the color of a variety of _Oriental Lily_ might have one allele for red petals `R`, and another allele for white petals `r`.

__Genotype:__
The combination of alleles an organism inherits for a specific gene (or multiple, depending on the context).

> __Ex.__ &nbsp; 
For the gene controlling petal color, one lily may have the genotype `Rr`.  
Another lily may have `RR`, and another still may have `rr`.

__Phenotype:__
The observable trait or characteristic expressed by an organism, as a result of its genotype and its interaction with the environment. 

> __Ex.__ &nbsp; 
A lily with the _genotype_ `Rr`, which codes for its color, may appear pink. The lily's _phenotype_ is then pink.

__Dominance:__
When one allele masks the effect of another it is _dominant_, and the masked allele is _recessive_.

> __Ex.__ &nbsp; If we have two alleles for eye color:
>
> - `B` for brown (dominant)
> - `b` for blue (recessive)
>
> Then, we have the following genotypes:
>
> - `BB` $\rightarrow$ brown
> - `Bb` $\rightarrow$ brown (`B` dominates)
> - `bB` $\rightarrow$ brown (`B` dominates)
> - `bb` $\rightarrow$ blue

__Incomplete Dominance:__
Neither allele is completely dominant. Instead, the phenotype is a blend of both alleles.

> The lily with the `Rr` genotype, that had pink petals instead of red or white, was an example of incomplete dominance.
>
> The conditions of the lily's environment may also affect its color.  
> For instance, the `Rr` lily may only express pink petals in warm climates, and appear red in cooler climates.

__Co-dominance:__
Both allele variants are fully expressed in the phenotype, without blending.

> __Ex.__ &nbsp; 
>
> - A chestnut horse (`HH`) has reddish-brown hair all over
> - A cremello horse (`hh`) has cream colored hair all over
>
> A chestnut and cremello may breed to produce a _palomino_.
> - A palomino horse (`Hh`) has the same body color as a chestnut horse, and the same mane and tail color as a cremello horse


__Homozygous__ – having two _identical_ alleles for a particular gene.

> `AA` or `aa`

__Heterozygous__ – having two _different_ alleles for a particular gene.

> `Aa`

__Overdominance:__
The heterozygous genotype has a trait _advantage_ over either homozygous genotype.

> __Ex.__ &nbsp; Sickle cell anemia
>
> - `AA` $\rightarrow$ No sickle cells, no malaria resistance
> - `Aa` $\rightarrow$ Some sickle cells, malaria resistance (best trait)
> - `aa` $\rightarrow$ Sickle cell disease, malaria resistance

__Underdominance:__
The heterozygous genotype has a _disadvantage_ compared to either homozygous genotype.

> __Ex.__ &nbsp; Let's say we have a population of fish where:
>
> - `BB` have bright, colorful scales making them attractive to mates.
> - `bb` have dull scales which aren't appealing to mates, but allow them to blend into their environment, avoiding predators.
>
> - Then we have `Bb` fish that get the worst of both worlds:
>
> 	- not bright enough to attract mates like the `BB` fish
>	- but not dull enough to avoid predators like the `bb` fish  
>   
> They are competing at a disadvantage in each category: _mating_ and _survival_.

__Epigenetics:__
The study of changes in gene expression that occur without alterations to the underlying DNA sequences. Such changes to an organism may result from its environment, lifestyle, aging, diet, chemical exposure, drug use, psychological states, or its own hormones.

> __Ex.__ &nbsp; 
Exposure to freezing temperatures triggers an epigenetic reaction in wolves, allowing them grow thicker fur coats to withstand the cold.

__Eugenics:__
A set of beliefs and practices aimed at improving the genetic quality of the human population by selective breeding or genetic manipulation.

---
<br>

__Autosome:__
One of the numbered, non-sex chromosomes. Humans have 22 pairs of autosomes and 1 pair of sex chromosomes.

__Sex Chromosome:__
A chromosome involved in determining the sex of an organism. Either `XX` (girls) or `XY` (boys) in humans.

__Homologous Chromosomes:__
A pair of chromosomes with the same structure, carrying genes for the same traits, but possibly different alleles. Humans have 23 pairs of homologous chromosomes because each chromosome in a pair, one from each parent, is of the same type.

__Sister Chromatids:__
Each of the two identical halves of a chromosome that have been replicated in preparation for cell division.

[__Mitosis:__](https://www.yourgenome.org/theme/what-is-mitosis/)
The process of cellular division by which a parent cell is split into two genetically identical daughter cells.

[__Meiosis:__](https://www.yourgenome.org/theme/what-is-meiosis/)
The process of cellular division by which gametes (sperm or egg cells) are created.

> In humans, 
>
> - A diploid _germ_ cell randomly splits aparts its chromosome pairs to produce two new sets of chromosomes.  
> 	Each chromosome from a pair can go to either set, so there are $2^{23}$ ways to form these two groups.
>
> - Some homologous chromosomes between these two groups may cross over sections of one of their chromatids
>
> - The two groups break off into two individual haploid cells.
>
> - The haploid cells each divide in two again by splitting their chromatids, producing a total of four unique haploid gametes for sexual reproduction.

__Crossover:__
The exchange of genetic material between homologous chromosome pairs during meiosis. This results in new combinations of alleles in gametes.

> `10|1101`  
> `01|1111`  
>
> Crossover yields:
>
> `10|1111`  
> `01|1101`

__Mutation:__
A change in the DNA sequence of a gene, which can result from errors in DNA replication, environmental factors, or random events. Mutations can be advantageous, neutral, or adverse. Mutations occur more frequently during meiosis compared to mitosis.

__Zygote:__
A single-celled, diploid organism formed by the fertilization of an egg (ovum) with a sperm cell. The zygote contains a full set of chromosome pairs, exactly half from each parent, forming the genetic blueprint for the entirety of the organism.


---
<br>

__Society:__
A _population_ interacting under a set of organized rules.

> __Ex.__ &nbsp;
A pack of wolves follows a strict social structure with an alpha pair leading the pack who determine hunting, resting and movement patterns.

__Population:__
A collection of organisms of the same species that live in the same area, at the same time, and are capable of breeding with each other.

> __Ex.__ &nbsp;
A population of flamingos inhabits the same saltwater lagoon during the mating season.

__Ecosystem:__
A geographical area of interacting organisms and their physical environment. This includes both living (biotic) and non-living (abiotic) entities, as well as the processes that exists between them to support life.

> __Ex.__ &nbsp;
A coral reef is an ecosystem where fish, coral, algae and other sea creatures interact with water currents, sunlight, minerals and each other to sustain life.

__Panmixia:__
Uniformly random mating within a breeding population.

> __Ex.__ &nbsp;
Many coral species release eggs and sperm into the ocean currents, leading to completely random fertilization.

__Gene Pool:__
The total stock of unique alleles in an interbreeding population.

> __Ex.__ &nbsp; 
In a large population of butterflies, the gene pool consists of the alleles for all their different colors, wing patterns, sizes, flight speed, lifespans, etc.

__Genetic Drift:__
The process by which allele frequencies in a population change due to random events. This is especially significant in small populations, and can lead to an allele's extinction over time.

> __Ex.__ &nbsp;
 A volcanic eruption isolates a small group of lizards on an island. Over generations, their allele frequencies differ significantly compared to the original population on the mainland.

__Evolution:__
The gradual process by which populations of organisms change over generations, driven by genetic variations that are passed on to offspring. Factors such as crossover, mutation, and _natural selection_ shape these changes, ultimately influencing the _gene pool_ of the entire population.

> __Ex.__ &nbsp; 
Over millions of years, species of finches evolved distinct beak shapes to help them access different food sources available in their respective environments:
>
> - Long, pointed beaks were useful for snatching insects.
> - Short, blunt beaks worked better for cracking seeds.

__Adaptation:__
A change by which an organism becomes better suited to survive its environment. Adaptation can be considered the localized and immediate version of evolution. It focuses on the actions or behaviour of an _individual_, whereas evolution is deals with the large-scale changes to the overall gene pool of a species.

> __Ex.__ &nbsp; 
In the extreme cold of the Antarctic isles, penguins have learned to huddle to conserve heat. Individual penguins take turns moving from the outer edges of the group to the center where its warm.
>
> This is an _adaptation_ because it's an immediate, behavioural strategy that penguins use to survive their harsh environment.

__Fitness:__
A quantitative score of an organism's ability to produce healthy offspring.

> __Ex.__ &nbsp;
Himalayan hares with white fur that blends more effectively into their snowy environment have higher fitness than hares with more detectable fur.

__Survival of the Fittest:__
The notion that organisms best adapted to their environment are more likely to survive, reproduce, and pass on their profitable traits to the next generation.

> __Ex.__ &nbsp;
In arid deserts, only lizards with the ability to conserve water through their thick skin and minimal sweating survive long enough to reproduce.

__Natural Selection:__
The process by which the organisms more fit for an environment survive while the weaker ones die off, resulting in the alleles of these organisms becoming more common in future generations. This is the core process that brings about evolution.

> __Ex.__ &nbsp; 
Giraffes with long necks could reach the leaves of tall acacia trees, giving them a survival advantage over giraffes with shorter necks. Over time, the alleles responsible for long necks became prominent in the giraffe gene pool.

__Selective Pressure:__
Any external factor that affects the survival and reproduction of organisms with a certain phenotype. Selective pressures can be environmental (such as climate), biological (like predators), or social (like mating preferences).

> __Ex.__ &nbsp; 
A selective pressure for faster running speed in gazelles occurs due to the presence of predators like cheetahs, leading to the evolution of faster gazelles over time.

__Competition:__
The interaction between species, or organisms within a species, where both require control over one or more resources that are limited in supply. Competition reduces the fitness all members involved, because they cannot all successfully coexist.

> __Ex.__ &nbsp;
In the African savanna, lions and hyenas compete for overlapping prey such as wildebeest.

__Carrying Capacity:__
The maximum population size of a species that can be sustained by a given environment due to constraints such as food, water and shelter. Exceeding the carrying capacity can lead to resource depletion, environmental degradation, and a population crash.

__Overcrowding:__
A condition where a population exceeds its carrying capacity, resulting in competition for resources, reduced health, population decline, and intensified selective pressure.

> __Ex.__ &nbsp;
A lake has a carrying capacity of 500 fish due to its limited algae and oxygen. If the population of fish exceeds this number, starvation and die-offs will occur from overcrowding. Natural selection may then favor fish who require less oxygen or smaller fish who depend on less algae.


## __GA-Specific Concepts__

__Elitism:__  
The existence of a dominating group in society that has superior qualities, such as intellect, wealth, strength, charisma, attractiveness, craftsmanship, or lineage.

In genetic algorithms, elitism is the practice of preserving the top-ranking individuals from one generation directly to the next generation without modification. This ensures the best solutions are not lost to selection or mutation, yielding a non-decreasing pattern of gene pool quality. However, relying too heavily on elitism can reduce genetic diversity and cause the algorithm to converge prematurely on suboptimal extrema.

__Premature Convergence:__  
An unwanted effect in GAs where the population becomes similar too quickly, reducing the algorithm's ability to explore new areas in the search-space, which may contain better solutions. Premature convergence stalls the evolutionary process on a local optimum, as it fails to find the best global genetic code. Diversity-promoting mechanisms such as selective mutation, incest-prevention, larger populations, segregating individuals with or _niching techniques_ can reduce the risk of this issue.

Premature convergence is biologically similar to inbreeding exposing genetic faults in offspring, or a lack of genetic diversity making a population less adaptable to a changing environment.

__Exploration vs. Exploitation:__  
Exploration refers to searching through the problem space for new and diverse solutions, whereas exploitation involves refining and optimizing solutions that are already known to be good. Finding the optimal balance between these two strategies is central to the success of a genetic algorithm.

> - Elitist selection and crossover contribute to __exploitation__, as they create offspring similar to the best individuals in society.
>
> - Panmictic selection and mutation contribute to __exploration__ of the problem space, as they promote diversity.

If the GA explores too heavily, then the population will never converge to an extrema (best fitness). However, without exploring enough the population will lack genetic diversity, opening the door to premature convergence. Similarly, over-exploiting the best members of society will kill off diversity, while making no effort to capitalize off of the elite will prevent local refinement. There is an apparent trade-off between the two strategies, and maximizing either one would prevent the algorithm from working—leading to either no convergence, or suboptimal convergence.

__Niching & Crowding:__  
These techniques are used to reduce the exploitation of existing solutions by maintaining diversity in a population. They are especially useful in multi-modal optimization problems, where you want to preserve multiple different solutions (or solution components) without letting one dominate the population. Both approaches aim to prevent premature convergence and explore various solution forms in parellel, though they achieve this in different ways.

__Crowding:__  
The situation where individuals in a species with similar traits compete directly with each other for limited resources or opportunities to reproduce. Crowding can lead to increased competition, aggression, reduced mating availability, and more selective pressure to evolve—which individuals will win the riches and prosper?

In GAs, crowding involves limiting the number of similar individuals that contribute to the next generation. This can be implemented by introducing a _crowding penalty_, which reduces the fitness of individuals with similar traits or gene patterns. The penality would favor the replacement of similar solutions with new offspring, reducing their numbers, thus reducing the penality (as crowding diminishes). 

Another approach involves finding solutions with the same (or similar) fitness values, and only choosing one of them for selection, reducing their mating ability.

Crowding focuses on eliminating redundant solutions by forcing competition among similar individuals, preventing any one archetype from dominating the population. Organisms with differing traits are designed to rule the population at all times. Crowding is best suited for cases where exploring a wide range of potential solutions is critical.

> __In Short__
> - Keeps only the strongest from groups of similar individuals
> - Maintains variety of unique solutions without a central convergence
> - Best for finding _all_ unique, strong solutions

__Niching:__  
In the wild, a niche refers to the specific role adapted to by a species (or individual) to avoid competition with others an ecosystem. For instance, consider a forest where tall trees compete for sunlight up in the canopy, while smaller shrubs thrive in the shade below. Both the tall trees and the shrubs occupy a different niche in the forest, allowing for coexistence without rivalry for sun exposure.

In GAs, niching is the strategy of maintaining distinct subpopulations (or _niches_) in a broader population, each focusing on different regions of the solution space. This technique preserves diversity by allowing multiple groups to evolve independently, leading to high-quality archetypes (local optima) for each group. Each such archetype then represents a possible solution to the problem with an additional constraint or specialization. Niching is particularly useful for problems with discrete forms or finite categories of solutions.

> __In Short__
> - Preserves subgroups with gradual variation
> - Encourages separate convergence for each niche
> - Best for problems with discrete forms

__Parameter Tuning:__  
Refers to adjusting the key variables and procedures of a genetic algorithm (such as population size, mutation rate, crossover method, etc.) to optimize its ability to converge to good solutions. This is about balancing exploration and exploitation, while also managing memory usage and execution time. A good understanding of the question being answered is imperative to choosing the best parameters for a specific GA. 

> __Population Size:__
>
> Generally, maintaining a larger population will find optimal solutions faster, as there are more options to choose from, but will require more memory to run. 
>
> __Random Events:__
>
> Increasing random events like mutation rate will lead to deeper exploration but can make it difficult for the GA to settle on a final destination. On the other hand, a low probability of random events will cost the algorithm more generations, and thus processing time, to find desired solutions.
>
> __Dynamic Parameters:__
>
> Chaning how the GA operates as generations advance may allow the algorithm to make best use of memory and efficiency at different stages of the population's development.
>
> For instance,  
> - __Early Stages:__ 
> 	- Start with a smaller population and focus on exploration.
> 	- Employ higher mutation rates and use niching to diversify the population
>
> - __Mid to Late Stages:__
> 	- Gradually increase population size and reduce mutation rate to converge towards promising archetypes.
>
> - __Endgame:__
> 	- Focus on heavy exploitation of the best members of society.
> 	- Pool all niches together and enforce elitism to consolidate the archetypes into the absolute best.
>
> __Tuning Strategies:__
> 
> - Test and tweak parameters on sample problems (with smaller solution spaces) to find optimal configurations before deploying the algorithm on larger data sets.
>
> - Train machine-learning models to predict optimal parameter settings across multple problems from the same domain.


## __Evolutionary Algorithms (EA)__

Genetic algorithms fall into a broader category known as _evolutionary algorithms_, where complex behaviour or problem-solving ability is built incrementally by adopting biological principles like natural selection.

Here are some examples of non-GA EAs:

- __Genetic Programming (GP)__
	- Like genetic algorithms but entire computer programs or mathematical equations (represented as trees) are evolved instead of strings.

- __Memetic Algorithms (MA)__
	- A hybrid of genetic algorithms and local search techniques.

	- In each generation, solutions are evolved using genetic algorithms, but once the children are created, they are further refined to improve smaller parts of the solution (local optimization).

- __Cellular Evolutionary Algorithm (cEA)__
	- Inspired by cellular automata (CA), cEAs model problems through a network of interacting agents (cells) that evolve and propagate across a grid.

	- Each agent represents a potential solution that is influenced by its local environment (neighbouring cells).

- __Differential Evolution (DE)__
	- An EA for solving real-valued optimization problems by iteratively improving vectors.

	- It does not require the problem to be differentiable, making it useful for problem spaces that are noisy, non-continuous, or dynamic problems.

- __Artificial Bee Colony (ABC)__
	- Inspired by the foraging behaviour of the honeybee swarm, an ABC simulates the search for optimal solutions by emulating the way bees explore food sources.

	- The ABC model consists of 3 groups of bees: employed bees, onlookers, and scouts.

	- __Employed bees__ search and refine existing solutions.

	- __Onlooker bees__ evaluate and select solutions based on information shared by the employed bees.

	- __Scout bees__ explore new areas of the solution space.

	- The locations of the food sources represent possible solutions, and the amount of nectar in them represents the fitness of that solution.


- __Particle Swarm Optimization (PSO)__
	- Inspired by the social behaviour of flocks of birds or schools of fish, PSO uses "particles" (solutions) that move through the search-space with a position and velocity.
	
	- Each particle's movement is influenced by its personal best-known position, but is also guided towards the best-known positions of the search-space, which can be updated other particles.

- __Cuckoo Search Algorithm (CSA)__
	- Inspired by the brood parasitism of cuckoos, who lay their eggs in the nests of other birds, tricking them into raising their cuckoo chicks, CSA solutions are modelled as "eggs laid in a nest."

	- The algorithm explores the search-space by replacing the least fit eggs with random cuckoo eggs.

	- The eggs are grouped in nests, each managed by a host bird who has a chance to discover a cuckoo egg. When this happens, the host may either discard the egg or abandon the nest and start a new one, leading to the creation of a new random solution.

- __Neuroevolution__
	- Uses the evolution lifecycle to optimize the structure and weights of artificial neural networks (ANNs).
