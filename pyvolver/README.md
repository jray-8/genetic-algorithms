## Pyevolver
A general-purpose template for quickly adapting a genetic algorithm to suit to your problem.

## How to Use:

- _factory functions_ will allow you to create new types of chromosomes and species
- chromosomes and organisms (of a species) will have all their parameters stored at the class level
- instances, then, could be created with no args
- objects will be lightweight and store only what sets them apart from others
- use the `dir()` and `help()` functions to learn more about any classes

#### 1. Create chromosome types for your species using the *factory functions*:

A new type which implements the Chromosome abstract class is returned.
```
create_symbolic_chromosome(
	_length: int 						= 10, 
	_gene_pool: tuple[Any] 				= ('0','1'),
	_crossover_type: CrossoverType 		= CrossoverType.SINGLE_PT,
	_mutation_rate: float 				= 0.05,
	_name: str							= None,
	_k_crosses: int						= 2,
	generate_func						= None,
	mutate_func							= None,
	str_func							= None,
	**kwargs
) -> SymbolicChromosome:
```
```
create_numeric_chromosome(
	_length: int 								= 1, 
	_gene_pool: tuple[Union[int, float]] 		= (-100, 100),
	_crossover_type: CrossoverType 				= CrossoverType.AX,
	_mutation_rate: float 						= 0.50,
	_name: str									= None,
	_k_crosses: int								= 2,
	_is_integer: bool							= False,
	_mutation_std_dev: float					= 0.5,
	generate_func 								= None,
	mutate_func									= None,
	str_func									= None,
	**kwargs
) -> NumericChromosome:
```
```
create_permutational_chromosome(
	_gene_pool: tuple[Any] 				= ('A', 'B', 'C', 'D', 'E'),
	_crossover_type: CrossoverType 		= CrossoverType.OX2,
	_mutation_rate: float 				= 0.25,
	_name: str							= None,
	generate_func 						= None,
	mutate_func							= None,
	str_func							= None,
	**kwargs
) -> PermutationalChromosome:
```

In general:
```
def create_chromosome(
	_length: 					int,
	_gene_pool: 				list[Union[str, int, float]],
	_crossover_type: 			CrossoverType,
	_mutation_rate: 			float,
	chromosome_type: 			type[Chromosome],
	# keyword args
	_name: 						str = None,
	_k_crosses: 				int = 2,
	_is_integer: 				bool = False,
	_mutation_std_dev: 			float = 0.5,
	generate_func				= None,
	mutate_func					= None,
	str_func					= None,
	**kwargs
):
'''
	Factory function to create a new Chromosome class which implements chromosome_type.

	## Must specify:
	- _length (int): 						number of genes in chromosome
	- _gene_pool (list[str, int, float]):	the domain each gene can be encoded by
	- _crossover_type (CrossoverType):		which method of crossover this chromosome uses
	- _mutation_rate (float):				probability of independent gene mutation events

	## May specify:
	- _name (str):							name this unique chromosome
	- _k_crosses (int):						number of points for k-point crossover

	For Numeric Chromosome:
	- _is_integer (bool):					domain is integers or real numbers
	- _mutation_std_dev (float):			std dev. of the normally distributed mutation delta

	## May override:
	- generate_func(self)					overrides self._generate()
	- mutate_func(self)						overrides self.mutate()
	- str_func(self)						overrides self.__str__()

Set custom class-level variables with **kwargs.
'''
```

Chromosome objects store only `self._genes`, which is a list of its genes.

Genes can be accessed directly from the instance.  
For example:
``` 
chrom = my_chromosome_type()
x = chrom[2] # the second gene

# iterate over genes in multiple ways:
for c in chrom:
	print(c)

for i in range(chrom.length): # access class-level length
	print(chrom[i])

print(chrom[:]) # show entire list of genes
```


#### 2. Create a new species using the chromosomes you've created:
```
create_species(
	name: 					str, 
	fitness_test: 			Callable[[Organism], Union[int, float]], 
	_chromosome_types: 		Union[list[type[Chromosome]], type], 
	_chromosome_count: 		int = 1, 
	_maximize: 				bool = True, 
	**kwargs: 				Any
) -> Organism:

class Organism(ABC):
''' A generic class to create a new species from.

	This species will be monoploid (only 1 unique set of chromosomes).

	Must implement (class-level):
	- chromosome_count (int):					number of unique chromosomes
	- chromosome_types (tuple[type]):			the class of each unique chromosome
	- maximize (bool):							maximization problem, else min

	Instance attributes:
	- genome (list[Chromosome]):				a list of chromosomes that together represent a unique solution
		- chromosome (Chromosome): 				an encoded list of genes
	- fitness (int, float): 					the health value of the individual
'''
```
A new species type which implements the Organism abstract class is returned.

Species instances have only:
- `self.fitness` - a read-only property describing the health of the organism
- `self.genome` - a list of the organism's chromosomes

Its chromosomes could also be accessed directly from the instance:
```
individual = my_species()

x = individual[0] # its first chromosome object

for chrom in individual:
	print(chrom)

for i in range(individual.chromosome_count):
	print(individual[i])
```


#### 3. Create a population of individuals from your species:
Use the `Population` class:
``` 
def __init__(self,
	species: 			Organism,
	pop_size: 			int = 8,
	max_capacity: 		int = None, 		# pop_size by default
	num_offspring: 		int = 2, 			# must be chosen directly (>= 1)
	birth_rate: 		int = None, 		# sustains intitial pop_size by default
	num_parents: 		int = None,
	ideal_fitness:		float = None,
	evolution_options: 	EvolveFlags = 0
) -> Population:
```

#### 4. Evolve!
```
my_population.evolve(n=100)
x = my_pop.solution # strongest Organism in society
```