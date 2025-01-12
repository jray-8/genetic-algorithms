from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TypeVar, Callable, Union, Any

from .chromosomes import Chromosome, CrossoverType, NumericChromosome

# Date: December 17, 2024
# Author: Jeffrey Ray
# 
# Description:
#	This module defines the Organism class, which holds the chromosomes and traits of a species at the class-level,
# 	and allows for the creation of individuals (instances) from the species.

__all__ = ['create_species', 'Organism']

T = TypeVar('T')

class Organism(ABC):
	''' A generic class to create a new species from.

		This species will be monoploid (only 1 unique set of chromosomes).

		Must implement (class-level):
		- `chromosome_count (int):`					number of unique chromosomes
		- `chromosome_types (tuple[type]):`			the class of each unique chromosome
		- `maximize (bool):`						maximization problem, else min.

		Instance attributes:
		- `genome (tuple[Chromosome]):`				a list of chromosomes that together represent a unique solution
		- `fitness (int, float):` 					the health value of the individual
	'''
	__slots__ = ('genome', '_fitness') 			# used by instances of the class
	
	chromosome_count: int						# the number of unique chromosomes an organism has (haploid number)
	chromosome_types: tuple[type[Chromosome]] 	# list of chromosome classes
	maximize: bool								# greater fitness score is healthier

	def __init__(self, genome: tuple[Chromosome] = None):
		if genome is not None:
			self.genome = genome 
		else:
			self._create_genome()
		try:
			self._fitness = self.evaluate_fitness()
		except Exception as e: # not convertible
			raise TypeError(f'evaluate_fitness should return type float or int, not {type(self._fitness)}') from e
	
	def clone(self):
		''' Creates a new Organism with identical chromosomes. '''
		return self.__class__(genome=self._copy_genome())
	
	def _copy_genome(self):
		''' Clones all of this individual's chromosomes. '''
		return tuple(chromosome.clone() for chromosome in self.genome)
	
	def _create_genome(self):
		''' Constructs a new genome for this individual. 
			- uses the chromosomes' generate method
		'''
		self.genome = tuple(chrom() for chrom in self.chromosome_types)

	@property
	def fitness(self):
		return self._fitness

	@abstractmethod
	def evaluate_fitness(self) -> Union[int, float]:
		''' Provide this individual a fitness score based on the health of the solution it represents. '''
		return
	
	def healthier(self, value: Union[Organism, int, float]):
		''' True if self is 'more fit' than some other organism (or fitness score). '''
		if isinstance(value, Organism):
			score = value._fitness
		else:
			try:
				score = float(value)
			except Exception as e:
				raise TypeError(f'value must be of type (Organism, int, or float), not {type(value).__name__}') from e
		# compare
		if self.maximize:
			return self._fitness > score
		else:
			return self._fitness < score

	#region For Procreation: 
	def _crossover(self: Organism, other: Organism, n: int) -> tuple[Organism]:
		''' Cross over the (haploid) chromosome set of this individual with another.
			- Returns a tuple of n new genomes
		'''
		# prepare storage for crossed chromosomes
		chromosomes = [[] for _ in range(self.chromosome_count)] # n of each type

		# perform crossover for each chromosome type
		for i in range(self.chromosome_count):
			chromosomes[i].extend(self.genome[i].crossover(other.genome[i], n))

		# create new organisms
		# zip(*chromosomes) produces a list of genomes (complete sets of chromosomes)
		zygotes = tuple(type(self)(genome=complete_set) for complete_set in zip(*chromosomes))
		return zygotes

	def mutate(self):
		''' Mutates this organism (in-place).'''
		for chrom in self.genome:
			chrom.mutate()

	def procreate(self, other: Organism, num_children: int) -> tuple[Organism]:
		''' Create n new children by crossing over the parents' chromosomes, 
		
		then mutating the resulting childrens' genes.
		'''
		cls = type(self)
		if type(other) is not cls:
			raise TypeError(f'parents must be from the same species, but received {cls} and {type(other)}')
		zygotes = self._crossover(other, num_children)
		for z in zygotes:
			z.mutate()
		return zygotes
	#endregion

	def __str__(self):
		''' String representation of this unique individual, showing: 
			- all encoded chromosomes
			- fitness score
		'''
		pad = 15
		fitness_tag = f'{"Fitness:":<{pad}}'
		s = ''
		for i, chrom in enumerate(self.genome):
			name = f'Chromosome {i+1}' if chrom.name is None else chrom.name
			name += ':'
			s += f'{name:<{pad}} {chrom}\n'
		s += f'{fitness_tag} {self._fitness:.6g}' # use up to 6 significant digits
		return s

	@classmethod
	def describe(cls):
		''' Prints a string displaying unique features of this Species. '''
		s = f'Species: {cls.__name__}\n'
		s += f'Number of chromosomes: {cls.chromosome_count}\n'
		if cls.maximize:
			s += 'Maximization problem - higher fitness levels are healthier.'
		else:
			s += 'Minimization problem - lower fitness levels are healthier.'
		for i, chrom_type in enumerate(cls.chromosome_types):
			s += '\n\n'
			name = f'Chromosome {i+1}' if chrom_type.name is None else chrom_type.name
			s += f'{name}\n'
			s += f'Number of genes: {chrom_type.length}\n'
			s += f'Gene pool: {chrom_type.gene_pool}\n'
			if issubclass(chrom_type, NumericChromosome):
				t = 'integers' if chrom_type.is_integer else 'real numbers'
				s += f'  Set: {t}\n'
				s += f'  Mutation std dev.: {chrom_type.mutation_std_dev:.2f}\n'
			s += f'Crossover type: {chrom_type.crossover_type.name}\n'
			if chrom_type.crossover_type == CrossoverType.MULTI_PT:
				s += f'  Cross points (k): {chrom_type.k_crosses}\n'
			s += f'Mutation rate: {chrom_type.mutation_rate:.2f}'
		print(s)

	# -- treat as container (for chromosomes)
	def __getitem__(self, i):
		return self.genome[i]
	
	def __setitem__(self, i, value):
		self.genome[i] = value

	def __len__(self):
		return len(self.genome)
	
	def __iter__(self):
		return iter(self.genome)
	
	def __contains__(self, value):
		return value in self.genome
	
	def __delitem__(self, i):
		raise AttributeError('cannot delete chromosomes')
	
	def __eq__(self, other: Organism):
		return all(self.genome[i] == other.genome[i] for i in range(self.chromosome_count))



# factory function to subclass Organism to a specific Species
def create_species(
	name: 					str, 
	fitness_test: 			Callable[[Organism], Union[int, float]], 
	_chromosome_types: 		Union[list[type[Chromosome]], type], 
	_chromosome_count: 		int = 1, 
	_maximize: 				bool = True, 
	**kwargs: 				Any
) -> Organism:
	''' Create a new species, an implementation of the Organism class.

		Specify all unique properies of this species directly.

		- `fitness_test(self)` must return a score based on the chromosomes in self.genome

		- `_chromosome_types` is a list of Chromosome classes,  
		or a single class that is used for all chromosomes

		- `_chromosome_count` is set automatically unless using a single class

		- `_maximize` indicates greater fitness is better

		Set custom class-level variables with **kwargs.
	'''

	# validate args -- convert to correct type if possible
	name = validate_species_name(name)
	validate_fitness_test(fitness_test)
	_chromosome_count = validate_chromosome_count(_chromosome_count)
	_chromosome_types, _chromosome_count = validate_chromosome_types(_chromosome_types, _chromosome_count)

	class Species(Organism):
		chromosome_count = _chromosome_count
		chromosome_types = tuple(_chromosome_types)
		maximize = _maximize
		evaluate_fitness = fitness_test

		def __str__(self):
			cls = self.__class__
			return f'Species: {cls.__name__}\n' + super().__str__()
		
	# Set arbitrary class-level variables
	for key, value in kwargs.items():
		setattr(Species, key, value)
		
	Species.__name__ = name
	return Species


## validators
def validate_species_name(name) -> str:
	try:
		string = str(name)
	except Exception as e:
		raise TypeError('failed to convert species name to type str') from e
	return string

def validate_fitness_test(fitness_test):
	if not callable(fitness_test):
		raise TypeError(f'fitness_test must be a function (Organism) -> fitness score, not {type(fitness_test).__name__}')
	
def validate_chromosome_count(chrom_count) -> int:
	try:
		chrom_int = int(chrom_count)
	except Exception as e:
		raise TypeError(f'_chromosome_count must be type int, not {type(chrom_count).__name__}') from e
	if chrom_int < 1:
		raise ValueError(f'Organisms must have at least 1 chromosome')
	return chrom_int

def validate_chromosome_types(chrom_types, chrom_count) -> tuple[list[type[Chromosome]], int]:
	''' Returns a valid tuple, (chromosome_types, chrom_count)'''

	type_error = TypeError('chromosome types must be subclasses of Chromosome')

	# received single chromosome type
	if isinstance(chrom_types, type):
		if not issubclass(chrom_types, Chromosome):
			raise type_error
		# create multiple references to the sole chromosome type
		chrom_types = [chrom_types] * chrom_count
		
	else: # treat as iterable
		for chrom in chrom_types:
			if not issubclass(chrom, Chromosome):
				raise type_error
		# set chromosome_count automatically
		chrom_count = len(chrom_types)
	return chrom_types, chrom_count
	