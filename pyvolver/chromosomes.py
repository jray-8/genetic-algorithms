from __future__ import annotations

import math
import random

from abc import ABC, abstractmethod
from enum import Enum, auto
from typing import TypeVar, Callable

from .helpers import restrict_bounds

# Date: December 24, 2024
# Author: Jeffrey Ray
# 
# Description: This module defines Chromosome types, their generation, crossover, and mutation

T = TypeVar('T')

class CrossoverType(Enum):
	''' Chromosomal crossover methods.

		For Symbolic / Numeric:
		- SINGLE 			- Single-Point Crossover
		- MULTI 			- Multi-Point Crossover
		- UX 				- Uniform Crossover
		- HUX 				- Half-Uniform Crossover
		- HDX				- Highly Disruptive Crossover

		For Numeric Only:
		- DC 				- Discrete Crossover
		- AX 				- Average Crossover
		- FC 				- Flat Crossover

		For Permutational Only:
		- OX1 				- Order-1 Crossover
		- OX2 				- Order-2 Crossover
		- PMX 				- Partially Mapped Crossover
		- CX 				- Cycle Crossover
	'''
	# symbolic
	SINGLE_PT = auto() # asdasd
	''' Single-Point Crossover '''
	MULTI_PT = auto()
	''' Multi-Point Crossover '''
	UX = auto()
	''' Uniform Crossover '''
	HUX = auto()
	''' Half-Uniform Crossover '''
	HDX = auto()
	''' Highly Disruptive Crossover '''

	# numeric
	DC = auto()
	''' Discrete Crossover '''
	AX = auto()
	''' Average Crossover '''
	FC = auto()
	''' Flat Crossover '''

	# permutational
	OX1 = auto()
	''' Order-1 Crossover '''
	OX2 = auto()
	''' Order-2 Crossover'''
	PMX = auto()
	''' Partially Mapped Crossover'''
	CX = auto()
	''' Cycle Crossover '''

	@classmethod
	def show_all(cls):
		''' Prints all available crossover choices. '''
		print(cls.__doc__)


class Chromosome(ABC):
	''' A base class for any type of chromosome. 

		Must implement:
		- `length (int)`
		- `gene_pool (tuple)`
		- `crossover_type (CrossoverType)`
		- `valid_crossovers (list[CrossoverType]):`		available crossover methods for subclasses
		- `mutation_rate (float)`

		May implement:
		- `name (str):`						name this unique chromosome
		- `k_crosses (int):`				number of points for k-point crossover
		
		Must override:
		- `_generate(self):` 				how the genes on this chromosome are created
		- `mutate(self):`					how mutation occurs
	'''
	length: 			int
	gene_pool: 			tuple
	crossover_type: 	CrossoverType
	valid_crossovers: 	list[CrossoverType]
	mutation_rate: 		float
	# -- optional --
	name: str = None
	k_crosses: int = 2
	# -- mutations --
	mutation_methods = [
		'_rotate_mutate',
		'_swap_mutate',
		'_reverse_mutate',
		'_scramble_mutate'
	]
	# -- crossover mappings --
	# each method must be: (other: Chromosome, n: int) -> list[Chromosome]
	crossover_map: dict[CrossoverType, str] = {
		CrossoverType.SINGLE_PT: '_single_point_crossover',
		CrossoverType.MULTI_PT: '_multi_point_crossover',
		CrossoverType.UX: '_uniform_random_crossover',
		CrossoverType.HUX: '_half_uniform_crossover',
		CrossoverType.HDX: '_highly_disruptive_crossover',
		CrossoverType.DC: '_uniform_random_crossover',
		CrossoverType.AX: '_average_crossover',
		CrossoverType.FC: '_flat_crossover',
		CrossoverType.OX1: '_order_1_crossover',
		CrossoverType.OX2: '_order_2_crossover',
		CrossoverType.PMX: '_partially_mapped_crossover',
		CrossoverType.CX: '_cycle_crossover'
	}

	# each chromosome instance will hold only a list of genes
	__slots__ = ('_genes')

	def __init__(self, genes: list=None):
		''' Create a new chromosome of this type. '''
		if genes is not None:
			self._genes = genes
		else:
			self._generate()
		
	def clone(self):
		''' Returns an identical copy of this chromosome object. '''
		return type(self)(self.length, self._genes.copy())
	
	def __repr__(self):
		name = self.name
		if name is None:
			name = type(self).__name__
		return f'{name}({self._genes})'
	
	def __str__(self):
		''' A string to display the genes on this chromosome. '''
		return ' '.join(str(gene) for gene in self._genes)

	@abstractmethod
	def _generate(self): # creates self._genes
		''' Initializes the genes on this chromosome. '''
		return
	
	@abstractmethod
	def mutate(self):
		''' Apply the chance to mutate any single gene on this chromosome.'''
		return
	
	def crossover(self, other: Chromosome, n: int) -> list[Chromosome]:
		''' Returns a list of (at least) n new chromosomes by performing crossover between the parents. '''
		# can assume type(self) == type(other)
		method_name = self.crossover_map.get(self.crossover_type)
		if method_name is None:
			raise ValueError(f'crossover_type: {self.crossover_type.name}, does not exist for {type(self).__name__}')
		crossover_method = getattr(self, method_name) # raises AttributeError if not found
		return crossover_method(other, n)
	
	# get the major type of this Chromosome class descendant
	@property
	def chromosome_type(self) -> type[Chromosome]:
		if isinstance(self, SymbolicChromosome):
			return SymbolicChromosome
		elif isinstance(self, NumericChromosome):
			return NumericChromosome
		elif isinstance(self, PermutationalChromosome):
			return PermutationalChromosome
		else:
			return 

	@classmethod
	def explain_crossover(cls, cross_type: CrossoverType):
		''' Prints information about the chosen crossover type. '''
		if not isinstance(cross_type, CrossoverType):
			raise TypeError(f'cross_type must be CrossoverType, not {type(cross_type)}')
		method_name = cls.crossover_map.get(cross_type)
		if method_name is None:
			raise ValueError(f'crossover: {cross_type.name}, not found in the crossover_map')
		# get method
		crossover_method = getattr(cls, method_name)
		if crossover_method.__doc__:
			print(crossover_method.__doc__)
		else:
			print(f'no documentation available for crossover: {CrossoverType.name}')
		
	#region mutation methods 
	# -- these methods randomize the genes on the chromosome (maintains counts of each)
	def _rotate_mutate(self):
		''' Rotates the permutation by a random number of positions: `1..(length-1)`. '''
		n = random.randrange(1, self.length)
		self._genes = self._genes[n:] + self._genes[:n]

	def _swap_mutate(self):
		''' Swap the place of 2 genes. '''
		i, j = random.sample(range(self.length), 2)
		self._genes[i], self._genes[j] = self._genes[j], self._genes[i]

	def _reverse_mutate(self):
		''' Reverses a section of DNA. '''
		i, j = random.sample(range(self.length), 2)
		if i > j:
			i, j = j, i
		self._genes[i:j+1] = reversed(self._genes[i:j+1])

	def _scramble_mutate(self):
		''' Shuffles a section of DNA. '''
		i, j = random.sample(range(self.length), 2)
		if i > j:
			i, j = j, i
		self._genes[i:j+1] = random.sample(self._genes[i:j+1], j-i+1)
	#endregion

	#region crossover methods 
	def _uniform_crossover(self, other: Chromosome, n: int, op: Callable[[T, T], T]) -> list[Chromosome]:
		''' Applies an operation on every pair of genes (from each chromosome) to obtain a new gene. '''
		chromosomes = []
		for _ in range(n):
			child = []
			for i in range(self.length):
				gene = op(self._genes[i], other._genes[i])
				child.append(gene)
			# create a new Chromosome obj of the same type as self
			chromosomes.append(type(self)(genes=child))
		return chromosomes
	
	def _k_point_crossover(self, other: Chromosome, n: int, k: int) -> list[Chromosome]:
		''' Crosses alternating sections of two chromosomes to form new ones. '''
		# each crossover yields 2 new chromosomes
		# -> we must run the crossover process at least n/2 times (rounded up)
		chromosomes = []
		repetitions = math.ceil(n / 2)
		for _ in range(repetitions):

			# determine the number of crossover points
			# there is a maximum of (self.length - 1) places to cut the chromosome
			num_points = min(k, self.length - 1)

			# locations to slice the chromosome
			# there will be (num_points + 1) sections [start, end) of genes
			start_points = [0] + random.sample(range(1, self.length), num_points)
			start_points.sort()
			end_points = start_points[1:] + [self.length] # added 'the end of the list of genes'

			# list of genes for two new chromosomes
			children = ([], [])

			j = 0 # alternate which parent chromosome the section comes from
			for start, end in zip(start_points, end_points):
				children[j].extend(self._genes[start:end])
				j = 1 - j # flip
				children[j].extend(other._genes[start:end])

			# add new child chromosomes
			chromosomes.append(type(self)(genes=children[0]))
			chromosomes.append(type(self)(genes=children[1]))
		return chromosomes
	
	def _single_point_crossover(self, other: Chromosome, n: int) -> list[Chromosome]:
		return self._k_point_crossover(other, n, 1)
	
	def _multi_point_crossover(self, other: Chromosome, n: int) -> list[Chromosome]:
		return self._k_point_crossover(other, n, self.k_crosses)
	
	def _uniform_random_crossover(self, other: Chromosome, n: int) -> list[Chromosome]:
		''' Creates a new chromosome by randomly selecting a parent chromosome for each gene. '''
		return self._uniform_crossover(other, n, lambda gene_1, gene_2: random.choice((gene_1, gene_2)))
	
	def _half_uniform_crossover(self, other: Chromosome, n: int) -> list[Chromosome]:
		''' Exactly half the genes that differ in the parents are swapped to produce 2 new children. '''
		chromosomes = []
		repetitions = math.ceil(n/2)
		for _ in range(repetitions):
			# get positions where parent genes differ
			diff_indices = [i for i in range(self.length) if self._genes[i] != other._genes[i]]
			# keep half of these positions, at random
			random.shuffle(diff_indices)
			diff_indices = diff_indices[:len(diff_indices) // 2]
			# operate on the selected half
			child_1 = self._genes[:]
			child_2 = other._genes[:]
			for i in diff_indices:
				child_1[i], child_2[i] = child_2[i], child_1[i]
			chromosomes.append(type(other)(genes=child_1))
			chromosomes.append(type(other)(genes=child_2))
		return chromosomes
	
	def _highly_disruptive_crossover(self, other: Chromosome, n: int) -> list[Chromosome]:
		''' Choose k random genes from parent 1 and swap them with another k random genes from parent 2 to produce offspring. '''
		chromosomes = []
		repetitions = math.ceil(n/2)
		for _ in range(repetitions):
			k = random.randrange(1, self.length) # number of genes to swap
			# choose 2 different sets of genes of length k
			indices_1 = random.sample(self.length, k)
			indices_2 = random.sample(self.length, k)
			# swap them
			child_1 = self._genes[:]
			child_2 = other._genes[:]
			for i in range(k):
				idx_1 = indices_1[i]
				idx_2 = indices_2[i]
				child_1[idx_1], child_2[idx_2] = child_2[idx_2], child_1[idx_1]
			chromosomes.append(type(other)(genes=child_1))
			chromosomes.append(type(other)(genes=child_2))
		return chromosomes
	
	def _average_crossover(self: NumericChromosome, other: NumericChromosome, n: int) -> list[NumericChromosome]:
		''' Creates a new chromosome by averaging the two parent chromosomes. '''
		def operation(gene_1, gene_2):
			x = (gene_1 + gene_2) / 2
			if self.is_integer:
				return round(x)
			return x
		return self._uniform_crossover(other, n, op=operation)
	
	def _flat_crossover(self: NumericChromosome, other: NumericChromosome, n: int) -> list[NumericChromosome]:
		''' Creates a new chromosome by selecting from a uniform distribution between each parent's genes. '''
		chromosomes = []
		for _ in range(n):
			data = []
			for i in range(self.length):
				if self.is_integer:
					gene = random.randint(self._genes[i], other._genes[i])
				else:
					gene = random.uniform(self._genes[i], other._genes[i])
				data.append(gene)
			chromosomes.append(type(self)(genes=data))
		return chromosomes
	
	def _order_1_crossover(self: PermutationalChromosome, other: PermutationalChromosome, n: int) -> list[PermutationalChromosome]:
		''' Randomly select a segment from one parent to be copied to the child.

			All missing indices will be filled out from the second parent, in the order in which they appear, 

			starting from the end of the segment and looping around.
		'''
		chromosomes = []
		parents = [self, other]
		repetitions = math.ceil(n / 2)
		for _ in range(repetitions):
			i, j = sorted(random.sample(range(self.length), 2)) # ensure i < j

			for p in range(2): # each chromosome will have a chance to be parent 1 and 2
				parent_1 = parents[p]
				parent_2 = parents[1 - p]

				child = [None] * self.length
				child[i:j+1] = parent_1._genes[i:j+1] # copied segment from parent_1
				# fill remaining positions from the parent_2, starting from (j+1) and looping around
				num_remaining = self.length - (j - i + 1)
				r = j + 1 # index of remaining gene from parent_2
				for k in range(num_remaining):
					idx = (j + 1 + k) % self.length
					gene = parent_2[r]
					while gene in child:
						r = (r + 1) % self.length
						gene = parent_2[r]
					child[idx] = gene
				# create new chromosome
				chromosomes.append(type(self)(genes=child))
		return chromosomes
	
	def _order_2_crossover(self: PermutationalChromosome, other: PermutationalChromosome, n: int) -> list[PermutationalChromosome]:
		''' 
			- randomly select genes from parent 2
			- find the corresponding genes in parent 1, and keep these positions blank in the child
			- copy all other genes from parent 1 directly to the child
			- copy the remaining genes from parent 2, in the order that they appear (left-to-right)
		'''
		chromosomes = []
		parents = [self, other]
		repetitions = math.ceil(n / 2)
		for _ in range(repetitions):
			# choose 1..(length-1) random positions
			positions = sorted(random.sample(range(self.length), random.randrange(1, self.length)))

			for p in range(2):
				parent_1 = parents[p]
				parent_2 = parents[1 - p]

				child = [None] * self.length

				# from parent 2
				selected_genes = [parent_2[i] for i in positions]
				j = 0 # index for next selected gene

				# copy all of parent 1 to child, except the `selected_genes`
				for i in range(self.length):
					gene = parent_1._genes[i]
					if gene not in selected_genes:
						child[i] = gene
					else: # fill remaining genes from the select in parent 2
						child[i] = selected_genes[j]
						j += 1
				# create new chromosome
				chromosomes.append(type(self)(genes=child))
		return chromosomes
	
	def _partially_mapped_crossover(self: PermutationalChromosome, other: PermutationalChromosome, n: int) -> list[PermutationalChromosome]:
		''' Copies a segment from one parent, and maps the genes at those locations to new ones in the other parent.
			- randomly select a segment from parent 1 to be copied to the child
			- for each gene in the same locations as the segment, but in parent 2:
				- if the gene has already been copied to the child, skip this gene
				- find what alternate gene it maps to (same location in parent 1)
				- find the location of that alternate gene back in parent 2
				- that location will be the destination of the gene in the child
				- if the location was already filled in the child, repeat the mapping process until a free location is found
			- all remaining missing genes are copied verbatim from parent 2
		'''
		chromosomes = []
		for _ in range(n):
			i, j = sorted(random.sample(range(self.length), 2))
			child = [None] * self.length
			child[i:j+1] = self._genes[i:j+1] # copy this segment

			copied_segment = set(child[i:j+1])

			for k in range(i, j+1):
				gene = other._genes[k]
				if gene not in copied_segment: # not from the parent 1 segment
					alt_gene = self._genes[k]
					loci = other._genes.index(alt_gene)
					while child[loci] is not None: # re-map (gene was mapped back onto segment)
						loci = other._genes.index(child[loci])
					child[loci] = gene
			# fill remaining child genes with other
			for k in range(self.length):
				if child[k] is None:
					child[k] = other._genes[k]
			chromosomes.append(type(self)(genes=child))
		return chromosomes
	
	def _cycle_crossover(self: PermutationalChromosome, other: PermutationalChromosome, n: int) -> list[PermutationalChromosome]:
		remaining = set(range(self.length))
		cycles = []
		# find all cycles
		while len(remaining) > 0: # find more cycles
			new_cycle = []
			i = remaining.pop() # random index
			start = self._genes[i] # cycle once we get back here
			gene = None
			while gene != start:
				gene = other._genes[i] # mapped to same loci in other
				i = self._genes.index(gene) # find in self
				new_cycle.append(i)
			cycles.append(new_cycle)
			# subtract from remaining
			remaining.difference_update(new_cycle)

		# create children
		parents = [self, other]
		# there will only be 2 unique children
		child_1 = [None] * self.length
		child_2 = [None] * self.length
		# give cycles from alternating parents to the child
		for j in range(len(cycles)):
			i = j % 2 # switch who parent 1 and 2 is each cycle
			parent_1 = parents[i]
			parent_2 = parents[1-i]
			for loc in new_cycle[j]: # add one full cycle
				child_1[loc] = parent_1[loc]
				child_2[loc] = parent_2[loc]

		# create chromosomes
		chrom_1 = type(self)(genes=child_1)
		chrom_2 = type(self)(genes=child_2)

		chromosomes = [chrom_1, chrom_2]
		repetitions = math.ceil(n / 2) - 2
		for _ in range(repetitions): # add clones of the 2 unique children
			chromosomes.append(chrom_1.clone())
			chromosomes.append(chrom_2.clone())
		return chromosomes
	
	#endregion

	# dunder methods
	# -- to treat Chromosome as a container
	def __getitem__(self, i):
		return self._genes[i]
	
	def __setitem__(self, i, value):
		self._genes[i] = value

	def __len__(self):
		return len(self._genes)
	
	def __iter__(self):
		return iter(self._genes)
	
	def __contains__(self, value):
		return value in self._genes
	
	def __delitem__(self, i):
		raise AttributeError('cannot delete genes')
	
	def __eq__(self, other: Chromosome):
		return self._genes == other._genes
	

class SymbolicChromosome(Chromosome):
	''' Genes can be any string found in the gene_pool. 
		- genes are generated randomly (by default)

		Must implement:
		- `length (int)`
		- `gene_pool (tuple):` 					possible symbols each gene can be encoded by
		- `crossover_type (CrossoverType):`		which method will be used
		- `mutation_rate (float):`				probability a gene is mutated

		Can implement:
		- `k_crosses (int):`					number of points for k-point crossover

		Can override:
		- `_generate(self)`
		- `mutate(self)`
		- `__str__(self)`
	'''
	valid_crossovers = (
		CrossoverType.SINGLE_PT, 
		CrossoverType.MULTI_PT, 
		CrossoverType.UX, 
		CrossoverType.HUX, 
		CrossoverType.HDX
	)

	# prevent instances from being created -- must be subclassed
	def __init__(self, genes=None):
		if type(self) is SymbolicChromosome:
			raise TypeError('SymbolicChromosome cannot be instantiated directly')
		super().__init__(genes)

	def _generate(self):
		self._genes = [random.choice(self.gene_pool) for _ in range(self.length)]

	def mutate(self):
		for i in range(self.length):
			if random.random() < self.mutation_rate:
				cur_gene = self._genes[i]
				new_gene = cur_gene
				while new_gene == cur_gene:
					new_gene = random.choice(self.gene_pool)
				self._genes[i] = new_gene


class NumericChromosome(Chromosome):
	''' Genes will be integers or floats from a specific domain
		described by the gene_pool: (min, max).
		- genes are randomly generated (by default)

		Must implement:
		- `length (int)`
		- `gene_pool (tuple[low, high]):` 		the lower and upper bound of each gene's domain
		- `crossover_type (CrossoverType):`		which method will be used
		- `mutation_rate (float):`				probability a gene is mutated
		- `is_integer (bool):`					use set of integers, or real numbers

		Can implement:
		- `k_crosses (int):`					number of points for k-point crossover
		- `mutation_std_dev (float):`			standard deviation of the normally distributed mutation delta

		Can override:
		- `_generate(self)`
		- `mutate(self)`
		- `__str__(self)`
	'''
	valid_crossovers = (
		CrossoverType.DC, 
		CrossoverType.AX, 
		CrossoverType.FC, 
		CrossoverType.SINGLE_PT, 
		CrossoverType.MULTI_PT, 
		CrossoverType.HUX,
		CrossoverType.HDX
	)
	is_integer: bool
	mutation_std_dev: float

	def __init__(self, genes=None):
		if type(self) is NumericChromosome:
			raise TypeError('NumericChromosome cannot be instantiated directly')
		super().__init__(genes)

	def __str__(self):
		if not self.is_integer: # real numbers
			return ' '.join(f'{gene:.6f}' for gene in self._genes)
		return super().__str__() # integers
	
	def _generate(self):
		lower, upper = self.gene_pool
		self._genes = [random.randint(lower, upper) if self.is_integer else random.uniform(lower, upper) for _ in range(self.length)]
	
	def mutate(self):
		lower, upper = self.gene_pool
		for i in range(self.length):
			if random.random() < self.mutation_rate:
				delta = random.gauss(0, self.mutation_std_dev)
				new_gene = restrict_bounds(self._genes[i] + delta, lower, upper)
				if self.is_integer:
					new_gene = round(new_gene)
				self._genes[i] = new_gene


class PermutationalChromosome(Chromosome):
	''' Genes will be unique integers from `0..len(gene_pool)`

		Every integer must be used exactly once in the chromosome (no duplicates).

		- genes are randomly permuted (by default)

		Must implement:
		- `length (int)`
		- `gene_pool (tuple):` 			the items we are ordering
		- `crossover_type (Enum):`		which method will be used
		- `mutation_rate (float):`		probability the whole chromosome is mutated

		Can override:
		- `_generate(self)`
		- `mutate(self)`
		- `__str__(self)`
	'''
	valid_crossovers = (
		CrossoverType.OX1, 
		CrossoverType.OX2, 
		CrossoverType.PMX, 
		CrossoverType.CX
	)
	
	def __init__(self, genes=None):
		if type(self) is PermutationalChromosome:
			raise TypeError('PermutationalChromosome cannot be instantiated directly')
		super().__init__(genes)

	def __str__(self):
		''' Display the ordering of the domain items represented by this chromosome. '''
		items = (str(self.gene_pool[self._genes[i]]) for i in range(self.length))
		return f'<{", ".join(items)}>'
	
	def _generate(self):
		''' Randomly generates a permutation of the integers from 0 to `self.length - 1`. '''
		self._genes = random.sample(range(self.length), self.length)
	
	def mutate(self):
		''' Apply a chance to randomly select a mutation to affect the whole chromosome:
			1. rotation mutation 	-- (cycles genes by a random amount)
			2. swap mutation 		-- (swap 2 genes)
			3. reverse mutation 	-- (reverses a section of DNA)
			4. scramble mutation 	-- (shuffles a section of DNA)
		'''
		if random.random() < self.mutation_rate:
			mutation = random.choice(self.mutation_methods)
			mutation = getattr(self, mutation)
			mutation()
