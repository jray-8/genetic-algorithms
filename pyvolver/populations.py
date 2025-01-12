
import math
from enum import Flag, auto
from typing import Union

from .species import Organism

# Date: December 17, 2024
# Author: Jeffrey Ray
# 
# Description:
# 	This module defines the Population class, which creates a society of individuals from a single Species.

__all__ = ['EvolveFlags', 'Population']

# options that can be set for a population to change how evolution occurs
class EvolveFlags(Flag):
	'''
		- MATE_WITH_OTHERS - half the parents will be randomly generated individuals (outsiders)
		- TWO_GROUP_MATING - Parents are divided into two groups, and mating occurs across the groups
	'''
	NONE = 0
	MATE_WITH_OTHERS = auto()
	TWO_GROUP_MATING = auto()


class Population:
	''' A society of individuals of the same `Species`.  
		Mating will occur between 2 unique individuals.

		Elitist by default. Strongest members keep living (until they are not the strongest).

		Attributes:
		- `members (list):` 				All individuals in the society.
		- `ideal_fitness (float):`			A fitness score which represents the desired solution.
		- `_evflags (EvolveFlags):` 		Bitmask of boolean options that affect the evolution process.

		Properties:
		- `species (Species):`				The class of all individuals in society.
		- `start_size (int):`				The number of individuals who established society.
		- `pop_size (int):`					The current population size.
		- `max_capacity (int):` 			The population size may never exceed this amount.
		- `generation (int):`				The number of evolution cycles since the initial population (gen 0).
		- `solution (Organism):`			The fittest member of society.

		Population Dynamics:
		- `birth_rate (int):` 				How many babies should be born each generation?
		- `num_parents (int):` 				The number of parents selected each generation.
		- `num_offspring (int):` 			The number of babies each mating pair will have.

		Changing a population dynamic will hold 2 variables static and recalculate the third.

		The variable being changed will be preserved, and the second held constant will be:

		Precedence (to keep constant):
		1. num_parents
		2. birth_rate		
		3. num_offspring
		
	'''
	MIN_CAP = 2 			# minimum capacity (to avoid extinction)
	MIN_PARENTS = 2 		# minimum number of parents allowed
	STR_INDENT = ' ' * 2	# indent for members info in __str__

	def __init__(self,
			species: 			type[Organism],
			pop_size: 			int = 8,
			max_capacity: 		int = None, 		# pop_size by default 
			birth_rate: 		int = None, 		# sustains intitial pop_size by default
			num_parents: 		int = None,
			num_offspring: 		int = None,
			ideal_fitness:		float = None,
			evolve_options: 	EvolveFlags = EvolveFlags.NONE
		):
		'''
			- max_capacity defaults to pop_size. Use `-1` to indicate no cap
			- birth_rate defaults to pop_size (to sustain initial population size)
			- num_parents defaults to half pop_size
		'''
		# instance variables
		self.members: list[Organism]
		self.ideal_fitness: float
		self._evflags: EvolveFlags
		self._species: type[Organism]
		self._generation: int
		self._start_size: int
		self._max_capacity: int			# pop_size by default
		self._birth_rate: int			# lower bound, initial pop_size by default
		self._num_parents: int			# upper bound
		self._num_offspring: int		# exact amount

		# set variables with implicit conversions, and validate
		self._set_species(species)
		self._set_start_size(pop_size)
		self._set_max_capacity(max_capacity)
		self._set_ideal_fitness(ideal_fitness)
		self._evflags = evolve_options # bitmask for evolution behaviour

		# form the initial society
		self.reset(pop_size)

		# set last
		self.set_population_dynamics(birth_rate, num_parents, num_offspring)

	#region Container 
	# -- treat as container (for Organisms)
	def __getitem__(self, i):
		return self.members[i]
	
	def __setitem__(self, i, value):
		self.members[i] = value

	def __len__(self):
		return len(self.members)
	
	def __iter__(self):
		return iter(self.members)
	
	def __contains__(self, value):
		return value in self.members
	
	def __delitem__(self, i):
		del self.members[i]
	#endregion

	#region Validators (and direct set) 
	def _set_species(self, species: type):
		''' Sets the population species for the first time. '''
		if not issubclass(species, Organism):
			if isinstance(species, type):
				type_name = species.__name__
			else:
				type_name = type(species).__name__
			raise TypeError(f'species must be a subclass of Organism, but is type {type_name}')
		self._species = species

	def _set_start_size(self, pop_size: int):
		''' Sets the size of the population for the first time. '''
		try:
			pop_size_int = int(pop_size)
		except Exception as e:
			raise TypeError(f'pop_size must be type int, not {type(pop_size).__name__}') from e
		if pop_size_int < Population.MIN_CAP:
			raise ValueError(f'pop_size must be greater than {Population.MIN_CAP}, or risk extinction')
		self._start_size = pop_size_int
		
	def _set_max_capacity(self, cap: int):
		''' Sets the upper-bound limit on the population size. '''
		if cap is None:
			self._max_capacity = self._start_size
			return
		try:
			cap_int = int(cap)
		except Exception as e:
			raise TypeError(f'max_capacity must be type int, not {type(cap).__name__}') from e
		# no limit indicator
		if cap_int < 0:
			self._max_capacity = None
			return
		elif cap_int < Population.MIN_CAP:
			raise ValueError(f'max_capacity must be at least {Population.MIN_CAP}, or risk extinction')
		self._max_capacity = cap_int

	def _set_ideal_fitness(self, fitness: float):
		''' Sets the fitness score of the `golden child` solution. '''
		if fitness is None:
			self.ideal_fitness = None
			return
		try:
			self.ideal_fitness = float(fitness)
		except Exception as e:
			raise TypeError(f'ideal_fitness must be type float, not {type(fitness).__name__}') from e
		
	def _set_birth_rate(self, n: int):
		''' Does not set any other population dynamic. '''
		try:
			n_int = int(n)
		except Exception as e:
			raise TypeError(f'birth_rate must be type int, not {type(n).__name__}') from e
		if n_int < 1:
			raise ValueError('birth_rate must be at least 1, to resupply the population')
		self._birth_rate = n_int

	def _set_num_parents(self, p: int):
		''' Does not set any other population dynamic. '''
		try:
			p_int = int(p)
		except Exception as e:
			raise TypeError(f'num_parents must be type int, not {type(p).__name__}') from e
		if p_int < Population.MIN_PARENTS:
			raise ValueError(f'num_parents must be at least {Population.MIN_PARENTS}, or risk extinction')
		self._num_parents = p_int

	def _set_num_offspring(self, b: int):
		''' Does not set any other population dynamic. '''
		try:
			b_int = int(b)
		except Exception as e:
			raise TypeError(f'num_offspring must be type int, not {type(b).__name__}') from e
		if b_int < 1:
			raise ValueError('mating couples must have at least 1 offspring')
		self._num_offspring = b_int
	#endregion

	#region Properties 
	@property
	def species(self):
		return self._species
	
	@property
	def pop_size(self):
		''' The number of members in the current population. '''
		return len(self.members)
	
	@property
	def start_size(self):
		''' The number of members who founded the population. '''
		return self._start_size
	
	@property
	def generation(self):
		return self._generation
	
	@property
	def solution(self):
		''' The strongest member in society. '''
		return self.members[0]
	
	# settable
	@property
	def max_capacity(self):
		''' Maximum size for the population. '''
		return self._max_capacity
	
	@max_capacity.setter
	def max_capacity(self, value):
		self._set_max_capacity(value)
	
	@property
	def birth_rate(self):
		''' Children born each generation. 

			Will update num_offspring, keeping num_parents constant.
		'''
		return self._birth_rate
	
	@birth_rate.setter
	def birth_rate(self, value):
		self._set_birth_rate(value)
		self._update_num_offspring()
	
	@property
	def num_parents(self):
		''' Number of parents each generation.

			Will update num_offspring, keeping the birth_rate constant. 
		'''
		return self._num_parents
	
	@num_parents.setter
	def num_parents(self, value):
		self._set_num_parents(value)
		self._update_num_offspring()

	@property
	def num_offspring(self):
		''' Babies made per couple. 

		 	Will update birth_rate, keeping num_parents constant.
		'''
		return self._num_offspring
		
	@num_offspring.setter
	def num_offspring(self, value):
		self._set_num_offspring(value)
		self._update_birth_rate()
	#endregion

	#region Population Dynamics 
	def _update_num_offspring(self):
		''' Sets the number of offspring of each couple based on birth_rate and num_parents. '''
		''' 
			Let's say there will be p parents per generation. 
			We enforce p >= 2.

			Then, let:

				n = birth_rate
				b = num_offspring
		
			If each selected parent mates with every other selected parent once,
			there will be (p choose 2) unique pairs.
			 
			If we want a total of n children born each generation,
			and we choose p parents per generation,
			then we can calculate the number of offspring, b, that each mating pair must produce:
				n = (p choose 2) * b
				n = [p! / ((p-2)! * 2!)] * b
				2n = [p! / (p-2)!] * b
				2n = [p(p-1)(p-2)! / (p-2)!] * b
				2n = p(p-1) * b
				b = 2n / p(p-1)
				
			For EvolveFlags.TWO_GROUPS_MATING:
				
			If, instead, we want 2 distinct groups of parents to mate across with,
			such that there are p/2 parents in each group, then:
				n = (p/2) * (p/2) * b
				n = (1/4) * (p^2) * b
				b = 4n/(p^2)

			Or, if p/2 is not an integer:
				n = floor(p/2) * ceil(p/2) * b
				b = n / (floor(p/2) * ceil(p/2))
				
			Remember, we must have a whole number of offspring, so we round b up to the nearest integer.
		'''
		if EvolveFlags.TWO_GROUP_MATING in self._evflags:
			self._num_offspring = math.ceil(self._birth_rate / ((self._num_parents // 2) * math.ceil(self._num_parents / 2)))
		else:
			self._num_offspring = math.ceil(2 * self._birth_rate / (self._num_parents * (self._num_parents - 1)))

	def _update_num_parents(self):
		''' Sets the number of parents based on birth rate and offspring count. '''
		''' 
			from _update_num_offspring():

				b = 2n / p(p-1)
				p(p-1) = 2n/b
				p^2 - p - 2n/b = 0

				Quadratic Formula:
					p = [-(-1) ± sqrt( (-1)^2 - 4(1)(-2n/b) )] / 2(1)
					p = [1 ± sqrt(1 + 8n/b)] / 2

					since p >= 2, we will choose the positive root
					p = [1 + sqrt(1 + 8n/b)] / 2

			For EvolveFlags.TWO_GROUPS_MATING:
			
				b = 4n/(p^2)
				p = sqrt(4n/b)
		'''
		if EvolveFlags.TWO_GROUP_MATING in self._evflags:
			self._num_parents = math.ceil(math.sqrt(self._birth_rate / self._num_offspring))
		else:
			self._num_parents = math.ceil((1 + math.sqrt(1 + (8 * self._birth_rate / self._num_offspring))) / 2)

	def _update_birth_rate(self):
		''' Sets the birth rate based on the number of parents and the offspring count. '''
		'''
			from _update_num_offspring():
				b = 2n / p(p-1)
				n = [b * p(p-1)] / 2

			For EvolveFlags.TWO_GROUPS_MATING:
				n = floor(p/2) * ceil(p/2) * b

			Here, we don't care if the birth rate is a whole number.
		'''
		if EvolveFlags.TWO_GROUP_MATING in self._evflags:
			self._birth_rate = self._num_offspring * (self._num_parents ** 2) / 4
		else:
			self._birth_rate = (self._num_offspring * self._num_parents * (self.num_offspring - 1)) / 2

	def set_population_dynamics(self, birth_rate: Union[int, None], num_parents: Union[int, None], num_offspring: Union[int, None]):
		''' Set and validate the 3 population dynamics.

			Will keep 2/3 of the variables constant, and calculates the third accordingly.

			Fills in missing data `None` with defaults.

			Precedence (to keep constant):
			1. num_parents
			2. birth_rate
			3. num_offspring
		'''
		# (offspring is not set) OR (everything is set)
		if num_offspring is None or all(x is not None for x in (birth_rate, num_parents, num_offspring)):
			# update offspring num
			if birth_rate is None:
				birth_rate = self._start_size
			if num_parents is None:
				num_parents = max(Population.MIN_PARENTS, self.pop_size // 2)
			self._set_birth_rate(birth_rate)
			self._set_num_parents(num_parents)
			self._update_num_offspring()
		# only birth rate and offspring are set
		elif birth_rate is not None:
			# update parents num
			self._set_birth_rate(birth_rate)
			self._set_num_offspring(num_offspring)
			self._update_num_parents()
		else: # offspring is set, birth rate is not set
			# update birth rate
			if num_parents is None:
				num_parents = max(Population.MIN_PARENTS, self.pop_size // 2)
			self._set_num_parents(num_parents)
			self._set_num_offspring(num_offspring)
			self._update_birth_rate()
		''' Cases for parameters: 

			None, None, None 	-> C **
			None, B, None 		-> C *
			A, None, None		-> C *
			A, B, None			-> C
			A, B, C				-> C

			A, None, C			-> B

			None, None, C 		-> A *
			None, B, C 			-> A
		'''
	#endregion

	def __str__(self):
		indent = Population.STR_INDENT
		s = f'--- Generation {self._generation} ---\n\n'
		for i, member in enumerate(self.members):
			s += f'Member {i+1})\n'
			s += indent + Organism.__str__(member).replace('\n', '\n' + indent) + '\n\n'
		return s.strip()
	
	def _rank_members(self):
		''' Sorts the members list from most fit to least fit. '''
		self.members.sort(reverse=self._species.maximize, key=lambda x: x.fitness)

	def _get_strongest(self, n):
		''' Returns the top n strongest members from society. '''
		return self.members[:n]

	def _select(self):
		''' Choose members of society to be this generation's parents,'''
		if EvolveFlags.MATE_WITH_OTHERS in self._evflags:
			others = [self._species() for _ in range(math.ceil(self._num_parents / 2))]
			return self._get_strongest(self._num_parents // 2) + others
		else:
			return self._get_strongest(self._num_parents)
		
	def breed(self):
		''' The chosen parents of this generation have at it. 
			- Returns a list of children 
		'''
		selected = self._select()
		n = len(selected)
		children = []
		if EvolveFlags.TWO_GROUP_MATING in self._evflags: # two groups
			for i in range(n // 2):
				for j in range(n // 2, n):
					children.extend(selected[i].procreate(selected[j], self._num_offspring))
		else: # one group
			for i in range(n):
				for j in range(i+1, n):
					children.extend(selected[i].procreate(selected[j], self._num_offspring))
		return children

	def evolve(self, n: int=10):
		''' Evolve the population for n more generations. '''
		for _ in range(n):
			# golden child -- perfect genetics
			if self.ideal_fitness is not None and abs(self.solution.fitness - self.ideal_fitness) < 1e-6:
				break
			# selection
			children = self.breed()
			# replacement
			self.members.extend(children)
			self._rank_members()
			self.members = self.members[:self._max_capacity] # keep only the strongest
			self._generation += 1

	#region Population Control 
	def reset(self, pop_size=None):
		''' Reinitializes the population from scratch, with the same initial conditions. '''
		if pop_size is not None:
			self._start_size = pop_size
		self.members = [self._species() for _ in range(self._start_size)]
		self._rank_members()
		self._generation = 0

	def spawn(self, n):
		''' Generate n new members to join the population. '''
		for _ in range(n):
			self.members.append(self._species())
		self._rank_members()

	def add_members(self, new_recruits: list):
		''' Add new recruits to the population. 
			- recruits must be of the same species as the population
		'''
		# join society if recruits are same species
		for person in new_recruits:
			if isinstance(person, self._species):
				self.members.append(person)
			else:
				raise TypeError(f'cannot add a {type(person).__name__} individual to a population of type {type(self).__name__}')
		self._rank_members()
			
	def replace_members(self, n: int, replacements: list):
		''' Replace the n weakest members of society. 
			- if new members are not supplied, random ones will be generated
		'''
		if not isinstance(n, int):
			raise TypeError('n must be an integer')
		if n < 1:
			raise ValueError('must replace at least 1 member')
		del self.members[-n:]
		if replacements:
			self.add_members(replacements)
		n -= len(replacements) # excess needed
		if n > 0:
			self.spawn(n)

	def genocide(self, threshold):
		''' Kill off all solutions with a fitness 'weaker' than the threshold. '''
		superior_population = [pure_blood for pure_blood in self.members if pure_blood.healthier(threshold)]
		del self.members # they were too weak...
		self.members = superior_population
	#endregion
