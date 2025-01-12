
from typing import Union, Any, Callable

from .helpers import restrict_bounds
from .chromosomes import (
	Chromosome,
	CrossoverType, 
	SymbolicChromosome, 
	NumericChromosome, 
	PermutationalChromosome
)

# Date: December 24, 2024
# Author: Jeffrey Ray
# 
# Description: This module defines factory functions to create new Chromosome types

__all__ = ['create_symbolic_chromosome', 'create_numeric_chromosome', 'create_permutational_chromosome']

## Chromosome factory
def _create_chromosome(
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
	''' Factory function to create a new NumericChromosome type.
	
		Must specify:
		- `_length`
		- `_gene_pool`
		- `_crossover_type`
		- `_mutation_rate`
		- `chromosome_type:` 				the parent Chromosome class

		Optional:
		- `_name`
		- `_k_crosses`

		For Numeric Chromosomes:
		- `_is_integer`
		- `_mutation_std_dev`

		May override:
		- `generate_func(self):`			sets self._genes list
		- `mutate_func(self):`				mutates self._genes
		- `str_func(self):`					returns string

		Set custom class-level variables with **kwargs.
	'''
	if not issubclass(chromosome_type, Chromosome):
		raise TypeError(f'chromosome_type must be a subclass of Chromosome')

	# check common args -- convert to appropriate type if possible
	_name = validate_chromosome_name(_name)
	_mutation_rate = validate_mutation_rate(_mutation_rate)
	validate_crossover_type(_crossover_type, chromosome_type)

	# for symbolic only
	if issubclass(chromosome_type, SymbolicChromosome):
		_gene_pool = validate_gene_pool(_gene_pool)
		_length = validate_length(_length, 3)
		
	# for numeric only
	elif issubclass(chromosome_type, NumericChromosome):
		_gene_pool = validate_numeric_gene_pool(_gene_pool, _is_integer)
		_length = validate_length(_length, 1)

	# for permutational only
	elif issubclass(chromosome_type, PermutationalChromosome):
		_gene_pool = validate_gene_pool(_gene_pool)
		_length = len(_gene_pool) # enforce: length must be exactly the number of items in an ordering
		validate_length(_length, 3)

	# k-point crossovers
	if _crossover_type == CrossoverType.MULTI_PT:
		_k_crosses = validate_k_crosses(_k_crosses, _length)

	class NewChrom(chromosome_type):
		name = _name
		length = _length
		gene_pool = _gene_pool
		crossover_type = _crossover_type
		mutation_rate = _mutation_rate
		k_crosses = _k_crosses
		if issubclass(chromosome_type, NumericChromosome):
			is_integer = _is_integer
			mutation_std_dev = _mutation_std_dev

	if _name is not None:
		NewChrom.__name__ = _name

	# Override methods
	method_overrides = {
		'_generate': generate_func,
		'mutate': mutate_func,
		'__str__': str_func
	}
	for method, override in method_overrides.items():
		if override is not None:
			setattr(NewChrom, method, override)

	# Set arbitrary class-level variables
	for key, value in kwargs.items():
		setattr(NewChrom, key, value)
		
	return NewChrom


# shortcuts with default args
def create_symbolic_chromosome(
	_length: int 						= 10, 
	_gene_pool: tuple[Any] 				= ('0','1'),
	_crossover_type: CrossoverType 		= CrossoverType.SINGLE_PT,
	_mutation_rate: float 				= 0.05,
	_name: str							= None,
	_k_crosses: int						= 2, # (default) 2-point crossover
	generate_func						= None,
	mutate_func							= None,
	str_func							= None,
	**kwargs
) -> type[SymbolicChromosome]:
	''' Factory function to create a new SymbolicChromosome type. 
	
		- `_length:` 					number of genes in chromosome
		- `_gene_pool:`					symbols any gene could be
		- `_mutation_rate:`				probability of independent gene mutation events		
		- `_k_crosses:`					number of points for k-point crossover

		May override:
		- `generate_func(self):`		sets self._genes list
		- `mutate_func(self):`			mutates self._genes
		- `str_func(self):`				returns string

		Set custom class-level variables with **kwargs.
	'''
	return _create_chromosome(
		_length 				= _length, 
		_gene_pool 				= _gene_pool, 
		_crossover_type 		= _crossover_type, 
		_mutation_rate 			= _mutation_rate, 
		chromosome_type 		= SymbolicChromosome,
		_name 					= _name,
		_k_crosses				= _k_crosses,
		generate_func 			= generate_func,
		mutate_func 			= mutate_func,
		str_func 				= str_func,
		kwargs 					= kwargs
	)

def create_numeric_chromosome(
	_length: int 								= 1, 
	_gene_pool: tuple[Union[int, float]] 		= (-100, 100),
	_crossover_type: CrossoverType 				= CrossoverType.AX,
	_mutation_rate: float 						= 0.50,
	_name: str									= None,
	_k_crosses: int								= 2,
	_is_integer: bool							= False, 	# (default) real numbers
	_mutation_std_dev: float					= 0.5, 		# (default) 95% of the time, a value within 0 ± 1
	generate_func 								= None,
	mutate_func									= None,
	str_func									= None,
	**kwargs
) -> type[NumericChromosome]:
	''' Factory function to create a new NumericChromosome type. 
	
		May specify:
		- `_length:` 					number of genes in chromosome
		- `_gene_pool:`					the lower and upper bound of each gene's domain
		- `_mutation_rate:`				probability of independent gene mutation events		
		- `_is_integer:`				domain is integers or real numbers
		- `_mutation_std_dev:`			of the normally distributed mutation delta
		- `_k_crosses:`					number of points for k-point crossover

		May override:
		- `generate_func(self):`		sets self._genes list
		- `mutate_func(self):`			mutates self._genes
		- `str_func(self):`				returns string

		Set custom class-level variables with **kwargs.
	'''
	return _create_chromosome(
		_length 				= _length, 
		_gene_pool 				= _gene_pool, 
		_crossover_type 		= _crossover_type, 
		_mutation_rate 			= _mutation_rate, 
		chromosome_type 		= NumericChromosome,
		_name 					= _name,
		_k_crosses				= _k_crosses,
		_is_integer 			= _is_integer,
		_mutation_std_dev		= _mutation_std_dev,
		generate_func 			= generate_func,
		mutate_func 			= mutate_func,
		str_func 				= str_func,
		kwargs 					= kwargs
	)

def create_permutational_chromosome(
	_gene_pool: tuple[Any] 				= ('A', 'B', 'C', 'D', 'E'),
	_crossover_type: CrossoverType 		= CrossoverType.OX2,
	_mutation_rate: float 				= 0.25,
	_name: str							= None,
	generate_func 						= None,
	mutate_func							= None,
	str_func							= None,
	**kwargs
) -> type[PermutationalChromosome]:
	''' Factory function to create a new PermutationalChromosome type. 
	
		May specify:
		- `_length:` 					number of genes in chromosome
		- `_gene_pool:`					the items we are ordering
		- `_mutation_rate:`				probability the whole chromosome is mutated

		May override:
		- `generate_func(self):`		sets self._genes list
		- `mutate_func(self):`			mutates self._genes
		- `str_func(self):`				returns string

		Set custom class-level variables with **kwargs.
	'''
	return _create_chromosome(
		_length 				= None, 
		_gene_pool 				= _gene_pool, 
		_crossover_type			= _crossover_type, 
		_mutation_rate 			= _mutation_rate, 
		chromosome_type 		= PermutationalChromosome,
		_name 					= _name,
		generate_func 			= generate_func,
		mutate_func 			= mutate_func,
		str_func 				= str_func,
		kwargs 					= kwargs
	)


## validators for factory function args
def validate_chromosome_name(name) -> Union[str, None]:
	if name is not None:
		try:
			string = str(name)
		except Exception as e:
			raise TypeError('failed to convert chromosome name to type str') from e
		return string
	return None

def validate_length(_length, min_size: int) -> int:
	try:
		int_length = int(_length)
	except Exception as e:
		raise TypeError(f'_length must be type int, not {type(_length).__name__}') from e
	if int_length < min_size:
		raise ValueError(f'length of chromosome must be at least {min_size}')
	return int_length

def validate_mutation_rate(_mutation_rate) -> float:
	''' Returns the bound-restricted mutation rate. '''
	try:
		p = float(_mutation_rate)
	except Exception:
		raise TypeError(f'_mutation_rate must be a float in [0,1]')
	return restrict_bounds(p, 0, 1)

def validate_iterable_gene_pool(_gene_pool):
	if not hasattr(_gene_pool, '__iter__'):
		raise TypeError(f'failed to iterate over _gene_pool of type {type(_gene_pool).__name__}')

def validate_gene_pool(_gene_pool) -> tuple:
	''' Any iterable allowed. '''
	validate_iterable_gene_pool(_gene_pool)
	domain = tuple(_gene_pool)
	if len(domain) < 1:
		raise ValueError('_gene_pool must be non-empty')
	return domain

def validate_string_gene_pool(_gene_pool) -> tuple:
	''' Convert an iterable into a tuple of strings. '''
	validate_iterable_gene_pool(_gene_pool)
	try:
		domain = [str(s) for s in _gene_pool]
	except Exception as e:
		raise TypeError('failed to convert all elements of _gene_pool to strings') from e
	if len(domain) < 1:
		raise ValueError('_gene_pool must be non-empty')
	return tuple(domain)

def validate_numeric_gene_pool(_gene_pool, use_ints) -> tuple:
	''' Convert an iterable into a 2-tuple of ints or floats. '''
	validate_iterable_gene_pool(_gene_pool)
	_type = int if use_ints else float
	_type_name = 'ints' if use_ints else 'floats'
	try:
		domain = [_type(i) for i in _gene_pool]
	except Exception as e:
		raise TypeError(f'failed to convert all elements of _gene_pool to {_type_name}') from e
	if len(domain) != 2:
		raise ValueError(f'_gene_pool must be of length 2: (lower_bound, upper_bound), not {domain}')
	return tuple(domain)
		
def validate_crossover_type(_crossover_type: CrossoverType, chrom_type: type[Chromosome]):
	if not isinstance(_crossover_type, CrossoverType):
		raise TypeError(f'_crossover_type must be CrossoverType, not {type(_crossover_type).__name__}')
	if _crossover_type not in chrom_type.valid_crossovers:
		acceptable_names = [c.name for c in chrom_type.valid_crossovers]
		raise ValueError(
			f'a {chrom_type.__name__} may only use the crossover types: {acceptable_names}\n' +
			f'but received: {_crossover_type.name}'
		)

def validate_k_crosses(_k_crosses, _length: int) -> int:
	try:
		cross_int = int(_k_crosses)
	except Exception as e:
		raise TypeError(f'_k_crosses must be type int, not {type(_k_crosses).__name__}') from e
	if cross_int < 1:
		raise ValueError(f'_k_crosses must be greater than 1')
	if cross_int > _length - 1:
		raise ValueError(f'a chromosome of size {_length} cannot have more than {_length - 1} cross points')
	return cross_int
