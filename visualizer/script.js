
// Global variables
let generation = 0; 			// Tracks the current generation
let isAnimating = false; 		// Tracks whether animation is running
let animationInterval = null; 	// Holds the animation interval ID
let evolveTime = null; 			// Time in milliseconds (per generation)

// Reference to elements
const speedSlider = document.getElementById('speed-slider'); 		// slider for animation speed
const mutationSlider = document.getElementById('mutation-slider')	// for mutation rate
const evolveButton = document.getElementById('evolve');
const pauseButton = document.getElementById('pause');

const fitnessInputs = {
	a: document.getElementById('fitness-a'),
	b: document.getElementById('fitness-b'),
	c: document.getElementById('fitness-c'),
	d: document.getElementById('fitness-d'),
	e: document.getElementById('fitness-e')
};

const cells = Array.from(document.querySelectorAll('.cell'));


// Class for the living creatures
class Creature {
	static length = 10; 		// length for all strings
	static mutationRate = 0.1; 	// 10% chance of mutation per character

	constructor(string = null) {
		this.string = string || Creature.generateString();
		this.fitness = 0; // Initialize with default fitness of 0
	}

	toString() {
		return this.string;
	}

	// Calculate this creature's fitness based on a mapping of chars to fitness values
	calculateFitness(fitnessMap) {
		this.fitness = Array.from(this.string).reduce((sum, char) => {
			return sum + (fitnessMap[char] || 0);
		}, 0);
	}

	static getRandomGene() {
		const genes = ['a', 'b', 'c', 'd', 'e', 'w', 'x', 'y', 'z'];
		return genes[Math.floor(Math.random() * genes.length)];
	}

	static generateString() {
		return Array.from({ length: Creature.length }, () => 
			Creature.getRandomGene()
		).join('');
	}

	// Use single-point crossover
	static crossover(parentA, parentB) {
		// At least segment of length 1 must be crossed
		const crossPoint = Math.floor(Math.random() * (Creature.length - 1)) + 1;
		return new Creature(
			parentA.string.slice(0, crossPoint) + parentB.string.slice(crossPoint)
		);
	}

	// Apply chance to randomly replace a gene with another for each gene in the string
	mutate() {
		this.string = Array.from(this.string).map(char => {
			if (Math.random() < Creature.mutationRate) {
				return Creature.getRandomGene();
			}
			return char;
		}).join('');
	}
}


// Our list of 8 individuals
const population = Array.from({ length: 8 }, () => new Creature());

// Initialize gene-fitness map
const fitnessMap = {a: 0, b: 0, c: 0, d: 0, e: 0};


// Add class colors to the genes in a string
function colorizeString(string) {
	return Array.from(string).map(char => {
		const className = `letter-${char}`;
		return `<span class="${className}">${char}</span>`;
	}).join('');
}

// Function to update fitness values from inputs
function updateFitnessMap() {
	for (const key in fitnessMap) {
		fitnessMap[key] = parseInt(fitnessInputs[key].value) || 0;
	}
}

// Updates fitness of all, and sorts the population in descending order of fitness
// -- based on current mapping
function rankPopulation() {
	population.forEach(individual => individual.calculateFitness(fitnessMap));
	population.sort((a, b) => b.fitness - a.fitness);
}

// Update the cells on the grid -- assumes a ranked population
function updateGrid() {

	// Update the top 8 best strings and display their fitness
	cells.forEach((cell, index) => {
		const {string, fitness} = population[index];

		cell.querySelector('.string').innerHTML = colorizeString(string);
		cell.querySelector('.fitness').textContent = `Fitness: ${fitness}`;
	});
}

// First re-rank the population, then update the grid
function rankAndUpdateGrid() {
	updateFitnessMap();
	rankPopulation();
	updateGrid();
}

// Generate the next generation of creatures and update our population
function evolve() {

	// Step 1: Select the parents (top-4 fit)
	const parents = population.slice(0, 4);

	// Step 2: Generate 8 new children using crossover and mutation
	const newGeneration = [];
	while (newGeneration.length < 8) {
		// Choose 2 random parents to procreate
		const parentA = parents[Math.floor(Math.random() * parents.length)];
		const parentB = parents[Math.floor(Math.random() * parents.length)];

		const child = Creature.crossover(parentA, parentB);
		child.mutate();
		child.calculateFitness(fitnessMap);
		newGeneration.push(child);
	}

	// Step 3: Combine old and new generations
	const combinedPool = [...population, ...newGeneration];
	combinedPool.sort((a, b) => b.fitness - a.fitness); // rank

	// Step 4: Keep the top 4 overall (elitism)
	const top4 = combinedPool.slice(0, 4);
	population.splice(0, population.length, ...top4); // now length 4

	// Step 5: Bring 4 additional members to the next generation from the group of children
	// -- there will be at least 4 not included in the population yet
	let childrenAdded = 0;
	for (const child of newGeneration){
		if (!population.includes(child)){
			population.push(child);
			++childrenAdded;
		}
		if (childrenAdded == 4) break;
	}

	// sort the added children to maintain ranking
	population.sort((a, b) => b.fitness - a.fitness);

	// Step 6: Update generation count and grid
	++generation;
	document.getElementById('generation').textContent = generation;
	updateGrid();
}


// Function to start animation
function startAnimation() {
	if (isAnimating) return; // Prevent multiple intervals
	isAnimating = true;

	animationInterval = setInterval(() => {
		evolve();
	}, evolveTime);
}

// Function to stop animation
function stopAnimation() {
	isAnimating = false;
	clearInterval(animationInterval);
}

function updateEvolveSpeed() {
	const evolveSpeed = parseFloat(speedSlider.value); // gen. per sec

	// Convert to frequency
	evolveTime = 1000 / evolveSpeed; // ms per gen.

	// Display generations per second
	document.getElementById('speed-display').textContent = evolveSpeed.toFixed(1);
}

function updateMutationRate() {
	const mutationRate = parseFloat(mutationSlider.value);
	// Fallback in case mutation rate was not set properly (0.1 default)
	const validMutationRate = (isNaN(mutationRate) || mutationRate < 0 || mutationRate > 1) ? 0.1 : mutationRate;
	Creature.mutationRate = validMutationRate;
	// Display to 2 decimal places
	document.getElementById('mutation-display').textContent = validMutationRate.toFixed(2);
}


// Event listeners to update fitness when inputs change
Object.values(fitnessInputs).forEach(input => {
	input.addEventListener('input', rankAndUpdateGrid);
});

// Event listeners for buttons
evolveButton.addEventListener('click', startAnimation);
pauseButton.addEventListener('click', stopAnimation);

// Event listener for speed slider
speedSlider.addEventListener('input', () => {
	updateEvolveSpeed();
	if (isAnimating) {
		stopAnimation();
		startAnimation(); // Restart animation with new speed
	}
});

mutationSlider.addEventListener('input', () =>{
	updateMutationRate();
});

// Keyboard shortcut for toggling evolve/pause
document.addEventListener('keyup', (event) => {
	if (event.code === 'Space') {
		event.preventDefault(); // Prevent default scrolling behavior

		const targetButton = isAnimating ? pauseButton : evolveButton;

		targetButton.classList.add('pressed');
		targetButton.click();

		// Remove 'pressed' after a short delay -- 0.15 seconds
		setTimeout(() => targetButton.classList.remove('pressed'), 150);
	}
});


// For initial render:
rankAndUpdateGrid(); // read fitness score, rank generation 0, then update cell information

// match the sliders to start
updateEvolveSpeed();
updateMutationRate();