// Replace links to .ipynb with links to .html
document.addEventListener("DOMContentLoaded", function () {
	// Select all anchor tags with an href attribute
	const links = document.querySelectorAll('a[href]');

	links.forEach(link => {
		const href = link.getAttribute('href');

		if (href && href.startsWith('./pyvolver')) {
			const updatedHref = 'https://github.com/jray-8/genetic-algorithms/tree/b38534c7031ae393b882e9faf5576596281ae832/pyvolver';
			link.setAttribute('href', updatedHref);
		}

		// Only modify links ending in .ipynb
		else if (href && href.endsWith('.ipynb')) {
			const updatedHref = href.replace(/\.ipynb$/, '.html');
			link.setAttribute('href', updatedHref);
		}
	});
});
