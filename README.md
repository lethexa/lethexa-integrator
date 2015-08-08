This is a very simple implementation of the Runge-Kutta-integation method.


Usage
-----
	var integrator = require('lethexa-rungekutta');

	var rungekutta = new integrator.RungeKutta( function(x, y) {
		return 3 * x*x * y;
	});

	var x0 = 1, y0 = 2, h = 0.1;

	var result = rungekutta.nextStep(x0, y0, h);

	console.log(result);

Test
----
Run 'npm install' and 'npm test' in the projects main directory


Contributors
------------

* lethexa 


