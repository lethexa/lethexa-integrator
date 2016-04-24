lethexa-integrator
==================

This is a very simple implementation of the integation methods Euler, Runge-Kutta first order and second order.
You can provide your own algebra to the integrator for example to integrate vectors with existing vector libraries.

For this purpose implement the algebra (add, mulScalar) in an object:

	var MathFunc = {
		add: function(x, y) {
			return x + y;
		},

		mulScalar: function(x, s) {
			return x * s;
		}
	};

	var euler = new integrator.Euler( function(x, y) {
                return 3 * x*x * y;
        }, MathFunc);
	...




Solve 1st order differential equation using Euler method
--------------------------------------------------------

        var integrator = require('lethexa-integrator');

        var euler = new integrator.Euler( function(x, y) {
                return 3 * x*x * y;
        });

        var input = {
		x: 1,
		y: 2
	}
	var h = 0.1;
        var result = euler.nextStep(input, h);
        console.log(result);


Solve 1st order differential equation using Runge-Kutta
-------------------------------------------------------

	var integrator = require('lethexa-integrator');

	var rungekutta1Order = new integrator.RungeKutta1Order( function(x, y) {
		return 3 * x*x * y;
	});

	var input = {
		x: 1,
		y: 2
	}
	var h = 0.1;
	var result = rungekutta1Order.nextStep(input, h);
	console.log(result);


Solve 2st order differential equation using Runge-Kutta
-------------------------------------------------------

	var integrator = require('lethexa-integrator');
	
	var rungekutta = new integrator.RungeKutta2Order( function(x, y, dy) {
		return 2*y + dy;
	});

	var input = {
		x: 0, 
		y: 3,
		dy: 0
	}
	var h = 0.1;
	var result = rungekutta.nextStep(input, h);
	console.log(result);


Test
----

Run 'npm install' and 'grunt test' in the projects main directory


License
-------

This library is published under MIT-license.

