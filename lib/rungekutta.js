
/**
 * Runge-Kutta integration scheme for first order differential equations
 * @class Euler
 * @constructor
 * @param func {function(x, y)} The function to integrate.
 * @param ops {Object} The operations add, mul, ...
 */
module.exports.Euler = function( func, ops ) {

	ops = ops || {};
	ops.add = ops.add || function(x, y) {
		return x + y;
	};

	ops.mulScalar = ops.mulScalar || function(x, c) {
		return x * c;
	};


	/**
	 * Calculates the next step via integrating
	 * @method nextStep
	 * @param x0 {Number} The first value
	 * @param y0 {Number} The first value
	 * @param h {Number} The stepsize
	 * @return {object} Returns an object in the following form 
	 * @example
	    {
	 	x1: <value>,
		y1: <value>
	    }
	 */
	this.nextStep = function( x0, y0, h ) {
		return {
			x1: x0 + h,
			y1: ops.add(y0, ops.mulScalar(func( x0, y0 ), h))
		};
	};
};


/**
 * Runge-Kutta integration scheme for first order differential equations
 * @class RungeKutta1Order
 * @constructor
 * @param func {function(x, y)} The function to integrate.
 * @param ops {Object} The operations add, mul, ...
 */
module.exports.RungeKutta1Order = function( func, ops ) {

	ops = ops || {};
	ops.add = ops.add || function(x, y) {
		return x + y;
	};

	ops.mulScalar = ops.mulScalar || function(x, c) {
		return x * c;
	};


	/**
	 * Calculates the next step via integrating
	 * @method nextStep
	 * @param x0 {Number} The first value
	 * @param y0 {Number} The first value
	 * @param h {Number} The stepsize
	 * @return {object} Returns an object in the following form 
	 * @example
	    {
	 	x1: <value>,
		y1: <value>
	    }
	 */
	this.nextStep = function( x0, y0, h ) {
		var hHalf = h * 0.5;

		var k1 = ops.mulScalar(func( x0, y0 ), h);
		var k2 = ops.mulScalar(func( x0 + hHalf, ops.add(y0, ops.mulScalar(k1, 0.5)) ), h);			
		var k3 = ops.mulScalar(func( x0 + hHalf, ops.add(y0, ops.mulScalar(k2, 0.5)) ), h);			
		var k4 = ops.mulScalar(func( x0 + h, ops.add(y0, k3) ), h);			

		var sumK = k1;
		sumK = ops.add(sumK, ops.mulScalar(k2, 2.0));
		sumK = ops.add(sumK, ops.mulScalar(k3, 2.0));
		sumK = ops.add(sumK, k4);

		return {
			x1: x0 + h,
			y1: ops.add(y0, ops.mulScalar(sumK, 1.0 / 6.0))
		};
	};
};



/**
 * Runge-Kutta integration scheme for second order differential equations
 * @class RungeKutta2Order
 * @constructor
 * @param func {function(x, y, dy)} The function to integrate.
 * @param ops {Object} The operations add, mul, ...
 */
module.exports.RungeKutta2Order = function( func, ops ) {

	ops = ops || {};
	ops.add = ops.add || function(x, y) {
		return x + y;
	};

	ops.mulScalar = ops.mulScalar || function(x, c) {
		return x * c;
	};


	/**
	 * Calculates the next step via integrating
	 * @method nextStep
	 * @param x0 {Number} The initial value
	 * @param y0 {Number} The initial value
	 * @param dy0 {Number} The initial value
	 * @param h {Number} The stepsize
	 * @return {object} Returns an object in the following form 
	 * @example
	    {
	 	x1: <value>,
		y1: <value>,
		dy1: <value>
	    }
	 */
	this.nextStep = function( x0, y0, dy0, h ) {
		var hHalf = h * 0.5;

		var k1 = ops.mulScalar(dy0, h);
		var m1 = ops.mulScalar(func(x0, y0, dy0), h);
		
		var k2 = ops.mulScalar(ops.add(dy0, ops.mulScalar(m1, 0.5)), h);			
		var m2 = ops.mulScalar(func( x0 + hHalf, ops.add(y0, ops.mulScalar(k1, 0.5)), ops.add(dy0, ops.mulScalar(m1, 0.5)) ), h);
		
		var k3 = ops.mulScalar(ops.add(dy0, ops.mulScalar(m2, 0.5)), h);			
		var m3 = ops.mulScalar(func( x0 + hHalf, ops.add(y0, ops.mulScalar(k2, 0.5)), ops.add(dy0, ops.mulScalar(m2, 0.5)) ), h);

		var k4 = ops.mulScalar(ops.add(dy0, m3), h);
		var m4 = ops.mulScalar(func( x0 + h, ops.add(y0, k3), ops.add(dy0, m3) ), h);		


		var sumK = k1;
		sumK = ops.add(sumK, ops.mulScalar(k2, 2.0));
		sumK = ops.add(sumK, ops.mulScalar(k3, 2.0));
		sumK = ops.add(sumK, k4);

		var sumM = m1;
		sumM = ops.add(sumM, ops.mulScalar(m2, 2.0));
		sumM = ops.add(sumM, ops.mulScalar(m3, 2.0));
		sumM = ops.add(sumM, m4);

		return {
			x1: x0 + h,
			y1: ops.add(y0, ops.mulScalar(sumK, 1.0 / 6.0)),
			dy1: ops.add(dy0, ops.mulScalar(sumM, 1.0 / 6.0))
		};
	};
};


