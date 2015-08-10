
/**
 * Runge-Kutta integration scheme for first order differential equations
 * @class RungeKutta1Ord
 * @constructor
 * @param func {function(x, y)} The function to integrate.
 * @param ops {Object} The operations add, mul, ...
 */
module.exports.RungeKutta1Order = function( func, ops ) {

	ops = ops || {};

	var add = function(x, y) {
		return x + y;
	};

	var mul = function(x, y) {
		return x * y;
	};

	var mulConst = function(x, c) {
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

		var k1 = mulConst(func( x0, y0 ), h);
		var k2 = mulConst(func( add(x0, hHalf), add(y0, mulConst(k1, 0.5)) ), h);			
		var k3 = mulConst(func( add(x0, hHalf), add(y0, mulConst(k2, 0.5)) ), h);			
		var k4 = mulConst(func( add(x0, h), add(y0, k3) ), h);			

		var sumK = k1;
		sumK = add(sumK, mulConst(k2, 2.0));
		sumK = add(sumK, mulConst(k3, 2.0));
		sumK = add(sumK, k4);

		return {
			x1: x0 + h,
			y1: y0 + sumK / 6.0
		};
	};
};



/**
 * Runge-Kutta integration scheme for second order differential equations
 * @class RungeKutta2Ord
 * @constructor
 * @param func {function(x, y, dy)} The function to integrate.
 * @param ops {Object} The operations add, mul, ...
 */
module.exports.RungeKutta2Order = function( func, ops ) {

	ops = ops || {};

	var add = function(x, y) {
		return x + y;
	};

	var mul = function(x, y) {
		return x * y;
	};

	var mulConst = function(x, c) {
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
		//var hHalf = h * 0.5;
		//var k1 = h * func( x0, y0 );
		//var k2 = h * func( add(x0, hHalf), add(y0, mul(k1, 0.5)) );			
		//var k3 = h * func( add(x0, hHalf), add(y0, mul(k2, 0.5)) );			
		//var k4 = h * func( add(x0, h), add(y0, k3) );			



		var hHalf = h * 0.5;

		var k1 = mulConst(dy0, h);
		var m1 = mulConst(func(x0, y0, dy0), h);
		
		var k2 = mulConst(dy0 + m1*0.5, h);			
		var m2 = mulConst(func( add(x0, hHalf), add(y0, mulConst(k1, 0.5)), add(dy0, mulConst(m1, 0.5)) ), h);
		
		var k3 = mulConst(dy0 + m2*0.5, h);			
		var m3 = mulConst(func( add(x0, hHalf), add(y0, mulConst(k2, 0.5)), add(dy0, mulConst(m2, 0.5)) ), h);

		var k4 = mulConst(dy0 + m3, h);
		var m4 = mulConst(func( add(x0, h), add(y0, k3), add(dy0, m3) ), h);		


		var sumK = k1;
		sumK = add(sumK, mulConst(k2, 2.0));
		sumK = add(sumK, mulConst(k3, 2.0));
		sumK = add(sumK, k4);

		var sumM = m1;
		sumM = add(sumM, mulConst(m2, 2.0));
		sumM = add(sumM, mulConst(m3, 2.0));
		sumM = add(sumM, m4);

		return {
			x1: x0 + h,
			y1: y0 + sumK / 6.0,
			dy1: dy0 + sumM / 6.0
		};
	};
};


