
/**
 * Runge-Kutta integration scheme for first order differential equations
 * @class RungeKutta1Ord
 * @constructor
 * @param func {function(x, y)} The function to integrate.
 */
module.exports.RungeKutta1Order = function( func ) {

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
		var k1 = func( x0, y0 );
		var k2 = func( add(x0, hHalf), add(y0, mul(hHalf, k1)) );			
		var k3 = func( add(x0, hHalf), add(y0, mul(hHalf, k2)) );			
		var k4 = func( add(x0, h), add(y0, mul(h, k3)) );			

		var sum = k1;
		sum = add(sum, mulConst(k2, 2.0));
		sum = add(sum, mulConst(k3, 2.0));
		sum = add(sum, k4);

		return {
			x1: x0 + h,
			y1: y0 + (h / 6.0) * (sum)
		};
	};
};



/**
 * Runge-Kutta integration scheme for second order differential equations
 * @class RungeKutta2Ord
 * @constructor
 * @param func {function(x, y, dy)} The function to integrate.
 */
module.exports.RungeKutta2Order = function( func ) {

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
	this.nextStep = function( x0, y0, dy0, h ) {
		var hHalf = h * 0.5;
		var k1 = h * dy0;
		var m1 = func(x0, y, dy0);
		
		var k2 = func( add(x0, hHalf), add(y0, mul(hHalf, k1)) );			
		var m2 = func( add(x0, hHalf), add(y0, mul(k1, 0.5), add(dy0, m1*0.5) );
		
		var k3 = func( add(x0, hHalf), add(y0, mul(hHalf, k2)) );			
		var m3 = func( add(x0, hHalf), add(y0, mul(k2, 0.5), add(dy0, m2*0.5) );

		var k4 = func( add(x0, h), add(y0, mul(h, k3)) );
		var m4 = func( add(x0, h), add(y0, k3), add(dy0, m3) );		


		var sumK = k1;
		sumK = add(sumK, mulConst(k2, 2.0));
		sumK = add(sumK, mulConst(k3, 2.0));
		sumK = add(sumK, k4);

		sumM = m1;
		sumM = add(sumM, mulConst(m2, 2.0));
		sumM = add(sumM, mulConst(m3, 2.0));
		sumM = add(sumM, m4);

		return {
			x1: x0 + h,
			y1: y0 + (h / 6.0) * sumK,
			dy1: dy0 + (h / 6.0) * sumM
		};
	};
};


