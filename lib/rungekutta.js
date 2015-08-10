
/**
 * Runge-Kutta integration scheme
 * @class RungeKutta1Ord
 * @constructor
 * @param func {function(x, y)} The function to integrate.
 */
module.exports.RungeKutta1Ord = function( func ) {

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
