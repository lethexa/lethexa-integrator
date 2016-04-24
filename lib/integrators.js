/* global exports */

(function (exports) {
  'use strict';

  var MathFunc = {
    add: function(x, y) {
      return x + y;
    },
    mulScalar: function(x, c) {
      return x * c;
    }
  };


  /**
   * Euler integration scheme for first order differential equations
   * @class Euler
   * @constructor
   * @param func {function(x, y)} The function to integrate.
   * @param ops {Object} The operations add, mul, ...
   */
  exports.Euler = function( func, ops ) {
    ops = ops || MathFunc;

    /**
     * Calculates the next step via integrating
     * @method nextStep
     * @param input {Object} The input values for x and y
     * @example
       {
         x: <value>,
         y: <value>
       }
     * @param h {Number} The stepsize
     * @return {object} Returns an object in the following form 
     * @example
       {
         x: <value>,
         y: <value>
       }
     */
    this.nextStep = function( input, h ) {
      var x0 = input.x;
      var y0 = input.y;
      return {
        x: x0 + h,
        y: ops.add(y0, ops.mulScalar(func(x0, y0), h))
      };
    };
  };


  /**
   * Runge-Kutta integration scheme for first order differential equations.
   * @class RungeKutta1Order
   * @constructor
   * @param func {function(x, y)} The function to integrate.
   * @param ops {Object} The operations add, mul, ...
   */
  exports.RungeKutta1Order = function( func, ops ) {
    ops = ops || MathFunc;

    /**
     * Calculates the next step via integrating
     * @method nextStep
     * @param input {Object} The input values for x and y
     * @example
       {
         x: <value>,
         y: <value>
       }
     * @param h {Number} The stepsize
     * @return {object} Returns an object in the following form 
     * @example
       {
         x: <value>,
         y: <value>
       }
     */
    this.nextStep = function( input, h ) {
      var x0 = input.x;
      var y0 = input.y;
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
        x: x0 + h,
        y: ops.add(y0, ops.mulScalar(sumK, 1.0 / 6.0))
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
  exports.RungeKutta2Order = function( func, ops ) {
    ops = ops || MathFunc;

    /**
     * Calculates the next step via integrating
     * @method nextStep
     * @param input {Object} The input values for x, y and dy
     * @example
       {
         x: <value>,
         y: <value>,
         dy: <value>
       }
     * @param h {Number} The stepsize
     * @return {object} Returns an object in the following form 
     * @example
       {
         x: <value>,
         y: <value>,
         dy: <value>
       }
     */
    this.nextStep = function( input, h ) {
      var x0 = input.x;
      var y0 = input.y;
      var dy0 = input.dy;
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
        x: x0 + h,
        y: ops.add(y0, ops.mulScalar(sumK, 1.0 / 6.0)),
        dy: ops.add(dy0, ops.mulScalar(sumM, 1.0 / 6.0))
      };
    };
  };

})(typeof exports === 'undefined' ? this.integrator = {} : exports);

