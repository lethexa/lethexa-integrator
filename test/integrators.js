var assert = require('assert');
// var integrator = require('../lib/integrators');
var integrator = require((process.env.APP_DIR_FOR_CODE_COVERAGE || '../lib/') + 'integrators');



describe('Euler', function () {
  describe('#nextStep()', function () {
    it('should return the correct result', function () {
      var rungekutta = new integrator.Euler( function(x, y) {
        return 3 * x*x * y;
      }); 
      var input = {
        x: 1, 
        y: 2
      }
      var h = 0.1;

      var result = rungekutta.nextStep(input, h);

      assert.equal(1.1, Math.round(result.x * 1000.0) / 1000.0);
      assert.equal(2.6, Math.round(result.y * 1000.0) / 1000.0);
    });
  });
});

describe('RungeKutta1Order', function () {
  describe('#nextStep()', function () {
    it('should return the correct result', function () {
      var rungekutta = new integrator.RungeKutta1Order( function(x, y) {
        return 3 * x*x * y;
      }); 
      var input = {
        x: 1,
        y: 2
      };
      var h = 0.1;

      var result = rungekutta.nextStep(input, h);

      assert.equal(1.1, Math.round(result.x * 1000.0) / 1000.0);
      assert.equal(2.785, Math.round(result.y * 1000.0) / 1000.0);
    });
  });
});

describe('RungeKutta2Order', function () {
  describe('#nextStep()', function () {
    it('should return the correct result', function () {
      var rungekutta = new integrator.RungeKutta2Order( function(x, y, dy) {
        return dy + 2.0 * y;
      }); 
      var input = {
        x: 0, 
        y: 3,
        dy: 0
      }
      var h = 0.1;

      var result = rungekutta.nextStep(input, h);

      assert.equal(0.1, Math.round(result.x * 1000.0) / 1000.0);
      assert.equal(3.031, Math.round(result.y * 1000.0) / 1000.0);
      assert.equal(0.633, Math.round(result.dy * 1000.0) / 1000.0);
    });
  });
});


