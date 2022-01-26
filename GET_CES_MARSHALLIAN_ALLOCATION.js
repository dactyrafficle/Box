
/*

  obj = {
    'delta':delta,
    'alpha':alpha,
    'beta':beta,
    'px':px,
    'py':py,
    'budget':budget
  }
  
*/

Box.prototype.GET_CES_MARSHALLIAN_ALLOCATION = function(obj) {

  let delta = obj.delta;
  let alpha = obj.alpha;
  let beta = obj.beta;
  let px = obj.px;
  let py = obj.py;
  let budget = obj.budget;

  let x = null;
  let y = null;
  let u = null;


  if (delta === 1) {
    
    // LINEAR
    if (alpha/px > beta/py) {
     x = budget / px;
     y = 0;
    }
    if (alpha/px < beta/py) {
     x = 0;
     y = budget / py;
    }
    if (alpha/px === beta/py) {
     x = 0;
     y = 0;
    }
    u = alpha*x + beta*y;
  } else if (delta === 0) {
    
    // COBB DOUGLAS
    x = alpha/(alpha + beta) * budget / px;
    y = beta/(alpha + beta) * budget / py;
    u = x**alpha*y**beta;
  } else if (delta < -100) {
    
    // LEONTIEFF
    x = budget / (px + alpha / beta * py);
    y = budget /(beta / alpha * px + py);
    u = (x/px > y/py ) ? (y/py) : (x/px);
  } else {
    
    // GENERAL CES : THE ELSE CASE
    let A = (alpha * py) / (beta * px);
    let B = A**(1 / (delta - 1));
    let B_INV = 1/ B;;
  
    x = budget / (px + py * B);
    y = budget / (px * B_INV + py);
    u = (alpha / delta) * x ** delta + (beta / delta) * y ** delta;
  }

  return {
    'x':x,
    'y':y,
    'u':u
  };

}