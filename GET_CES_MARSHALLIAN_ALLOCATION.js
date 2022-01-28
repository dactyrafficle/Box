
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
  let delta_inv = 1/delta;
  let alpha = obj.alpha;
  let beta = obj.beta;
  let px = obj.px;
  let py = obj.py;
  let budget = obj.budget;
  
  let A = (alpha * py) / (beta * px);
  let B = (1 / (delta - 1));
  let C = A**B;
  let C_INV = 1/ C;
  
  let x = null;
  let y = null;
  let u = null;

  let output = {
    'delta':delta,
    'alpha':alpha,
    'beta':beta,
    'px':px,
    'py':py,
    'budget':budget,
    'A':A,
    'B':B,
    'C':C
  }
  
  // SET THE TYPE
  output.type = "ELSE";
  output.delta = delta;
  if (delta === 1) {output.type = "LINEAR"}
  if (delta === 0) {output.type = "LOG"}
  if (delta < -100) {output.type = "LEONTIEFF"}

  if (output.type === "LINEAR") {
    
    if (alpha/px > beta/py) {
     x = budget / px;
     y = 0;
    }
    if (alpha/px < beta/py) {
     x = 0;
     y = budget / py;
    }
    if (alpha/px === beta/py) {
     x = 0.5 * budget / px;
     y = 0.5 * budget / py;
    }
    u = alpha*x + beta*y;

    output.expenditure = px*x + py*y; 
    output.u = u;
    output.x = x;
    output.y = y;
    return output;
  }
 
  if (output.type === "LOG") {
    
    x = alpha/(alpha + beta) * budget / px;
    y = beta/(alpha + beta) * budget / py;
    u = x**alpha*y**beta;
    
    // i should find a way to add both log utilities, because sometimes we use cobb douglas, other times the log version
    
    output.expenditure = px*x + py*y; 
    output.u = u;
    output.x = x;
    output.y = y;
    return output;
  } 
  
  if (output.type === "LEONTIEFF") {
    
    x = budget / (px + py);
    y = budget / (px + py);
    u = x;
    
    output.expenditure = px*x + py*y; 
    output.u = u;
    output.x = x;
    output.y = y;
    return output;
  } 
  

  // GENERAL CES : THE ELSE CASE
  x = budget / (px + py * C);
  y = budget / (px * C_INV + py);
  u = (alpha * x ** delta + beta * y ** delta)**delta_inv;

  output.expenditure = px*x + py*y; 
  output.u = u;
  output.x = x;
  output.y = y;
  return output;

}