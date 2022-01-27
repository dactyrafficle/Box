/*

 u = (alpha/delta)*x**delta + (beta/delta)*y**delta

 b.LINE_WIDTH(2);
 b.STROKE_STYLE('#fc0a');
 b.SHOW_INDIFFERENCE_CURVE({
  'delta':delta,    // must
  'alpha':alpha,    // must
  'beta':beta,      // must
  'u':u,            // SUFFICIENT, BUT NOT A GOOD CHOICE
  'x':x,            // PREFERABLE
  'y':y,            // PREFERABLE
  'budget':budget,
  'px':px,
  'py':py
 }

  FUNCTIONAL DEPENDENCIES
  1. this.CONNECT_VALUES(arr)
  
  ** verifiy when delta = -145
  ** verify when delta = 0.90 to 1.00
  ** when does it hit 0 ?
  ** how arr is sorted

*/


// CES UTILITY
Box.prototype.SHOW_CES_INDIFFERENCE_CURVE = function(obj) {
  
  // must have delta, alpha and beta
  // and we must either have utility, or have what we need to make utility
 
  // case 1 : u
  // case 2 : u=f(x,y)
  // case 3 : x=f(M,px) : y=f(M,py) : u=f(x,y)
  
  
  // I CAN ACTUALLY SET DELTA = -INFINITY HERE, SO I SHOULD DO THAT, IF DELTA IS BIG ENOUGH
  
  
  let output = {
    'delta':null,
    'alpha':null,
    'beta':null,
    'u':null,
    'x':null,            
    'y':null,
    'budget':null,
    'px':null,
    'py':null,
    'type':null
  }

  let alpha, alpha_inv, beta, beta_inv, delta, delta_inv;
  let u;
  let x, y;

  alpha = obj.alpha;
  alpha_inv = 1/alpha;
  beta = obj.beta;
  beta_inv = 1/beta;
  
  delta = obj.delta;
  delta_inv = (1/delta || null);
  
  // SET THE TYPE
  output.type = "ELSE";
  output.delta = delta;
  if (delta === 1) {output.type = "LINEAR"}
  if (delta === 0) {output.type = "LOG"}
  if (delta === -Infinity) {output.type = "LEONTIEFF"}


  // FIRST CASE :
  if (obj.u) {
    u = obj.u;
  }

  // SECOND CASE : u = u(x, y)
  if (!obj.u && (obj.x && obj.y)) {
    x = obj.x;
    y = obj.y;
    
    if (output.type === "LINEAR") {u = alpha*x + beta*y};
    if (output.type === "LOG") {u = x**alpha*y**beta};
    if (output.type === "LEONTIEFF") {u = 1};
    if (output.type === "ELSE") {u = (alpha*x**delta + beta*y**delta)**delta_inv};
    
    output.u = u;
  }

  /*
  // THIRD CASE : u = u(budget, px, py)
  let budget, px, py;
  if (!obj.u && (!obj.x && !obj.y) && (obj.budget && obj.px && obj.py)) {
    // console.log('CASE 3');
    budget = obj.budget;
    px = obj.px;
    py = obj.py;
    x = alpha*budget/px;
    y = beta*budget/py;
    u = x**alpha*y**beta;
  }
  */
  
    
 
  // LINEAR
  if (output.type === "LINEAR") {
    this.CONNECT_VALUES([{'x':0,'y':u/beta},{'x':u/alpha,'y':0}]);
    return output;
  }
  
  // LEONTIEFF
  if (output.type === "LEONTIEFF") {

    return output;
  }
  
  // FOR LOG AND ELSE, WE LOOP
  
  // LOG
  if (output.type === "LOG") {
    
    let dx = this.data.range.x.max / 50;
    let dy = this.data.range.y.max / 50;
    let x_index = 0, y_index = 0;
    let arr = [];

    for (let i = 0; i <= 50; i++) {

      arr[2*i+0] = {'x':x_index,'y':(u/x_index**alpha)**beta_inv};
      arr[2*i+1] = {'x':(u/x_index**beta)**alpha_inv,'y':y_index};
      x_index += dx;
      y_index += dy;
    }
    
    arr.sort(function(a,b) {
      return a.x - b.x;
    });
    this.CONNECT_VALUES(arr);
    return output;
  }

  



  // THIS IS THE ELSE CASE

  let dx = this.data.range.x.max / 50;
  let dy = this.data.range.y.max / 50;
  
  let x_index = 0, y_index = 0;
  let arr = [];

  // RELEVANT IF delta < 0
  let x_min = 0; //((delta/alpha)*u)**delta_inv; ************************************
  let y_min = 0; // ((delta/beta)*u)**delta_inv;

  for (let i = 0; i <= 50; i++) {
  
    // THE ELSE CASES
    if (delta > 0 || (delta < 0 && x_index >= x_min)) {
      arr.push({
        'x':x_index,
        'y':((u**delta - alpha*x_index**delta)/beta)**delta_inv
      });
    }

    if (delta > 0 || (delta < 0 && y_index >= y_min)) {
      arr.push({
        'x':((u**delta - beta*y_index**delta)/alpha)**delta_inv,
        'y':y_index
      });
    }

    x_index += dx;
    y_index += dy;
    
  } // closing for loop

  

  // SORT ARR : IF it's almost leontieff, ie. delta = -140, then i need to sort by x, from lowest to smallest
  // then if 2 points have the same x, from largest to smallest
  arr.sort(function(a,b) {
    return a.x - b.x;
  });
  // console.log(arr);
  
  // once its sorted, it might make sense to do a while loop to remove the points at the beginning and end that are not on the grid
  
  this.CONNECT_VALUES(arr);
  
  return output;
}

