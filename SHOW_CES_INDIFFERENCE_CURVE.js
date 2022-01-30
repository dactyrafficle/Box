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
    'type':null,
    'procedure':null,
    'x_c':null,
    'y_c':null
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
  if (delta < -100) {output.type = "LEONTIEFF"}

  output.alpha = alpha;
  output.beta = beta;



  // FIRST CASE :
  if (obj.u) {
    output.procedure = 'FIRST : UTILITY SUPPLIED';
    u = obj.u;
    output.u = u;
  }

  // SECOND CASE : u = u(x, y)
  if (!obj.u && (obj.x && obj.y)) {
    output.procedure = 'SECOND : UTILITY CALCULATED';

    x = obj.x;
    y = obj.y;
    
    output.x = x;
    output.y = y;
    
    if (output.type === "LINEAR") {u = alpha*x + beta*y};
    if (output.type === "LOG") {u = x**alpha*y**beta};
    if (output.type === "LEONTIEFF") {u = 1};
    if (output.type === "ELSE") {u = (alpha*x**delta + beta*y**delta)**delta_inv};
    
    output.u = u;
  }

  // THE DOMAIN AND RANGE
  output.x_c = null;
  output.y_c = null;
  if (delta !== 0) {
    output.x_c = (u**delta / alpha)**(1/delta);
    output.y_c = (u**delta / beta)**(1/delta);
  }
  

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

  for (let i = 0; i <= 50; i++) {
  
    // THE ELSE CASES **
    if (u**delta >= alpha*x_index**delta) {
      arr.push({
        'x':x_index,
        'y':((u**delta - alpha*x_index**delta)/beta)**delta_inv
      });
    }
    if (u**delta >= beta*y_index**delta) {
      arr.push({
        'x':((u**delta - beta*y_index**delta)/alpha)**delta_inv,
        'y':y_index
      });
    }

    x_index += dx;
    y_index += dy;
    
  } // closing for loop

  // ** WE USE X AND Y TO GET U
  // WHEN WE LOOP OVER X, WE TRY TO FIND THE Y THAT MAKES THAT U POSSIBLE
  // BUT WHEN DELTA < 0, THERE MAY NOT EXIST SUCH A Y
  // AND IF WE DONT ACCOUNT FOR THAT, WELL GET WACKY VALUES PLOTTED ON THE GRAPH
  
  // WHEN 0 < DELTA < 1, WE CAN FIND AN X OR Y TO MEET U GIVEN AND OTHER Y OR X, RESPECTIVELY 

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

