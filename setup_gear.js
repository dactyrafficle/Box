(function() {

  let b, c;
  let g = [];
  
  window.addEventListener('load', function() {

    b = new Box();
    c = b.RETURN_CANVAS();
    container_special.appendChild(c);
    
    b.CANVAS_SIZE(500, 500);    // this is the number of pixels
    b.RANGE_X(0, 12);          // set the range in x
    b.RANGE_Y(0, 12);          // set the range in y 

    b.ADD_CLICK();
    b.ADD_MOUSEMOVE();

    b.CLEAR_CANVAS();
    
    // SHOW GRIDLINES 
    b.LINE_WIDTH(1);
    b.STROKE_STYLE('#ddd');
    b.SHOW_GRID_X();
    b.SHOW_GRID_Y();
   
    g[0] = new Gear({
      'val':{'x':6,'y':6},
      'n':23,
      'r_inner':3
    });
    b.GEAR(g[0]);
    
    g[1] = new Gear({
      'val':{'x':1,'y':1},
      'n':23,
      'r_inner':5,
      'r_outer':null,
      'theta':null,
      'd_theta':null,
      'omega':null,
      'line_width':2,
      'stroke_style':'#58da',
      'fill_style':null
    });
    b.GEAR(g[1]);
    
    window.setInterval(function(){
      b.CLEAR_CANVAS();
      
      // SHOW GRIDLINES 
      b.LINE_WIDTH(1);
      b.STROKE_STYLE('#ddd');
      b.SHOW_GRID_X();
      b.SHOW_GRID_Y();

      b.GEAR(g[0]);
      
      g[1].update();
      b.GEAR(g[1]);
      
    }, 1000/30);
    
    


  }); // CLOSING window.onload
})(); // CLOSING anon