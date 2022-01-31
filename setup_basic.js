(function() {

  let b, c;
  window.addEventListener('load', function() {

    b = new Box();
    c = b.RETURN_CANVAS();
    container_basic.appendChild(c);
    
    b.CANVAS_SIZE(500, 500);    // this is the number of pixels
    b.RANGE_X(-1, 11);          // set the range in x
    b.RANGE_Y(-1, 11);          // set the range in y 

    b.ADD_CLICK();
    b.ADD_MOUSEMOVE();

    b.CLEAR_CANVAS();

    b.LINE_WIDTH(1);
    b.STROKE_STYLE('#ddd');
    b.SHOW_GRID_X();
    b.SHOW_GRID_Y();
    
    b.LINE_WIDTH(2);
    b.STROKE_STYLE('#999');
    b.SHOW_AXES();
    
    // DRAW POINTS
    b.FILL_STYLE('#fc0a');
    b.RADIUS(2);
    b.SHOW_VALUE({'x':3,'y':4.5});
    b.POINT({'x':4,'y':4.0});
    b.POINT({'x':5,'y':3.5});
    
    // DRAW A LINE
    b.LINE_WIDTH(3);
    b.STROKE_STYLE('#9fdfbf');
    b.CONNECT_VALUES([
      {'x':5,'y':9},
      {'x':9,'y':8},
      {'x':9,'y':7}
    ]);
    
    // DRAW A RECTANGLE
    b.LINE_WIDTH(3);
    b.STROKE_STYLE('#99b3e6');
    b.FILL_STYLE('#d6e0f555');
    b.RECT([
      {'x':7,'y':3},
      {'x':9,'y':7},
      {'x':5,'y':9},
      {'x':3,'y':5}
    ]);
    
    // DRAW A SHAPE
    b.LINE_WIDTH(2);
    b.STROKE_STYLE('#fbfa');
    b.FILL_STYLE('#fcf5');
    b.SHAPE([
      {'x':1,'y':8},
      {'x':2,'y':10},
      {'x':4,'y':9},
      {'x':2,'y':5}
    ]);

  }); // CLOSING window.onload
})(); // CLOSING anon