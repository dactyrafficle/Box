(function() {

  let b, c;
  window.addEventListener('load', function() {
  
    myinput = document.getElementById('myinput');

    b = new Box();
    c = b.RETURN_CANVAS();
    container_ces_idc.appendChild(c);
    
    b.CANVAS_SIZE(500, 500);    // this is the number of pixels
    b.RANGE_X(-1, 11);          // set the range in x
    b.RANGE_Y(-1, 11);          // set the range in y 

    b.CLEAR_CANVAS();

    b.LINE_WIDTH(1);
    b.STROKE_STYLE('#ddd');
    b.SHOW_GRID_X();
    b.SHOW_GRID_Y();
    
    b.LINE_WIDTH(2);
    b.STROKE_STYLE('#999');
    b.SHOW_AXES();

    b.LINE_WIDTH(2);
    b.STROKE_STYLE('#5d5d');
    b.SHOW_CES_INDIFFERENCE_CURVE({
      'delta':0.5,
      'alpha':0.5,
      'beta':0.5,
      'u':null,
      'x':b.data.range.x.avg,
      'y':b.data.range.y.avg
    });
    
    b.LINE_WIDTH(2);
    b.STROKE_STYLE('#fc0a');
    b.SHOW_CES_INDIFFERENCE_CURVE({
      'delta':0,
      'alpha':0.5,
      'beta':0.5,
      'u':null,
      'x':b.data.range.x.avg,
      'y':b.data.range.y.avg
    });

    b.LINE_WIDTH(2);
    b.STROKE_STYLE('#adc2eb');
    b.SHOW_CES_INDIFFERENCE_CURVE({
      'delta':-1,
      'alpha':0.5,
      'beta':0.5,
      'u':null,
      'x':b.data.range.x.avg,
      'y':b.data.range.y.avg
    });
    
    b.LINE_WIDTH(2);
    b.STROKE_STYLE('#f93a');
    b.SHOW_CES_INDIFFERENCE_CURVE({
      'delta':-10,
      'alpha':0.5,
      'beta':0.5,
      'u':null,
      'x':b.data.range.x.avg,
      'y':b.data.range.y.avg
    });

  }); // CLOSING window.onload
})(); // CLOSING anon