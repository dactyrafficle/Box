
/*

let b = new Box();
let c = b.RETURN_CANVAS();
container.appendChild(c);

b.CANVAS_SIZE(500, 500);     // this is the number of pixels 
b.RANGE_X(-1, 11);           // set the range in x 
b.RANGE_Y(-1, 11);           // set the range in y 

b.ADD_CLICK();       // updates b.current.click.value and .pixel 
b.ADD_MOUSEMOVE();   // updates b.current.mousemove.value and .pixel 

*/

function Box() {

  this.c = document.createElement('canvas');
  this.ctx = this.c.getContext('2d');
  
  // STRUCTURE
  this.data = {
    'dimension':{
      'w':100,
      'h':100
    },
    'zoom':{
      'x':0,
      'y':0
    },
    'translate':{
      'x':0,
      'y':0
    },
    'range':{
      'x':{
        'min':0,
        'max':0,
        'avg':0,
        'span':0
      },
      'y':{
        'min':0,
        'max':0,
        'avg':0,
        'span':0
      }
    }
  };

  this.current = {
    'click':{
      'value':null,
      'pixel':null
    },
    'mousemove':{
      'value':null,
      'pixel':null
    }
  }
     
  // DEFAULTS
  this.RANGE_X(0, 100);
  this.RANGE_Y(0, 100);
  this.CANVAS_SIZE(100, 100);

  this.ctx.radius = 3;
}



Box.prototype.RETURN_CANVAS = function() {
  return this.c;
}

Box.prototype.CLEAR_CANVAS = function() {
  this.ctx.fillStyle = '#fff';
  this.ctx.beginPath();
  this.ctx.rect(0, 0, this.data.dimension.w, this.data.dimension.h);
  this.ctx.fill();
}

Box.prototype.ADD_CLICK = function() {
  this.c.addEventListener('click', function(e) {
    let pixel = {'x':e.offsetX,'y':e.offsetY};
    let value = this.PIXEL2VALUE(pixel);
    this.current.click.value = value;
    this.current.click.pixel = pixel;
  }.bind(this));
}

Box.prototype.ADD_MOUSEMOVE = function() {
  this.c.addEventListener('mousemove', function(e) {
    let pixel = {'x':e.offsetX,'y':e.offsetY};
    let value = this.PIXEL2VALUE(pixel);
    this.current.mousemove.value = value;
    this.current.mousemove.pixel = pixel;
  }.bind(this));
}

Box.prototype.RANGE_X = function(min, max) {
  this.data.range.x.min = min;
  this.data.range.x.max = max;
  this.data.range.x.avg = (max + min) / 2;
  this.data.range.x.span = max - min;
  this.data.zoom.x = this.data.dimension.w / this.data.range.x.span;
  this.data.translate.x = -this.data.range.x.min;
}

Box.prototype.RANGE_Y = function(min, max) {
  this.data.range.y.min = min;
  this.data.range.y.max = max;
  this.data.range.y.avg = (max + min) / 2;
  this.data.range.y.span = max - min;
  this.data.zoom.y = this.data.dimension.h / this.data.range.y.span;
  this.data.translate.y = -this.data.range.y.min;
}

Box.prototype.CANVAS_SIZE = function(w, h) {
  this.data.dimension.w = w;
  this.data.dimension.h = h; 
  this.c.width = this.data.dimension.w;
  this.c.height = this.data.dimension.h;
 
  // RESET ZOOM AND XLATE BY REAPPLYING THE RANGES
  this.RANGE_X(this.data.range.x.min, this.data.range.x.max);
  this.RANGE_Y(this.data.range.y.min, this.data.range.y.max);
}

Box.prototype.VALUE2PIXEL = function(val) {  // val : (0,0) is bottom-left (Cartesian)
 return {
  'x':(val.x+this.data.translate.x)*this.data.zoom.x,
  'y':this.data.dimension.h - (val.y+this.data.translate.y)*this.data.zoom.y
 }
}

Box.prototype.PIXEL2VALUE = function(pixel) { // pixel : (0,0) is top-left (standard computer/matrix grid)
 return {
  'x':(pixel.x/this.data.zoom.x)-this.data.translate.x,
  'y':(this.data.dimension.h-pixel.y)/this.data.zoom.y-this.data.translate.y
 }
}

Box.prototype.STROKE_STYLE = function(x) {
  this.ctx.strokeStyle = x; 
}

Box.prototype.FILL_STYLE = function(x) {
  this.ctx.fillStyle = x; 
}

Box.prototype.LINE_WIDTH = function(x) {
  this.ctx.lineWidth = x; 
}
Box.prototype.RADIUS = function(rx) {
  this.ctx.radius = rx;
}

Box.prototype.SHOW_GRID_X = function() {

  let dx = 1;
  let x_start = Math.floor(this.data.range.x.min/dx)*dx;
  
  let y0 = this.data.range.y.min;
  let y1 = this.data.range.y.max;

  for (let x = x_start; x <= this.data.range.x.max; x += dx) {

    let v0 = {'x':x,'y':y0};
    let v1 = {'x':x,'y':y1};

    let p0 = this.VALUE2PIXEL(v0);
    let p1 = this.VALUE2PIXEL(v1);
    
    this.ctx.beginPath();
    this.ctx.moveTo(p0.x, p0.y);
    this.ctx.lineTo(p1.x, p1.y);
    this.ctx.stroke();   
  }
};


Box.prototype.SHOW_GRID_Y = function() {

  let dy = 1;
  let y_start = Math.floor(this.data.range.y.min/dy)*dy;
  
  let x0 = this.data.range.x.min;
  let x1 = this.data.range.x.max;

  for (let y = y_start; y <= this.data.range.y.max; y += dy) {

    let v0 = {'x':x0,'y':y};
    let v1 = {'x':x1,'y':y};

    let p0 = this.VALUE2PIXEL(v0);
    let p1 = this.VALUE2PIXEL(v1);
    
    this.ctx.beginPath();
    this.ctx.moveTo(p0.x, p0.y);
    this.ctx.lineTo(p1.x, p1.y);
    this.ctx.stroke();   
  }
};

Box.prototype.SHOW_AXES = function() {
  
  this.CONNECT_VALUES([
    {'x':0,'y':this.data.range.y.min},
    {'x':0,'y':this.data.range.y.max}
  ]);
  
  this.CONNECT_VALUES([
    {'x':this.data.range.x.min,'y':0},
    {'x':this.data.range.x.max,'y':0}
  ]);
  
}


// val = {'x':x,'y':y}
Box.prototype.VALUE_IN_RANGE = function(val) {

  let x = val.x;
  let y = val.y;
  let x_min = this.data.range.x.min;
  let x_max = this.data.range.x.max; 
  let y_min = this.data.range.y.min;
  let y_max = this.data.range.y.max;
  
  if (x >= x_min && x <= x_max && y >= y_min && y <= y_max) {
    return true;
  } else {
    return false;
  }

}

// val = {'x':x,'y':y}
// this.RADIUS(rx)
Box.prototype.SHOW_VALUE = function(val) {

 let p = this.VALUE2PIXEL(val);
 let r = this.ctx.radius;

 this.ctx.beginPath();
 this.ctx.arc(p.x, p.y, r*2, 0, 2*Math.PI);
 this.ctx.fill();
}
Box.prototype.POINT = function(val) {

 let p = this.VALUE2PIXEL(val);
 let r = this.ctx.radius;

 this.ctx.beginPath();
 this.ctx.arc(p.x, p.y, r*2, 0, 2*Math.PI);
 this.ctx.fill();
}

// let vals = [val, val,...]
// b.LINE_WIDTH(1);
// b.STROKE_STYLE('#ddd');

Box.prototype.CONNECT_VALUES = function(vals) {

  for (let i = 0; i < vals.length-1; i++) {

    let pixel_0 = this.VALUE2PIXEL(vals[i+0]);
    let pixel_1 = this.VALUE2PIXEL(vals[i+1]);

    this.ctx.beginPath();
    this.ctx.moveTo(pixel_0.x, pixel_0.y);
    this.ctx.lineTo(pixel_1.x, pixel_1.y);
    this.ctx.stroke();
  }
}
Box.prototype.CONNECT_POINTS = function(vals) {

  for (let i = 0; i < vals.length-1; i++) {

    let pixel_0 = this.VALUE2PIXEL(vals[i+0]);
    let pixel_1 = this.VALUE2PIXEL(vals[i+1]);

    this.ctx.beginPath();
    this.ctx.moveTo(pixel_0.x, pixel_0.y);
    this.ctx.lineTo(pixel_1.x, pixel_1.y);
    this.ctx.stroke();
  }
}
// let vals = [val, val,...]
Box.prototype.RECT = function(vals) {

  let pixels = [];
  pixels[0] = this.VALUE2PIXEL(vals[0]);
  pixels[1] = this.VALUE2PIXEL(vals[1]);
  pixels[2] = this.VALUE2PIXEL(vals[2]);
  pixels[3] = this.VALUE2PIXEL(vals[3]);

  this.ctx.beginPath();
  this.ctx.moveTo(pixels[0].x, pixels[0].y);
  this.ctx.lineTo(pixels[1].x, pixels[1].y);
  this.ctx.lineTo(pixels[2].x, pixels[2].y);
  this.ctx.lineTo(pixels[3].x, pixels[3].y);
  this.ctx.lineTo(pixels[0].x, pixels[0].y);
  this.ctx.stroke();
  this.ctx.fill();

}






/*

Box.prototype.SHOW_GRID_Y = function(dy, color_string) {

 if (!arguments[0]) {
  dy = 1;
 }
 let zoom = this.data.zoom.x;

 let alpha_1 = 1;
 let alpha_10 = -1+(1/200)*zoom;
 let color_string_1 = 'rgba(208, 208, 208,' + alpha_1 + ')';
 let color_string_10 = 'rgba(224, 224, 224,' + alpha_10 + ')';
 
 let y_start = Math.floor(this.data.range.y.min/dy)*dy;
 
 let i = 0; // BC FLOATING PT NUMBERS MAKE TESTING y%dy TOUGH
 for (let y = y_start; y < this.data.range.y.max; y += dy/10) {
  if (i%10===0) {
   let val0 = {'x':this.data.range.x.min,'y':y};
   let val1 = {'x':this.data.range.x.max,'y':y};
   this.CONNECTVALUES(val0, val1, color_string_1);
  }

   let val0 = {'x':this.data.range.x.min,'y':y};
   let val1 = {'x':this.data.range.x.max,'y':y};
  this.CONNECTVALUES(val0, val1, color_string_10);
  i++;
 }
};
*/
Box.prototype.SHOW_FLOATING_LOG_Y_AXIS = function(n) {

 let sh = this.data.dimension.h/n;
 let sw = this.data.dimension.w/n;
 
 let n1 = 2;
 let n2 = n-n1;
 
 let p0 = {'x':sw,'y':sh*n1};
 let p1 = {'x':sw,'y':sh*n2};

 let v0 = this.PIXEL2VALUE(p0);
 let v1 = this.PIXEL2VALUE(p1);
 
 this.CONNECTVALUES(v0, v1, '#333', 0.5);
 
 let dsw = 5;
 for (let i = n1; i <= n2; i++) {
  let p0 = {'x':sw-dsw,'y':sh*i};
  let p1 = {'x':sw+dsw,'y':sh*i};
  let v0 = this.PIXEL2VALUE(p0);
  let v1 = this.PIXEL2VALUE(p1);
  this.ctx.fillStyle = '#333';
  this.ctx.textAlign = 'left';
  this.ctx.textBaseline = 'middle';
  this.ctx.fillText((v0.y).toFixed(2), p0.x+2*dsw, p0.y);
  this.ctx.stroke();
  this.CONNECTVALUES(v0, v1, '#333', 0.5); 
 }
 
 this.ctx.textAlign = 'center';
 this.ctx.textBaseline = 'top';
 this.ctx.fillText((this.data.label.y), sw, sh*(n2+0.25));
};
Box.prototype.SHOW_FLOATING_Y_AXIS = function(n) {

 let sh = this.data.dimension.h/n;
 let sw = this.data.dimension.w/n;
 
 let n1 = 2;
 let n2 = n-n1;
 
 let p0 = {'x':sw,'y':sh*n1};
 let p1 = {'x':sw,'y':sh*n2};

 let v0 = this.PIXEL2VALUE(p0);
 let v1 = this.PIXEL2VALUE(p1);
 
 this.CONNECTVALUES(v0, v1, '#333', 0.5);
 
 let dsw = 5;
 for (let i = n1; i <= n2; i++) {
  let p0 = {'x':sw-dsw,'y':sh*i};
  let p1 = {'x':sw+dsw,'y':sh*i};
  let v0 = this.PIXEL2VALUE(p0);
  let v1 = this.PIXEL2VALUE(p1);
  this.ctx.fillStyle = '#333';
  this.ctx.textAlign = 'left';
  this.ctx.textBaseline = 'middle';
  this.ctx.fillText((v0.y).toFixed(0), p0.x+2*dsw, p0.y);
  this.ctx.stroke();
  this.CONNECTVALUES(v0, v1, '#333', 0.5); 
 }
 
 this.ctx.textAlign = 'center';
 this.ctx.textBaseline = 'top';
 this.ctx.fillText((this.data.label.y), sw, sh*(n2+0.25));
};
Box.prototype.SHOW_FLOATING_LOG_X_AXIS = function(obj) {

  
  let x_0 = this.data.dimension.w/5;
  let x_1 = this.data.dimension.w - x_0;
  let y = this.data.dimension.h*4.25/5;
  
  let val_0 = this.PIXEL2VALUE({'x':x_0,'y':y});
  let val_1 = this.PIXEL2VALUE({'x':x_1,'y':y});
  
  this.CONNECT_VALUES({
    'vals':[val_0, val_1]
  })



  
};

Box.prototype.SHOW_FLOATING_X_AXIS = function(n, y_val) {

 //let n = 9;
 let sh;
 if (arguments[1] !== null) {
   let v3 = {'x':0,'y':y_val};
   //console.log(v3);
   let p3 = this.VALUE2PIXEL(v3);
   //console.log(p3);
   sh = p3.y;
 } else {
  sh = this.data.dimension.h - this.data.dimension.h/n;
 }
 let sw = this.data.dimension.w/n;
 
 let n1 = 1;
 let n2 = n-n1;
 
 let p0 = {'x':sw*n1,'y':sh};
 let p1 = {'x':sw*n2,'y':sh};

 let v0 = this.PIXEL2VALUE(p0);
 let v1 = this.PIXEL2VALUE(p1);
 
 this.CONNECTVALUES(v0, v1, '#333', 0.5);

 let dsh = 5;
 for (let i = n1; i <= n2; i++) {
  let p0 = {'x':sw*i,'y':sh+dsh};
  let p1 = {'x':sw*i,'y':sh-dsh};
  let v0 = this.PIXEL2VALUE(p0);
  let v1 = this.PIXEL2VALUE(p1);
  this.ctx.fillStyle = '#333';
  this.ctx.textAlign = 'center';
  this.ctx.textBaseline = 'top';
  this.ctx.fillText(v0.x.toFixed(0), p0.x, p0.y+1*dsh);
  this.ctx.stroke();
  this.CONNECTVALUES(v0, v1, '#333', 0.5); 
 }

 this.ctx.textAlign = 'right';
 this.ctx.textBaseline = 'middle';
 this.ctx.fillText((this.data.label.x).toUpperCase(), sw*(n1-0.25), sh);

};

Box.prototype.showAxes = function(fontSize) {

 this.ctx.strokeStyle = '#333';
 
 this.ctx.beginPath();
 this.ctx.moveTo((0+this.data.translate.x)*this.data.zoom.x, 0);
 this.ctx.lineTo((0+this.data.translate.x)*this.data.zoom.x, this.data.dimension.h);
 this.ctx.stroke();
 
 // X AXIS LABEL
 this.ctx.fillStyle = '#333';
 this.ctx.font = (fontSize || this.data.dimension.w/100*3.5) + 'px Monospace';
 this.ctx.textAlign = "right";
 this.ctx.fillText(this.data.label.y, this.data.translate.x*this.data.zoom.x-this.data.dimension.w/100, this.data.translate.x*this.data.zoom.x);

 // X AXIS
 this.ctx.beginPath();
 this.ctx.moveTo(0, this.data.dimension.h-(0+this.data.translate.y)*this.data.zoom.y);
 this.ctx.lineTo(this.data.dimension.w, this.data.dimension.h-(0+this.data.translate.y)*this.data.zoom.y);
 this.ctx.stroke();
 
 // X AXIS LABEL
 this.ctx.fillStyle = '#333';
 this.ctx.font = (fontSize || this.data.dimension.w/100*3.5) + 'px Monospace';
 this.ctx.fillText(this.data.label.x, this.data.dimension.w-(0+this.data.translate.x)*this.data.zoom.x, this.data.dimension.h-(0+this.data.translate.y)*this.data.zoom.y + this.data.dimension.w/100*3.5);
 
 if (this.data.range.x.min > 0 || this.data.range.y.min > 0) {
    
   let s0 = 20;
   let s1 = 30;
   let c = {'x':s0, 'y':this.data.dimension.h - s0};

   let o = this.PIXEL2VALUE(c);
   console.log(o);
   console.log(this.data.range);
   
   let cx = {'x':s0+s1, 'y':this.data.dimension.h - s0};
   let ox = this.PIXEL2VALUE(cx);
   
   let cy = {'x':s0, 'y':this.data.dimension.h - (s0 + s1)};
   
    this.ctx.lineWidth = 1;
    this.ctx.strokeStyle = '#333';
    this.ctx.beginPath();
    this.ctx.moveTo(cy.x, cy.y);
    this.ctx.lineTo(c.x, c.y);
    this.ctx.lineTo(cx.x, cx.y);
    this.ctx.stroke();
    
    this.ctx.textAlign = 'start';
    this.ctx.textBaseline = 'top';
    this.ctx.fillText(this.data.label.x, c.x, c.y+2);
    
    
 }
 
}

Box.prototype.SHOWVALUE = function(val, colorstring, rx) {

 let pixel = this.VALUE2PIXEL(val);
 this.ctx.fillStyle = colorstring;
 this.ctx.beginPath();
 this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
 this.ctx.fill();

}

/*
  let obj = {
    'val':val,
    'color_string':'#999',
    'rx':3
  }
*/
Box.prototype.DRAW_POINT = function(obj) {

  let color_string = '#999';
  let rx = 3;
  if (obj.color_string) {color_string = obj.color_string}
  if (obj.rx) {rx = obj.rx}

  let pixel = this.VALUE2PIXEL(obj.val);
  this.ctx.fillStyle = color_string;
  this.ctx.beginPath();
  this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
  this.ctx.fill();
}

/*
let obj = {
  'val':val,
  'slope':slope,
  'color_string':'#999',
  'line_width':2
}
*/
Box.prototype.DRAW_LINE = function(obj) {


  let line_width = obj.line_width;
  let color_string = obj.color_string;

  // console.log(val);
  let val0 = {
    'x':this.data.range.x.min,
    'y':obj.val.y - (obj.val.x - this.data.range.x.min)*obj.slope
  };
  // console.log(val0);
  let val1 = {
    'x':this.data.range.x.max,
    'y':val0.y + this.data.range.x.span*obj.slope
  };
  
 let pixel0 = this.VALUE2PIXEL(val0);
 let pixel1 = this.VALUE2PIXEL(val1);

 this.ctx.lineWidth = line_width;
 
 this.ctx.strokeStyle = color_string;
 this.ctx.fillStyle = color_string;
 this.ctx.beginPath();
 this.ctx.moveTo(pixel0.x, pixel0.y);
 this.ctx.lineTo(pixel1.x, pixel1.y);
 this.ctx.stroke();

}



Box.prototype.CONNECTVALUES = function(val0, val1, color_string, line_width) {

 let pixel0 = this.VALUE2PIXEL(val0);
 let pixel1 = this.VALUE2PIXEL(val1);

 this.ctx.lineWidth = (line_width || 1);
 
 this.ctx.strokeStyle = color_string;
 this.ctx.fillStyle = color_string;
 this.ctx.beginPath();
 this.ctx.moveTo(pixel0.x, pixel0.y);
 this.ctx.lineTo(pixel1.x, pixel1.y);
 this.ctx.stroke();
}

/*

arr = [{},{}]

*/
Box.prototype.CONNECTVALUES2 = function(arr, color_string, line_width) {

  for (let i = 0; i < arr.length-1; i++) {

    let pixel0 = this.VALUE2PIXEL(arr[i+0]);
    let pixel1 = this.VALUE2PIXEL(arr[i+1]);

    this.ctx.lineWidth = (line_width || 1);
    this.ctx.strokeStyle = color_string;
    this.ctx.beginPath();
    this.ctx.moveTo(pixel0.x, pixel0.y);
    this.ctx.lineTo(pixel1.x, pixel1.y);
    this.ctx.stroke();

  }
}

Box.prototype.RECT_OUTLINE = function(val, w, h, color_string, line_width) {
 let pixel = this.VALUE2PIXEL(val);
 this.ctx.lineWidth = (line_width || 1);
 this.ctx.strokeStyle = (color_string || '#333');
 this.ctx.beginPath();
 this.ctx.rect(pixel.x, pixel.y, w*this.data.zoom.x, h*this.data.zoom.y);
 this.ctx.stroke();
}
Box.prototype.RECT_SOLID = function(val, w, h, color_string, line_width) {
 let pixel = this.VALUE2PIXEL(val);
 this.ctx.fillStyle = (color_string || '#fff');
 this.ctx.fillRect(pixel.x, pixel.y, w*this.data.zoom.x, h*this.data.zoom.y);
}
Box.prototype.TEXT = function(str, val, color_string, font_size, font_family) {
 let pixel = this.VALUE2PIXEL(val);
 //console.log(pixel);
 font_size = (font_size || 15);
 font_family = (font_family || 'Arial');
 this.ctx.font = font_size + 'px ' + font_family;
 this.ctx.strokeStyle = (color_string || '#fff');
 this.ctx.fillText(str, pixel.x, pixel.y);
}


function abc(x, arr) {

 if (arr[0] === 'get') {
   if (arr[1] === 'x') {
     return x;
   } else {
     return eval(arr[1]);        // why do ppl hate eval so much? why should i not use it?
   }
 }
 if (arr%arr===0) {
   return arr;
 }
 
 if (arr[0] === '**') {
   return abc(x, arr[1])**abc(x, arr[2]);
 }
 
 if (arr[0] === '/') {
   return abc(x, arr[1]) / abc(x, arr[2]);
 }
 if (arr[0] === '*') {
   return abc(x, arr[1]) * abc(x, arr[2]);
 }
 if (arr[0] === '-') {
   return abc(x, arr[1]) - abc(x, arr[2]);
 }
 
}
Box.prototype.DRAW_HISTOGRAM = function(obj) {
  
  // obj.data is what i need
  let n_bins = obj.k;
  let bin_width = this.data.range.x.span / n_bins;
  console.log(obj);
  for (let i = 0; i < n_bins; i++) {
    //console.log('this');
    this.RECT_OUTLINE({'x':i*bin_width,'y':0}, bin_width, -obj.data[i].count, '#aaa', 1); 
  }
  
}





/*
  let obj = {
    'val':val,
    'color_string':'#000',
    'rx':1
  }
*/
Box.prototype.DRAW_VALUE = function(obj) {
 let pixel = this.VALUE2PIXEL(obj.val);
 this.ctx.fillStyle = obj.colorstring;
 this.ctx.beginPath();
 this.ctx.arc(pixel.x, pixel.y, obj.rx, 0, 2*Math.PI);
 this.ctx.fill();
}





/*
  obj = {
    'alpha':[alpha0, alpha1],
    'budget':[budget0, budget1],
    'px':[px0, px1],
    'py':[py0, py1],
    'marshallian':true,
    'hicksian':true,
    'log':false,
    'color_string':'#000',
    'line_width':2
    'rx':4
  }
*/

Box.prototype.DRAW_DEMAND_CURVE = function(obj) {

  let line_width = 2;
  let rx = 4;
  
  let line_color = {
    'marshallian':{
      'initial':{
        'x':'#ffe6b3', // CASE 0
        'y':'#ffe6b3'  // CASE 1
      },
      'final':{
        'x':'#d1e0e0',  // CASE 2
        'y':'#d1e0e0'   // CASE 3
      }
    },
    'hicksian':{
      'initial':{
        'x':'#f937',  // CASE 4 : INITIAL HICKSIAN X
        'y':'#ffc266'   // CASE 5 : INITIAL HICKSIAN Y
      },
      'final':{
        'x':'#c2d1f0',  // CASE 6 : FINAL HICKSIAN X
        'y':'#c2d1f0'   // CASE 7
      }
    }
  }


  // CORE PARAMETERS

  let alpha_0 = (obj.alpha) ? (obj.alpha[0]) : (null); // BE CAREFUL WITH NULL, NULL**NULL = 1
  let alpha_1 = (obj.alpha) ? (obj.alpha[1]) : (null);
  
  let beta_0 = (obj.beta) ? (obj.beta[0]) : (null);
  let beta_1 = (obj.beta) ? (obj.beta[1]) : (null);
  
  let budget_0 = (obj.budget) ? (obj.budget[0]) : (null);
  let budget_1 = (obj.budget) ? (obj.budget[1]) : (null);

  let px_0 = (obj.px) ? (obj.px[0]) : (null);
  let px_1 = (obj.px) ? (obj.px[1]) : (null);
  
  let py_0 = (obj.py) ? (obj.py[0]) : (null);
  let py_1 = (obj.py) ? (obj.py[1]) : (null);


  // THE RESULTING ALLOCATIONS

  // MARSHALLIAN
  let x_0 = ((alpha_0 * budget_0 / px_0) || null); // BE CAREFUL WITH NULL, NULL**NULL = 1
  let y_0 = ((beta_0 * budget_0 / py_0) || null);
  let u_0 = (x_0 === null || y_0 === null) ? (null) : (x_0**alpha_0*y_0**beta_0);

  let x_1 = ((alpha_1 * budget_1 / px_1) || null);
  let y_1 = ((beta_1 * budget_1 / py_1) || null);
  let u_1 = (x_1 === null || y_1 === null) ? (null) : (x_1**alpha_1*y_1**beta_1);
  
  // HICKSIAN; COMPENSATED
  let xc_0 = (u_0*(py_0/px_0*alpha_0/beta_0)**beta_0);
  let yc_0 = (u_0*(px_0/py_0*beta_0/alpha_0)**alpha_0);
  
  let xc_1 = (!(u_0*(py_1/px_1*alpha_1/beta_1)) || !beta_1) ? (null) : (u_0*(py_1/px_1*alpha_1/beta_1)**beta_1);
  let yc_1 = (!(u_0*(px_1/py_1*beta_1/alpha_1)) || !alpha_1) ? (null) : (u_0*(px_1/py_1*beta_1/alpha_1)**alpha_1);

  
  // PREP FOR LOOP

  let dx_pixel = 2; // if dx_pixel = 2, then we calculate x-y every 2nd pixel
  let dx = dx_pixel * (this.data.range.x.span / this.data.dimension.w );
  
  let temp = new Array(8);
  temp[0] = {'x':0,'y':alpha_0 * budget_0 / (dx/1000)};   // INITIAL MARSHALLIAN X
  temp[1] = {'x':0,'y':beta_0 * budget_0 / (dx/1000)};    // INITIAL MARSHALLIAN Y
  
  temp[2] = {'x':0,'y':alpha_1 * budget_1 / (dx/1000)};   // INITIAL HICKSIAN X
  temp[3] = {'x':0,'y':beta_1 * budget_1 / (dx/1000)};    // INITIAL HICKSIAN Y
 
  temp[4] = {'x':0,'y':u_0*(py_0/(dx/1000)*alpha_0/beta_0)**beta_0}; // MARSHALLIAN RESPONSE X
  temp[5] = {'x':0,'y':u_0*(px_0/(dx/1000)*beta_0/alpha_0)**alpha_0}; // MARSHALLIAN RESPONSE Y
  
  temp[6] = {'x':0,'y':u_0*(py_1/(dx/1000)*alpha_1/beta_1)**beta_1}; // HICKSIAN RESPONSE X
  temp[7] = {'x':0,'y':u_0*(px_1/(dx/1000)*beta_1/alpha_1)**alpha_1}; // HICKSIAN RESPONSE Y

  let initial_marshallian_x = (obj.x[0] && obj.marshallian);
  let initial_marshallian_y = (obj.y[0] && obj.marshallian);
  let final_marshallian_x = (obj.x[1] && obj.marshallian);
  let final_marshallian_y = (obj.y[1] && obj.marshallian);
  let initial_hicksian_x = (obj.x[0] && obj.hicksian);
  let initial_hicksian_y = (obj.y[0] && obj.hicksian);
  let final_hicksian_x = (obj.x[1] && obj.hicksian);
  let final_hicksian_y = (obj.y[1] && obj.hicksian);


  // DRAW IN LOOP
  for (let x = dx; x < this.data.range.x.max+dx; x += dx) {

    // CASE 6 : FINAL HICKSIAN X
    if (final_hicksian_x) {
      
      /*
      let price = x;
      let qty = u_0*(py_1/x*alpha_1/beta_1)**beta_1;
      
      if (obj.log) {
        let ln_price = Math.log(price);
        let ln_qty = Math.log(qty);
      }
      */
      
      if (obj.log) {
        
        let lnpx_1 = Math.log(px_1);
        let lnxc_1 = Math.log(u_0) + beta_1 * Math.log(py_1*alpha_1/beta_1) - beta_1 * lnpx_1;
        
        this.DRAW_LINE({
          'val':{'x':lnpx_1,'y':lnxc_1},
          'slope':-beta_1,
          'color_string':line_color.hicksian.final.x,
          'line_width':line_width
        });
        
        this.DRAW_POINT({
          'val':{'x':lnpx_1,'y':lnxc_1},
          'color_string':line_color.hicksian.final.x,
          'rx':rx
        });
      
      } else {
        
        let price = x;
        let qty = u_0*(py_1/x*alpha_1/beta_1)**beta_1;
        let val = {'x':price,'y':qty};
        
        this.CONNECT_VALUES({
          'vals':[temp[6], val],
          'color_string':line_color.hicksian.final.x,
          'line_width':line_color
        });
        
        this.DRAW_POINT({
          'val':{'x':px_1,'y':xc_1},
          'color_string':line_color.hicksian.final.x,
          'rx':rx
        });
        
        temp[6] = val;
        
      }
    }
    
    // CASE 7 : FINAL HICKSIAN Y
    if (final_hicksian_y) {
      
      if (obj.log) {
        
        let lnpy_1 = Math.log(py_1);
        let lnyc_1 = Math.log(u_0) + alpha_1 * Math.log(px_1*beta_1/alpha_1) - alpha_1 * lnpy_1;
        
        this.DRAW_LINE({
          'val':{'x':lnpy_1,'y':lnyc_1},
          'slope':-alpha_1,
          'color_string':line_color.hicksian.final.y,
          'line_width':line_width
        });
        
        this.DRAW_POINT({
          'val':{'x':lnpy_1,'y':lnyc_1},
          'color_string':line_color.hicksian.final.y,
          'rx':rx
        });
      
      } else {
      
        let y = u_0*(px_1/x*beta_1/alpha_1)**alpha_1;
        let val = {'x':x,'y':y};
        this.CONNECTVALUES(temp[7], val, line_color.hicksian.final.y, line_width);
        temp[7] = val;
        
        let pixel = this.VALUE2PIXEL({'x':py_1,'y':yc_1});
        this.ctx.fillStyle = line_color.hicksian.final.y;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
        this.ctx.fill();
      
      }
    }


    // CASE 4
    if (initial_hicksian_x) {
      
      if (obj.log) {
        
        let lnpx_0 = Math.log(px_0);
        let lnxc_0 = Math.log(u_0) + beta_0 * Math.log(py_0*alpha_0/beta_0) - beta_0 * lnpx_0;
        
        let pixel = this.VALUE2PIXEL({'x':lnpx_0,'y':lnxc_0});
        this.ctx.fillStyle = line_color.hicksian.initial.x;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
        this.ctx.fill();
        
        this.DRAW_LINE({
          'val':{'x':lnpx_0,'y':lnxc_0},
          'slope':-beta_0,
          'color_string':line_color.hicksian.initial.x,
          'line_width':line_width
        });
        
      } else {
      
        let y = u_0*(py_0/x*alpha_0/beta_0)**beta_0;
        let val = {'x':x,'y':y};
        this.CONNECTVALUES(temp[4], val, line_color.hicksian.initial.x, line_width);
        temp[4] = val;
        
        let pixel = this.VALUE2PIXEL({'x':px_0,'y':u_0*(py_0/px_0*alpha_0/beta_0)**beta_0});
        this.ctx.fillStyle = line_color.hicksian.initial.x;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI); // THIS POINT WILL ACTUALLY NOT BE SEEN
        this.ctx.fill();
        
      }
    }

    // CASE 5
    if (initial_hicksian_y) {
      
      if (obj.log) {
        
        let lnpy_0 = Math.log(py_0);
        let lnyc_0 = Math.log(u_0) + alpha_0 * Math.log(px_0*beta_0/alpha_0) - alpha_0 * lnpy_0;
        
        this.DRAW_LINE({
          'val':{'x':lnpy_0,'y':lnyc_0},
          'slope':-alpha_0,
          'color_string':line_color.hicksian.initial.y,
          'line_width':line_width
        });
        
        this.DRAW_POINT({
          'val':{'x':lnpy_0,'y':lnyc_0},
          'color_string':line_color.hicksian.initial.y,
          'rx':rx
        });
      
      } else {
      
        let y = u_0*(px_0/x*beta_0/alpha_0)**alpha_0;
        let val = {'x':x,'y':y};
        this.CONNECTVALUES(temp[5], val, line_color.hicksian.initial.y, line_width);
        temp[5] = val;
        
        let pixel = this.VALUE2PIXEL({'x':py_0,'y':u_0*(px_0/py_0*beta_0/alpha_0)**alpha_0});
        this.ctx.fillStyle = line_color.hicksian.initial.y;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI); // THIS POINT WILL ACTUALLY NOT BE SEEN
        this.ctx.fill();
      
      }
    }

    // CASE 3 : FINAL MARSHALLIAN Y
    if (final_marshallian_y) {
      
      if (obj.log) {
        
        let lnpy_1 = Math.log(py_1);
        let lny_1 = Math.log(beta_1*budget_1) - lnpy_1;
        
        this.DRAW_LINE({
          'val':{'x':lnpy_1,'y':lny_1},
          'slope':-1,
          'color_string':line_color.marshallian.final.y,
          'line_width':line_width
        });
        
        this.DRAW_POINT({
          'val':{'x':lnpy_1,'y':lny_1},
          'color_string':line_color.marshallian.final.y,
          'rx':rx
        });
      
      } else {
      
        let y = beta_1 * budget_1 / x;
        let val = {'x':x,'y':y};
        this.CONNECTVALUES(temp[3], val, line_color.marshallian.final.y, line_width);
        temp[3] = val;
        
        let pixel = this.VALUE2PIXEL({'x':py_1,'y':(beta_1 * budget_1 / py_1)});
        this.ctx.fillStyle = line_color.marshallian.final.y;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
        this.ctx.fill();
      
      }
    }

    // CASE 2 : FINAL MARSHALLIAN X
    if (final_marshallian_x) {
      
      if (obj.log) {
        
        // need x, but we have px, alpha, M
        
        let lnpx_1 = Math.log(px_1);
        let lnx_1 = Math.log(alpha_1*budget_1) - lnpx_1;
        
        let pixel = this.VALUE2PIXEL({'x':lnpx_1,'y':lnx_1});
        this.ctx.fillStyle = line_color.marshallian.final.x;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
        this.ctx.fill();
        
        this.DRAW_LINE({
          'val':{'x':lnpx_1,'y':lnx_1},
          'slope':-1,
          'color_string':line_color.marshallian.final.x,
          'line_width':line_width
        });
        
      } else {
      
        let y = alpha_1 * budget_1 / x;
        let val = {'x':x,'y':y};
        this.CONNECTVALUES(temp[2], val, line_color.marshallian.final.x, line_width);
        temp[2] = val;
        
        let pixel = this.VALUE2PIXEL({'x':px_1,'y':(alpha_1 * budget_1 / px_1)});
        this.ctx.fillStyle = line_color.marshallian.final.x;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
        this.ctx.fill();
      
      }
    }
    
    // CASE 1 : INITIAL MARSHALLIAN Y
    if (initial_marshallian_y) {
      
      if (obj.log) {
        
        let lnpy_0 = Math.log(py_0);
        let lny_0 = Math.log(beta_0*budget_0) - lnpy_0;
        
        this.DRAW_LINE({
          'val':{'x':lnpy_0,'y':lny_0},
          'slope':-1,
          'color_string':line_color.marshallian.initial.y,
          'line_width':line_width
        });
        
        this.DRAW_POINT({
          'val':{'x':lnpy_0,'y':lny_0},
          'color_string':line_color.marshallian.initial.y,
          'rx':rx
        });
      
      } else {
      
        let y = beta_0 * budget_0 / x;
        let val = {'x':x,'y':y};
        this.CONNECTVALUES(temp[1], val, line_color.marshallian.initial.y, line_width);
        temp[1] = val;
        
        let pixel = this.VALUE2PIXEL({'x':py_0,'y':(beta_0 * budget_0 / py_0)});
        this.ctx.fillStyle = line_color.marshallian.initial.y;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
        this.ctx.fill();
      
      }
    }
    
    // CASE 0 : INITIAL MARSHALLIAN X
    if (initial_marshallian_x) {
      
      if (obj.log) {
        
        // need x, but we have px, alpha, M
        
        let lnpx_0 = Math.log(px_0);
        let lnx_0 = Math.log(alpha_0*budget_0) - lnpx_0;
        
        let pixel = this.VALUE2PIXEL({'x':lnpx_0,'y':lnx_0});
        this.ctx.fillStyle = line_color.marshallian.initial.x;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
        this.ctx.fill();

        this.DRAW_LINE({
          'val':{'x':lnpx_0,'y':lnx_0},
          'slope':-1,
          'color_string':line_color.marshallian.initial.x,
          'line_width':line_width
        });

      } else {
      
        let y = alpha_0 * budget_0 / x;
        let val = {'x':x,'y':y};
        this.CONNECTVALUES(temp[0], val, line_color.marshallian.initial.x, line_width);
        temp[0] = val;
        
        let pixel = this.VALUE2PIXEL({'x':px_0,'y':(alpha_0 * budget_0 / px_0)});
        this.ctx.fillStyle = line_color.marshallian.initial.x;
        this.ctx.beginPath();
        this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
        this.ctx.fill();
      
      }
    }

  } // closing the loop

  return {
    'initial':{
      'marshallian':{'x':x_0,'y':y_0,'u':u_0},
      'hicksian':{'x':xc_0,'y':yc_0,'u':u_0}
    },
    'final':{
      'marshallian':{'x':x_1,'y':y_1,'u':u_1},
      'hicksian':{'x':xc_1,'y':yc_1,'u':u_0}
    }
  };
  
}




/*
  let obj = {
    'budget':budget,
    'px':px,
    'py':py
  }
*/
Box.prototype.DRAW_BUDGET_LINE = function(obj) {

  let budget = obj.budget;
  let px = obj.px;
  let py = obj.py;

  let line_width = 2;
  let color_string = '#000';

  if (obj.color_string) {color_string = obj.color_string};
  if (obj.line_width) {line_width = obj.line_width};
  
  let x_int = {'x':budget/px,'y':0};
  let y_int = {'x':0,'y':budget/py};

  let pixel0 = this.VALUE2PIXEL(x_int);
  let pixel1 = this.VALUE2PIXEL(y_int);

  this.ctx.lineWidth = line_width;
 
  this.ctx.strokeStyle = color_string;
  this.ctx.fillStyle = color_string;
  this.ctx.beginPath();
  this.ctx.moveTo(pixel0.x, pixel0.y);
  this.ctx.lineTo(pixel1.x, pixel1.y);
  this.ctx.stroke();
};

/*
let obj = {
  'points':[
    {'x':0,'y':0},
    {'x':0,'y':0},
    {'x':0,'y':0}
  ]
}
*/
Box.prototype.RESCALE_BASED_ON_CM = function(obj) {

  console.log(obj);

  let cx = 0;
  let cy = 0;
  
  let x_min = +(9**5);
  let x_max = -(9**5);
  let y_min = +(9**5);
  let y_max = -(9**5);

  for (let i = 0; i < obj.points.length; i++) {
    cx += obj.points[i].x;
    cy += obj.points[i].y; 
    
    if (obj.points[i].x > x_max) {
      x_max = obj.points[i].x;
    }
    if (obj.points[i].x < x_min) {
      x_min = obj.points[i].x;
    }
    if (obj.points[i].y > y_max) {
      y_max = obj.points[i].y;
    }
    if (obj.points[i].y < y_min) {
      y_min = obj.points[i].y;
    }

  }

  let dx = x_max - x_min;
  let dy = y_max - y_min;

  this.RANGE_X(0, x_max*5/4);
  this.RANGE_Y(0, y_max*5/4);

  return {
    'cx':cx,
    'cy':cy,
    'dx':dx,
    'dy':dy,
    'x_min':x_min,
    'x_max':x_max,
    'y_min':y_min,
    'y_max':y_max
  }

}

/*
let obj = {
  'points':[
    {'x':0,'y':0},
    {'x':0,'y':0},
    {'x':0,'y':0}
  ]
}
*/
Box.prototype.RESCALE_BASED_ON_CENTROID = function(obj) {

  console.log(obj);
  let a = obj.points[0];
  let b = obj.points[1];
  let c = obj.points[2];

  let cx = (a.x + b.x + c.x)/3;
  let cy = (a.y + b.y + c.y)/3;
  
  let ra = ((a.x - cx)**2 + (a.y - cy)**2)**0.5;
  let rb = ((b.x - cx)**2 + (b.y - cy)**2)**0.5;
  let rc = ((c.x - cx)**2 + (c.y - cy)**2)**0.5;

  let r = 0;
  if (ra > rb) {
    r = ra;
  } else {
    r = rb;
  }
  if (rc > r) {
    r = rc;
  }
  
  if (a.x == b.x && b.x == c.x &&  a.y == b.y && b.y == c.y) {
    r = 0.1;
  }
  
  let z = r*3;
  if (z < 0.5) {
    z = 0.5;
  }

 // RESCALE
 this.rangex(cx-z, cx+z);
 this.rangey(cy-z, cy+z);
  
}