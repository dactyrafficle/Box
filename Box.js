
// let box = new Box();
// let c = box.returnCanvas();
// container.appendChild(c);  // append it to a container div element or sthing

// box.dimension(500, 500); // this is the number of pixels

// you may need to reset these on some other condition
// box.rangex(5, 50); // set the range in x
// box.rangey(3, 130); // set the range in y 


// box.SHOWGRIDY(5);

// box.SHOWVALUE({'x':5,'y':5}, '#fco', 2);



 
function Box(label_x, label_y) {
 this.c = document.createElement('canvas');
 this.ctx = this.c.getContext('2d');
 
 this.data = {
  'dimension':{'w':100,'h':100},
  'zoom':{'x':1,'y':1},
  'translate':{'x':0,'y':0},
  'range':{
    'x':{
      'min':0,
      'max':100
    },
    'y':{
      'min':0,
      'max':100
    }
  },
  'label':{'x':label_x,'y':label_y}
 };
 
 this.maps = [
  //['**','get',0.5],
  //['/',23.5,['**','get',0.5]]
 ]
}

Box.prototype.VAL2PIXEL = function(val) {  // val : (0,0) is bottom-left (Cartesian)
 return {
  'x':(val.x+this.data.translate.x)*this.data.zoom.x,
  'y':this.data.dimension.h - (val.y+this.data.translate.y)*this.data.zoom.y
 }
}

Box.prototype.PIXEL2VAL = function(pixel) { // pixel : (0,0) is top-left (standard computer/matrix grid)
 return {
  'x':(pixel.x/this.data.zoom.x)-this.data.translate.x,
  'y':(this.data.dimension.h-pixel.y)/this.data.zoom.y-this.data.translate.y
 }
}

// centerOnValue(val) // takes the value, and places it in the center of the canvas, scaled so that it goes from (0,0) to (cx, cy)

Box.prototype.recenter = function(val) {
  let spanx = this.data.range.x.span;
  let spany = this.data.range.y.span;
  this.rangex(val.x-spanx/2, val.x+spanx/2);
  this.rangey(val.y-spany/2, val.y+spany/2);
  //console.log(this);
  // BUT I NEED TO RECALC EVERYTHING
  //this.clear();
  //this.SHOW_GRID_X(1);
  //this.draw();
}

Box.prototype.returnCanvas = function() {
 this.c.addEventListener('click', function(e) {
  let pixel = {'x':e.offsetX,'y':e.offsetY};
  let val = this.PIXEL2VAL(pixel);
  //this.recenter(val);
 }.bind(this)); // makes this.c inherit 'this' from Box
 return this.c;
}




Box.prototype.border = function(x) {
 this.c.style.border = x;
}

Box.prototype.dimension = function(w, h) {
 this.data.dimension.w = w;
 this.data.dimension.h = h; 
 this.c.width = this.data.dimension.w;
 this.c.height = this.data.dimension.h;
}

Box.prototype.zoom = function(zx, zy) {
 this.data.zoom.x = zx;
 this.data.zoom.y = zy;
}

Box.prototype.translate = function(x, y) {
 this.data.translate.x = x;
 this.data.translate.y = y; 
}

Box.prototype.rangex = function(min, max) {
 this.data.range.x.min = min;
 this.data.range.x.max = max;
 this.data.range.x.span = max-min;

 this.data.zoom.x = this.data.dimension.w / this.data.range.x.span;
 this.data.translate.x = -this.data.range.x.min;
}

Box.prototype.rangey = function(min, max) {
 this.data.range.y.min = min;
 this.data.range.y.max = max;
 this.data.range.y.span = max-min;

 this.data.zoom.y = this.data.dimension.h / this.data.range.y.span;
 this.data.translate.y = -this.data.range.y.min;
}



Box.prototype.CLEAR_CANVAS = function() {
 this.ctx.fillStyle = '#fff';
 this.ctx.beginPath();
 this.ctx.rect(0, 0, this.data.dimension.w, this.data.dimension.h);
 this.ctx.fill();
}

Box.prototype.SHOWGRIDX = function(dx) {

 if (!arguments[0]) {
   dx = 50* this.data.range.x.span / this.data.dimension.w; // default is 25 px
 }

 for (let x = -dx; x > this.data.range.x.min; x -= dx) {
   let val0 = {'x':x,'y':this.data.range.y.min};
   let val1 = {'x':x,'y':this.data.range.y.max};
   this.CONNECTVALUES(val0, val1, '#ddd');
 }

 for (let x = 0; x < this.data.range.x.max; x += dx) {
   let val0 = {'x':x,'y':this.data.range.y.min};
   let val1 = {'x':x,'y':this.data.range.y.max};
   this.CONNECTVALUES(val0, val1, '#ddd');
 }


};

Box.prototype.SHOWGRIDY = function(dy) {

 if (!arguments[0]) {
   dy = 50* this.data.range.y.span / this.data.dimension.h; // default is 25 px
 }

 for (let y = -dy; y > this.data.range.y.min; y -= dy) {
   let val0 = {'x':this.data.range.x.min,'y':y};
   let val1 = {'x':this.data.range.x.max,'y':y};
   this.CONNECTVALUES(val0, val1, '#ddd');
 }

 for (let y = 0; y < this.data.range.y.max; y += dy) {
   let val0 = {'x':this.data.range.x.min,'y':y};
   let val1 = {'x':this.data.range.x.max,'y':y};
   this.CONNECTVALUES(val0, val1, '#ddd');
 }


};

Box.prototype.SHOW_GRID_X = function(dx, color_string) {

 if (!arguments[0]) {
  dx = 1;
 }
 let zoom = this.data.zoom.x;

 let alpha_1 = 1;
 let alpha_10 = -1+(1/200)*zoom;
 let color_string_1 = 'rgba(208, 208, 208,' + alpha_1 + ')';
 let color_string_10 = 'rgba(224, 224, 224,' + alpha_10 + ')';
 
 let x_start = Math.floor(this.data.range.x.min/dx)*dx;
 
 let i = 0; // BC FLOATING PT NUMBERS MAKE TESTING x%dx TOUGH
 for (let x = x_start; x < this.data.range.x.max; x += dx/10) {
  if (i%10===0) {
   let val0 = {'x':x,'y':this.data.range.y.min};
   let val1 = {'x':x,'y':this.data.range.y.max};
   this.CONNECTVALUES(val0, val1, color_string_1);
  }

  let val0 = {'x':x,'y':this.data.range.y.min};
  let val1 = {'x':x,'y':this.data.range.y.max};
  this.CONNECTVALUES(val0, val1, color_string_10);
  i++;
 }
};
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
Box.prototype.SHOW_FLOATING_LOG_Y_AXIS = function(n) {

 let sh = this.data.dimension.h/n;
 let sw = this.data.dimension.w/n;
 
 let n1 = 2;
 let n2 = n-n1;
 
 let p0 = {'x':sw,'y':sh*n1};
 let p1 = {'x':sw,'y':sh*n2};

 let v0 = this.PIXEL2VAL(p0);
 let v1 = this.PIXEL2VAL(p1);
 
 this.CONNECTVALUES(v0, v1, '#333', 0.5);
 
 let dsw = 5;
 for (let i = n1; i <= n2; i++) {
  let p0 = {'x':sw-dsw,'y':sh*i};
  let p1 = {'x':sw+dsw,'y':sh*i};
  let v0 = this.PIXEL2VAL(p0);
  let v1 = this.PIXEL2VAL(p1);
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

 let v0 = this.PIXEL2VAL(p0);
 let v1 = this.PIXEL2VAL(p1);
 
 this.CONNECTVALUES(v0, v1, '#333', 0.5);
 
 let dsw = 5;
 for (let i = n1; i <= n2; i++) {
  let p0 = {'x':sw-dsw,'y':sh*i};
  let p1 = {'x':sw+dsw,'y':sh*i};
  let v0 = this.PIXEL2VAL(p0);
  let v1 = this.PIXEL2VAL(p1);
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
Box.prototype.SHOW_FLOATING_LOG_X_AXIS = function(n, y_val) {

  let sh;
  if (arguments[1] && arguments[1] !== null) {
    let v3 = {'x':0,'y':y_val};
    let p3 = this.VAL2PIXEL(v3);
    sh = p3.y;
  } else {
    sh = this.data.dimension.h*(n-1)/n;
  }
  let sw = this.data.dimension.w/n;
  //console.log(sh);
 
 let n1 = 2;
 let n2 = n-n1;
 
 let p0 = {'x':sw*n1,'y':sh};
 let p1 = {'x':sw*n2,'y':sh};

 let v0 = this.PIXEL2VAL(p0);
 let v1 = this.PIXEL2VAL(p1);
 
 this.CONNECTVALUES(v0, v1, '#333', 0.5);

 let dsh = 5;
 for (let i = n1; i <= n2; i++) {
  let p0 = {'x':sw*i,'y':sh+dsh};
  let p1 = {'x':sw*i,'y':sh-dsh};
  let v0 = this.PIXEL2VAL(p0);
  let v1 = this.PIXEL2VAL(p1);
  this.ctx.fillStyle = '#333';
  this.ctx.textAlign = 'center';
  this.ctx.textBaseline = 'top';
  this.ctx.fillText((v0.x).toFixed(2), p0.x, p0.y+1*dsh);
  this.ctx.stroke();
  this.CONNECTVALUES(v0, v1, '#333', 0.5); 
 }

 this.ctx.textAlign = 'right';
 this.ctx.textBaseline = 'middle';
 this.ctx.fillText((this.data.label.x), sw*(n1-0.25), sh);

};

Box.prototype.SHOW_FLOATING_X_AXIS = function(n, y_val) {

 //let n = 9;
 let sh;
 if (arguments[1] !== null) {
   let v3 = {'x':0,'y':y_val};
   //console.log(v3);
   let p3 = this.VAL2PIXEL(v3);
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

 let v0 = this.PIXEL2VAL(p0);
 let v1 = this.PIXEL2VAL(p1);
 
 this.CONNECTVALUES(v0, v1, '#333', 0.5);

 let dsh = 5;
 for (let i = n1; i <= n2; i++) {
  let p0 = {'x':sw*i,'y':sh+dsh};
  let p1 = {'x':sw*i,'y':sh-dsh};
  let v0 = this.PIXEL2VAL(p0);
  let v1 = this.PIXEL2VAL(p1);
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

   let o = this.PIXEL2VAL(c);
   console.log(o);
   console.log(this.data.range);
   
   let cx = {'x':s0+s1, 'y':this.data.dimension.h - s0};
   let ox = this.PIXEL2VAL(cx);
   
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

 let pixel = this.VAL2PIXEL(val);

 this.ctx.fillStyle = colorstring;
 this.ctx.beginPath();
 this.ctx.arc(pixel.x, pixel.y, rx, 0, 2*Math.PI);
 this.ctx.fill();

}

Box.prototype.CONNECTVALUES = function(val0, val1, color_string, line_width) {

 let pixel0 = this.VAL2PIXEL(val0);
 let pixel1 = this.VAL2PIXEL(val1);

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

    let pixel0 = this.VAL2PIXEL(arr[i+0]);
    let pixel1 = this.VAL2PIXEL(arr[i+1]);

    this.ctx.lineWidth = (line_width || 1);
    this.ctx.strokeStyle = color_string;
    this.ctx.beginPath();
    this.ctx.moveTo(pixel0.x, pixel0.y);
    this.ctx.lineTo(pixel1.x, pixel1.y);
    this.ctx.stroke();

  }
}

Box.prototype.RECT_OUTLINE = function(val, w, h, color_string, line_width) {
 let pixel = this.VAL2PIXEL(val);
 this.ctx.lineWidth = (line_width || 1);
 this.ctx.strokeStyle = (color_string || '#333');
 this.ctx.beginPath();
 this.ctx.rect(pixel.x, pixel.y, w*this.data.zoom.x, h*this.data.zoom.y);
 this.ctx.stroke();
}
Box.prototype.RECT_SOLID = function(val, w, h, color_string, line_width) {
 let pixel = this.VAL2PIXEL(val);
 this.ctx.fillStyle = (color_string || '#fff');
 this.ctx.fillRect(pixel.x, pixel.y, w*this.data.zoom.x, h*this.data.zoom.y);
}
Box.prototype.TEXT = function(str, val, color_string, font_size, font_family) {
 let pixel = this.VAL2PIXEL(val);
 //console.log(pixel);
 font_size = (font_size || 15);
 font_family = (font_family || 'Arial');
 this.ctx.font = font_size + 'px ' + font_family;
 this.ctx.strokeStyle = (color_string || '#fff');
 this.ctx.fillText(str, pixel.x, pixel.y);
}

Box.prototype.draw = function() {
  
  this.maps.forEach(function(a, b, c) {

    let arr = a;
    //console.log(arr);
      /*
   // what is the operation?
   let operation = (function() {
     if (a[0] === 'get') {
       return function(a) {
         return a;
       }
     }
     if (a[0] === '**') {
       return function(a,b) {
         return a**b;
       }
     }
    if (a[0] === '+') {  /// its not arguments, i need arguments of 
      return function(a,b) {
        return a + b;
      }
    };
    if (a[0] === '*') {
      return function(a,b) {
        return a*b;
      }
    };
   })();
   console.log(operation);
   */
   
   for (let i = this.data.range.x.min; i < this.data.range.x.max; i++) {
     let val = {'x':i,'y':abc(i,arr)};
     this.SHOWVALUE(val, '#fc0', 2);
   }
  }.bind(this));

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


/*
 obj = {
  'alpha':alpha,
  'beta':beta,
  'm',m,
  'x':x,
  'y':y  
 }
*/
Box.prototype.DRAW_ISOQUANT = function(obj, color_string, line_width, inverted) {
 // if you have x, y, alpha, beta : you have utility to start with
 // equate the marshallian demands alpha*M/px with beta*M/py on M
 // y = x*px*beta/alpha
 // if you have x, y use U = x**alpha*y**beta
 
 // this.CONNECTVALUES(v0, v1, '#333', 0.5);
 
  let a = obj.x**obj.alpha;
  let b = obj.y**obj.beta;
  let u = a*b;
  let beta_inv = 1/obj.beta;

  let dx = 0.5;
  let dy = 0.25;
  let temp = {'x':this.data.range.x.min, 'y':(u/this.data.range.x.min**obj.alpha)**beta_inv};
  
  for (x = this.data.range.x.min+dx; x < this.data.range.x.max+dx; x += dx) {

    let y = (u/x**obj.alpha)**beta_inv;
    let val;
    
    if (Math.abs(y - temp.y) > dy) {
      
      if (inverted) {
        val = {'x':x, 'y':y};
        this.CONNECTVALUES({'x':this.data.range.x.max-temp.x, 'y':this.data.range.y.max-temp.y}, {'x':this.data.range.x.max-x, 'y':this.data.range.y.max-y}, color_string, line_width);
        temp = val;
      } else {
        val = {'x':x, 'y': y};
        this.CONNECTVALUES(temp, val, color_string, line_width);
        temp = val;
      }
    } 
  }



/*
 for (x = this.data.range.x.min; x < this.data.range.x.max+dx; x += dx) {



    // add all the points
    let x0 = x;
    let x1 = x+dx;
    let y0 = (u/x0**obj.alpha)**beta_inv;
    let y1 = (u/x1**obj.alpha)**beta_inv;
    let val0 = {'x':x0, 'y': y0};
    let val1 = {'x':x1, 'y': y1};
  
   if ((y0 >= this.data.range.y.min && y0 <= this.data.range.y.max) || (y1 >= this.data.range.y.min && y1 <= this.data.range.y.max)) {
   
    // AT LEAST ONE POINT WORKS
    // let p = [val0, val1];
    // console.log(p);
    if (inverted) {
      val0 = {'x':(this.data.range.x.max-x0), 'y':(this.data.range.y.max-y0)};
      val1 = {'x':(this.data.range.x.max-x1), 'y':(this.data.range.y.max-y1)};
    }
    
    this.CONNECTVALUES(val0, val1, color_string, line_width);
   }
 }
*/
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

Box.prototype.RESCALE_BASED_ON_CM = function(alpha, x_scale_min, y_scale_min, VERTICES) {

  let cx = alpha*VERTICES[0].x + (1-alpha)*VERTICES[1].x;
  let cy = alpha*VERTICES[0].y + (1-alpha)*VERTICES[1].y;
  
  let dx = Math.abs(VERTICES[0].x - VERTICES[1].x);
  let dy = Math.abs(VERTICES[0].y - VERTICES[1].y);
  
 // if (VERTICES[0].x == VERTICES[1].x && VERTICES[0].y == VERTICES[1].y) {
    this.rangex(cx-x_scale_min/2, cx+x_scale_min/2);
    this.rangey(cy-y_scale_min/2, cy+y_scale_min/2);
    //return;
  //}
  
  if (dx > 0.25*this.data.range.x.span) {
    this.rangex(cx-dx*2, cx+dx*2);
  }
  if (dy > 0.25*this.data.range.y.span) {
    this.rangey(cy-dy*2, cy+dy*2);
  }


  
}