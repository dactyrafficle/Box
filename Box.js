
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
  console.log(this);
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



Box.prototype.clear = function() {
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

Box.prototype.SHOW_GRID_X = function(dx) {

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

Box.prototype.SHOW_FLOATING_LOG_Y_AXIS = function() {

 let n = 9;
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
  this.ctx.fillText((Math.E**v0.y).toFixed(2), p0.x+2*dsw, p0.y);
  this.ctx.stroke();
  this.CONNECTVALUES(v0, v1, '#333', 0.5); 
 }
 
 this.ctx.textAlign = 'center';
 this.ctx.textBaseline = 'top';
 this.ctx.fillText((this.data.label.y).toUpperCase(), sw, sh*(n2+0.25));
};

Box.prototype.SHOW_FLOATING_LOG_X_AXIS = function() {

 let n = 9;
 let sh = this.data.dimension.h - this.data.dimension.h/n;
 let sw = this.data.dimension.w/n;
 
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
  this.ctx.fillText((Math.E**v0.x).toFixed(2), p0.x, p0.y+1*dsh);
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
 this.ctx.beginPath();
 this.ctx.moveTo(pixel0.x, pixel0.y);
 this.ctx.lineTo(pixel1.x, pixel1.y);
 this.ctx.stroke();
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
