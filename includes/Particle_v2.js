
// GET_DISTANCE(other);
// CONSTRAIN_X(a=-999, b=999);
// CONSTRAIN_Y(a=-999, b=999);
// UPDATE_PARTICLE_POSITION();

/*

// tethered is embedded in update particle position
// particle1.pivot_point = particle2
// particle1.pivot_length = L
// if it is tethered, it works off of torque, alpha, omega, theta
UPDATE_TETHERED_POSITION(dt_=1)
UPDATE_PARTICLE_POSITION(dt_=1)

in each case, we must set all forces particle.F.x = 0, and particle.F.y = 0

  // ADD_REPULSION_BASED_ON_PARTICLE_MASS(other_, distance_ = "r2", multiplier_=1)
  // ADD_ATTRACTION_BASED_ON_PARTICLE_MASS(other_, distance_="r2", multiplier_=1)

*/



class Particle_system {
  
  // particles. tethered.
  // vector to theta (function i added here)
  // distance matrix. a matrix of the distance of every particle relative to every other particle
  
  constructor(obj) {
    
    this.x_min = (obj.x_min || 0);
    this.x_max = (obj.x_max || 100);
    this.y_min = (obj.y_min || 0);
    this.y_max = (obj.y_max || 100);
    
    this.gravitational_constant = 0;
    
    
    
    this.particle_max_speed = (obj.particle_max_speed || 5);
    
    this.particle_arr = [];
  };  // CLOSING CONSTRUCTOR
  
  GENERATE_PARTICLES(n_=10, r_ = 1, v_=1) {
    for (let i = 0; i < n_; i++) {
      let x_ = this.x_min+Math.random()*(this.x_max-this.x_min);
      let y_ = this.y_min+Math.random()*(this.y_max-this.y_min);
      let vx_ = -v_ + Math.random()*2*v_;
      let vy_ = -v_ + Math.random()*2*v_;
      this.particle_arr.push(new Particle({'x':x_,'y':y_,'vx':vx_,'vy':vy_,'particle_radius':r_})); 
    }
  };
  
}

class Particle {

  // obj = {x, y}
  constructor(obj) {
    
    this.particle_name = (obj.particle_name || 'Particle');
    this.particle_charge = (obj.particle_charge || 0);
    this.particle_mass = (obj.particle_mass || 1);
    this.particle_radius = (obj.particle_radius || 1);
    this.particle_radius_of_neutrality = (obj.particle_radius_of_neutrality || 0); 

    // IS THE PARTICLE TETHERED TO A POINT ? (IE. ATTACHED TO A POINT BY A RIGID STRING)
    this.IS_TETHERED = false;
    
    this.pivot_point = null;
    this.pivot_length = null;
    
    //this.pivot.x = null;
    //this.pivot.y = null;

    this.stroke_style = (obj.stroke_style || '#f000');
    this.fill_style = (obj.fill_style || '#f000');

    // POSITION VECTOR
    this.p = {};
    this.p.x = (obj.x || 0);
    this.p.y = (obj.y || 0);
    
    // VELOCITY VECTOR
    this.v = {};
    this.v.x = (obj.vx || 0);
    this.v.y = (obj.vy || 0);
    this.max_particle_speed = (obj.max_particle_speed || 9999);  // MAX VELOCITY
    this.minimum_damping_speed = (obj.minimum_damping_speed || 0);
    this.damping_factor = (obj.damping_factor || 0.9);

    this.a = {
      'x':0,
      'y':0
    };
    
    
    // if tethered
    this.L = 1;
    this.theta = 0;
    this.omega = 0;
    this.alpha = 0;
    
    this.F = {
      'x':0,
      'y':0
    };

  };  // CLOSING CONSTRUCTOR
  
  
  
  // UPDATE ACCUMULATED FORCES
  // 
  ADD_ATTRACTION_BASED_ON_PARTICLE_MASS(other_, distance_="r2", multiplier_=1) {
    
    let dx_ = this.p.x-other_.p.x;
    let dy_ = this.p.y-other_.p.y;
    
    let rho_ = ((dx_)**2+(dy_)**2)**(0.5); 
    
    // displacement unit vector    
    let nx_ = dx_/rho_;
    let ny_ = dy_/rho_;
    
    let R2_inv = 1/(rho_**2);
    let m1 = this.particle_mass;
    let m2 = other_.particle_mass;
    
    // F = m1*m2/(r**2)    
    let F = multiplier_*m1*m2*R2_inv;
       
    if (distance_ == "r2") {
      F = multiplier_*m1*m2*R2_inv;
    }
    if (distance_ == "r") {
      F = multiplier_*m1*m2/rho_;
    }
       
    let Fx = nx_*F;
    let Fy = ny_*F;
        
    this.F.x += -Fx;
    this.F.y += -Fy;
    other_.F.x += Fx;
    other_.F.y += Fy;
        
  }; // closing


  // UPDATE ACCUMULATED FORCES
  // 
  ADD_REPULSION_BASED_ON_PARTICLE_MASS(other_, distance_ = "r2", multiplier_=1) {
    
    let dx_ = this.p.x-other_.p.x;
    let dy_ = this.p.y-other_.p.y;
    
    let rho_ = ((dx_)**2+(dy_)**2)**(0.5); 
    
    // displacement unit vector    
    let nx_ = dx_/rho_;
    let ny_ = dy_/rho_;
    
    let R2_inv = 1/(rho_**2);
    let m1 = this.particle_mass;
    let m2 = other_.particle_mass;
    
    // F = m1*m2/(r**2)    
    let F = multiplier_*m1*m2*R2_inv;
    
    if (distance_ == "r2") {
      F = multiplier_*m1*m2*R2_inv;
    }
    if (distance_ == "r") {
      F = multiplier_*m1*m2/rho_;
    }    
    
    
    let Fx = nx_*F;
    let Fy = ny_*F;
        
    this.F.x += Fx;
    this.F.y += Fy;
    other_.F.x += -Fx;
    other_.F.y += -Fy;
        
  }; // closing
  
  ADD_PARTICLE_TO_PARTICLE_FORCE(obj_) {
    
    // default
    /*
     line of action along displacement vector : p1 to p2 (position vectors)
     proportional to inverse square of displacement
     repel
     
     obj_.other
     obj_.based_on = "particle_mass"
     obj_.multiplier = 1
    */
    
    
    /*
    obj. how to distinguish particle position ? as opposed to velocity
    obj.other = particle;
    obj.based_on = "particle_mass" or "particle_charge"
    obj.mutual = false;
    obj.range = "r2"; (distance, distance_sq) etc
    obj.force_multiplier = 1;
    obj.direction = (attract || repel || other)
    other_, distance_ = "r2", multiplier_=1
    
    */
    
    // DISPLACEMENT VECTOR 
    let dx_ = this.p.x-obj_.other_.p.x;
    let dy_ = this.p.y-obj_.other_.p.y;
    let rho_ = ((dx_)**2+(dy_)**2)**(0.5); 
    
    // DISPLACEMENT UNIT VECTOR    
    let nx_ = dx_/rho_;
    let ny_ = dy_/rho_;
    
    let R2_inv = 1/(rho_**2);
    
    let m1 = 0;
    let m2 = 0;
    if (obj_.based_on == "particle_mass") {
      m1 = this.particle_mass;
      m2 = obj_.other.particle_mass;
    }
    
    
    
    
    // F = m1*m2/(r**2)    
    let F = multiplier_*m1*m2*R2_inv;
    
    if (distance_ == "r2") {
      F = multiplier_*m1*m2*R2_inv;
    }
    if (distance_ == "r") {
      F = multiplier_*m1*m2/rho_;
    }    
    
    
    let Fx = nx_*F;
    let Fy = ny_*F;
    
    
    
    
    this.F.x += Fx;
    this.F.y += Fy;
    
    if (obj_.mutual) {
      other_.F.x += -Fx;
      other_.F.y += -Fy;
    }
        
  }; // closing
  
  UPDATE_TETHERED_POSITION(dt_=1) {

    // IF THERE IS NO PIVOT POINT nor PIVOT LENGTH, THEN RETURN
    if (!this.pivot_point || !this.pivot_length) {
      return;
    }
    // p2 is anchored to p1 at a distance of L
    
    // let vector R be the displacement from the pivot to the particle
    let dx_ = this.p.x - this.pivot_point.p.x;
    let dy_ = this.p.y - this.pivot_point.p.y;
    
    // INITIAL MAGNITUDE OF R
    let rho_ = (dx_**2+dy_**2)**(0.5);
    
    // UNIT VECTOR OF R
    let nx_ = dx_/rho_;
    let ny_ = dy_/rho_;
  
    // REPOSITION particle.p.x and particle.p.y to reflect that MAG(R) needs to be equal to L

    // reposition p2 along the new vector A, at a distance of r
    dx_ = nx_*this.pivot_length;
    dy_ = ny_*this.pivot_length;
    
    // RECALCULAE THE MAGNITUDE OF R
    rho_ = (dx_**2+dy_**2)**(0.5);
    
    let R = {
      "x":dx_,
      "y":dy_,
      "rho":rho_,
      "nx":nx_,
      "ny":ny_
    }
    
    this.p.x = this.pivot_point.p.x + dx_; 
    this.p.y = this.pivot_point.p.y + dy_;

    // THE TANGENTIAL COMPONENTS
    
    // now we need to calculate ANY vector that is perpendicular to (nx,ny)
    // the idea is to use this vector as a reference to which we can compare the tangential component of F
    let nx_perp = nx_;
    let ny_perp = -ny_;
    
    
    // DIRECTION OF TORQUE (use the cross-product)
    
    // torque = r*F*sin(theta)
    
    let tau = R.x*this.F.y - R.y*this.F.x; // torque
    
    // IN DEGREES
    this.theta = vector_to_degrees_theta(R.x, R.y);
 

    // from the total forces accumulated on p2, take the component perpendicular to A
    // the dot product between some vector B and (nx_,ny_) = 0 if they are perpendicular
    // the vector projection of F onto (nx_, ny_)

    // (F dot n)(n dot n) n
    let A_ = (this.F.x*nx_+this.F.y*ny_)/(nx_*nx_+ny_*ny_);
    let ux = A_*nx_;
    let uy = A_*ny_;
    
    let Fx_perp = this.F.x - ux;
    let Fy_perp = this.F.y - uy;
    
    let F_tangent = {};
    F_tangent.x = this.F.x - ux;
    F_tangent.y = this.F.y - uy;
    F_tangent.mag = ((F_tangent.x)**2+(F_tangent.y)**2)**0.5;
   
    // this perpendicular force exerts angular acceleration, imposing an angular velocity, and changing the angle between p1 and p2
    // this.alpha = F_tangent.mag / this.particle_mass;
    
    // ALPHA IS ANGULAR ACCELERATION
    this.alpha = tau / this.particle_mass;
    
    // OMEGA IS ANGULAR VELOCITY IN DEGREES PER SECOND
    this.omega += this.alpha * dt_;
    
    // IN DEGREES
    this.theta += this.omega * dt_;
    
    this.p.x = this.pivot_point.p.x + this.pivot_length * Math.cos(this.theta*Math.PI/180);
    this.p.y = this.pivot_point.p.y + this.pivot_length * Math.sin(this.theta*Math.PI/180);


    function vector_to_degrees_theta(x_, y_) {
      
      let R_ = (x_**2+y_**2)**0.5;
      
      let dx_ = Math.abs(x_);
      let dy_ = Math.abs(y_);
      let a_ = Math.asin(dy_ / R_)*180/Math.PI
      
      if (x_ > 0 && y_== 0) {
        return 0;
      }
      
      if (x_ > 0 && y_ > 0) {
        return a_;
      }
      
      if (x_ == 0 && y_ > 0) {
        return 90;
      }
      
      if (x_ < 0 && y_ > 0) {
        return 180 - a_;
      }

      if (x_ < 0 && y_ == 0) {
        return 180;
      }
      
      if (x_ < 0 && y_ < 0) {
        return 180 + a_;
      }

      if (x_ == 0 && y_ < 0) {
        return 270;
      }

      if (x_ > 0 && y_ < 0) {
        return 360 - a_;
      }
      
      return 1;
      
    }; // closing vector_to_degrees_theta
    

 
  }; // closing UPDATE 
  
  
  GET_DISTANCE(other) {
    
    let obj = {};
    
    obj.dx = this.p.x - other.p.x;
    obj.dy = this.p.y - other.p.y;
    obj.r = (obj.dx**2 + obj.dy**2)**0.5;    // SCALAR MAGNITUDE

    // NORMALIZE
    obj.nx = obj.dx/obj.r;
    obj.ny = obj.dy/obj.r;
    
    return obj;
    
  }; // closing GET_DISTANCE

  GET_DISTANCE2(other) {
  // STEP-3. GET THE FORCE VECTOR
  let F = {};
  F.mag = ( F_electric * this.particle_charge * other.particle_charge ) / dr**2;

  // STEP-3b. THE INFLUENCE OF GRAVITY is always A FORCE OF ATTRACTION
  F.mag -= ( F_gravity * this.particle_mass * other.particle_mass ) / dr**2;
  
  // STEP-5. IF THE DISTANCE IS LESS THEN THEIR COMBINED RADII
  if ( dr < (this.particle_radius_of_neutrality + other.particle_radius_of_neutrality) ) {
    // F.mag = maximum_force * this.particle_charge * other.particle_charge;
    F.mag = 0;
  }

  F.x = n.x * F.mag;
  F.y = n.y * F.mag;
  F.r = (F.x**2 + F.y**2)**0.5;
  
  return F;

  }; // closing GET_DISTANCE2
  
  // RETURNS A VECTOR
  GET_FORCE_AB(other) {

    // STEP-1. GET THE DISTANCE
    let dx = this.p.x - other.p.x;
    let dy = this.p.y - other.p.y;
    
    let dr = (dx**2 + dy**2)**0.5;  // magnitude
    
    // STEP-2. THEN NORMALIZE THE DISTANCE VECTOR.
    let n = {};
    n.x = dx/dr;
    n.y = dy/dr;
    n.r = (n.x**2 + n.y**2)**0.5;  // normalized. nr=1

    // STEP-3. GET THE FORCE VECTOR
    let F = {};
    F.mag = ( F_electric * this.particle_charge * other.particle_charge ) / dr**2;

    // STEP-3b. THE INFLUENCE OF GRAVITY is always A FORCE OF ATTRACTION
    F.mag -= ( F_gravity * this.particle_mass * other.particle_mass ) / dr**2;
    
    // STEP-5. IF THE DISTANCE IS LESS THEN THEIR COMBINED RADII
    if ( dr < (this.particle_radius_of_neutrality + other.particle_radius_of_neutrality) ) {
      // F.mag = maximum_force * this.particle_charge * other.particle_charge;
      F.mag = 0; 
    }
    
    if (this.particle_name === 'Positive terminal' && other.particle_name === 'Electron') {
      if ( dr < this.particle_radius ) {
        n_current_positive_terminal_electrons++;
      }

      
    }


    F.x = n.x * F.mag;
    F.y = n.y * F.mag;
    F.r = (F.x**2 + F.y**2)**0.5;
    
    return F;
  };
  
  // APPLY THE ACCUMULATED FORCES
  // dt_=1 means that the net force will be allowed to accelerate the particle for 1 second
  // meaning our new velocity will be velocity += acceleration * dt_
  // and the particles position will be subjected to the new velocity for a duration of dt_
  // so dt_ is quite important
  UPDATE_PARTICLE_POSITION(dt_=1) {
    
    
    if (this.pivot_point && this.pivot_length) {
      this.UPDATE_TETHERED_POSITION(dt_);
      return;
    }
    
    if (!this.pivot) {
      // then continue
    }
    
    // impulse. i need a variable like dt.
    // because implicit in this fn right now is the idea that the force acts on the body for 1 whole second
    // obviously it does not
 
    // IMPULSE DURATION DOES NOT IMPACT ACCELERATION
    this.a.x = this.F.x / this.particle_mass;
    this.a.y = this.F.y / this.particle_mass;
    
    // BUT THAT ACCELERATION IS ONLY APPLIED FOR THE DURATION OF THE IMPULSE
    // SO IT IMPACT VELOCITY
    this.v.x += this.a.x * dt_;
    this.v.y += this.a.y * dt_;
    
    this.v.mag = (this.v.x**2 + this.v.y**2)**0.5;
    
    // TEXT ITS MAX VELOCITY
    if (this.v.mag > this.max_particle_speed) {
      this.v.x *= this.max_particle_speed / this.v.mag;
      this.v.y *= this.max_particle_speed / this.v.mag;
      this.v.mag = this.max_particle_speed;
    }

    // AND THAT VELOCITY ONLY INFLUENCES POSITION FOR THAT DURATION ALSO
    this.p.x += this.v.x * dt_;
    this.p.y += this.v.y * dt_;
   
    
    // DAMPING THE VELOCITY
    if (this.v.mag > this.minimum_damping_speed) {
      this.v.x *= this.damping_factor;
      this.v.y *= this.damping_factor;
    }
    
    this.F.x = 0;
    this.F.y = 0;
    
  }; // CLOSING UPDATE_PARTICLE_POSITION
  
  // CONSTRAIN THE PARTICLE IN THE X DIRECTION
  CONSTRAIN_X(a_=-999, b_=999) {

    if (this.p.x < (a_+this.particle_radius)) {
      this.p.x = (a_ + this.particle_radius) + ((a_ + this.particle_radius)-this.p.x);
      this.v.x *= -1;
      return;
    }
    if (this.p.x > (b_-this.particle_radius) ) {
      this.p.x = (b_ - this.particle_radius) - (this.p.x - (b_ - this.particle_radius));
      this.v.x *= -1;
      return;
    }
  };
  
  LOOP_X(a_=-999, b_=999) {

    if (this.p.x <= (a_-this.particle_radius)) {
      this.p.x = (b_+this.particle_radius);
      // this.v.x *= -1;
      return;
    }
    if (this.p.x >= (b_+this.particle_radius) ) {
      this.p.x = (a_+this.particle_radius);
      // this.v.x *= -1;
      return;
    }
  };
  
  CONSTRAIN_Y(a_=-999, b_=999) {

    if (this.p.y < (a_+this.particle_radius)) {
      this.p.y = (a_ + this.particle_radius) + ((a_ + this.particle_radius)-this.p.y);
      this.v.y *= -1;
      return;
    }
    if (this.p.y > (b_-this.particle_radius) ) {
      this.p.y = (b_ - this.particle_radius) - (this.p.y - (b_ - this.particle_radius));
      this.v.y *= -1;
      return;
    }
  };
  
  LOOP_Y(a_=-999, b_=999) {

    if (this.p.y <= (a_-this.particle_radius)) {
      this.p.y = (b_+this.particle_radius);
      return;
    }
    if (this.p.y >= (b_+this.particle_radius) ) {
      this.p.y = (a_-this.particle_radius);
      return;
    }
  };
  
  
  

  
  COLLIDE_WITH_LINE_SEGMENT(arr_, vx_, vy_) {
  // accepts an array of 2 objects
  
    // this will be hard without rotations
    // i really want the line segment to be flat
  
    // the point is P
    // the line segment is AB
    
    // let u be the vector made by joining P and any point along AB (lets say AP or AB, for example)
    // let u1 = Proj(AP onto AB)
    // let u2 = u - u1
    // mag(u2) is therefore the shortest distance between P and AB
  
    let A = arr_[0];
    let B = arr_[1];
    let d = {
      'x':B.x-A.x,
      'y':B.y-A.y
    }
    let P = this.p;
    let u = {
      'x':P.x-A.x,
      'y':P.y-A.y
    }
    
    // the projection of u onto d
    let u1 = {
      'x':1/(d.x**2+d.y**2) * (d.x*u.x + d.y*u.y) * d.x,
      'y':1/(d.x**2+d.y**2) * (d.x*u.x + d.y*u.y) * d.y
    }
    
    // the perpendicular
    // vector u2 represents the shortest distance between P and AB
    let u2 = {
      'x':u.x - u1.x,
      'y':u.y - u1.y      
    }
    
    let mag_u2 = (u2.x**2 + u2.y**2)**0.5;
    // console.log(mag_u2);
    

    
    if ( (mag_u2 <= this.particle_radius) && (this.p.x+this.particle_radius) >= A.x && (this.p.x-this.particle_radius) <= B.x) {
      
      // console.log(Math.random());
      // we need the relative velocity of the circle relative to the line segment
      
      let v_rel = {
        'x':this.v.x - vx_,
        'y':this.v.y - vy_
      };
      
      /* 
       we want to reposition the circle to a new point Q such that
        Q.x = P.x - t*v_rel.x
        Q.y = P.y - t*v_rel.y
        
      we need the projection of the relative velocity onto the perpendicular
      
      we need to move the center of the circle along v_rel
      and Proj(v_rel onto u2) tells us that by moving mag(v_rel) in the v_rel direction
      we move Proj() in the u2 direction
      

      */
      
      // Proj(v_rel onto u2)
      let w = {
        'x':1/(u2.x**2+u2.y**2) * (v_rel.x*u2.x + v_rel.y*u2.y) * u2.x,
        'y':1/(u2.x**2+u2.y**2) * (v_rel.x*u2.x + v_rel.y*u2.y) * u2.y
      }
      
      // the magnitude of v_rel
      let mag_w = (w.x**2+w.y**2)**0.5;
      let wn = {
        'x':w.x/mag_w,
        'y':w.y/mag_w
      }
      
      // we need to know mag(w) so that we know how many units of w we need to move
      // which will also be the number of units of v_rel we need to subtract to get our proper position
      
      let z = this.particle_radius - mag_u2; // the distance we need to cover in the direction of u2
      let c = z / mag_w;  // the number of units of vector w we need, which is also how many v_rels we need
      
      this.p.x -= c*v_rel.x;
      this.p.y -= c*v_rel.y;
      
      
      // segment AB trades momentum along its perpendicular
      // so P gives up Prog(v_rel onto perp) in one direction
      // but gains it in the other direction
      

      
      
      this.v.x = this.v.x - w.x + (-w.x);
      this.v.y = this.v.y - w.y + (-w.y);
      
      this.p.x += c*this.v.x;
      this.p.y += c*this.v.y;
      
      return;

/*
      // the unit vector of u2 (the perpendicular)
      let nx = u2.x / mag_u2;
      let ny = u2.y / mag_u2;
      
      // now we need the projection of the relative velocity onto the perpendicular u2
      // projection of v_rel onto u2

*/    

    }  // closing collision if statement

  };
  
  
  COLLIDE(other) {
    
    let obj = this.GET_DISTANCE(other);

    // IF THE CIRCLES OVERLAP, THEN THERE WAS A COLLISION
    if (obj.r < this.particle_radius + other.particle_radius) {
      
      // BACKTRACK THE PARTICLE POSITIONS BY 1 FRAME
      let Apx = this.p.x - this.v.x;
      let Avx = this.v.x;
      let Apy = this.p.y - this.v.y;
      let Avy = this.v.y;
      
      let Bpx = other.p.x - other.v.x;
      let Bvx = other.v.x;
      let Bpy = other.p.y - other.v.y;
      let Bvy = other.v.y;

      let alpha = Apx - Bpx;
      let beta = Avx - Bvx;
      let gamma = Apy - Bpy;
      let delta = Avy - Bvy;
      let rho = this.particle_radius + other.particle_radius;
      
      let a_ = beta**2 + delta**2;
      let b_ = 2*(alpha*beta+gamma*delta);
      let c_ = alpha**2+gamma**2-rho**2;
      
      let t1_ = -0.5*(b_/a_) + 0.5*((b_**2-4*a_*c_)**0.5)/a_;
      let t2_ = -0.5*(b_/a_) - 0.5*((b_**2-4*a_*c_)**0.5)/a_;
      let t_ = 0;
      
      // t will be the positive value closest to 0, while less than 1
      
      let t1_is_good = true;
      let t2_is_good = true;
      
      // CANNOT BE NEGATIVE
      if (t1_ < 0) {t1_is_good = false}
      if (t2_ < 0) {t2_is_good = false}
      
      // MUST BE LESS THAN 1
      if (t1_ > 1) {t1_is_good = false}
      if (t2_ > 1) {t2_is_good = false}
      
      // IF ONE OR THE OTHER IS GOOD
      if (t1_is_good && !t2_is_good) {
        t_ = t1_;
      }  
      
      if (!t1_is_good && t2_is_good) {
        t_ = t2_;
      }  
      
      // IF THEY ARE BOTH BETWEEN 0 AND 1, THEN TAKE THE SMALLER OF THE 2
      if (t1_is_good && t2_is_good) {
        if (t1_ < t2_) {t_ = t1_}
        if (t2_ < t1_) {t_ = t2_}
      }
      
      // NOT SURE WHAT TO DO IN THIS CASE
      if (!t1_is_good && !t2_is_good) {
        t_ = 0;
      }
      
      return {
        't1':t1_,
        't2':t2_,
        't':t_,
        't_desc':'this tells us when the collision happened, as a ratio of the current frame'
      }

    }
        
  };  // CLOSING COLLIDE()
  
  
  dies() {
    console.log(this.particle_name + ' was unable to obtain even one single apple.');
    console.log('As a consequence, ' + this.particle_name + ' has died.');
    console.log('His status has been changed from alive to dead.');
    console.log('And the cause of death is written as starvation.');
  };
};

class Electron extends Particle {

  constructor(obj) {
    super(obj);
    this.particle_name = 'Electron';
    this.particle_charge = -1;
    this.particle_mass = 1;
    this.particle_radius = 1;
    this.particle_radius_of_neutrality = 0;
  };
};


class Proton extends Particle {

  constructor(obj) {
    super(obj);
    this.particle_name = 'Proton';
    this.particle_charge = 1;
    this.particle_mass = 1000;
    this.particle_radius = 1;
    this.particle_radius_of_neutrality = 0;
  };
};

class Nucleus extends Particle {

  constructor(obj) {
    super(obj);
    this.particle_name = 'Nucleus';
    this.n_protons = (obj.n_protons || 3);
  };
  
  
  

};

class Copper extends Particle {

  constructor(obj) {
    super(obj);
    this.particle_name = 'Copper';
    this.particle_charge = 17;
    this.particle_mass = 63546;
    this.particle_radius = 5;
    this.particle_radius_of_neutrality = 20;
  };
};

class Positive_terminal extends Particle {

  constructor(obj) {
    super(obj);
    this.particle_name = 'Positive terminal';
    this.particle_charge = (obj.particle_charge || 55);
    this.particle_mass = 1;
    this.particle_radius = 20;
    this.particle_radius_of_neutrality = 7;
  };

};


class Negative_terminal extends Particle {

  constructor(obj) {
    super(obj);
    this.particle_name = 'Negative terminal';
    this.particle_charge = (obj.particle_charge || -55);
    this.particle_mass = 1;
    this.particle_radius = 10;
    this.particle_radius_of_neutrality = 10;
  };

};

Particle.prototype.DISPLAY_PARTICLE = function(b) {

  // PROTONS
  if (this.particle_name === 'Proton') {

    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius,
      'line_width':1,
      'stroke_style':'#f007',
      'fill_style':'#ffff'
    });
    return;
  }

  // COPPER
  if (this.particle_name === 'Copper') {

    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius*2,
      'line_width':1,
      'stroke_style':'#fc00',
      'fill_style':'#fc02'
    });
    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius*3.5,
      'line_width':1,
      'stroke_style':'#fc00',
      'fill_style':'#fc02'
    });
    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius*5,
      'line_width':0.5,
      'stroke_style':'#fc0f',
      'fill_style':'#fc01'
    });
    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius*6.5,
      'line_width':1,
      'stroke_style':'#fc00',
      'fill_style':'#ffcc0007'
    });
    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius_of_neutrality,
      'line_width':1,
      'stroke_style':'#fc0f',
      'fill_style':'#fc00'
    });
    b.FILL_STYLE('#fc0f');
    b.RADIUS(1);
    b.POINT(this.p);

    return;
  }
  
  // ELECTRONS
  if (this.particle_name === 'Electron') {
    b.FILL_STYLE('#47c7');
    b.RADIUS(this.particle_radius);
    b.POINT(this.p);
    return;
  }

  // POSITIVE TERMINAL
  if (this.particle_name === 'Positive terminal') {
    
    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius,
      'line_width':1,
      'stroke_style':'#f007',
      'fill_style':'#f003'
    });
    
    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius*1.25,
      'line_width':1,
      'stroke_style':'#f000',
      'fill_style':'#f002'
    });
    
    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius*1.5,
      'line_width':1,
      'stroke_style':'#f000',
      'fill_style':'#f001'
    });
    
    b.FILL_STYLE('#f00f');
    b.RADIUS(this.particle_radius);
    // b.POINT(this.p);
    return;
  }
 
  // NEGATIVE TERMINAL
  if (this.particle_name === 'Negative terminal') {
    
    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius*2,
      'line_width':1,
      'stroke_style':'#47ce',
      'fill_style':'#47c3'
    });
    
    b.FILL_STYLE('#47ce');
    b.RADIUS(this.particle_radius);
    // b.POINT(this.p);
    return;
  }
  
  // ELSE
  if (this.particle_name === 'Particle') {

    b.CIRCLE({
      'val':this.p,
      'r':this.particle_radius,
      'line_width':1,
      'stroke_style':'#f00a',
      'fill_style':'#fff0'
    });
    return;
  }

};