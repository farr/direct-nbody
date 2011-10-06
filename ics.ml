(*  ics.ml: Generate initial conditions from various distributions.
    Copyright (C) 2006--2008 Will M. Farr <farr@mit.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*)

open Body;;
open Base;;

module type IC = 
  sig
    type b

    (** Orbital elements for binary description. *)
    type orbit_elements = {
      m : float; (** Mean anomaly. *)
      a : float; (** Semi-major axis. *)
      e : float; (** Eccentricity. *)
      i : float; (** Inclination. *)
      capom : float; (** Longitude of ascending node. *)
      omega : float (** Argument of periapsis. *)
    }

    (** Construct a constant-density spherical distribution with the
        given number of bodies at zero kinetic energy. *)
    val make_cold_spherical : int -> b array

    (** Consturct a plummer model with the given number of bodies. *)
    val make_plummer : int -> b array

    (** Read in the state from the given file that is the format
        produced by NBODY2. *)
    val from_nbody2_file : string -> int -> b array

    (** Make a virialized, constant-density sphere with velocities
        chosen uniformly in a sphere about the velocity-space origin. *)
    val make_hot_spherical : int -> b array

    (** Like {!Ics.IC.make_hot_spherical}, but with the given virial
        ratio (KE/PE).  *)
    val make_spherical_out_of_equilibrium : int -> float -> b array

    (** Construct a new system that represents the given system, but
        in "standard" units: Mtot = 1, E = -1/4, virial ratio =
        0.5. *)
    val rescale_to_standard_units : b array -> b array

    (** Adjust the coordinate frame of the system to the
        center-of-mass frame. *)
    val adjust_frame : b array -> b array

    (** [add_mass_spectrum gen_mass bs] Construct a new system with
        masses chosen according to the [gen_mass] function.  *)
    val add_mass_spectrum : (unit -> float) -> b array -> b array

    (** [generate_binary m1 m2 elts] generates a binary system with the
        given masses and orbital elements, in the center-of-mass frame. *)
    val generate_binary : float -> float -> orbit_elements -> b array

    (** [shift_system bs r0] moves the origin of the coordinate system
        of [bs] to [r0]. *)
    val shift_system : b array -> float array -> b array

    (** [combine_systems bs1 bs2] combines the systems of bodies into
        a single system. *)
    val combine_systems : b array -> b array -> b array

  end

module Make(B : BODY) : (IC with type b = B.b)  = 
struct
  module A = Analysis.Make(B)
  module E = Energy.Make(B)

  type b = B.b

  type orbit_elements = {
    m : float;
    a : float;
    e : float;
    i : float;
    capom : float;
    omega : float
  }

  let make_body = B.make

  let adjust_frame bs = 
    let p_tot = A.total_momentum bs and 
	com = A.center_of_mass bs and 
	m_tot = A.total_mass bs in 
      Array.map
	(fun b -> 
	  let m = B.m b and q = B.q b and p = B.p b in 
          let qnew = Array.make 3 0.0 and 
              pnew = Array.make 3 0.0 in
	    for i = 0 to 2 do 
	      qnew.(i) <- q.(i) -. com.(i);
	      pnew.(i) <- p.(i) -. p_tot.(i)*.m/.m_tot
	    done;
            B.make (B.t b) m qnew pnew)
	bs

  let random_between a b = 
    let dx = b -. a in 
      a +. Random.float dx

  let random_from_dist ?(xmin = 0.0) ?(xmax = 1.0) ?(ymin = 0.0) ?(ymax = 1.0) dist = 
    let rec loop x y = 
      if y < dist x then 
	x
      else 
	loop (random_between xmin xmax) (random_between ymin ymax) in 
      loop (random_between xmin xmax) (random_between ymin ymax)

  let random_vector r = 
    let pi = 4.0 *. atan 1.0 in 
    let cos_theta = random_between (-1.0) 1.0 and 
	phi = random_between 0.0 (2.0 *. pi) in 
    let theta = acos cos_theta in 
    let sin_theta = sin theta in 
      [| r *. (cos phi) *. sin_theta;
	 r *. (sin phi) *. sin_theta;
	 r *. cos_theta |]
	
  let make_cold_spherical n = 
    let m = 1.0 /. (float_of_int n) in 
    let bs = 
      Array.init n
	(fun i -> 
	  make_body 0.0 m
	    (random_vector 
	       (random_from_dist ~xmax: (12.0/.5.0) 
		  ~ymax: (144.0/.25.0) square))
	    (Array.make 3 0.0)) in 
      adjust_frame bs
	
  let make_plummer n = 
    let m = 1.0 /. (float_of_int n) and 
	pi = 4.0 *. atan 1.0 in 
    let sf = 16.0 /. (3.0 *. pi) in 
    let bs = 
      Array.init n 
	(fun i -> 
	  let r = 1.0 /. 
	    (sqrt ((random_between 0.0 1.0)**(-2.0/.3.0) -. 1.0)) and 
	      x = random_from_dist ~ymax:0.1 
	    (fun x -> (square x)*.(1.0 -. (square x))**3.5) in 
	  let q = random_vector (r /. sf) and 
	      p = random_vector 
	    ((sqrt sf) *. m *. x *. (sqrt 2.0) *. 
		(1.0 +. (square r))**(-0.25)) in 
	    make_body 0.0 m q p) in 
      adjust_frame bs

  let make_hot_spherical n =
    let v0 = sqrt (21.0/.10.0) and
	m = 1.0/.(float_of_int n) in 
    let bs = 
      Array.init n
	(fun i -> 
	  let r = random_from_dist square and
	      v = random_between 0.0 v0 in 
	    make_body 0.0 m (random_vector r) (random_vector (m *. v))) in 
      adjust_frame bs

  let from_nbody2_file name n = 
    let in_buffer = Scanf.Scanning.from_file name in 
      Array.init n (fun i -> 
	let m = Scanf.bscanf in_buffer " %g " (fun x -> x) in 
	let q = Array.init 3 (fun i -> 
	  Scanf.bscanf in_buffer " %g " (fun x -> x)) in 
	let p = Array.init 3 (fun i -> 
	  Scanf.bscanf in_buffer " %g " (fun v -> m*.v)) in 
	  make_body 0.0 m q p)

    (* [q] is desired KE/PE ratio (i.e. q = 0.5 for virial equilibrium). *)
    let make_spherical_out_of_equilibrium n q =
      assert (q >= 0.0 && q < 1.0);
      let bs = make_hot_spherical n and 
	  a = sqrt(q/.(1.0-.q)) and
	  b = 2.0*.(1.0 -. q) in 
      Array.iter 
	(fun body -> 
	  let q = B.q body and p = B.p body in 
	  for i = 0 to 2 do 
	    p.(i) <- p.(i) *. a;
	    q.(i) <- q.(i) *. b
	  done)
	bs;
      bs

    let rescale_to_standard_units bs = 
      let mtot = A.total_mass bs in 
      let bs = 
        Array.map 
          (fun b -> 
            let t = B.t b and q = B.q b and p = B.p b and m = B.m b in 
            let mnew = m /. mtot in 
            let pnew = Array.make 3 0.0 in 
              for i = 0 to 2 do 
                pnew.(i) <- p.(i) *. mnew /. m
              done;
              B.make t mnew q pnew)
          bs in
      let e = E.energy bs in
        assert(e < 0.0);
        let efactor = ((-0.25) /. e) in 
          Array.map 
            (fun b -> 
              let q = B.q b and 
                  p = B.p b in 
              let qnew = Array.make 3 0.0 and 
                  pnew = Array.make 3 0.0 in 
                for i = 0 to 2 do 
                  qnew.(i) <- q.(i) /. efactor;
                  pnew.(i) <- p.(i) *. (sqrt efactor)
                done;
                B.make (B.t b) (B.m b) qnew pnew)
            bs

    let rescale_mass b mnew = 
      let m = B.m b and 
          p = B.p b in 
      let pnew = Array.make 3 0.0 in 
        for i = 0 to 2 do 
          pnew.(i) <- p.(i) *. mnew/.m
        done;
        B.make (B.t b) mnew (B.q b) pnew

    let add_mass_spectrum gen_mass bs = 
      Array.map (fun b -> rescale_mass b (gen_mass ())) bs

    let gen_salpeter mmin mmax = 
      let ind = -2.3 in
      let ind1 = ind +. 1.0 in 
        ((mmax**ind1 -. mmin**ind1)*.(Random.float 1.0) +. mmin**ind1)**(1.0 /. ind1)
          
    let kepler_eqn e manom eanom =
      manom -. eanom +. e*.(sin eanom)

    let bisect_solve epsabs f xmin xmax = 
      let rec bisect_loop fxmin fxmax xmin xmax =
        assert (fxmin *. fxmax <= 0.0);
        let dx = xmax -. xmin in 
          if dx <= epsabs then 
            ~-.fxmin*.dx /. (fxmax -. fxmin)
          else
            let xmid = 0.5*.(xmin +. xmax) in 
            let fxmid = f xmid in 
              if fxmid *. fxmin <= 0.0 then 
                bisect_loop fxmin fxmid xmin xmid
              else
                bisect_loop fxmid fxmax xmid xmax in 
        bisect_loop (f xmin) (f xmax) xmin xmax

    let solve_kepler e m = 
      let f ecc = kepler_eqn e m ecc in 
        bisect_solve 1e-12 f (m +. e) (m -. e)

    let rotate_z v theta = 
      let c = cos theta and 
          s = sin theta in 
        match v with 
          | [|x; y; z|] -> 
            [|c*.x -. s*.y; c*.y +. s*.x; z|]
          | _ -> raise (Invalid_argument "rotate_z: expected three-vector")

    let rotate_x v theta = 
      let c = cos theta and 
          s = sin theta in 
        match v with 
          | [|x; y; z|] -> 
            [|x; c*.y -. s*.z; c*.z +. s*.y|]
          | _ -> raise (Invalid_argument "rotate_x: expected three-vector")

    let xy_to_orbit r v i capom om = 
      let rotate v = 
        rotate_z (rotate_x (rotate_z v om) capom) i in 
        (rotate r, rotate v)

    let orbit_rv_xy m a e = 
      let ecc = solve_kepler e m in 
      let se = sin ecc and 
          ce = cos ecc and 
          a2 = a*.a and 
          l = sqrt (1.0 -. e*.e) in
      let x = a*.(ce -. e) and
          y = a*.l*.se in 
      let r = sqrt (x*.x +. y*.y) in 
      let vx = ~-.a2*.se/.r and 
          vy = a2*.l*.ce/.r in 
        ([|x; y; 0.0|], [|vx; vy; 0.0|])

    let generate_binary m1 m2 elts = 
      let (rxy, vxy) = orbit_rv_xy elts.m elts.a elts.e in 
      let (r, v1) = xy_to_orbit rxy vxy elts.i elts.capom elts.omega in 
      let a = elts.a in 
      let a3 = a*.a*.a in 
      let n = sqrt ((m1+.m2)/.a3) in 
      let v = Array.map (fun vx -> vx*.n) v1 in 
      let p1 = Array.map (fun v -> ~-.m1*.m2*.v/.(m1+.m2)) v and 
          p2 = Array.map (fun v -> m1*.m2*.v/.(m1+.m2)) v in 
      let r1 = Array.map (fun r -> ~-.m2*.r/.(m1+.m2)) r and 
          r2 = Array.map (fun r -> m1*.r/.(m1+.m2)) r in 
        [|make_body 0.0 m1 r1 p1;
          make_body 0.0 m2 r2 p2|]

    let shift_system bs r0 = 
      Array.map 
        (fun b -> 
          let t = B.t b and 
              r = B.q b and 
              p = B.p b and 
              m = B.m b in 
          let r' = Array.mapi (fun i rx -> rx +. r0.(i)) r in 
            make_body t m r' p)
        bs

    let combine_systems : (b array -> b array -> b array) = Array.append 

    let pi = 4.0*.(atan 1.0)
    let twopi = 2.0*.pi

    let random_elements amin amax = 
      let lamax = log amax and 
          lamin = log amin in 
      let la = random_between lamin lamax in 
      let a = exp la and 
          m = random_between 0.0 twopi and 
          capom = random_between 0.0 twopi and 
          omega = random_between 0.0 twopi and 
          cos_i = random_between (-1.0) 1.0 and
          e2 = random_between 0.0 1.0 in 
      let e = sqrt e2 and 
          i = acos cos_i in 
        {m = m; a = a; e = e; i = i; capom = capom; omega = omega}        
  end
