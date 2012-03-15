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

(** Initial conditions of various types. *)

open Body;;
open Base;;

(** Output type of the {!Ics.Make} functor. *)
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

    (** [random_from_dist ?xmin ?xmax ?ymin ?ymax f] draws a random
        number from the distribution given by [f].  The number will be
        in the range [xmin] to [xmax]; [ymin] to [ymax] should bound
        the magnitude of [f] over this range of arguments. *)
    val random_from_dist : ?xmin : float -> ?xmax : float -> ?ymin : float -> 
      ?ymax : float -> (float -> float) -> float

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

    (** [make_flat_sis n] produces a cold distribution of stars with a
        flattened singular-isothermal-sphere profile: rho(r) = 1 if r
        < 1, rho(r) = 1/r^2 if 1 < r < 90/8, rho(r) = 0 otherwise. *)
    val make_flat_sis : int -> b array

    (** [make_hubble n] produces a cold distribution of stars with a
        Hubble density profile: rho(r) = 1/(1+r^2) for r < 40/pi. *)
    val make_hubble : int -> b array

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

    (** [elements_to_rv m1 m2 elts] returns [(r1,v1,r2,v2)] for the
        binary with masses [(m1,m2)] and the given orbital
        elements. *)
    val elements_to_rv : float -> float -> orbit_elements -> 
      (float array * float array * float array * float array)

    (** [generate_binary m1 m2 elts] generates a binary system with the
        given masses and orbital elements, in the center-of-mass frame. *)
    val generate_binary : float -> float -> orbit_elements -> b array

    (** [shift_system bs r0 v0] moves the origin of the coordinate
        system of [bs] to [r0], and the velocity reference frame to
        [v0]. *)
    val shift_system : b array -> float array -> float array -> b array

    (** [combine_systems bs1 bs2] combines the systems of bodies into
        a single system. *)
    val combine_systems : b array -> b array -> b array

    (** [random_elements amin amax] Produce a random set of orbital
        elements with sma log-distributed between amin and amax,
        eccentricity thermally distributed. *)
    val random_elements : float -> float -> orbit_elements

    (** [king_r_and_m_samples w0] returns (rs, ms, ws), radii, masses,
        and potential values for a King model with central potential
        w0.  The units are such that the central density is 1, the
        core radius is 1, and G = 1. *)
    val king_r_and_m_samples : float -> (float array) * (float array) * (float array)

    (** [make_king w0 n] returns a king model in the COM frame with N
        bodies and central concentration w0.  The units are chosen so
        that G = 1, M = 1, and the central radius, r_c = 1.  This should
        be {b nearly} standard units. *)
    val make_king : float -> int -> b array

    (** Returns the analytic rho^2 dM weighted radius for a king model
        with the given central potential. *)
    val king_analytic_density_squared_radius : float -> float

    (** Returns the fraction of the total mass contained within the
        {!Ics.IC.king_analytic_density_squared_radius}. *)
    val king_analytic_density_squared_core_fraction : float -> float

    (** [make_from_cmc_snapshot file] returns [(bs, v, gradV)].  Where
        [bs] are the black holes in the cluster, [v] is the potential
        from the other stars, and [gradV] computes the gradient of v.
        The units are such that the total mass of black holes = 1, the
        half-mass radius of the black holes = 1, and G = 1.  The black
        holes will be translated into the center of momentum frame. *)
    val make_from_cmc_snapshot : string -> b array * (b -> float) * (b -> float array)
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

  let pi = 3.1415926535897932385

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

  let make_flat_sis n = 
    let m = 1.0 /. (float_of_int n) in 
    let bs = 
      Array.init n 
        (fun i -> 
          let mfrac = Random.float 1.0 in 
          let mfrac1 = 4.0 /. 127.0 in
          let r = 
            if mfrac <= mfrac1 then 
              (* Inside the core. *)
              (mfrac /. mfrac1)**(1.0/.3.0)
            else
              (mfrac /. mfrac1 +. 2.0)/.3.0 in
            make_body 0.0 m (random_vector r) (Array.make 3 0.0)) in
      adjust_frame bs

  let solve f xmin xmax = 
    let rec loop xmin xmax fmin fmax = 
      if xmax -. xmin < 1e-8 then 
        (fmax*.xmin -. fmin*.xmax) /. (fmax -. fmin)
      else 
        let xmid = 0.5*.(xmin +. xmax) in 
        let fmid = f xmid in 
          if fmid *. fmin < 0.0 then 
            loop xmin xmid fmin fmid
          else
            loop xmid xmax fmid fmax in
      loop xmin xmax (f xmin) (f xmax)

  let make_hubble n = 
    let m = 1.0 /. (float_of_int n) in 
    let bs = Array.init n (fun i -> 
      let mfrac = Random.float 1.0 in 
      let f r = pi*.((atan r) -. r) /. (-40.0 +. pi*.(atan (40.0/.pi))) -. mfrac in 
      let r = solve f 0.0 (40.0/.pi) in 
        make_body 0.0 m (random_vector r) (Array.make 3 0.0)) in 
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
            (fxmax*.xmin -. fxmin*.xmax)/.(fxmax -. fxmin)
          else
            let xmid = 0.5*.(xmin +. xmax) in 
            let fxmid = f xmid in 
              if fxmid *. fxmin <= 0.0 then 
                bisect_loop fxmin fxmid xmin xmid
              else
                bisect_loop fxmid fxmax xmid xmax in 
        bisect_loop (f xmin) (f xmax) xmin xmax

    let solve_kepler e m = 
      if e = 0.0 then 
        m
      else
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
        rotate_z (rotate_x (rotate_z v om) i) capom in 
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

    let elements_to_rv m1 m2 elts = 
      let (rxy, vxy) = orbit_rv_xy elts.m elts.a elts.e in 
      let (r, v1) = xy_to_orbit rxy vxy elts.i elts.capom elts.omega in 
      let a = elts.a in 
      let a3 = a*.a*.a in 
      let n = sqrt ((m1+.m2)/.a3) in 
      let v = Array.map (fun vx -> vx*.n) v1 in 
      let v1 = Array.map (fun v -> ~-.m2*.v/.(m1+.m2)) v and 
          v2 = Array.map (fun v -> m1*.v/.(m1+.m2)) v in 
      let r1 = Array.map (fun r -> ~-.m2*.r/.(m1+.m2)) r and 
          r2 = Array.map (fun r -> m1*.r/.(m1+.m2)) r in 
        (r1,v1,r2,v2)

    let generate_binary m1 m2 elts = 
      let (r1,v1,r2,v2) = elements_to_rv m1 m2 elts in 
      let p1 = Array.map (fun x -> x*.m1) v1 and 
          p2 = Array.map (fun x -> x*.m2) v2 in 
        [|make_body 0.0 m1 r1 p1;
          make_body 0.0 m2 r2 p2|]

    let shift_system bs r0 v0 = 
      Array.map 
        (fun b -> 
          let t = B.t b and 
              r = B.q b and 
              p = B.p b and 
              m = B.m b in 
          let r' = Array.mapi (fun i rx -> rx +. r0.(i)) r and 
              p' = Array.mapi (fun i px -> px +. m*.v0.(i)) p in
            make_body t m r' p')
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

    (* This is an analytic formula for psi(W) from King 1966. *)
    let psi w = 
      let sqrt_w = sqrt w in 
      let sqrt_pi = 1.7724538509055160273 in 
        0.25*.(3.0*.(exp w)*.sqrt_pi*.(Gsl_sf.erf sqrt_w) -. 2.0*.sqrt_w*.(3.0+.2.0*.w))

    let rho_normalized w0 w = 
      (psi w) /. (psi w0)

    (* RHS for King model equations when R <~ 1. *)
    let small_R_RHS w0 w ys dydws = 
      let x = ys.(0) and y = ys.(1) in 
      let rho = rho_normalized w0 w in
      let y2 = y*.y in 
      let y3 = y2*.y in 
        dydws.(0) <- y; (* dX/dW = Y *)
        dydws.(1) <- 1.0/.x*.(9.0/.4.0*.rho*.y3 +. 3.0/.2.0*.y2); (* dY/dW = 1/X(9/4 rho y^3 + 3/2 y^2) *)
        dydws.(2) <- 6.2831853071795864769*.y*.rho*.(sqrt x) (* dM/dW = 2 pi rho X^(1/2) Y *)

    (* RHS for King model equations when R >~ 1. *) 
    let large_R_RHS w0 w ys dydws = 
      let p = ys.(0) and q = ys.(1) in 
      let p2 = p*.p and q2 = q*.q in 
      let p4 = p2*.p2 and q3 = q*.q2 in 
      let rho = rho_normalized w0 w in 
        dydws.(0) <- q; (* dP/dW = Q *)
        dydws.(1) <- 9.0*.rho/.p4*.q3; (* dQ/dW = 9 rho Q^3/P^4 *)
        dydws.(2) <-  (-12.566370614359172954)*.rho*.q/.p4 (* dM/dW = - 4 pi rho Q / P^4 *)

    let xy_state_to_pq_state ys = 
      let x = ys.(0) and y = ys.(1) in 
      let x12 = sqrt x in 
        ys.(0) <- 1.0 /. x12;
        ys.(1) <- ~-.y /. (2.0*.x*.x12)

    let king_r_and_m_samples w0 = 
      let eps = 1e-8 in 
      let h0 = min ((sqrt eps) /. 100.0) (w0 /. 2.0) in 
      let small_rhs w ys dydws = small_R_RHS w0 w ys dydws and 
          large_rhs w ys dydws = large_R_RHS w0 w ys dydws in 
      let small_system = Gsl_odeiv.make_system small_rhs 3 and 
          large_system = Gsl_odeiv.make_system large_rhs 3 in 
      let step = Gsl_odeiv.make_step Gsl_odeiv.RK8PD ~dim:3 in
      let control = Gsl_odeiv.make_control_y_new ~eps_abs:eps ~eps_rel:eps in 
      let evolve = Gsl_odeiv.make_evolve 3 in 
      let w = w0 -. h0 in 
      let ys = [|2.0/.3.0*.h0; -2.0/.3.0; 0.0|] in 
      let get_r y system = 
        if system == small_system then 
          sqrt y.(0) 
        else
          1.0 /. y.(0) in 
      let rec loop system h w ys rs ms ws = 
        let (w_new, h_new) = Gsl_odeiv.evolve_apply evolve control step system ~t:w ~t1:0.0 ~h:h ~y:ys in 
          if w_new <= 0.0 then 
            (* Done; Normalize the masses, and return. *)
            ((Array.of_list (List.rev ((get_r ys system) :: rs))),
             (Array.of_list (List.rev (ys.(2) :: ms))),
             (Array.of_list (List.rev (w_new :: ws))))
          else
            let r_new = get_r ys system and 
                r_old = List.hd rs in 
              if r_new > 1.0 && r_old <= 1.0 then begin
                (* Swap X = R^2, P = 1/R. *)
                xy_state_to_pq_state ys;
                Gsl_odeiv.step_reset step;
                Gsl_odeiv.evolve_reset evolve;
                loop large_system h_new w_new ys (r_new :: rs) (ys.(2) :: ms) (w_new :: ws)
              end else begin
                (* Continue with the same system. *)
                loop system h_new w_new ys (r_new :: rs) (ys.(2) :: ms) (w_new :: ws)
              end in 
        loop small_system (~-.h0) w ys [get_r ys small_system] [0.0] [w]
              
    let king_v_cumulative jv jve = 
      let jv2 = jv*.jv and 
          jve2 = jve*.jve in 
      let jv3 = jv2*.jv and 
          jve3 = jve2*.jve in 
      let ejve2 = exp jve2 in 
        (6.0*.(exp (jve2 -. jv2))*.jv +. 4.0*.jv3 -. 5.3173615527165480819*.ejve2*.(Gsl_sf.erf jv)) /. 
          (6.0*.jve +. 4.0*.jve3 -. 5.3173615527165480819*.ejve2*.(Gsl_sf.erf jve))

    let binary_search (x : float) (xs : float array) = 
      let n = Array.length xs in 
        assert(x >= xs.(0));
        assert(x <= xs.(n-1));
      let rec loop imin imax = 
        if imax - imin <= 1 then 
          (imin, imax)
        else
          let imid = imin + (imax - imin) / 2 in 
            if x < xs.(imid) then 
              loop imin imid
            else
              loop imid imax in 
        loop 0 (n-1)

    let r_from_m m rs ms = 
      let (imin, imax) = binary_search m ms in 
      let frac = (ms.(imax) -. m)/.(ms.(imax) -. ms.(imin)) in 
        rs.(imin) +. frac*.(rs.(imax) -. rs.(imin))

    let w_from_m m ws ms = 
      let (imin, imax) = binary_search m ms in 
      let frac = (ms.(imax) -. m)/.(ms.(imax) -. ms.(imin)) in 
        ws.(imin) +. frac*.(ws.(imax) -. ws.(imin))

    let j mmax = 0.59841342060214901691*.(sqrt mmax)

    let draw_king_body mbody ms rs_interp ws_interp = 
      let n = Array.length ms in 
      let mmax = ms.(n-1) in
      let m = Random.float mmax in 
      let r = Gsl_interp.eval rs_interp m and 
          w = Gsl_interp.eval ws_interp m in
      let j = j mmax in 
      let jve = sqrt w in
      let vfrac = Random.float 1.0 in
      let jv = bisect_solve 1e-8 (fun jv -> (king_v_cumulative jv jve) -. vfrac) 0.0 jve in
      let v = jv /. j in 
        B.make 0.0 mbody (random_vector r) (random_vector (mbody*.v))

    let make_king w0 n = 
      let mbody = 1.0 /. (float_of_int n) and 
          (rs,ms,ws) = king_r_and_m_samples w0 in 
      let rs_interp = Gsl_interp.make_interp Gsl_interp.AKIMA ms rs and 
          ws_interp = Gsl_interp.make_interp Gsl_interp.AKIMA ms ws in
      let bs = Array.init n (fun _ -> draw_king_body mbody ms rs_interp ws_interp) in 
        adjust_frame bs

    (* Want to compute the density-squared-weighted king radius
       (integral done over mass). *)
    let king_analytic_density_squared_radius w0 = 
      let (rs, ms, ws) = king_r_and_m_samples w0 in 
      let mmax = ms.(Array.length ms - 1) in
      let numerator_integrand = Array.mapi (fun i w -> let rho = rho_normalized w0 w in rs.(i)*.rho*.rho) ws and 
          denominator_integrand = Array.map (fun w -> let rho = rho_normalized w0 w in rho*.rho) ws in 
      let num_interp = Gsl_interp.make_interp Gsl_interp.AKIMA ms numerator_integrand and 
          denom_interp = Gsl_interp.make_interp Gsl_interp.AKIMA ms denominator_integrand in 
      let num = Gsl_interp.eval_integ num_interp 0.0 mmax and 
          den = Gsl_interp.eval_integ denom_interp 0.0 mmax in 
        num /. den

    let king_analytic_density_squared_core_fraction w0 = 
      let (rs,ms,ws) = king_r_and_m_samples w0 in 
      let rc = king_analytic_density_squared_radius w0 in 
      let m_interp = Gsl_interp.make_interp Gsl_interp.AKIMA rs ms in 
        (Gsl_interp.eval m_interp rc) /. ms.(Array.length ms - 1)

    let mi = 1
    let ri = 2
    let vri = 3
    let vti = 4 
    let sti = 14
    let phii = 23

    let read_cmc chan = 
      let lexbuf = Lexing.from_channel chan in 
        Cmc_parser.main Cmc_lexer.tokenize lexbuf

    let split_black_holes cmcs = 
      let bhs = ref [] and 
          others = ref [] in 
        for i = 0 to Array.length cmcs - 1 do
          let star = cmcs.(i) in 
            if star.(sti) = 14.0 then 
              bhs := star :: !bhs
            else
              others := star :: !others
        done;
        let bhs = Array.of_list (List.rev !bhs) and 
            others = Array.of_list (List.rev !others) in 
          Array.fast_sort (fun o1 o2 -> if o1.(ri) < o2.(ri) then -1 else if o1.(ri) > o2.(ri) then 1 else 0) others;
          (bhs, others)
        
    let cmc_total_mass cmcs = 
      let mtot = ref 0.0 in 
        for i = 0 to Array.length cmcs - 1 do 
          mtot := !mtot +. cmcs.(i).(mi)
        done;
        !mtot +. 0.0

    let cmc_half_mass_radius cmcs = 
      let mtot = cmc_total_mass cmcs in 
      let n = Array.length cmcs in
      let rec loop m i = 
        if i >= n then 
          raise (Failure "cmc_half_mass_radius")
        else
          let mi = cmcs.(i).(mi) in 
            if m +. mi > 0.5*.mtot then 
              cmcs.(i).(ri)
            else
              loop (m+.mi) (i+1) in 
        loop 0.0 0

    let cmc_to_new_units m r cmcs = 
      for i = 0 to Array.length cmcs - 1 do 
        let star = cmcs.(i) in 
          star.(mi) <- star.(mi) /. m;
          star.(ri) <- star.(ri) /. r;
          star.(vri) <- star.(vri) /. r;
          star.(vti) <- star.(vti) /. r
      done

    let cmc_to_body cmc = 
      let r = cmc.(ri) and 
          m = cmc.(mi) and 
          vr = cmc.(vri) and 
          vt = cmc.(vti) in 
      let rhat = random_vector 1.0 in
      let t = cross rhat (random_vector 1.0) in 
      let tnorm = norm t in 
      let r = Array.map (fun x -> r *. x) rhat and 
          p = Array.mapi (fun i x -> m*.(vr *. x +. vt *. t.(i) /. tnorm)) rhat in 
        B.make 0.0 m r p

    let cmcs_to_potential cmcs = 
      let n = Array.length cmcs in 
      let rs = Array.make (n+1) 0.0 and 
          fs = Array.make (n+1) 0.0 and 
          menc = ref 0.0 in
        for i = 0 to n - 1 do 
          let star = cmcs.(i) in 
            rs.(i+1) <- star.(ri);
            menc := !menc +. star.(mi);
            fs.(i+1) <- !menc /. (rs.(i+1)*.rs.(i+1));
            if rs.(i+1) <= rs.(i) then rs.(i+1) <- rs.(i)*.(1.0+.1e-8)
        done;
        let mmax = cmcs.(n-1).(mi) and 
            rmax = cmcs.(n-1).(ri) in
        let f_of_r = Gsl_interp.make_interp Gsl_interp.AKIMA rs fs in
          fun b -> 
            let r = norm (B.q b) and 
                m = B.m b in 
              if r > rmax then 
                ~-.m*.mmax/.r
              else
                m*.(~-.(Gsl_interp.eval_integ f_of_r r rmax) -. mmax/.rmax)

    let cmcs_to_gradV cmcs = 
      let n = Array.length cmcs in
      let rs = Array.make (n+1) 0.0 and 
          mencs = Array.make (n+1) 0.0 in 
      let mtot = ref 0.0 in 
        for i = 0 to n - 1 do 
          let star = cmcs.(i) in 
            mtot := !mtot +. star.(mi);
            rs.(i+1) <- star.(ri);
            if rs.(i+1) <= rs.(i) then rs.(i+1) <- (1.0 +. 1e-8)*.rs.(i);
            mencs.(i+1) <- !mtot
        done;
        let rmax = cmcs.(n-1).(ri) in 
        let mmax = cmcs.(n-1).(mi) in 
        let menc_of_r = Gsl_interp.make_interp Gsl_interp.AKIMA rs mencs in 
          fun b -> 
            let q = B.q b in
            let br = norm q in 
            let menc = if br > rmax then mmax else Gsl_interp.eval menc_of_r br in 
            let fnorm = (B.m b)*.menc/.(br*.br) in 
              Array.map (fun x -> fnorm*.x/.br) q (* Grad_v is -F, so prop. to rhat. *)

    let make_from_cmc_snapshot file = 
      let inp = open_in file in 
      let cmcs = read_cmc inp in 
        close_in inp;
        let (bhs, others) = split_black_holes cmcs in 
        let mscale = cmc_total_mass bhs and 
            rscale = cmc_half_mass_radius bhs in 
          cmc_to_new_units mscale rscale bhs;
          cmc_to_new_units mscale rscale others;
          let pot = cmcs_to_potential others in 
          let force = cmcs_to_gradV others in
          let bs = Array.map cmc_to_body bhs in 
            (adjust_frame bs), pot, force
end
