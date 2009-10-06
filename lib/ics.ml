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
    val make_cold_spherical : int -> b array
    val make_plummer : int -> b array
    val from_nbody2_file : string -> int -> b array
    val make_hot_spherical : int -> b array
    val make_spherical_out_of_equilibrium : int -> float -> b array
    val rescale_to_standard_units : b array -> b array
  end

module Make(B : BODY) : (IC with type b = B.b)  = 
  struct
    module A = Analysis.Make(B)
    module E = Energy.Make(B)

    type b = B.b

    let make_body = B.make

    let adjust_frame bs = 
      let p_tot = A.total_momentum bs and 
	  com = A.center_of_mass bs and 
	  m_tot = A.total_mass bs in 
      Array.iter 
	(fun b -> 
	  let m = B.m b and q = B.q b and p = B.p b in 
	  for i = 0 to 2 do 
	    q.(i) <- q.(i) -. com.(i);
	    p.(i) <- p.(i) -. p_tot.(i)*.m/.m_tot
	  done)
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
      adjust_frame bs;
      bs
	
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
      adjust_frame bs;
      bs

    let make_hot_spherical n =
      let v0 = sqrt (21.0/.10.0) and
	  m = 1.0/.(float_of_int n) in 
      let bs = 
	Array.init n
	  (fun i -> 
	    let r = random_from_dist square and
		v = random_between 0.0 v0 in 
	    make_body 0.0 m (random_vector r) (random_vector (m *. v))) in 
      adjust_frame bs;
      bs      

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
      let ke = E.total_kinetic_energy bs and 
          pe = E.total_potential_energy bs in 
      let p_factor = sqrt (0.25 /. ke) and 
          r_factor = 1.0/.((-0.5)/.pe) in 
      Array.map 
        (fun b -> 
          let t = B.t b and 
              m = B.m b and 
              q = B.q b and 
              p = B.p b in 
          B.make t m (Array.map (( *. ) r_factor) q) (Array.map (( *. ) p_factor) p))
        bs             
  end
