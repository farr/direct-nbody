(*  analysis.ml: Code for analyzing systems of bodies.
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

open Base

(** Output type of {!Analysis.Make}. *)
module type ANALYSIS = sig
  (** Bodies. *)
  type b

  (** Embedded energy calculating module. *)
  module E : Energy.ENERGY with type b = b

  (** [nearest_neighbors n bs b] returns the [n] nearest neighbors to
      [b] out of the bodies [bs]. *)
  val nearest_neighbors : int -> b array -> b -> b array

  (** [density_estimator n bs b] returns the
      [n]th-nearest-neighbor-based density estimator for the local
      density at the position of [b].  See S. Casertano, P. Hut.  Core
      radius and density measurements in N-body experiments:
      connections with theoretical and observational definitions.
      Astrophysical Journal, vol. 298, p. 80, 1985. *)
  val density_estimator : int -> b array -> b -> float

  (** [density_squared_estimator n bs b] returns an unbiased estimator
      for the density squared in the system [bs] at the position of
      [b].  See Casertano and Hut. *)
  val density_squared_estimator : int -> b array -> b -> float

  (** [density_center rhos bs] returns the location of the
      density-weighted center of the system [bs].  [rhos] is an array
      of containing the density estimate at the location of each [b].
      See Casertano and Hut
      ({!Analysis.ANALYSIS.density_estimator}).*)
  val density_center : float array -> b array -> float array

  (** [density_radius n bs] returns the density-weighted radius
      estimate for the core radius (using [n] nearest neighbors to
      estimate the density).  See Casertano and Hut
      ({!Analysis.ANALYSIS.density_estimator}).*)
  val density_radius : int -> b array -> float

  (** [density_density n bs] returns the density-weighted core density
      estimate (using [n] nearest neighbors to estimate the density).
      See Casertano and Hut ({!Analysis.ANALYSIS.density_estimator}).*)
  val density_density : int -> b array -> float

  (** [core_crossing_time n bs] uses the [n]th nearest neighbor
      density estimator to compute a core density for the system [bs],
      and then returns the core crossing time. *)
  val core_crossing_time : int -> b array -> float

  (** Returns the total momentum of the given system. *)
  val total_momentum : b array -> float array

  (** Returns the total mass of the given system. *)
  val total_mass : b array -> float

  (** Returns the center of mass of the given system. *)
  val center_of_mass : b array -> float array

  (** Returns the total angular momentum of the system about the
      origin. *)
  val total_angular_momentum : b array -> float array

  (** Given a filename that stores a sequence of marshalled [b
      array]s, load the given states and return them in order. *)
  val load_states : string -> b array array 

  (** [load_state n filename] loads the [n]th state from the given
      file, which should contain the concatenation of marshalled [b array]s. *)
  val load_state : int -> string -> b array

  (** Returns the system time. *)
  val time : b array -> float

  (** [boxcar_average n ys] returns an array of length equal to [ys]
      where each value is an average of the [n] nearest values in [ys]
      (this smoothing operation is called a width-[n] boxcar
      average). *)
  val boxcar_average : int -> float array -> float array

  val zip : 'a array -> 'b array -> ('a * 'b) array

  val unzip : ('a * 'b) array -> ('a array) * ('b array)

  (** [with_random_state state thunk] applies [thunk] with the random
      state set to [state] for the extent of the application,
      returning the random state to its initial value when [thunk]
      finishes or raises an exception.  The result is the return value
      of [thunk].*)
  val with_random_state : Random.State.t -> (unit -> 'a) -> 'a

  (** [perturbed_bodies scale bs] generates a new system where each
      position and momentum have been adjusted by adding a random
      increment with typical size [scale].*)
  val perturbed_bodies : float -> b array -> b array
    
  (** [polynomial_interpolation xs ys] returns a function that
      computes the value of the interpolating polynomial that passes
      through each pair of [(x,y)].*)
  val polynomial_interpolation : float array -> float array -> (float -> float)

  (** Computes the binding energy of the binary composed of the given
      bodies (which can be positive, if the bodies are unbound relative to
      each other). *)
  val binary_energy : b -> b -> float

  (** Returns the binary angular momentum, in the orbit frame. *)
  val binary_angular_momentum : b -> b -> float array

  (** Returns a list of all binaries in the given system and their
      binding energy. *)
  val binaries : b array -> ((b * b) * float) list

  (** Returns the most tightly bound binary in the system. *)
  val tightest_binary : b array -> b * b 

  (** Returns the SMA for bound orbits. *)
  val semi_major_axis : b -> b -> float

  (** Returns the eccentricity of an orbit. *)
  val eccentricity : b -> b -> float

  (** Returns the inclination of an orbit. *)
  val inclination : b -> b -> float

  (** Returns T such that ke = 3/2 N T. *)
  val body_temperature : b array -> float

  (** [bodies_to_el_phase_space bs] returns an array containing the
      energy, specific energy, angular momentum, and specific angular
      momentum of the corresponding member of [bs]. *)
  val bodies_to_el_phase_space : b array -> float array array

  (** [lagrange_radius ?origin bs frac] returns the radius inside
      which a fraction of at least [frac] of the total mass of the
      system resides.  [frac] must be between 0.0 (inclusive) and 1.0
      (exclusive).  If [origin] is given, the radius is with respect
      to it; otherwise the radius is with respect to the center of
      mass of the system. *)
  val lagrange_radius : ?origin : float array -> b array -> float -> float

  (** Given a system, produce a body with the same total mass, center of
      mass, and total momentum. *)
  val summary_body : b array -> b
end

module Make (B : Body.BODY) : ANALYSIS with type b = B.b = struct
  module E = Energy.Make(B)

  type b = B.b

  let pi = 4.0*.atan 1.0

  let all_equal compare xs start endd = 
    let x0 = xs.(start) in 
    let rec ae_loop i = 
      if i >= endd then 
        true
      else if compare x0 xs.(i) = 0 then 
        ae_loop (i+1)
      else
        false in 
      ae_loop (start+1)
        
  let find_nthf ?(copy = true) compare nth xs = 
    if nth < 0 || nth >= Array.length xs then raise (Invalid_argument "find_nthf: nth outside array bounds");
    let xs = if copy then Array.copy xs else xs in 
    let rec find_nth_loop start nth endd = 
      if endd - start = 1 then 
        if nth <> 0 then 
          raise (Failure "find_nthf: nth index not in array bounds")
        else
          xs.(start)
      else if all_equal compare xs start endd then 
        xs.(start)
      else begin
        let part = xs.(start + (Random.int (endd - start))) in 
        let rec swap_loop low high = 
          let new_low = let rec new_low_loop l = 
                          if l >= endd then l else if compare xs.(l) part <= 0 then new_low_loop (l+1) else l in new_low_loop low and
              new_high = let rec new_high_loop h = 
                           if h < start then h else if compare xs.(h) part > 0 then new_high_loop (h-1) else h in new_high_loop high in
            if new_low > new_high then new_low else if new_low >= endd then endd else if new_high < start then new_low else begin
              let tmp = xs.(new_low) in 
                xs.(new_low) <- xs.(new_high);
                xs.(new_high) <- tmp;
                swap_loop new_low new_high
            end in
        let ilow = swap_loop start (endd - 1) in
          if nth < (ilow - start) then 
            find_nth_loop start nth ilow
          else
            find_nth_loop ilow (nth - (ilow - start)) endd
      end in 
      find_nth_loop 0 nth (Array.length xs)

  let sorted_array_insert comp arr obj = 
    let n = Array.length arr in 
      if comp obj arr.(n-1) < 0 then 
        let rec loop = function 
          | 0 -> 
            if comp obj arr.(0) < 0 then begin
              arr.(1) <- arr.(0);
              arr.(0) <- obj
            end else
              arr.(1) <- obj
          | i -> 
            if comp obj arr.(i) < 0 then begin
              arr.(i+1) <- arr.(i);
              loop (i - 1)
            end else
              arr.(i+1) <- obj in 
          loop (n - 2)
            
  let nearest_neighbors n bs b = 
    let q = B.q b in 
    let distance b1 b2 = 
      let q1 = B.q b1 and q2 = B.q b2 in 
        if q1 == q then 
          1
        else if q2 == q then 
          -1
        else 
          let d1 = Base.distance_squared q q1 and 
              d2 = Base.distance_squared q q2 in 
            if d1 > d2 then 
              1
            else if d1 < d2 then 
              -1
            else 
              0 in 
    let candidates = Array.sub bs 0 n in 
      Array.fast_sort distance candidates;
      for i = n to Array.length bs - 1 do 
        sorted_array_insert distance candidates bs.(i)
      done;
      candidates

  let fold_left2 fn start a1 a2 = 
    let n = Array.length a1 in 
      assert( n = Array.length a2); 
      let value = ref (fn start a1.(0) a2.(0)) in 
        for i = 1 to n - 1 do
          value := fn !value a1.(i) a2.(i)
        done;
        !value      

  let density_estimator n bs b = 
    let neighbors = nearest_neighbors n bs b in 
    let r = Base.distance (B.q b) (B.q neighbors.(n-1)) in 
    let mtot = (Array.fold_left (fun mt b -> mt +. (B.m b)) 0.0 neighbors) -.
      (B.m neighbors.(n-1)) in 
      3.0*.mtot/.(4.0*.pi*.(Base.cube r)) 

  let density_squared_estimator n bs b = 
    let rho = density_estimator n bs b and 
        nf = float_of_int n in 
      (nf -. 2.0) /. (nf -. 1.0) *. rho *. rho

  let density_center rhos bs = 
    let rhotot = Array.fold_left (+.) 0.0 rhos in 
      fold_left2 
        (fun rd b r -> 
          let q = B.q b in 
            for i = 0 to 2 do 
              rd.(i) <- rd.(i) +. r/.rhotot*.q.(i)
            done;
            rd)
        (Array.make 3 0.0)
        bs 
        rhos 
        
  let density_radius n bs = 
    let rhos = Array.map (density_estimator n bs) bs in 
    let dc = density_center rhos bs in 
    let rtot = Array.fold_left (+.) 0.0 rhos in 
      fold_left2 
        (fun dr b r -> 
          let d = Base.distance dc (B.q b) in 
            dr +. d*.r/.rtot)
        0.0
        bs
        rhos

  let density_density n bs = 
    let rhos = Array.map (density_estimator n bs) bs in 
    let rho_tot = Array.fold_left (+.) 0.0 rhos in 
      Array.fold_left 
        (fun dd r -> 
          dd +. (Base.square r)/.rho_tot)
        0.0 
        rhos      

  (* We have that the circular velocity at radius r_c for an object of
     mass m_c is v_circ = sqrt(m_c/r_c).  This implices a
     crossing-timescale of t_cr = r_c/v_circ = r_c^(3/2)/m_c^(1/2).
     But, if m_c = 4/3*pi*rho_c*r_c^3, then we have t_cr =
     1/(4/3*pi*rho_c)^(1/2). *)
  let core_crossing_time n bs = 
    let rhoc = density_density n bs in 
    1.0/.(sqrt ((4.0/.3.0)*.pi*.rhoc))

  let total_momentum bs = 
    Array.fold_left 
      (fun pt b -> 
        let p = B.p b in 
        for i = 0 to 2 do 
          pt.(i) <- pt.(i) +. p.(i)
        done;
        pt)
      (Array.make 3 0.0)
      bs

  let total_mass bs = 
    Array.fold_left 
      (fun mt b ->
        let m = B.m b in 
        mt +. m)
      0.0
      bs

  let center_of_mass bs = 
    let mt = total_mass bs in 
    Array.fold_left 
      (fun cm b -> 
        let m = B.m b and q = B.q b in 
        for i = 0 to 2 do 
          cm.(i) <- cm.(i) +. q.(i)*.(m/.mt)
        done;
        cm)
      (Array.make 3 0.0)
      bs

  let total_angular_momentum bs = 
    Array.fold_left 
      (fun l b ->
        let lb = E.angular_momentum b in
        for i = 0 to 2 do 
          l.(i) <- l.(i) +. lb.(i)
        done;
        l)
      (Array.make 3 0.0)
      bs

  let load_states file = 
    let input = open_in file in 
    let bs = ref [] in 
    let rec loop () = 
      try 
        let new_bs = (input_value input : B.b array) in 
        bs := new_bs :: !bs;
        loop ()
      with 
      | _ -> close_in input; Array.of_list (List.rev !bs) in
    loop ()

  let load_state n file = 
    let input = open_in file in 
    let rec loop i = 
      try 
        let bs = (input_value input : B.b array) in 
        if i = n then bs else loop (i + 1) 
      with x -> close_in input; raise x in 
    loop 0

  let time bs = 
    B.t bs.(0)

  let interpolating_diffs xs ys = 
    let n = Array.length xs in 
    let diffs = Array.make n 0.0 in 
    diffs.(0) <- (ys.(1) -. ys.(0)) /. (xs.(1) -. xs.(0));
    for i = 1 to n - 2 do
      diffs.(i) <- (ys.(i+1) -. ys.(i-1)) /. (2.0 *. (xs.(i+1) -. xs.(i-1)))
    done;
    diffs.(n-1) <- (ys.(n-1) -. ys.(n-2)) /. (xs.(n-1) -. xs.(n-2));
    diffs

  let rough_rhos bs = 
    let com = center_of_mass bs in
    let distance b1 b2 = 
      let r1 = Base.distance (B.q b1) com and 
          r2 = Base.distance (B.q b2) com in 
      compare r1 r2 in
    let bs = Array.copy bs in 
    Array.fast_sort distance bs;
    let ms = 
      Array.of_list (List.rev 
                       (Array.fold_left
                          (fun ms b -> 
                            match ms with 
                            | [] -> [B.m b]
                            | m :: rest -> (m +. (B.m b)) :: ms)
                          [] 
                          bs)) in                      
    let rs = Array.map (fun b -> Base.distance (B.q b) com) bs in 
    let diffs = interpolating_diffs rs ms in 
    let pi = 4.0 *. atan 1.0 in
    Array.mapi (fun i dmdr -> 
      let r = rs.(i) in 
      (r, dmdr/.(4.0*.pi*.(Base.square r)))) diffs

  let merge_states (states : B.b array list) = 
    Array.concat states 

  let even n = n mod 2 = 0
  let odd n = not (even n)

  let boxcar_average n ys = 
    let nys = Array.length ys in 
    if n > nys then raise 
        (Invalid_argument ("cannot boxcar average over" ^ 
         " larger interval than array length"));
    let range = n / 2 in 
    Array.mapi (fun i y -> 
      let lower_i = max (i - range) 0 and 
          upper_i = min (nys - 1) (i + range) in 
      let delta = float_of_int (upper_i - lower_i + 1) in 
      let sum = ref 0.0 in 
      for i = lower_i to upper_i do 
        sum := !sum +. ys.(i)
      done;
      !sum /. delta)
      ys

  let zip a b = 
    Array.mapi 
      (fun i a -> 
        let b = b.(i) in 
        (a,b))
      a

  let unzip a = 
    (Array.map fst a, Array.map snd a)

  let boxcar_average_file n file = 
    let inbuf = Scanf.Scanning.from_file file in 
    let rec loop xs ys = 
      if Scanf.Scanning.end_of_input inbuf then 
        zip (Array.of_list (List.rev xs)) 
	  (boxcar_average n (Array.of_list (List.rev ys)))
      else
        let x, y = Scanf.bscanf inbuf " %g %g " (fun x y -> (x,y)) in 
        loop (x :: xs) (y :: ys) in 
    loop [] []

  let with_random_state s thunk = 
    let s_old = Random.get_state () in 
    try 
      Random.set_state s;
      let result = thunk () in 
      Random.set_state s_old;
      result
    with x -> Random.set_state s_old; raise x

  let perturbed_bodies scale bs = 
    let new_bs = Array.map B.copy bs in 
    Array.iter 
      (fun b -> 
        let q = B.q b and p = B.p b in 
        for i = 0 to 2 do 
          q.(i) <- q.(i) +. scale*.(Random.float 1.0 -. 0.5);
          p.(i) <- p.(i) +. scale*.(Random.float 1.0 -. 0.5)
        done)
      new_bs;
    new_bs

  let polynomial_interpolation xs ys = 
    let xs = Array.copy xs and 
        ys = Array.copy ys and 
        n = Array.length ys in 
    assert (n = Array.length xs);
    fun x -> 
      let storage = Array.copy ys in 
      for range = 1 to n - 1 do 
        for i = 0 to n - (1 + range) do 
          let upper = i + range in 
          storage.(i) <- ((x -. xs.(upper)) *. storage.(i) +. 
                            ((xs.(i) -. x) *. storage.(i+1))) /. 
            (xs.(i) -. xs.(upper))
        done
      done;
      storage.(0)

  let map2 fn a b = 
    let n = Array.length a in 
    assert (Array.length b = n);
    let result = Array.make n (fn a.(0) b.(0)) in 
    for i = 1 to n - 1 do 
      result.(i) <- fn a.(i) b.(i)
    done;
    result

  let relative_error a b = 
    abs_float(2.0 *. (a -. b) /. ((abs_float a) +. (abs_float b)))

  let richardson_extrapolate nmax tol h0 err_inc f = 
    let fs = Array.make nmax 0.0 and 
        hs = Array.make nmax 0.0 in 
    let rec loop i h est err = 
      if i >= nmax then 
        if err > tol then 
          raise (Failure "could not coverge")
        else 
          est
      else begin
        hs.(i) <- h;
        fs.(i) <- f h;
        if i >= 2 then begin
          let f0i = (polynomial_interpolation
                       (Array.sub hs 0 i) 
                       (Array.sub fs 0 i)) 0.0 and
              f0im1 = (polynomial_interpolation
                         (Array.sub hs 0 (i - 1))
                         (Array.sub fs 0 (i - 1))) 0.0 in 
          let new_err = relative_error f0i f0im1 in 
          if new_err > err_inc *. err then 
            if err < tol then 
              est
            else
              raise (Failure "error not within tolerance")
          else if new_err < err then 
            loop (i + 1) (h /. 2.0) f0i new_err
          else
            loop (i + 1) (h /. 2.0) est err
        end else
          loop (i + 1) (h /. 2.0) est err 
      end in 
    loop 0 h0 (0.0/.0.0) infinity
      
  
  let sum = Array.fold_left (+.) 0.0 

  let average a = (sum a) /. (float_of_int (Array.length a))

  let flatten arr = 
    Array.of_list
      (List.flatten
         (Array.to_list
            (Array.map
               Array.to_list
               arr)))

  let delta_q_delta_p bs rbs = 
    let qs, ps = 
      (unzip 
         (map2 
            (fun b rb -> 
              (map2 (-.) (B.q b) (B.q rb),
               map2 (-.) (B.p b) (B.p rb)))
            bs
            rbs)) in 
    (flatten qs, flatten ps)

  let poincaire_integral_invariant r b0 b1 = 
    let (deltaq0, deltap0) = delta_q_delta_p b0 r and 
        (deltaq1, deltap1) = delta_q_delta_p b1 r in 
    (sum (map2 ( *. ) deltaq0 deltap1)) -. 
      (sum (map2 ( *. ) deltap0 deltaq1))

  let estimate_at_zero h f = 
    let fm = f (~-. h) and 
        fm2 = f (h /. (-2.0)) and 
        fp2 = f (h /. 2.0) and 
        fp = f h in 
    (1.0/.6.0*.(4.0*.fm2 -. fm +. 4.0*.fp2 -. fp), 
     abs_float(2.0/.3.0*.(fm -. fm2 +. fp -. fp2)))

  let delta_pushforward_omega h0 map bs = 
    let s = Random.get_state () in 
    let f h = 
      with_random_state s (fun () ->
        let bs0 = perturbed_bodies h bs and 
            bs1 = perturbed_bodies h bs in 
        abs_float (1.0 -. 
                     (poincaire_integral_invariant 
                        (map bs) (map bs0) (map bs1))/.
                   (poincaire_integral_invariant bs bs0 bs1))) in 
    estimate_at_zero h0 f

  let explore_pushforward_space map sfs hs bs = 
    List.map 
      (fun sf -> 
        List.fold_left
          (fun ((log_sf_min, log_dI_min, rel_err_min) as a)
              ((log_sf, log_dI, rel_err) as b) -> 
                if rel_err < rel_err_min then 
                  b
                else 
                  a)
          (0.0, 0.0, infinity)
          (List.map (fun h -> 
            let (x, err) = delta_pushforward_omega h (map sf) bs in 
            (sf**0.2, x, abs_float(err/.x)))
             hs))
      sfs           
  
  let binary_energy b1 b2 = 
    let vrel = Array.make 3 0.0 in 
    let m1 = B.m b1 and 
        p1 = B.p b1 and 
        p2 = B.p b2 and
        m2 = B.m b2 in 
    for i = 0 to 2 do 
      vrel.(i) <- p2.(i)/.m2 -. p1.(i)/.m1
    done;
    let mu = m1*.m2/.(m1+.m2) in 
    let ke = 0.5 *. mu *. (dot vrel vrel) and 
        pe = E.potential_energy b1 b2 in 
    ke +. pe

  let binary_angular_momentum b1 b2 = 
    let m1 = B.m b1 and m2 = B.m b2 and 
        r1 = B.q b1 and r2 = B.q b2 and 
        p1 = B.p b1 and p2 = B.p b2 in 
    let mu = m1*.m2/.(m1+.m2) in
    let v1 = Array.map (fun x -> x /. m1) p1 and 
        v2 = Array.map (fun x -> x /. m2) p2 in 
    let r = Array.mapi (fun i x -> x -. r1.(i)) r2 and 
        v = Array.mapi (fun i x -> x -. v1.(i)) v2 in 
    let lspec = cross r v in 
      Array.map (fun x -> x *.mu) lspec

  (* Returns a list of ((b1,b2), e_rel). *)
  let binaries bs = 
    let binaries = ref [] in
    for i = 0 to Array.length bs - 1 do 
      let bi = bs.(i) in 
      for j = i+1 to Array.length bs - 1 do 
        let bj = bs.(j) in 
        let e = binary_energy bi bj in 
        if e < 0.0 then 
          binaries := ((bi,bj), e) :: !binaries
      done
    done;
    !binaries

  let tightest_binary bs = 
    let bins = binaries bs in 
    match bins with 
    | [] -> raise (Failure "no binaries in system")
    | bin :: bins -> 
        fst
          (List.fold_left
             (fun ((b1b2_best, e_best) as best) ((b1b2, e) as elt) -> 
               if abs_float e > abs_float e_best then 
                 elt
               else
                 best)
             bin
             bins)

  let semi_major_axis b1 b2 = 
    let ebin = binary_energy b1 b2 and 
        m1 = B.m b1 and 
        m2 = B.m b2 in 
      if ebin > 0.0 then 
        raise (Failure "semi_major_axis: no SMA for unbound orbit")
      else
        ~-.m1*.m2/.(2.0*.ebin)

  let eccentricity b1 b2 =
    let l = norm (binary_angular_momentum b1 b2) and 
        e = binary_energy b1 b2 in
    let m1 = B.m b1 and m2 = B.m b2 in 
    let mu = m1*.m2/.(m1+.m2) and 
        alpha = m1*.m2 in 
      sqrt (1.0 +. 2.0*.e*.l*.l/.(mu*.alpha*.alpha))

  let inclination b1 b2 = 
    let l = binary_angular_momentum b1 b2 in 
      acos (l.(2) /. (norm l))      

  let body_temperature bs = 
    let n = float_of_int (Array.length bs) and 
        e = E.total_kinetic_energy bs in 
    2.0*.e/.(3.0*.n)

  let bodies_to_el_phase_space bs = 
    Array.mapi 
      (fun i b -> 
        let lb = Base.norm (E.angular_momentum b) and 
            ke = E.kinetic_energy b and 
            pe = ref 0.0 in 
          for j = 0 to i - 1 do 
            pe := !pe +. E.potential_energy b bs.(j)
          done;
          for j = i + 1 to Array.length bs - 1 do 
            pe := !pe +. E.potential_energy b bs.(j)
          done;
          let etot = ke +. !pe in 
          let m = B.m b in
            [| etot; etot /. m; lb; lb /. m |]) 
      bs

  let lagrange_radius ?origin bs frac = 
    assert(0.0 <= frac && frac < 1.0);
    let origin = match origin with 
      | Some(o) -> o
      | None -> center_of_mass bs in 
    let compare b1 b2 = 
      let q1 = B.q b1 and 
          q2 = B.q b2 in 
      let d1 = Base.distance_squared q1 origin and 
          d2 = Base.distance_squared q2 origin in 
        Pervasives.compare d1 d2 in 
    let bs = Array.copy bs in
      Array.fast_sort compare bs;
      let mtot = total_mass bs in 
      let rec lagrange_radius_loop i m = 
        if i >= Array.length bs then raise (Failure "lagrange_radius") else
          let m = m +. B.m bs.(i) in 
            if m /. mtot >= frac then begin
              Base.distance (B.q bs.(i)) origin
            end else
              lagrange_radius_loop (i+1) m in
        lagrange_radius_loop 0 0.0

  let summary_body bs = 
    B.make (B.t bs.(0)) (total_mass bs) (center_of_mass bs) (total_momentum bs)
end
