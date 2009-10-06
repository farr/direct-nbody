(*  dFloat_advancer.ml: Automitac differentiated version of the advancer code.
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

open DFloat;;

module DFloatBase = struct 
  let nsteps = ref 0

  let eps2 = ref (C 0.0)

  let square x = x*.x

  let cube x = x*.x*.x

  let set_eps eps = 
    eps2 := square (C eps)

  let norm_squared v = 
    let sum = ref (C 0.0) in 
    for i = 0 to 2 do 
      sum := !sum +. square v.(i)
    done;
    !sum

  let norm v = sqrt (norm_squared v)

  let distance_squared x y = 
    let sum = ref (C 0.0) in 
    for i = 0 to 2 do 
      sum := !sum +. square (x.(i) -. y.(i))
    done;
    !sum

  let distance x y = sqrt (distance_squared x y)

  let dot x y = 
    let sum = ref (C 0.0) in 
    for i = 0 to 2 do 
      sum := !sum +. x.(i)*.y.(i)
    done;
    !sum

  let grad_v m1 q1 m2 q2 gv = 
    let r2 = distance_squared q1 q2 in 
    let re = r2 +. !eps2 in 
    let r3 = re*.(sqrt re) in 
    let factor = m1*.m2/.r3 in 
    for i = 0 to 2 do 
      gv.(i) <- factor*.(q1.(i) -. q2.(i))
    done

  let grad_v_dgrad_v m1 q1 p1 m2 q2 p2 gv dgv = 
    for i = 0 to 2 do
      gv.(i) <- q1.(i) -. q2.(i);
      dgv.(i) <- p1.(i)/.m1 -. p2.(i)/.m2
    done;
    let r2 = norm_squared gv and 
	rdv = dot gv dgv in 
    let re = r2 +. !eps2 in 
    let r3 = (sqrt re)*.re in 
    let r5 = r3*.re in 
    let factorg = m1*.m2/.r3 in
    let factord = (C 3.0)*.rdv*.m1*.m2/.r5 in 
    for i = 0 to 2 do 
      dgv.(i) <- factorg*.dgv.(i) -. factord*.gv.(i);
      gv.(i) <- factorg*.gv.(i)
    done

(* let relative_error x t = 
   abs_float ((x -. t) /. t);; *)

end

module A = Advancer.A

type db = {
    id : int;
    mutable t : diff;
    m : diff;
    q : diff array;
    p : diff array;
    mutable tau : diff;
    mutable h : diff;
    f : diff array;
    df : diff array;
    ddf : diff array;
    q1 : diff array;
    q2 : diff array;
    dVdq0 : diff array;
    dVdq1 : diff array;
    dVdq2 : diff array;
  }

type b = A.b

let promote_arr a = 
  Array.map (fun x -> C x) a

let promote_body 
    {A.id = id; A.t = t; A.m = m; A.q = q;
     A.p = p; A.tau = tau; A.h = h; A.f = f; A.df = df; A.ddf = ddf;
     A.q1 = q1; A.q2 = q2;
     A.dVdq0 = dVdq0; A.dVdq1 = dVdq1; A.dVdq2 = dVdq2} = 
  {id = id;
   t = C(t);
   m = C(m);
   q = promote_arr q;
   p = promote_arr p;
   tau = C(tau);
   h = C(h);
   f = promote_arr f;
   df = promote_arr df;
   ddf = promote_arr ddf;
   q1 = promote_arr q1;
   q2 = promote_arr q2;
   dVdq0 = promote_arr dVdq0;
   dVdq1 = promote_arr dVdq1;
   dVdq2 = promote_arr dVdq2}

let gen_id =
  let id = ref 0 in 
  fun () -> 
    let res = !id in 
    incr id;
    res

let t b = b.t
let m b = b.m
let dq b = b.q
let dp b = b.p
let q {A.q = q} = q
let p {A.p = p} = p
let next_t b = b.t +. b.h

let body_copy {id = id; t = t; m = m; q = q; p = p; tau = tau; f = f;
	       df = df; ddf = ddf} = 
  {id = id;
   t = t;
   m = m;
   q = Array.copy q;
   p = Array.copy p;
   tau = tau;
   h = tau;
   f = Array.copy f;
   df = Array.copy df;
   ddf = Array.copy ddf;
   q1 = Array.make 3 (C 0.0);
   q2 = Array.make 3 (C 0.0);
   dVdq0 = Array.make 3 (C 0.0);
   dVdq1 = Array.make 3 (C 0.0);
   dVdq2 = Array.make 3 (C 0.0);}

let make_body t m q p = 
  {id = gen_id ();
   t = t;
   m = m;
   q = q;
   p = p;
   tau = (C (-1.0));
   h = (C (-1.0));
   f = Array.make 3 (C 0.0);
   df = Array.make 3 (C 0.0);
   ddf = Array.make 3 (C 0.0);
   q1 = Array.make 3 (C 0.0);
   q2 = Array.make 3 (C 0.0);
   dVdq0 = Array.make 3 (C 0.0);
   dVdq1 = Array.make 3 (C 0.0);
   dVdq2 = Array.make 3 (C 0.0);}

let make = make_body
let copy = body_copy

let kinetic_energy {m = m; p = p} = 
  (DFloatBase.norm_squared p)/. ((C 2.0) *. m)

let potential_energy {m = m1; q = q1} {m = m2; q = q2} = 
  let r = DFloatBase.distance q1 q2 in 
  zero-.m1*.m2/.r

let energy bs = 
  let eg = ref (C 0.0) in 
  for i = 0 to Pervasives.(-) (Array.length bs) 1 do 
    let b = bs.(i) in 
    eg := !eg +. (kinetic_energy b);
    for j = Pervasives.(+) i 1 to Pervasives.(-) (Array.length bs) 1 do 
      eg := !eg +. (potential_energy b bs.(j))
    done
  done;
  !eg

let hermite_predict {m = m; q = q; p = p; f = f; df = df; ddf = ddf; q1 = q1; q2 = q2} dt = 
  let dt2 = dt /. (C 2.0) in 
  for i = 0 to 2 do 
    (* Changed 24.0 in ddf term to 12.0, per reviewer note. *)
    q1.(i) <- q.(i) +. p.(i)/.m*.dt2 +. 
	f.(i)/.((C 2.0)*.m)*.(DFloatBase.square dt2) +. 
	df.(i)/.((C 6.0)*.m)*.(DFloatBase.cube dt2) +. 
	ddf.(i)/.((C 12.0)*.m)*.(dt2*.(DFloatBase.cube dt2));
    q2.(i) <- q.(i) +. p.(i)/.m*.dt +. 
	f.(i)/.((C 2.0)*.m)*.(DFloatBase.square dt) +. 
	df.(i)/.((C 6.0)*.m)*.(DFloatBase.cube dt) +. 
	ddf.(i)/.((C 24.0)*.m)*.(dt*.(DFloatBase.cube dt))
  done

let q0_coeff dt h = 
  ((h -. ((C 2.0)*.dt))*.(h -. dt))/.(DFloatBase.square h)

let q1_coeff dt h = 
  ((C 4.0) *. (h -. dt) *. dt)/.(DFloatBase.square h)

let q2_coeff dt h = 
  ((((C 2.0)*.dt) -. h)*.dt)/.(DFloatBase.square h)

let predict_position {t = tb; h = h; q = q0; q1 = q1; q2 = q2} t q_pred = 
  let dt = t -. tb in 
  let q0c = q0_coeff dt h and 
      q1c = q1_coeff dt h and 
      q2c = q2_coeff dt h in 
  for i = 0 to 2 do 
    q_pred.(i) <- q0c*.q0.(i) +. q1c*.q1.(i) +. q2c*.q2.(i)
  done

let predict_momentum {t = tb; h = h; m = m; q = q0; q1 = q1; q2 = q2} t p_pred = 
  let dt = t -. tb in 
  let q0c = (m *. (((C 4.0)*.dt) -. ((C 3.0)*.h)))/.(DFloatBase.square h) and 
      q1c = m*.(C 4.0)*.(h -. (C 2.0)*.dt)/.(DFloatBase.square h) and 
      q2c = m*.((C 4.0)*.dt -. h) /. (DFloatBase.square h) in 
  for i = 0 to 2 do 
    p_pred.(i) <- q0c*.q0.(i) +. q1c*.q1.(i) +. q2c*.q2.(i)
  done

let compute_dVs bs i = 
  let b = bs.(i) and 
      gv = Array.make 3 (C 0.0) and 
      q_pred = Array.make 3 (C 0.0) in 
  let {h = h; t = t; m = m; q = q0; q1 = q1; q2 = q2; dVdq0 = dVdq0; dVdq1 = dVdq1; dVdq2 = dVdq2} = b in 
  let wt = h /. (C 6.0) in 
  let wtm = (C 4.0)*.wt and 
      ho2 = h /. (C 2.0) in
  for j = Pervasives.(+) i 1 to Pervasives.(-) (Array.length bs) 1 do 
    let {t = t2; h = h2; m = m2; dVdq0 = dVdq02; dVdq1 = dVdq12; dVdq2 = dVdq22} as b2 = bs.(j) in 
    if h = h2 then
      begin
	DFloatBase.grad_v m q0 m2 b2.q gv;
	for i = 0 to 2 do 
	  dVdq0.(i) <- dVdq0.(i) +. wt*.gv.(i);
	  dVdq02.(i) <- dVdq02.(i) -. wt*.gv.(i)
	done;
	DFloatBase.grad_v m q1 m2 b2.q1 gv;
	for i = 0 to 2 do 
	  dVdq1.(i) <- dVdq1.(i) +. wtm*.gv.(i);
	  dVdq12.(i) <- dVdq12.(i) -. wtm*.gv.(i)
	done;
	DFloatBase.grad_v m q2 m2 b2.q2 gv;
	for i = 0 to 2 do 
	  dVdq2.(i) <- dVdq2.(i) +. wt*.gv.(i);
	  dVdq22.(i) <- dVdq22.(i) -. wt*.gv.(i)
	done
      end 
    else
      begin
	let dt = t -. t2 in 
	let wt0 = wt *. (q0_coeff dt h2) and 
	    wt1 = wt *. (q1_coeff dt h2) and 
	    wt2 = wt *. (q2_coeff dt h2) in 
	predict_position b2 t q_pred;
	DFloatBase.grad_v m q0 m2 q_pred gv;
	for i = 0 to 2 do
	  dVdq0.(i) <- dVdq0.(i) +. wt *. gv.(i);
	  dVdq02.(i) <- dVdq02.(i) -. wt0 *. gv.(i);
	  dVdq12.(i) <- dVdq12.(i) -. wt1 *. gv.(i);
	  dVdq22.(i) <- dVdq22.(i) -. wt2 *. gv.(i);
	done;
	let t = t +. ho2 in
	let dt = t -. t2 in 
	let wt0 = wtm *. (q0_coeff dt h2) and 
	    wt1 = wtm *. (q1_coeff dt h2) and 
	    wt2 = wtm *. (q2_coeff dt h2) in 
	predict_position b2 t q_pred;
	DFloatBase.grad_v m q1 m2 q_pred gv;
	for i = 0 to 2 do
	  dVdq1.(i) <- dVdq1.(i) +. wtm *. gv.(i);
	  dVdq02.(i) <- dVdq02.(i) -. wt0 *. gv.(i);
	  dVdq12.(i) <- dVdq12.(i) -. wt1 *. gv.(i);
	  dVdq22.(i) <- dVdq22.(i) -. wt2 *. gv.(i);
	done;
	let t = t +. ho2 in 
	let dt = t -. t2 in 
	let wt0 = wt *. (q0_coeff dt h2) and 
	    wt1 = wt *. (q1_coeff dt h2) and 
	    wt2 = wt *. (q2_coeff dt h2) in 
	predict_position b2 t q_pred;
	DFloatBase.grad_v m q2 m2 q_pred gv;
	for i = 0 to 2 do
	  dVdq2.(i) <- dVdq2.(i) +. wt *. gv.(i);
	  dVdq02.(i) <- dVdq02.(i) -. wt0 *. gv.(i);
	  dVdq12.(i) <- dVdq12.(i) -. wt1 *. gv.(i);
	  dVdq22.(i) <- dVdq22.(i) -. wt2 *. gv.(i);
	done
      end
  done

let fill_whole_array (a : diff array) obj = 
  Array.fill a 0 (Array.length a) obj

let h_const h b = 
  b.tau <- h

let h_adaptive sf ({m = m; q = q; q2 = q2; f = f; h = h} as b) = 
  let r = DFloatBase.distance q q2 in
  let tau = if r = (C 0.0) then 
    h*.((C 2.0)**(C 0.2))
  else
    let f = DFloatBase.norm f in 
    h*.(sf*.h*.h*.f/.
	(m*.r))**(C 0.2) in 
  b.tau <- tau

let advance_body b sf = 
  let {h = h; tau = tau; m = m; q = q; p = p; q2 = q2; dVdq0 = dVdq0; dVdq1 = dVdq1; dVdq2 = dVdq2; f = f; df = df; ddf = ddf} = b in 
  for i = 0 to 2 do 
    q.(i) <- q.(i) +. p.(i)/.m*.h -. (h/.((C 2.0)*.m))*.((C 2.0)*.dVdq0.(i) +. dVdq1.(i));
    p.(i) <- p.(i) -. (dVdq0.(i) +. dVdq1.(i) +. dVdq2.(i))
  done;
  b.t <- b.t +. h;
  for i = 0 to 2 do 
    f.(i) <- (C (-6.0))*.dVdq2.(i)/.h;
    df.(i) <- (C (-6.0))*.(dVdq0.(i) -. dVdq1.(i) +. (C 3.0)*.dVdq2.(i))/.(DFloatBase.square h);
    ddf.(i) <- (C (-12.0))/.(DFloatBase.cube h)*.((C 2.0)*.dVdq0.(i) -. dVdq1.(i) +. (C 2.0)*.dVdq2.(i))
  done;
  h_adaptive sf b;
  fill_whole_array dVdq0 (C 0.0);
  fill_whole_array dVdq1 (C 0.0);
  fill_whole_array dVdq2 (C 0.0);
  incr DFloatBase.nsteps

let can_advance b dt = 
  dt < b.tau

let tau_lt {tau = h1} {tau = h2} = 
  if h1 < h2 then 
    -1
  else if h1 > h2 then
    1
  else
    0

let sort_range a cmp low high = 
  let temp = Array.sub a low high in 
  Array.fast_sort cmp temp;
  for i = low to Pervasives.(-) high 1 do 
    a.(i) <- temp.(Pervasives.(-) i low)
  done

let rec advance_bodies bs i dt sf = 
  let b = bs.(i) in 
  if can_advance b dt then 
    begin
      b.h <- dt;
      hermite_predict b dt;
      compute_dVs bs i;
      if Pervasives.(>) i 0 then 
        advance_bodies bs (Pervasives.(-) i 1) dt sf;
      advance_body b sf
    end
  else
    begin
      advance_bodies bs i (dt /. (C 2.0)) sf;
      sort_range bs tau_lt 0 (Pervasives.(+) i 1);
      advance_bodies bs i (dt /. (C 2.0)) sf
    end
      
let set_up_bs bs sf = 
  let gv = Array.make 3 (C 0.0) and 
      dgv = Array.make 3 (C 0.0) in 
  if bs.(0).tau < (C 0.0) then 
    Array.iteri 
      (fun i b -> 
	let q = b.q and m = b.m and p = b.p and 
	    f = b.f and df = b.df in 
	Array.iteri 
	  (fun j b2 -> 
	    let q2 = b2.q and m2 = b2.m and p2 = b2.p in 
	    if Pervasives.(=) i j then 
	      ()
	    else
	      begin
		DFloatBase.grad_v_dgrad_v m q p m2 q2 p2 gv dgv;
		for i = 0 to 2 do
		  f.(i) <- f.(i) -. gv.(i);
		  df.(i) <- df.(i) -. dgv.(i)
		done
	      end)
	  bs;
	b.tau <- (C 6.0)*.sf*.(DFloatBase.norm f)/.(DFloatBase.norm df))
      bs;
  Array.iter (fun b -> hermite_predict b b.tau) bs;
  Array.sort tau_lt bs

let sorted_advance bs dt sf = 
  set_up_bs bs sf;
  advance_bodies bs (Pervasives.(-) (Array.length bs) 1) dt sf;
  bs

let advance bs dt sf = 
  let bs = Array.map body_copy bs in 
  let new_bs = sorted_advance bs dt sf in 
  Array.fast_sort (fun {id = id1} {id = id2} -> Pervasives.compare id1 id2) new_bs;
  new_bs
    
