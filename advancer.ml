(*  advancer.ml: Actual integrator code.
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

open Base;;

module type ADVANCER = 
  sig
    type b
    val advance : b array -> float -> float -> b array
  end
      
module type ADWBODY = 
  sig
    type b
    val advance : b array -> float -> float -> b array
    val t : b -> float
    val m : b -> float
    val q : b -> float array
    val p : b -> float array
    val next_t : b -> float
    val make : float -> float -> float array -> float array -> b
    val copy : b -> b
  end
      
module A = 
  struct
    type body = {
	id : int;
	mutable t : float;
	m : float;
	q : float array;
	p : float array;
	mutable tau : float;
	mutable h : float;
	f : float array;
	df : float array;
	ddf : float array;
	q1 : float array;
	q2 : float array;
	dVdq0 : float array;
	dVdq1 : float array;
	dVdq2 : float array;
        tps : float array;
        qps : float array array;
        coeffs : float array array;
      }

    type b = body

    let gen_id =
      let id = ref 0 in 
      fun () -> 
	let res = !id in 
	incr id;
	res

    let t b = b.t
    let m b = b.m
    let q b = b.q
    let p b = b.p
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
       q1 = Array.make 3 0.0;
       q2 = Array.make 3 0.0;
       dVdq0 = Array.make 3 0.0;
       dVdq1 = Array.make 3 0.0;
       dVdq2 = Array.make 3 0.0;
       tps = Array.make 3 nan;
       qps = Array.init 3 (fun _ -> Array.make 3 nan);
       coeffs = Array.init 3 (fun _ -> Array.make 3 nan)
     }

    let make_body t m q p = 
      {id = gen_id ();
       t = t;
       m = m;
       q = q;
       p = p;
       tau = -1.0;
       h = -1.0;
       f = Array.make 3 0.0;
       df = Array.make 3 0.0;
       ddf = Array.make 3 0.0;
       q1 = Array.make 3 0.0;
       q2 = Array.make 3 0.0;
       dVdq0 = Array.make 3 0.0;
       dVdq1 = Array.make 3 0.0;
       dVdq2 = Array.make 3 0.0;
       tps = Array.make 3 0.0;
       qps = Array.init 3 (fun _ -> Array.make 3 nan);
       coeffs = Array.init 3 (fun _ -> Array.make 3 nan)
     }

    let make = make_body
    let copy = body_copy

    module Eg = Energy.Make(struct
      type b = body
      let t b = b.t
      let m b = b.m
      let q b = b.q 
      let p b = b.p
      let make = make
      let copy = copy
    end)

    let kinetic_energy {m = m; p = p} = 
      (Base.norm_squared p)/. (2.0 *. m)

    let potential_energy {m = m1; q = q1} {m = m2; q = q2} = 
      let r = Base.distance q1 q2 in 
      ~-.m1*.m2/.r

    let energy bs = 
      let eg = ref 0.0 in 
      for i = 0 to Array.length bs - 1 do 
	let b = bs.(i) in 
	eg := !eg +. (kinetic_energy b);
	for j = i + 1 to Array.length bs - 1 do 
	  eg := !eg +. (potential_energy b bs.(j))
	done
      done;
      !eg

    let hermite_predict {m = m; q = q; p = p; f = f; df = df; ddf = ddf; q1 = q1; q2 = q2} dt = 
      let dt2 = dt /. 2.0 in 
      for i = 0 to 2 do 
	(* Changed the 24.0 to 12.0 in q1 dt2^4 term.  This changes
	the angular momentum conservation performance to 6th order in
	N = 2 system for a single step *)
	q1.(i) <- q.(i) +. p.(i)/.m*.dt2 +. 
	    f.(i)/.(2.0*.m)*.(Base.square dt2) +. 
	    df.(i)/.(6.0*.m)*.(Base.cube dt2) +. 
	    ddf.(i)/.(12.0*.m)*.(dt2*.(Base.cube dt2));
	q2.(i) <- q.(i) +. p.(i)/.m*.dt +. 
	    f.(i)/.(2.0*.m)*.(Base.square dt) +. 
	    df.(i)/.(6.0*.m)*.(Base.cube dt) +. 
	    ddf.(i)/.(24.0*.m)*.(dt*.(Base.cube dt))
      done

    let q0_coeff dt h h2r = 
      ((h -. (2.0*.dt))*.(h -. dt))*.h2r

    let q1_coeff dt h h2r = 
      (4.0 *. (h -. dt) *. dt)*.h2r

    let q2_coeff dt h h2r = 
      (((2.0*.dt) -. h)*.dt)*.h2r

    let predict_position i {t = tb; h = h; q = q0; q1 = q1; q2 = q2; tps = tps; qps = qps; coeffs = coeffs} t q_pred = 
      if tps.(i) = t then 
        qps.(i)
      else begin
        let dt = t -. tb in
        let h2r = 1.0/.(h*.h) in 
        let q0c = q0_coeff dt h h2r and 
	    q1c = q1_coeff dt h h2r and 
	    q2c = q2_coeff dt h h2r in 
        for i = 0 to 2 do 
	  q_pred.(i) <- q0c*.q0.(i) +. q1c*.q1.(i) +. q2c*.q2.(i)
        done;
        tps.(i) <- t;
        Array.blit q_pred 0 qps.(i) 0 3;
        coeffs.(i).(0) <- q0c;
        coeffs.(i).(1) <- q1c;
        coeffs.(i).(2) <- q2c;
        q_pred
      end

    let compute_dVs bs i = 
      let b = bs.(i) and 
	  gv = Array.make 3 0.0 and 
	  q_pred = Array.make 3 0.0 in 
      let {h = h; t = t; m = m; q = q0; q1 = q1; q2 = q2; dVdq0 = dVdq0; dVdq1 = dVdq1; dVdq2 = dVdq2} = b in 
      let wt = h *. 0.16666666666666666666 in 
      let wtm = 4.0*.wt and 
	  ho2 = h *. 0.5 in
      for j = i + 1 to Array.length bs - 1 do 
	let {t = t2; h = h2; m = m2; dVdq0 = dVdq02; dVdq1 = dVdq12; dVdq2 = dVdq22; tps = tps; coeffs = coeffs} as b2 = bs.(j) in 
	if h = h2 then
	  begin
	    Base.grad_v m q0 m2 b2.q gv;
	    for i = 0 to 2 do 
	      dVdq0.(i) <- dVdq0.(i) +. wt*.gv.(i);
	      dVdq02.(i) <- dVdq02.(i) -. wt*.gv.(i);
	    done;
	    Base.grad_v m q1 m2 b2.q1 gv;
	    for i = 0 to 2 do 
	      dVdq1.(i) <- dVdq1.(i) +. wtm*.gv.(i);
	      dVdq12.(i) <- dVdq12.(i) -. wtm*.gv.(i)
	    done;
	    Base.grad_v m q2 m2 b2.q2 gv;
	    for i = 0 to 2 do 
	      dVdq2.(i) <- dVdq2.(i) +. wt*.gv.(i);
	      dVdq22.(i) <- dVdq22.(i) -. wt*.gv.(i)
	    done
	  end 
	else
	  begin
            let qp = predict_position 0 b2 t q_pred in 
            let c = coeffs.(0) in 
	    let wt0 = wt *. c.(0) and 
		wt1 = wt *. c.(1) and 
		wt2 = wt *. c.(2) in 
	    Base.grad_v m q0 m2 qp gv;
	    for i = 0 to 2 do
	      dVdq0.(i) <- dVdq0.(i) +. wt *. gv.(i);
	      dVdq02.(i) <- dVdq02.(i) -. wt0 *. gv.(i);
	      dVdq12.(i) <- dVdq12.(i) -. wt1 *. gv.(i);
	      dVdq22.(i) <- dVdq22.(i) -. wt2 *. gv.(i);
	    done;
	    let t = t +. ho2 in
            let qp = predict_position 1 b2 t q_pred in 
            let c = coeffs.(1) in 
	    let wt0 = wtm *. c.(0) and 
		wt1 = wtm *. c.(1) and 
		wt2 = wtm *. c.(2) in 
	    Base.grad_v m q1 m2 qp gv;
	    for i = 0 to 2 do
	      dVdq1.(i) <- dVdq1.(i) +. wtm *. gv.(i);
	      dVdq02.(i) <- dVdq02.(i) -. wt0 *. gv.(i);
	      dVdq12.(i) <- dVdq12.(i) -. wt1 *. gv.(i);
	      dVdq22.(i) <- dVdq22.(i) -. wt2 *. gv.(i);
	    done;
	    let t = t +. ho2 in 
	    let qp = predict_position 2 b2 t q_pred in 
            let c = coeffs.(2) in 
	    let wt0 = wt *. c.(0) and 
		wt1 = wt *. c.(1) and 
		wt2 = wt *. c.(2) in 
	    Base.grad_v m q2 m2 qp gv;
	    for i = 0 to 2 do
	      dVdq2.(i) <- dVdq2.(i) +. wt *. gv.(i);
	      dVdq02.(i) <- dVdq02.(i) -. wt0 *. gv.(i);
	      dVdq12.(i) <- dVdq12.(i) -. wt1 *. gv.(i);
	      dVdq22.(i) <- dVdq22.(i) -. wt2 *. gv.(i);
	    done
	  end
      done

    let fill_whole_array (a : float array) obj = 
      Array.fill a 0 (Array.length a) obj

    let h_const h b = 
      b.tau <- h

    let assign_timestep b dt = 
      if dt < 100.0*.epsilon_float*.b.t then 
        raise (Failure "too small timestep in advancer")
      else
        b.tau <- dt

    let h_adaptive sf ({h = h; q2 = q2; q = q; f = f; m = m} as b) = 
      let dq = Base.distance q q2 in (* dq ~ h^5 *) 
      if dq = 0.0 then 
        assign_timestep b (h*.2.01)
      else
        let f_dq = 0.5*.h*.h*.(Base.norm f)/.m in (* f_dq ~ h^2 *) 
        let ratio = f_dq /. dq in 
        assign_timestep b (h*.(sf*.ratio)**(1.0/.5.0))

    let advance_body b sf = 
      let {h = h; tau = tau; m = m; q = q; p = p; q1 = q1; q2 = q2; dVdq0 = dVdq0; dVdq1 = dVdq1; dVdq2 = dVdq2; f = f; df = df; ddf = ddf} = b in 
      for i = 0 to 2 do 
	q.(i) <- q.(i) +. p.(i)/.m*.h -. 
	    (h/.(2.0*.m))*.(2.0*.dVdq0.(i) +. dVdq1.(i));
	p.(i) <- p.(i) -. (dVdq0.(i) +. dVdq1.(i) +. dVdq2.(i))
      done;
      b.t <- b.t +. h;
      for i = 0 to 2 do 
	f.(i) <- (-6.0)*.dVdq2.(i)/.h;
	df.(i) <- (-6.0)*.(dVdq0.(i) -. dVdq1.(i) +. 3.0*.dVdq2.(i))/.
	  (Base.square h);
	ddf.(i) <- (-12.0)/.(Base.cube h)*.
	    (2.0*.dVdq0.(i) -. dVdq1.(i) +. 2.0*.dVdq2.(i))
      done;
      h_adaptive sf b;
      fill_whole_array dVdq0 0.0;
      fill_whole_array dVdq1 0.0;
      fill_whole_array dVdq2 0.0;
      incr Base.nsteps

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
      for i = low to high - 1 do 
	a.(i) <- temp.(i - low)
      done

    let rec advance_bodies bs i dt sf = 
      let b = bs.(i) in 
      if can_advance b dt then begin
        b.h <- dt;
        hermite_predict b dt;
        compute_dVs bs i;
        if i > 0 then 
          advance_bodies bs (i - 1) dt sf;
        advance_body b sf
      end else begin
        advance_bodies bs i (dt /. 2.0) sf;
        sort_range bs tau_lt 0 (i+1);
        advance_bodies bs i (dt /. 2.0) sf
      end
	  
    let set_up_bs bs sf = 
      let gv = Array.make 3 0.0 and 
	  dgv = Array.make 3 0.0 in 
      if bs.(0).tau < 0.0 then 
	Array.iteri 
	  (fun i b -> 
	    let q = b.q and m = b.m and p = b.p and 
		f = b.f and df = b.df in 
	    Array.iteri 
	      (fun j b2 -> 
		let q2 = b2.q and m2 = b2.m and p2 = b2.p in 
		if i = j then 
		  ()
		else
		  begin
		    Base.grad_v_dgrad_v m q p m2 q2 p2 gv dgv;
		    for i = 0 to 2 do
		      f.(i) <- f.(i) -. gv.(i);
		      df.(i) <- df.(i) -. dgv.(i)
		    done
		  end)
	      bs;
	    b.tau <- 6.0*.sf*.(Base.norm f)/.(Base.norm df))
	  bs;
      Array.iter (fun b -> hermite_predict b b.tau) bs;
      Array.sort tau_lt bs

    let sorted_advance bs dt sf = 
      set_up_bs bs sf;
      advance_bodies bs (Array.length bs - 1) dt sf;
      bs

    let advance bs dt sf = 
      let bs = Array.map body_copy bs in 
      let new_bs = sorted_advance bs dt sf in 
      Array.fast_sort (fun {id = id1} {id = id2} -> compare id1 id2) new_bs;
      new_bs
  end
