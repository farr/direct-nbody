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

module type ADVANCER = sig
  type b
  val advance : b array -> float -> float -> b array
end

module A = 
  struct
    type body = {
        id : int;
        mutable t : float;
        m : float;
        q : float array;
        p : float array;
        mutable h : float; (* Always matches q0, q1, q2. *)
        mutable hmax : float; (* hmax < 0 implies initial step *)
        mutable t0 : float; (* Always matches up with q0,q1,q2. *)
        dVdq0 : float array;
        dVdq1 : float array;
        dVdq2 : float array;
        q0 : float array;
        q1 : float array;
        q2 : float array;
        cs : float array;
        mutable output_counter : int
      }

    type b = body

    type id = int

    let follow_flag = ref false
    let nthin = ref 1

    let gen_id =
      let id = ref 0 in 
      fun () -> 
	let res = !id in 
	incr id;
	res

    let id b = b.id
    let t b = b.t
    let m b = b.m
    let q b = b.q
    let p b = b.p
    let next_t b = b.t +. b.h

    let body_copy b = 
      {id = b.id;
       t = b.t;
       m = b.m;
       q = Array.copy b.q;
       p = Array.copy b.p;
       t0 = b.t0;
       h = b.h;
       hmax = b.hmax;
       q0 = Array.copy b.q0;
       q1 = Array.copy b.q1;
       q2 = Array.copy b.q2;
       dVdq0 = Array.copy b.dVdq0;
       dVdq1 = Array.copy b.dVdq1;
       dVdq2 = Array.copy b.dVdq2;
       cs = Array.copy b.cs;
       output_counter = b.output_counter
     }

    let print chan {id = id; t = t; m = m; q = q; p = p} = 
      Printf.fprintf chan "--- !!Particle\n";
      Printf.fprintf chan "id: %d\n" id;
      Printf.fprintf chan "t: %g\n" t;
      Printf.fprintf chan "r:\n";
      Printf.fprintf chan "  - %g\n  - %g\n  - %g\n" q.(0) q.(1) q.(2);
      Printf.fprintf chan "v:\n";
      Printf.fprintf chan "  - %g\n  - %g\n  - %g\n" (p.(0)/.m) (p.(1)/.m) (p.(2)/.m);
      Printf.fprintf chan "m: %g\n" m

    let make_body_id id t m q p = 
      {id = id;
       t = t;
       m = m;
       q = Array.copy q;
       p = Array.copy p;
       t0 = nan;
       h = nan;
       hmax = nan;
       q0 = Array.make 3 nan;
       q1 = Array.make 3 nan;
       q2 = Array.make 3 nan;
       dVdq0 = Array.make 3 nan;
       dVdq1 = Array.make 3 nan;
       dVdq2 = Array.make 3 nan;
       cs = Array.make 3 nan;
       output_counter = 0
     }

    let read chan = 
      let id = Scanf.bscanf chan " --- !!Particle id: %d " (fun x -> x) in 
      let t = Scanf.bscanf chan " t: %g " (fun t -> t) in
      let q = Scanf.bscanf chan " r: - %g - %g - %g " (fun x y z -> [|x; y; z|]) in
      let v = Scanf.bscanf chan " v: - %g - %g - %g " (fun vx vy vz -> [|vx; vy; vz|]) in
      let m = Scanf.bscanf chan " m: %g " (fun m -> m) in
        make_body_id id t m q (Array.map (fun v -> v *. m) v)

    let make t m q p = make_body_id (gen_id ()) t m q p
    let copy = body_copy

    module Eg = Energy.Make(struct
      type b = body
      type id = int
      let id b = b.id
      let t b = b.t
      let m b = b.m
      let q b = b.q 
      let p b = b.p
      let make = make
      let copy = copy
      let read = read
      let print = print
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

    let predict_q012 ({m = m; q = q; p = p; 
                       dVdq0 = dVdq0; dVdq1 = dVdq1; dVdq2 = dVdq2; 
                       h = h; q0 = q0; q1 = q1; q2 = q2} as b)
        hnew = 
      let hnew2 = hnew*.hnew and 
          h2 = h*.h in 
      let hnew3 = hnew2*.hnew and 
          h3 = h2*.h in 
      for i = 0 to 2 do 
        q0.(i) <- q.(i);
        q1.(i) <- q0.(i) +. hnew*.p.(i)/.(2.0*.m) 
            -. hnew3*.(h+.hnew)/.(8.0*.h3*.m)*.dVdq0.(i) 
            +. hnew3*.(2.0*.h+.hnew)/.(16.0*.h3*.m)*.dVdq1.(i)
            -. hnew2*.(6.0*.h2 +. 3.0*.h*.hnew +. hnew2)/.(8.0*.h3*.m)*.dVdq2.(i);
        q2.(i) <- q0.(i) +. hnew*.p.(i)/.m
            -. hnew3*.(h+.hnew)/.(h3*.m)*.dVdq0.(i)
            +. hnew3*.(2.0*.h+.hnew)/.(2.0*.h3*.m)*.dVdq1.(i)
            -. hnew2*.(3.0*.h2 +. 3.0*.h*.hnew +. hnew2)/.(h3*.m)*.dVdq2.(i)
      done;
      b.h <- hnew;
      b.t0 <- b.t

    let clear_before_step b = 
      Array.fill b.dVdq0 0 3 0.0;
      Array.fill b.dVdq1 0 3 0.0;
      Array.fill b.dVdq2 0 3 0.0;
      Array.fill b.cs 0 3 0.0;
      b.cs.(0) <- 1.0

    let predict_position b t = 
      if b.t = t then 
        ()
      else begin
        let dt = t -. b.t0 and 
            h = b.h in 
        let h2 = h*.h in 
        let cs = b.cs in 
        cs.(0) <- ((2.0*.dt*.dt -. 3.0*.dt*.h +. h2)/.h2);
        cs.(1) <- ((4.0*.dt*.(h -. dt))/.h2);
        cs.(2) <- (dt*.(2.0*.dt-.h)/.h2);
        let c0 = cs.(0) and 
            c1 = cs.(1) and 
            c2 = cs.(2) in 
        let q = b.q and 
            q0 = b.q0 and 
            q1 = b.q1 and 
            q2 = b.q2 in 
        for i = 0 to 2 do 
          q.(i) <- c0*.q0.(i) +. c1*.q1.(i) +. c2*.q2.(i)
        done;
        b.t <- t
      end

    let interact_bodies bs imin imax t vC = 
      let n = Array.length bs in 
      let gv = Array.make 3 0.0 in 
      for i = imin to n - 1 do 
        predict_position bs.(i) t
      done;
      for i = imin to imax - 1 do 
        let b = bs.(i) in 
        let q = b.q and m = b.m and cs = b.cs in 
        let c0 = vC*.cs.(0) and c1 = vC*.cs.(1) and c2 = vC*.cs.(2) in 
        let dVdq0 = b.dVdq0 and dVdq1 = b.dVdq1 and dVdq2 = b.dVdq2 in 
        for j = i+1 to n-1 do 
          let b2 = bs.(j) in 
          let dVdq02 = b2.dVdq0 and dVdq12 = b2.dVdq1 and dVdq22 = b2.dVdq2 in
          let cs2 = b2.cs in 
          let c02 = vC*.cs2.(0) and 
              c12 = vC*.cs2.(1) and 
              c22 = vC*.cs2.(2) in 
          Base.grad_v m q b2.m b2.q gv;
          for i = 0 to 2 do 
            let gvi = gv.(i) in 
            dVdq0.(i) <- dVdq0.(i) +. c0*.gvi;
            dVdq1.(i) <- dVdq1.(i) +. c1*.gvi;
            dVdq2.(i) <- dVdq2.(i) +. c2*.gvi;
            dVdq02.(i) <- dVdq02.(i) -. c02*.gvi;
            dVdq12.(i) <- dVdq12.(i) -. c12*.gvi;
            dVdq22.(i) <- dVdq22.(i) -. c22*.gvi
          done
        done
      done
            
    let finish_body b t = 
      let q = b.q and m = b.m in 
      let p = b.p and q0 = b.q0 and h = b.h in 
      let dVdq0 = b.dVdq0 and dVdq1 = b.dVdq1 and dVdq2 = b.dVdq2 in 
      let cp = h /. m in 
      let cV0 = cp and cV1 = 0.5*.cp in 
      for i = 0 to 2 do 
        let v0i = dVdq0.(i) and 
            v1i = dVdq1.(i) and 
            v2i = dVdq2.(i) and
            pi = p.(i) in
        q.(i) <- q0.(i) +. cp*.pi -. cV0*.v0i -. cV1*.v1i;
        p.(i) <- pi -. v0i -. v1i -. v2i
      done;
      b.t <- t

    let nsmall_tsteps = ref 0
    let print_small_timesteps = ref 1

    let assign_timestep b dt = 
      if dt < 100.0*.epsilon_float*.(abs_float b.t) +. 100.0*.epsilon_float then begin
        incr nsmall_tsteps;
        if !nsmall_tsteps = !print_small_timesteps then begin
          Printf.eprintf "Warning: small timestep in assign_timestep (number %d): %g\n" !nsmall_tsteps dt;
          print_small_timesteps := !print_small_timesteps*2
        end;
        b.hmax <- 100.0*.epsilon_float*.b.t +. 100.0*.epsilon_float
      end else
        b.hmax <- dt

    (* We try below to keep the predictor position error a fixed
       fraction of the local acceleration distance.  This recipe comes
       from Makino, Optimal Order and Time-Step Criterion for
       Aarseth-type N-body Integrators, ApJ 369:200--212, 1991. *)
    let h_adaptive sf ({h = h; q2 = q2; q = q; dVdq2 = dVdq2; m = m} as b) = 
      let dq = Base.distance q q2 in (* dq ~ h^5 *) 
      if dq = 0.0 then 
        assign_timestep b (h*.2.01)
      else begin
        let a = (Base.norm dVdq2)*.6.0/.h/.m in  
        let a_dq = 0.5*.h*.h*.a in (* Should scale as h^2 *)
        let ratio = a_dq /. dq in (* should scale as h^-3 *)
        assign_timestep b (h*.(sf*.ratio)**0.3333333) 
      end

    let hmax_lt {hmax = h1} {hmax = h2} = 
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

    let rec advance_range bs imax dt sf = 
      if bs.(imax-1).hmax < dt then begin
        advance_range bs imax (dt/.2.0) sf;
        advance_range bs imax (dt/.2.0) sf
      end else begin
        let imin = let rec loop i = 
          if i < 0 then 
            0
          else if dt < bs.(i).hmax then 
            loop (i-1)
          else
            i+1 in 
        loop (imax - 1) in 
        for i = imin to imax - 1 do 
          predict_q012 bs.(i) dt;
          clear_before_step bs.(i)
        done;
        let t0 = bs.(imin).t in 
        interact_bodies bs imin imax t0 (dt/.6.0);
        if imin > 0 then advance_range bs imin (dt/.2.0) sf;
        interact_bodies bs imin imax (t0 +. 0.5*.dt) (2.0*.dt/.3.0);
        if imin > 0 then advance_range bs imin (dt/.2.0) sf;
        interact_bodies bs imin imax (t0 +. dt) (dt/.6.0);
        for i = imin to imax - 1 do 
          finish_body bs.(i) (t0+.dt);
          h_adaptive sf bs.(i);
          if !follow_flag then begin
            if bs.(i).output_counter mod !nthin = 0 then 
              print stdout bs.(i);
            bs.(i).output_counter <- bs.(i).output_counter + 1
          end
        done;  
        sort_range bs hmax_lt 0 imax
      end

    (* The idea here is to create a fictional previous step (with
       length h = 1.0) so that the predictor has something to work off
       of for the true first step.  We also try to derive a timescale
       from a comparison between the force and the force derivative at
       the starting time: dt = safety_factor*F/Fdot. *)
    let initialize_before_first_step bs sf = 
      let n = Array.length bs and 
          gv = Array.make 3 0.0 and 
          dgv = Array.make 3 0.0 and 
          fdot = Array.make 3 0.0 and 
          h = 1.0 in 
      let c0 = h /. 6.0 and c1 = 2.0*.h/.3.0 in 
      let c2 = c0 in 
      for i = 0 to n - 1 do 
        clear_before_step bs.(i);
        bs.(i).h <- h;
      done;
      for i = 0 to n - 1 do 
        let {m = m; 
             q = q; 
             p = p;
             dVdq0 = dVdq0; 
             dVdq1 = dVdq1;
             dVdq2 = dVdq2} = bs.(i) in 
        for j = i+1 to n - 1 do 
          let {m = mj;
               q = qj;
               p = pj;
               dVdq0 = dVdq0j;
               dVdq1 = dVdq1j;
               dVdq2 = dVdq2j} = bs.(j) in 
          Base.grad_v_dgrad_v m q p mj qj pj gv dgv;
          for k = 0 to 2 do 
            dVdq0.(k) <- dVdq0.(k) +. c0*.(gv.(k) -. h*.dgv.(k));
            dVdq1.(k) <- dVdq1.(k) +. c1*.(gv.(k) -. 0.5*.h*.dgv.(k));
            dVdq2.(k) <- dVdq2.(k) +. c2*.gv.(k);
            dVdq0j.(k) <- dVdq0j.(k) -. c0*.(gv.(k) -. h*.dgv.(k));
            dVdq1j.(k) <- dVdq1j.(k) -. c1*.(gv.(k) -. 0.5*.h*.dgv.(k));
            dVdq2j.(k) <- dVdq2j.(k) -. c2*.gv.(k)
          done
        done
      done;
      for i = 0 to n - 1 do 
        let {dVdq0 = dVdq0;
             dVdq1 = dVdq1;
             dVdq2 = dVdq2;
             h = h} as b = bs.(i) in 
        for j = 0 to 2 do 
          (* Not exactly fdot, but proportional to it *)
          fdot.(j) <- dVdq0.(j) -. dVdq1.(j) +. 3.0*.dVdq2.(j)
        done;
        assign_timestep b
          (1e-2*.sf**0.33333*.h*.
            (Base.norm dVdq2)/.(Base.norm fdot)) (* Timescale ~ f/fdot. *)
      done

    let vec_any pred v = 
      let n = Array.length v in 
      let rec loop i = 
        if i >= n then 
          false
        else
          (pred v.(i)) || loop (i+1) in 
      loop 0

    let first_step b = 
      match classify_float b.h with 
      | FP_normal -> false
      | _ -> true

    let advance bs dt sf = 
      let bs = Array.map body_copy bs in 
      if vec_any first_step bs then 
        initialize_before_first_step bs sf;
      Array.fast_sort hmax_lt bs;
      advance_range bs (Array.length bs) dt sf;
      bs
  end
