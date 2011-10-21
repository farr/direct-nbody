(** Various leapfrog integrators. *)

open Base

type b = {
  id : int;
  mutable t : float;
  mutable dt_max : float;
  m : float;
  q : float array;
  v : float array
}

type id = int

let id b = b.id
let t b = b.t
let m b = b.m
let q b = b.q
let p b = Array.map (fun x -> x *. b.m) b.v

let gen_id = 
  let id = ref 0 in 
    fun () -> 
      incr id;
      !id

let make t m q p = 
  {id = gen_id();
   t = 0.0;
   dt_max = 0.0;
   m = m;
   q = Array.copy q;
   v = Array.map (fun x -> x /. m) p}

let copy b = 
  {id = b.id;
   t = b.t;
   dt_max = b.dt_max;
   m = b.m;
   q = Array.copy b.q;
   v = Array.copy b.v}

let print chan {id = id; m = m; t = t; q = q; v = v; dt_max = dtm} = 
  Printf.fprintf chan "--- !Particle\n";
  Printf.fprintf chan "id: %d\n" id;
  Printf.fprintf chan "t: %g\n" t;
  Printf.fprintf chan "r:\n";
  Printf.fprintf chan "  - %g\n  - %g\n  - %g\n" q.(0) q.(1) q.(2);
  Printf.fprintf chan "v:\n";
  Printf.fprintf chan "  - %g\n  - %g\n  - %g\n" v.(0) v.(1) v.(2);
  Printf.fprintf chan "m: %g\n" m;
  Printf.fprintf chan "t_max: %g\n" (t +. dtm)

let read chan = 
  let id = Scanf.bscanf chan " --- !Particle id: %d " (fun x -> x) in 
  let t = Scanf.bscanf chan " t: %g " (fun t -> t) in
  let q = Scanf.bscanf chan " r: - %g - %g - %g " (fun x y z -> [|x; y; z|]) in
  let v = Scanf.bscanf chan " v: - %g - %g - %g " (fun vx vy vz -> [|vx; vy; vz|]) in
  let m = Scanf.bscanf chan " m: %g " (fun m -> m) in 
  let tm = Scanf.bscanf chan " t_max: %g " (fun dtm -> dtm) in
    {id = id; m = m; t = t; q = q; v = v; dt_max = (tm -. t)}

let drift dt ({t = t; q = q; v = v} as b) = 
  b.t <- t +. dt;
  q.(0) <- q.(0) +. dt*.v.(0);
  q.(1) <- q.(1) +. dt*.v.(1);
  q.(2) <- q.(2) +. dt*.v.(2)

let kick 
    eta dt 
    ({m = m1; q = q1; v = v1} as b1)
    ({m = m2; q = q2; v = v2} as b2) =
      let dx = q2.(0) -. q1.(0) and 
          dy = q2.(1) -. q1.(1) and 
          dz = q2.(2) -. q1.(2) in 
      let dvx = v2.(0) -. v1.(0) and 
          dvy = v2.(1) -. v1.(1) and 
          dvz = v2.(2) -. v1.(2) in
      let r2 = dx*.dx +. dy*.dy +. dz*.dz and
          vsq = dvx*.dvx +. dvy*.dvy +. dvz*.dvz and
          vdr = dx*.dvx +. dy*.dvy +. dz*.dvz in
      let r = sqrt r2 in 
      let r3 = r2*.r in 
      let c = dt /. r3 in 
      let c1 = m2 *. c and 
          c2 = m1 *. c in
        v1.(0) <- v1.(0) +. c1*.dx;
        v1.(1) <- v1.(1) +. c1*.dy;
        v1.(2) <- v1.(2) +. c1*.dz;
        v2.(0) <- v2.(0) -. c2*.dx;
        v2.(1) <- v2.(1) -. c2*.dy;
        v2.(2) <- v2.(2) -. c2*.dz;
        let mtot = m1 +. m2 in 
        let tff = eta *. (sqrt (r3 /. mtot)) and 
            tc = eta *. (sqrt (r2 /. vsq)) in 
        let scaled_vdr = vdr /. r2 in 
        let dtff = 1.5 *. scaled_vdr *. tff and
            dtc = scaled_vdr *. tc *. (1.0 +. mtot /. (vsq *. r)) in
        let tff = abs_float (tff /. (1.0 -. 0.5 *. dtff)) and
            tc = abs_float (tc /. (1.0 -. 0.5 *. dtc)) in
        let tscale = if tff < tc then tff else tc in
          if b1.dt_max > tscale then b1.dt_max <- tscale;
          if b2.dt_max > tscale then b2.dt_max <- tscale

let sort_sub bs iend = 
  let subbs = Array.sub bs 0 iend in 
    Array.fast_sort (fun b1 b2 -> compare b1.dt_max b2.dt_max) subbs;
    for i = 0 to iend - 1 do 
      bs.(i) <- subbs.(i)
    done

let rec find_can_advance_index bs iend dt eta = 
  if iend = 0 || bs.(iend-1).dt_max < dt then 
    iend
  else 
    find_can_advance_index bs (iend-1) dt eta

let rec advance_subsystem dt eta bs iend = 
  if iend <= 0 then 
    ()
  else begin
    let ibegin = find_can_advance_index bs iend dt eta and 
        dt2 = dt /. 2.0 in 
      for i = ibegin to iend - 1 do 
        bs.(i).dt_max <- infinity
      done;
      advance_subsystem dt2 eta bs ibegin;
      for i = ibegin to iend - 1 do
        drift dt2 bs.(i)
      done;
      for i = ibegin to iend - 1 do 
        for j = 0 to i - 1 do 
          kick eta dt bs.(i) bs.(j)
        done
      done;
      for i = ibegin to iend - 1 do 
        drift dt2 bs.(i)
      done;
      advance_subsystem dt2 eta bs ibegin;
      sort_sub bs iend
  end

let advance bs dt eta = 
  let bs = Array.map copy bs in 
    (* "Kick" the entire system with dt = 0.0 to set up the correct
       free-fall and collision timescales. *)
    for i = 0 to Array.length bs - 1 do 
      bs.(i).dt_max <- infinity
    done;
    for i = 0 to Array.length bs - 1 do 
      for j = i+1 to Array.length bs - 1 do 
        kick eta 0.0 bs.(i) bs.(j)
      done
    done;
    Array.fast_sort (fun b1 b2 -> compare b1.dt_max b2.dt_max) bs;
    advance_subsystem dt eta bs (Array.length bs);
    bs

let advance_const_dt bs dt n = 
  let bs = Array.map copy bs in
    for i = 0 to n - 1 do 
      advance_subsystem dt infinity bs (Array.length bs)
    done;
    bs
