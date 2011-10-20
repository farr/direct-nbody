(** Various leapfrog integrators. *)

open Base

type b = {
  id : int;
  mutable t : float;
  mutable tff2_min : float;
  mutable tc2_min : float;
  mutable dt_max : float;
  m : float;
  q : float array;
  v : float array
}

type id = int

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
   tff2_min = infinity;
   tc2_min = infinity;
   dt_max = 0.0;
   m = m;
   q = Array.copy q;
   v = Array.map (fun x -> x /. m) p}

let copy b = 
  {id = b.id;
   t = b.t;
   tff2_min = b.tff2_min;
   tc2_min = b.tc2_min;
   dt_max = b.dt_max;
   m = b.m;
   q = Array.copy b.q;
   v = Array.copy b.v}

let print chan {id = id; m = m; t = t; q = q; v = v} = 
  Printf.fprintf chan "--- !Particle\n";
  Printf.fprintf chan "id: %d\n" id;
  Printf.fprintf chan "t: %g\n" t;
  Printf.fprintf chan "r:\n";
  Printf.fprintf chan "  - %g\n  - %g\n  - %g\n" q.(0) q.(1) q.(2);
  Printf.fprintf chan "v:\n";
  Printf.fprintf chan "  - %g\n  - %g\n  - %g\n" v.(0) v.(1) v.(2);
  Printf.fprintf chan "m: %g\n" m

let read chan = 
  let id = Scanf.bscanf chan " --- !Particle id: %d " (fun x -> x) in 
  let t = Scanf.bscanf chan " t: %g " (fun t -> t) in
  let q = Scanf.bscanf chan " r: - %g - %g - %g " (fun x y z -> [|x; y; z|]) in
  let v = Scanf.bscanf chan " v: - %g - %g - %g " (fun vx vy vz -> [|vx; vy; vz|]) in
  let m = Scanf.bscanf chan " m: %g " (fun m -> m) in
    {id = id; m = m; t = t; q = q; v = v; tff2_min = infinity; tc2_min = infinity; dt_max = 0.0}

let drift dt ({t = t; q = q; v = v} as b) = 
  b.t <- t +. dt;
  q.(0) <- q.(0) +. dt*.v.(0);
  q.(1) <- q.(1) +. dt*.v.(1);
  q.(2) <- q.(2) +. dt*.v.(2)

let kick 
    dt 
    ({m = m1; q = q1; v = v1; tff2_min = tff2_min1; tc2_min = tc2_min1} as b1)
    ({m = m2; q = q2; v = v2; tff2_min = tff2_min2; tc2_min = tc2_min2} as b2) =
      let dx = q2.(0) -. q1.(0) and 
          dy = q2.(1) -. q1.(1) and 
          dz = q2.(2) -. q2.(2) in 
      let dvx = v2.(0) -. v1.(0) and 
          dvy = v2.(1) -. v1.(1) and 
          dvz = v2.(2) -. v1.(2) in
      let r2 = dx*.dx +. dy*.dy +. dz*.dz and
          vsq = dvx*.dvx +. dvy*.dvy +. dvz*.dvz and
          vdr = dx*.dvx +. dy*.dvy +. dz*.dvz in
      let r = sqrt r2 in 
      let r3 = r2*.r in 
        v1.(0) <- v1.(0) +. m2*.dt*.dx /. r3;
        v1.(1) <- v1.(1) +. m2*.dt*.dy /. r3;
        v1.(2) <- v1.(2) +. m2*.dt*.dz /. r3;
        v2.(0) <- v2.(0) -. m1*.dt*.dx /. r3;
        v2.(1) <- v2.(1) -. m1*.dt*.dy /. r3;
        v2.(2) <- v2.(2) -. m1*.dt*.dz /. r3;
        let tff2 = r3 /. (m1 +. m2) and 
            tc2 = r2 /. vsq in 
        let dtff2 = 3.0 *. vdr *. tff2 /. r2 and 
            dtc2 = 2.0 *. vdr *. tc2 *. (1.0 +. (m1 +. m2) /. (vsq *. r)) /. r2 in
        let tff2 = tff2 /. (1.0 -. 0.5 *. dtff2) and 
            tc2 = tc2 /. (1.0 -. 0.5 *. dtc2) in 
          if tff2 < tff2_min1 then b1.tff2_min <- tff2;
          if tc2 < tc2_min1 then b1.tc2_min <- tc2;
          if tff2 < tff2_min2 then b2.tff2_min <- tff2;
          if tc2 < tc2_min2 then b2.tc2_min <- tc2

let set_dt b = 
  b.dt_max <- sqrt (min b.tff2_min b.tc2_min);
  b.tff2_min <- infinity;
  b.tc2_min <- infinity

let sort_sub bs iend = 
  let subbs = Array.sub bs 0 iend in 
    Array.fast_sort (fun b1 b2 -> compare b1.dt_max b2.dt_max) subbs;
    for i = 0 to iend - 1 do 
      bs.(i) <- subbs.(i)
    done

let rec find_can_advance_index bs iend dt = 
  if iend = 0 || bs.(iend-1).dt_max < dt then 
    iend
  else 
    find_can_advance_index bs (iend-1) dt

let rec advance_subsystem dt eta bs iend = 
  if iend <= 0 then 
    ()
  else begin
    let ibegin = find_can_advance_index bs iend dt and 
        dt2 = dt /. 2.0 in 
      advance_subsystem dt2 eta bs ibegin;
      for i = ibegin to iend - 1 do
        drift dt2 bs.(i)
      done;
      for i = ibegin to iend - 1 do 
        for j = 0 to ibegin - 1 do 
          kick dt bs.(i) bs.(j)
        done
      done;
      for i = ibegin to iend - 1 do 
        drift dt2 bs.(i)
      done;
      advance_subsystem dt2 eta bs ibegin;
      for i = ibegin to iend - 1 do 
        set_dt bs.(i)
      done;
      sort_sub bs iend
  end

let advance bs dt eta = 
  let bs = Array.map copy bs in 
    (* "Kick" the entire system with dt = 0.0 to set up the correct
       free-fall and collision timescales. *)
    for i = 0 to Array.length bs - 1 do 
      for j = i+1 to Array.length bs - 1 do 
        kick 0.0 bs.(i) bs.(j)
      done
    done;
    Array.iter set_dt bs;
    advance_subsystem dt eta bs (Array.length bs)
