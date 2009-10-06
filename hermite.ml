(* hermite.ml: Hermite code advancer (used for comparison with
   advancer.ml code).

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

type b = {
    id : int;
    t : float;
    m : float;
    q : float array;
    p : float array;
    h : float;
    f : float array;
    df : float array
  }

let new_id = 
  let counter = ref 0 in 
  fun () -> 
    incr counter;
    !counter

let t {t = t} = t
let m {m = m} = m
let q {q = q} = q
let p {p = p} = p

let make t m q p = 
  {id = new_id (); t = t; m = m; q = Array.copy q; 
   p = Array.copy p; h = (-1.0); f = Array.make 3 0.0; 
   df = Array.make 3 0.0}

let copy b = 
  {b with q = Array.copy b.q;
   p = Array.copy b.p;
   f = Array.copy b.f;
   df = Array.copy b.df}

let predict_position_into qp 
    {t = t0; m = m; q = q; p = p; f = f; df = df} t = 
  let dt = t -. t0 in 
  let pc = dt /. m in
  let fc = pc *. dt *. 0.5 in
  let dfc = fc *. dt /. 3.0 in 
  for i = 0 to 2 do 
    qp.(i) <- q.(i) +. pc*.p.(i) +. fc*.f.(i) +. dfc*.df.(i)
  done;
  qp
    
let predict_momentum_into pp
    {t = t0; p = p; f = f; df = df} t = 
  let dt = t -. t0 in 
  let fc = dt in
  let dfc = fc *. dt *. 0.5 in 
  for i = 0 to 2 do 
    pp.(i) <- p.(i) +. fc*.f.(i) +. dfc*.df.(i)
  done;
  pp

let next_time {t = t; h = h} = 
  t +. h

(* Returns 0 on NaN; otherwise a total order on the next time for
  potential evaluations of these bodies. *)
let compare_next_time b1 b2 =
  let nt1 = next_time b1 and
      nt2 = next_time b2 in 
  if nt1 < nt2 then -1 else if nt1 > nt2 then 1 else 0

let compare_id {id = id1} {id = id2} = 
  if id1 < id2 then -1 else if id1 > id2 then 1 else 0

let timescale eta m q qp f h =
  let r = Base.distance q qp and 
      f = Base.norm f in 
  if r = 0.0 then begin
    h*.(2.0**0.25)
  end else
    h*.(eta*.(Base.square h)*.f/.(m*.r))**0.25

let advance_body eta 
    ({t = t; m = m; q = q; p = p; f = f; df = df; h = h} as b) bs = 
  let nt = next_time b in
  let qp = predict_position_into (Array.make 3 0.0) b nt and 
      pp = predict_momentum_into (Array.make 3 0.0) b nt in 
  let fp = Array.make 3 0.0 and 
      dfp = Array.make 3 0.0 and 
      gv = Array.make 3 0.0 and 
      dgv = Array.make 3 0.0 and
      qp2 = Array.make 3 0.0 and 
      pp2 = Array.make 3 0.0 in 
  List.iter 
    (fun ({m = m2} as b2) -> 
      let qp2 = predict_position_into qp2 b2 nt and 
	  pp2 = predict_momentum_into pp2 b2 nt in 
      Base.grad_v_dgrad_v m qp pp m2 qp2 pp2 gv dgv;
      for i = 0 to 2 do 
	fp.(i) <- fp.(i) -. gv.(i);
	dfp.(i) <- dfp.(i) -. dgv.(i)
      done)
    bs;
  let p_new = Array.make 3 0.0 and 
      q_new = Array.make 3 0.0 in 
  for i = 0 to 2 do 
    p_new.(i) <- p.(i) +. 0.5*.h*.(f.(i) +. fp.(i)) +.
	(Base.square h)/.12.0*.(df.(i) -. dfp.(i));
    q_new.(i) <- q.(i) +. 0.5*.h/.m*.(p.(i) +. p_new.(i)) +.
	(Base.square h)/.(12.0*.m)*.(f.(i) -. fp.(i))
  done;
  incr Base.nsteps;
  {b with t = t +. h;
   q = q_new;
   p = p_new;
   f = fp;
   df = dfp;
   h = timescale eta m q_new qp fp h}

let finish_bs eta stop_time bs = 
  let rec loop done_bs = function 
    | [] -> 
	let bs = Array.of_list done_bs in 
	Array.fast_sort compare_id bs;
	bs
    | b :: bs -> 
	loop
	  ((advance_body eta {b with h = stop_time -. b.t} (bs @ done_bs)) 
	   :: done_bs)
	  bs in 
  loop [] bs

let start_timescale eta {m = m; f = f; df = df} = 
  6.0*.eta*.(Base.norm f)/.(m*.(Base.norm df))

let start_bs eta bs = 
  let bs = List.map 
      (fun b -> {b with f = Array.make 3 0.0; df = Array.make 3 0.0})
      bs and
      gv = Array.make 3 0.0 and 
      dgv = Array.make 3 0.0 in 
  let rec loop done_bs = function 
    | [] -> done_bs 
    | ({m = m; q = q; p = p; f = f; df = df} as b) :: bs -> 
	List.iter 
	  (fun {m = m2; q = q2; p = p2; f = f2; df = df2} -> 
	    Base.grad_v_dgrad_v m q p m2 q2 p2 gv dgv;
	    for i = 0 to 2 do 
	      f.(i) <- f.(i) -. gv.(i);
	      f2.(i) <- f2.(i) +. gv.(i);
	      df.(i) <- df.(i) -. dgv.(i);
	      df2.(i) <- df2.(i) +. dgv.(i)
	    done)
	  bs;
	loop (b :: done_bs) bs in 
  ignore(loop [] bs);
  List.map (fun b -> {b with h = start_timescale eta b}) bs

let rec advance_bodies eta stop_time bs = 
  if next_time (List.hd bs) > stop_time then 
    finish_bs eta stop_time bs
  else
    let b = List.hd bs and 
	bs = List.tl bs in 
    let new_b = advance_body eta b bs in 
    advance_bodies eta stop_time (List.merge compare_next_time [new_b] bs)

let advance bs dt sf = 
  let bs = Array.to_list bs in 
  let bs = 
    if (List.hd bs).h < 0.0 then 
      start_bs sf bs
    else 
      bs in 
  let bs = List.fast_sort compare_next_time bs in 
  let stop_time = (List.hd bs).t +. dt in 
  advance_bodies sf stop_time bs
