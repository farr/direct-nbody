(*  omelyan_advancer.ml: Integrator which uses the Omelyan technique for estimating the potential gradient in the Chin and Chen algorithm.
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

module OA = struct 
  type b = {
      m : float;
      q : float array;
      p : float array;
      q_star : float array;
      mutable t : float;
    }

  let t {t = t} = t
  let m {m = m} = m
  let q {q = q} = q
  let p {p = p} = p

  let make t m q p = 
    {m = m; q = Array.copy q; p = Array.copy p; q_star = Array.make 3 0.0;
     t = 0.0}
  let copy {m = m; q = q; p = p; q_star = q_star; t = t} = 
    {m = m; q = Array.copy q; p = Array.copy p; q_star = Array.copy q_star;
     t = t}

  let rec iter_pairwise op = function 
    | [] -> ()
    | [b] -> ()
    | b :: bs -> op b bs; iter_pairwise op bs

  let b_operate h bs = 
    let f = h /. 6.0 and 
	gv = Array.make 3 0.0 in 
    let op {m = m; q = q; p = p} bs = 
      List.iter 
	(fun {m = m2; q = q2; p = p2} -> 
	  Base.grad_v m q m2 q2 gv;
	  for i = 0 to 2 do 
	    p.(i) <- p.(i) -. f*.gv.(i);
	    p2.(i) <- p2.(i) +. f*.gv.(i)
	  done)
	bs in 
    iter_pairwise op bs

  let a_operate h bs = 
    List.iter 
      (fun {m = m; q = q; p = p} -> 
	let f = h /. (2.0*.m) in 
	for i = 0 to 2 do 
	  q.(i) <- q.(i) +. p.(i)*.f
	done)
      bs

  let c_operate h bs = 
    let pf = 2.0*.h/.3.0 and 
	gv = Array.make 3 0.0 in 
    let op_q_star {m = m; q = q; q_star = q_star} bs = 
      let f = (Base.square h)/.(24.0*.m) in 
      List.iter
	(fun {m = m2; q = q2; q_star = q_star2} -> 
	  let f2 = (Base.square h)/.(24.0*.m2) in 
	  Base.grad_v m q m2 q2 gv;
	  for i = 0 to 2 do 
	    q_star.(i) <- q_star.(i) -. f*.gv.(i);
	    q_star2.(i) <- q_star2.(i) +. f2*.gv.(i)
	  done)
	bs in
    let op_p {m = m; q_star = qs; p = p} bs = 
      List.iter
	(fun {m = m2; q_star = qs2; p = p2} -> 
	  Base.grad_v m qs m2 qs2 gv;
	  for i = 0 to 2 do 
	    p.(i) <- p.(i) -. pf*.gv.(i);
	    p2.(i) <- p2.(i) +. pf*.gv.(i)
	  done)
	bs in 
    List.iter 
      (fun {q = q; q_star = q_star} -> Array.blit q 0 q_star 0 3)
      bs;
    iter_pairwise op_q_star bs;
    iter_pairwise op_p bs

  let advance h bs = 
    let bs = Array.to_list (Array.map copy bs) in 
    b_operate h bs;
    a_operate h bs;
    c_operate h bs;
    a_operate h bs;
    b_operate h bs;
    List.iter (fun b -> b.t <- b.t +. h) bs;
    Array.of_list bs
end
