(*  base.ml: Basic gravity definitions.
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

let nsteps = ref 0;;

let eps2 = ref 0.0;;

let square x = x*.x;;

let cube x = x*.x*.x;;

let set_eps eps = 
  eps2 := square eps;;

let norm_squared v = 
  let sum = ref 0.0 in 
  for i = 0 to 2 do 
    sum := !sum +. square v.(i)
  done;
  !sum +. 0.0;;

let norm v = sqrt (norm_squared v);;

let distance_squared x y = 
  let sum = ref 0.0 in 
  for i = 0 to 2 do 
    sum := !sum +. square (x.(i) -. y.(i))
  done;
  !sum +. 0.0;;

let distance x y = sqrt (distance_squared x y);;

let dot x y = 
  let sum = ref 0.0 in 
  for i = 0 to 2 do 
    sum := !sum +. x.(i)*.y.(i)
  done;
  !sum +. 0.0;;

let v m1 q1 m2 q2 = 
  let r2 = distance_squared q1 q2 in 
  let r = sqrt(r2 +. !eps2) in 
  ~-.m1*.m2/.r

let grad_v m1 q1 m2 q2 gv = 
  let r2 = distance_squared q1 q2 in 
  let re = r2 +. !eps2 in 
  let r3 = re*.(sqrt re) in 
  let factor = m1*.m2/.r3 in 
  for i = 0 to 2 do 
    gv.(i) <- factor*.(q1.(i) -. q2.(i))
  done;;

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
  let factord = 3.0*.rdv*.m1*.m2/.r5 in 
  for i = 0 to 2 do 
    dgv.(i) <- factorg*.dgv.(i) -. factord*.gv.(i);
    gv.(i) <- factorg*.gv.(i)
  done;;

let relative_error x t = 
  abs_float ((x -. t) /. t);;
