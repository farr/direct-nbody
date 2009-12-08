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

(** Basic utilities and gravitational computations. *)

let nsteps = ref 0;;
(** Stores the number of timesteps taken. *)

let eps2 = ref 0.0;;
(** Square of the softening parameter. *)

let square x = x*.x;;
(** Utility function. *)

let cube x = x*.x*.x;;
(** Utility function. *)

let set_eps eps = 
  eps2 := square eps;;
(** Allows to set the square of the softening parameter. *)

let norm_squared v = 
  let x = v.(0) and
      y = v.(1) and 
      z = v.(2) in
  x*.x +. y*.y +. z*.z
(** Square of the Euclidean norm in 3-D. *)

let norm v = sqrt (norm_squared v);;
(** Euclidean norm in 3-D. *)

let distance_squared x y = 
  let dx = x.(0) -. y.(0) and 
      dy = x.(1) -. y.(1) and
      dz = x.(2) -. y.(2) in 
  dx*.dx +. dy*.dy +. dz*.dz
(** Euclidean distance squared in 3-D. *)

let distance x y = sqrt (distance_squared x y);;
(** Euclidean distance in 3-D. *)

let dot x y = 
  x.(0)*.y.(0) +. x.(1)*.y.(1) +. x.(2)*.y.(2)
(** Dot product in 3-D. *)

let v m1 q1 m2 q2 = 
  let r2 = distance_squared q1 q2 in 
  let r = sqrt(r2 +. !eps2) in 
  ~-.m1*.m2/.r
(** The gravitational potential (including softening) between two
    bodies: [v m1 q1 m2 q2]. 
    
    @param m1 The mass of body 1.
    @param q1 The position of body 1.
    @param m2 The mass of body 2.
    @param q2 The position of body 2.
 *)

let grad_v m1 q1 m2 q2 gv = 
  let r2 = distance_squared q1 q2 in 
  let re = r2 +. !eps2 in 
  let r3 = re*.(sqrt re) in 
  let factor = m1*.m2/.r3 in 
  gv.(0) <- factor*.(q1.(0) -. q2.(0));
  gv.(1) <- factor*.(q1.(1) -. q2.(1));
  gv.(2) <- factor*.(q1.(2) -. q2.(2))
(** [grad_v m1 q1 m2 q2 gv] stores the gradient of the two body
    potential with respect to [q1] in [gv].

    @param m1 The mass of body 1.
    @param q1 The position of body 1.
    @param m2 The mass of body 2.
    @param q2 The position of body 2.
    @param gv The vector into which the gradient should be stored. 
*)

let grad_v_dgrad_v m1 q1 p1 m2 q2 p2 gv dgv = 
  gv.(0) <- q1.(0) -. q2.(0);
  dgv.(0) <- p1.(0)/.m1 -. p2.(0)/.m2;
  gv.(1) <- q1.(1) -. q2.(1);
  dgv.(1) <- p1.(1)/.m1 -. p2.(1)/.m2;
  gv.(2) <- q1.(2) -. q2.(2);
  dgv.(2) <- p1.(2)/.m1 -. p2.(2)/.m2;
  let r2 = norm_squared gv and 
      rdv = dot gv dgv in 
  let re = r2 +. !eps2 in 
  let r3 = (sqrt re)*.re in 
  let r5 = r3*.re in 
  let factorg = m1*.m2/.r3 in
  let factord = 3.0*.rdv*.m1*.m2/.r5 in 
  dgv.(0) <- factorg*.dgv.(0) -. factord*.gv.(0);
  gv.(0) <- factorg*.gv.(0);
  dgv.(1) <- factorg*.dgv.(1) -. factord*.gv.(1);
  gv.(1) <- factorg*.gv.(1);
  dgv.(2) <- factorg*.dgv.(2) -. factord*.gv.(2);
  gv.(2) <- factorg*.gv.(2)
(** [grad_v_dgrad_v m1 q1 p1 m2 q2 p2 gv dgv] computes both the
    gradient and its time derivative, storing into [gv] and [dgv].
    The arguments are mass, position, and momentum.  See
    {!Base.grad_v}. *)

let relative_error x t = 
  abs_float ((x -. t) /. t);;
(** [relative_error x true] computes the relative error on [x]
    relative to [true]. *)
