(*  dFloat_pushforward.ml: Computes the jacobian of an evolution mapping.
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

module type In = sig
  type b
  type db 

  val promote_body : b -> db

  val q : b -> float array
  val p : b -> float array

  val dq : db -> DFloat.diff array
  val dp : db -> DFloat.diff array

  val advance : db array -> DFloat.diff -> DFloat.diff -> db array
end

module Make(I : In) = struct
  open DFloat;;

  type b = I.b
  type db = I.db

  let advance_wrapper dt sf bs qparr = 
    let pbs = Array.map I.promote_body bs in 
    let n = Array.length bs in 
    Array.iteri 
      (fun i b ->
	let q = I.dq b and 
	    p = I.dp b in 
	let qi = Pervasives.( * ) 3 i and 
	    pi = Pervasives.(+) (Pervasives.( * ) 3 n) (Pervasives.( * ) 3 i) in 
	Array.blit qparr qi q 0 3;
	Array.blit qparr pi p 0 3)
      pbs;
    let rbs = I.advance pbs (C dt) (C sf) in
    let resarr = Array.make (Pervasives.( * ) 6 n) (C 0.0) in 
    Array.iteri 
      (fun i b ->
	let q = I.dq b and 
	    p = I.dp b in 
	let qi = Pervasives.( * ) 3 i and 
	    pi = Pervasives.(+) (Pervasives.( * ) 3 n) (Pervasives.( * ) 3 i) in 
	Array.blit q 0 resarr qi 3;
	Array.blit p 0 resarr pi 3)
      rbs;
    resarr

  let jacobian_advance dt sf bs = 
    let n = Array.length bs in 
    let qparr = Array.make (Pervasives.( * ) 6 n) 0.0 in 
    Array.iteri 
      (fun i b ->
	let q = I.q b and 
	    p = I.p b in 
	let qi = Pervasives.( * ) 3 i and 
	    pi = Pervasives.(+) (Pervasives.( * ) n 3) (Pervasives.( * ) i 3) in 
	Array.blit q 0 qparr qi 3;
	Array.blit p 0 qparr pi 3)
      bs;
    lower_jacobian (jacobian (advance_wrapper dt sf bs)) qparr

  let dimensions m = 
    (Array.length m, Array.length m.(0))

  let ( *|| ) m1 m2 = 
    let (n1,k1) = dimensions m1 and 
	(n2,k2) = dimensions m2 in 
    assert(Pervasives.(=) k1 n2);
    let result = Array.make_matrix n1 k2 0.0 in 
    for i = 0 to Pervasives.(-) n1 1 do 
      for k = 0 to Pervasives.(-) k2 1 do 
	let sum = ref 0.0 in 
	for j = 0 to Pervasives.(-) k1 1 do 
	  sum := Pervasives.(+.) !sum (Pervasives.( *. ) m1.(i).(j) m2.(j).(k))
	done;
	result.(i).(k) <- !sum
      done
    done;
    result

  let ( *| ) mat v = 
    let n = Array.length v and 
	m = Array.length mat in 
    let result = Array.make m 0.0 in 
    for i = 0 to Pervasives.(-) m 1 do 
      let mi = mat.(i) and
	  sum = ref 0.0 in 
      for j = 0 to Pervasives.(-) n 1 do 
	sum := Pervasives.(+.) !sum (Pervasives.( *. ) mi.(j) v.(j))
      done;
      result.(i) <- !sum
    done;
    result

  let dot v1 v2 = 
    let n = Array.length v1 in 
    assert (Pervasives.(=) (Array.length v2) n);
    let sum = ref 0.0 in 
    for i = 0 to Pervasives.(-) n 1 do 
      sum := Pervasives.(+.) !sum (Pervasives.( *. ) v1.(i) v2.(i))
    done;
    !sum
      
  let transpose mat = 
    let m = Array.length mat and 
	n = Array.length mat.(0) in 
    let result = Array.make_matrix n m 0.0 in 
    for i = 0 to Pervasives.(-) m 1 do 
      let mi = mat.(i) in 
      for j = 0 to Pervasives.(-) n 1 do 
	result.(j).(i) <- mi.(j)
      done
    done;
    result

  let symplectic n = 
    let result = Array.make_matrix n n 0.0 in 
    let no2 = Pervasives.(/) n 2 in 
    for i = 0 to Pervasives.(-) no2 1 do
      result.(i).(Pervasives.(+) no2 i) <- (-1.0);
      result.(Pervasives.(+) no2 i).(i) <- 1.0
    done;
    result
      
  let omega_pushforward dt sf bs = 
    let j = jacobian_advance dt sf bs in 
    let n = Array.length j in 
    let no2 = Pervasives.(/) n 2 in 
    let jtsj = (transpose j) *|| (symplectic n) *|| j in 
    let sum = ref 0.0 in 
    for i = 0 to Pervasives.(-) no2 1 do 
      sum := Pervasives.(+.) !sum jtsj.(Pervasives.(+) i no2).(i)
    done;
    Pervasives.(/.) !sum (float_of_int no2) 

end
