(*  energy.ml: Computes the energy of n-body systems.
    Copyright (C) <year>  <name of author>

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
open Body;;

module type ENERGY = 
  sig
    type b 
    val kinetic_energy : b -> float
    val potential_energy : b -> b -> float
    val total_kinetic_energy : b array -> float
    val total_potential_energy : b array -> float
    val energy : b array -> float
  end

module Make(B : BODY) : (ENERGY with type b = B.b) = 
  struct 
    type b = B.b

    let kinetic_energy b = 
      (norm_squared (B.p b)) /. (2.0 *. (B.m b))

    let potential_energy b1 b2 = 
      v (B.m b1) (B.q b1) (B.m b2) (B.q b2)

    let total_kinetic_energy bs = 
      Array.fold_left (fun e b -> e +. (kinetic_energy b)) 0.0 bs 

    let total_potential_energy bs = 
      let pe = ref 0.0 in 
      for i = 0 to Array.length bs - 1 do 
	for j = i + 1 to Array.length bs - 1 do
	  pe := !pe +. (potential_energy bs.(i) bs.(j))
	done
      done;
      !pe

    let energy bs = 
      let ke = total_kinetic_energy bs in 
      let pe = total_potential_energy bs in 
      ke +. pe
  end
