(*  Output PSDF tracks for simple systems.
    Copyright (C) 2011-2012 Will M. Farr <wmfarr@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. *)

module A = Advancer.A
module Ic = Ics.Make(A)
module E = Energy.Make(A)

let nbs = ref 16
let dt = ref 0.1
let eta = ref 1e-5
let de = ref 1e-3
let nthin = ref 1

let options = 
  [("-n", Arg.Set_int nbs,
    "N number of bodies");
   ("-nthin", Arg.Set_int nthin,
    "N number of steps between output");
   ("-dt", Arg.Set_float dt,
    "DT timestep for output");
   ("-eta", Arg.Set_float eta,
    "ETA safety factor for timestep adjustment");
   ("-de", Arg.Set_float de,
    "DE maximum relative error in energy tolerated")]

let print_array v = 
  Array.iter (fun x -> Printf.printf "%20.15g " x) v

let _ = 
  Arg.parse 
    (Arg.align options)
    (fun _ -> ())
    "follow_tracks.native OPTIONS ...";
  Random.self_init ();
  A.follow_flag := true;
  A.nthin := !nthin;
  let bs = Ic.rescale_to_standard_units (Ic.adjust_frame (Ic.make_plummer !nbs)) in 
  let bs = Array.mapi (fun i b -> {b with A.id = i}) bs in
  let filter bs = 
    Array.iter (A.print stdout) bs;
    Printf.printf "\n";
    bs in
  let e0 = E.energy bs in
    Printf.eprintf "Initial energy = %g\n%!" e0;
  let rec loop bs = 
    let new_bs = A.advance bs !dt !eta in
    let e = E.energy new_bs in 
      Printf.eprintf "Advanced to %g, e = %g\n%!" (A.t new_bs.(0)) e;
    let delta_e = abs_float ((e -. e0) /. e0) in 
      if delta_e > !de then begin
        Printf.eprintf "Energy error limit exceeded!";
        exit 1
      end else begin
        loop new_bs
      end in
    loop (filter bs)
          
    
