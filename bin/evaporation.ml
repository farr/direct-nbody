(*  Evaporation of clusters.
    Copyright (C) 2012 Will M. Farr <wmfarr@gmail.com>

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

let n = ref 10
let sf = ref 1e-5
let dt = ref 1.0
let de = ref 1e-3

let options = 
  [("-n", Arg.Set_int n, "N number of bodies");
   ("-sf", Arg.Set_float sf, "SF timestep safety factor");
   ("-dt", Arg.Set_float dt, "DT diagnostic output timestep");
   ("-de", Arg.Set_float de, "DE maximum energy error")]

let _ = 
  Random.self_init ();
  Arg.parse 
    (Arg.align options)
    (fun _ -> ())
    "evaporation.native [OPTIONS ...]";
  let bs = Ic.rescale_to_standard_units (Ic.adjust_frame (Ic.make_plummer !n)) in
    A.follow_flag := true;
    A.nthin := 100;
  let filter bs = 
    let e = E.energy bs in 
    let delta_e = abs_float ((e +. 0.25) /. 0.25) in 
      Printf.eprintf "Time = %g, Energy = %g, de/E = %g\n%!"
        bs.(0).A.t e delta_e;
      if delta_e > !de then begin
        Printf.eprintf "Terminating simulation due to energy error!\n%!";
        exit 1
      end;
      Array.iter (A.print stdout) bs;
      flush stdout;
      bs in 
  let rec loop bs = 
    let new_bs = filter (A.advance bs !dt !sf) in 
      loop new_bs in 
    loop (filter bs)
