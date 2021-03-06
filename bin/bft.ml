(*  Formation time of the first hard binary.
    Copyright (C) 2009-2012 Will M. Farr <wmfarr@gmail.com>

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

(** bft: Binary formation time. 

    When run, prints to stdout the time (in n-body units) of the first
    formation of a hard binary in the system.  Prints various messages
    to stderr during the evolution. 

*)

module B = Advancer.A
module Ic = Ics.Make(B)
module A = Analysis.Make(B)
module E = A.E
module Ad = Advancer.A
module I = Integrator.Make(B)(Ad)

let n = ref 25
let nrep = ref 100
let dt = ref 0.1
let tmax = ref 250.0
let errmax = ref 1e-2
let ktfac = ref 5.0
let sf = ref 2e-3

let _ = Random.self_init ()

let _ = 
  Arg.parse 
    [("-n", Arg.Set_int n, 
      (Printf.sprintf "number of bodies to evolve (default %d)" !n));
     ("-nrep", Arg.Set_int nrep, 
      (Printf.sprintf "number of separate evolutions (default %d)" !nrep));
     ("-dt", Arg.Set_float dt, 
      (Printf.sprintf "time interval to check for hard binaries (default %g)" !dt));
     ("-tmax", Arg.Set_float tmax, 
      (Printf.sprintf "maximum evolution time (default %g)" !tmax));
     ("-errmax", Arg.Set_float errmax, 
      (Printf.sprintf "maximum relative energy error (default %g)" !errmax));
     ("-ktfac", Arg.Set_float ktfac, 
      (Printf.sprintf "number of kT for 'hard' binary (default %g)" !ktfac));
     ("-sf", Arg.Set_float sf, 
      (Printf.sprintf "timestep safety factor (default %g)" !sf));
     ("-seed", Arg.Int (fun s -> Random.init s), "RNG seed (default self_init)")]
    (fun _ -> ())
    "bft [OPTIONS ...]"

(* Arguments are (time, temp). *)
exception Hard_binary of float * float

let _ = 
  Base.set_eps 0.0;
  for i = 0 to !nrep - 1 do 
    Printf.eprintf "Beginning rep %d\n" i;
    let bs0 = Ic.rescale_to_standard_units (Ic.make_plummer !n) in 
    let e0 = E.energy bs0 in 
    try
      let new_bs = 
        I.constant_sf_integrate
          ~ci:!dt
          ~max_eg_err:!errmax
          ~sf:!sf
          ~filter:(fun bs -> 
            let temp = A.body_temperature bs in 
            let (b1,b2) = A.tightest_binary bs in 
            let tight_e = A.binary_energy b1 b2 in 
            if (abs_float tight_e) > (abs_float 5.0*.temp) then 
              raise (Hard_binary ((B.t bs.(0)), tight_e))
            else
              bs)
          bs0 !tmax in 
      (* Integrated to tmax *)
      let e1 = E.energy new_bs in 
      Printf.eprintf "Integrated system to %g, rel error = %g, no hard binary.\n"
        (B.t new_bs.(0)) (abs_float ((e0 -. e1)/.e0));
      Printf.printf "%g\n%!" infinity
    with 
    | Hard_binary(time, _) -> 
        Printf.eprintf "Found binary at time %g\n" time;
        Printf.printf "%g\n%!" time
    | I.Eg_err(bs) -> 
        Printf.eprintf "Aborting due to energy error\n";
        Printf.printf "%g\n%!" nan
    | _ -> 
        Printf.eprintf "Aborting computation due unspecified error\n";
        Printf.printf "%g\n%!" nan
  done
