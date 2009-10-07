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
let errmax = ref 1e-6
let ktfac = ref 5.0

let _ = 
  Arg.parse 
    [("-n", Arg.Set_int n, "number of bodies to evolve");
     ("-nrep", Arg.Set_int nrep, "number of separate evolutions");
     ("-dt", Arg.Set_float dt, "time interval to check for hard binaries");
     ("-tmax", Arg.Set_float tmax, "maximum evolution time");
     ("-errmax", Arg.Set_float errmax, "maximum relative energy error");
     ("-ktfac", Arg.Set_float ktfac, "number of kT for 'hard' binary")]
    (fun _ -> ())
    "bft [OPTIONS ...]"

(* Arguments are (time, temp). *)
exception Hard_binary of float * float

let _ = Random.self_init ()

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
          ~sf:1e-8
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
      Printf.printf "%g\n" infinity
    with 
    | Hard_binary(time, _) -> 
        Printf.eprintf "Found binary at time %g\n" time;
        Printf.printf "%g\n" time
    | I.Eg_err(bs) -> 
        Printf.eprintf "Aborting due to energy error";
        Printf.printf "%g\n" nan
    | _ -> 
        Printf.eprintf "Aborting computation due unspecified error"
  done
