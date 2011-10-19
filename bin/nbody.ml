(** A program to integrate an N-body system forward in time, dumping
    diagnostics and system states. *)

let summary_file = ref "summary.dat"
let state_file = ref "states.dat"

let n = ref 1024
let eps = ref 0.0

let dt = ref 1.0

let eemax = ref 1e-3

let eta = ref 1e-4

let mmin = ref 0.1
let mmax = ref 100.0
let alpha = ref (-2.35)

let nden = ref 6

let opts = 
  Arg.align
    [("-n", Arg.Set_int n, "N number of bodies");
     ("-eps", Arg.Set_float eps, "EPS softening parameter");
     ("-summary-file", Arg.Set_string summary_file, "FILE summary output filename");
     ("-state-file", Arg.Set_string state_file, "FILE state snapshot filename");
     ("-dt", Arg.Set_float dt, "DT timestep between outputs");
     ("-max-de", Arg.Set_float eemax, "DE maximum energy error");
     ("-eta", Arg.Set_float eta, "ETA timestep safety factor");
     ("-seed", Arg.Int (fun x -> Random.init x), "SEED RNG seed");
     ("-mmin", Arg.Set_float mmin, "MMIN minimum stellar mass");
     ("-mmax", Arg.Set_float mmax, "MMAX maximum stellar mass");
     ("-alpha", Arg.Set_float alpha, "ALPHA power-law slope of the stellar mass function");
     ("-nden", Arg.Set_int nden, "N number of neighbors to average over for density estimators")]

module A = Advancer.A
module Ic = Ics.Make(A)
module An = Analysis.Make(A)
module E = An.E 

let draw_power_law alpha min max = 
  let big = if alpha < 0.0 then min**alpha else max**alpha in 
    Ic.random_from_dist ~xmin:min ~xmax:max ~ymin:0.0 ~ymax:big (fun x -> x**alpha)

let summarize bs = 
    let sum = open_out_gen [Open_append; Open_creat; Open_text] 0o644 !summary_file and 
        state = open_out_gen [Open_append; Open_creat; Open_binary] 0o644 !state_file in 
      try 
        let t = An.time bs and 
            twall = Sys.time () and
            rc = An.density_radius !nden bs and 
            rhoc = An.density_density !nden bs and 
            lrs = Array.init 9 (fun i -> An.lagrange_radius bs ((float_of_int i)/.10.0 +. 0.1)) and 
            e = E.energy bs in
        let de = abs_float ((e +. 0.25) /. 0.25) in
          Printf.fprintf sum "%g %g %g %g %g %g " t rc rhoc e de twall;
          for i = 0 to 7 do 
            Printf.fprintf sum "%g " lrs.(i)
          done;
          Printf.fprintf sum "%g\n" lrs.(8);
          close_out sum;
          Marshal.to_channel state bs [];
          close_out state;
          bs
      with 
        | x -> 
          close_out sum;
          close_out state;
          raise x
      
    

let _ = 
  Random.self_init ();
  Arg.parse opts (fun _ -> ()) "nbody.native OPT ...";
  let bs0 = 
    Ic.rescale_to_standard_units 
      (Ic.add_mass_spectrum (fun () -> draw_power_law !alpha !mmin !mmax) 
         (Ic.adjust_frame 
            (Ic.make_hot_spherical !n))) in 
  let rec loop bs = 
    let new_bs = A.advance bs !dt !eta in 
    let e = E.energy new_bs in 
    let de = abs_float ((e +. 0.25) /. 0.25) in 
      if de > !eemax then begin
        Printf.fprintf stderr "Maximum energy error exceeded at time %g, exiting.\n" (An.time new_bs);
        ignore(summarize new_bs);
        exit 1
      end else 
        loop (summarize new_bs) in 
    loop (summarize bs0)
                                    
  
