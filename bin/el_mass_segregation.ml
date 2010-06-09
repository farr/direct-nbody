open Printf

module Adv = Advancer.A
module A = Analysis.Make(Adv)
module Ic = Ics.Make(Adv)

let n = ref 3000
let dtcheck = ref 1.0
let tevolve = ref 1000.0
let states_file = ref "states.dat"
let cpfile = ref "checkpoint.dat"
let sumfile = ref "summary.dat"
let heavyfrac = ref 0.01
let heavyfac = ref 10.0
let sf = ref 1e-4
let specific = ref true

let options = 
  [("-seed", Arg.Int (fun i -> Random.init i), "set the random seed");
   ("-n", Arg.Set_int n, 
    sprintf "number of bodies (default %d)" !n);
   ("-dtcheck", Arg.Set_float dtcheck,
    sprintf "time between checkpoint and data output (default %g)" !dtcheck);
   ("-tevolve", Arg.Set_float tevolve,
    sprintf "total evolution time (default %g)" !tevolve);
   ("-statesfile", Arg.Set_string states_file,
    sprintf "name of state file (default %s)" !states_file);
   ("-cpfile", Arg.Set_string cpfile,
    sprintf "name of checkpoint file (default %s)" !cpfile);
   ("-sumfile", Arg.Set_string sumfile,
    sprintf "name of summary file (default %s)" !sumfile);
   ("-heavyfrac", Arg.Set_float heavyfrac,
    sprintf "fraction of bodies that are heavy (default %g)" !heavyfrac);
   ("-heavyfac", Arg.Set_float heavyfac,
    sprintf "factor of increase in mass for heavy bodies (default %g)" !heavyfac);
   ("-sf", Arg.Set_float sf,
    sprintf "safety factor for integrator timestep (default %g)" !sf);
   ("-specific", Arg.Set specific,
    sprintf "output E and L in specific units (default %b)" !specific)
  ]

(* Should happen before we parse args *)
let _ = 
  let devr = open_in_bin "/dev/random" in 
  let state = Array.init 55 (fun _ -> input_binary_int devr) in 
    Random.full_init state;
    close_in devr

let _ = Base.set_eps 0.0

let summary_init () = 
  let out = open_out !sumfile in 
    fprintf out "# t E L R0.1 R0.25 R0.5 R0.75 R0.9\n";
    flush out;
    close_out out

let summary bs = 
  let out = open_out_gen [Open_wronly; Open_append; Open_creat; Open_text] 0o644 !sumfile in
    fprintf out "%g %g %g %g %g %g %g %g\n"
      (A.time bs)
      (A.E.energy bs)
      (Base.norm (A.total_angular_momentum bs))
      (A.lagrange_radius bs 0.1)
      (A.lagrange_radius bs 0.25)
      (A.lagrange_radius bs 0.5)
      (A.lagrange_radius bs 0.75)
      (A.lagrange_radius bs 0.9);
    flush out;
    close_out out;
    bs
    
let checkpoint bs = 
  let cp_tmpname = !cpfile ^ ".tmp" in
  let cp_tmp = open_out_bin cp_tmpname in 
    Marshal.to_channel cp_tmp bs [];
    flush cp_tmp;
    close_out cp_tmp;
    Unix.rename cp_tmpname !cpfile;
    bs

let write_states bs = 
  let out = open_out_gen [Open_wronly; Open_append; Open_creat; Open_binary] 0o644 !states_file in
    Marshal.to_channel out bs [];
    flush out;
    close_out out;
    bs

let write_el bs = 
  let el = A.bodies_to_el_phase_space bs in
  let fname = sprintf "el-%g.dat" (A.time bs) in
  let out = open_out fname in 
    Array.iter 
      (fun el -> 
        fprintf out "%g %g %g %g\n" el.(0) el.(1) el.(2) el.(3))
      el;
    close_out out;
    let light_fname = sprintf "el-light-%g.dat" (A.time bs) and 
        heavy_fname = sprintf "el-heavy-%g.dat" (A.time bs) in
    let light = open_out light_fname and 
        heavy = open_out heavy_fname in 
    let mhigh = 
      Array.fold_left (fun mh b -> if mh > Adv.m b then mh else Adv.m b) neg_infinity bs in 
      Array.iteri 
        (fun i el -> 
          if Adv.m bs.(i) = mhigh then 
            fprintf heavy "%g %g %g %g\n" el.(0) el.(1) el.(2) el.(3)
          else
            fprintf light "%g %g %g %g\n" el.(0) el.(1) el.(2) el.(3))
        el;
      close_out light;
      close_out heavy;
      bs
      

let filter bs = 
  checkpoint (write_states (write_el (summary bs)))

let _ = 
  Arg.parse options (fun _ -> ()) "el_mass_segregation.native OPTIONS";
  let bs = 
    Ic.rescale_to_standard_units 
      (Ic.adjust_frame 
         (Ic.add_mass_spectrum (fun _ -> if Random.float 1.0 < !heavyfrac then !heavyfac else 1.0)
            (Ic.make_plummer !n))) in 
  let _ = filter bs in
  let rec advance_loop bs = 
    if A.time bs > !tevolve then 
      ()
    else
      let new_bs = Adv.advance bs !dtcheck !sf in 
        advance_loop (filter new_bs) in 
    advance_loop bs
