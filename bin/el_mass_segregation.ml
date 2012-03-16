(*  Follow mass segregation in the E-L phase space.
    Copyright (C) 2010-2012 Will M. Farr <wmfarr@gmail.com>

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
let massspec = ref false
let alpha = ref (-2.3) (* Salpeter IMF. *)
let mmin = ref 0.1
let mmax = ref 10.0
let sf = ref 1e-4

let options = 
  Arg.align
    [("-seed", Arg.Int (fun i -> Random.init i), "seed set the random seed");
     ("-n", Arg.Set_int n, 
      sprintf "n number of bodies (default %d)" !n);
     ("-dtcheck", Arg.Set_float dtcheck,
      sprintf "dt time between checkpoint and data output (default %g)" !dtcheck);
     ("-tevolve", Arg.Set_float tevolve,
      sprintf "t total evolution time (default %g)" !tevolve);
     ("-statesfile", Arg.Set_string states_file,
      sprintf "file name of state file (default %s)" !states_file);
     ("-cpfile", Arg.Set_string cpfile,
      sprintf "file name of checkpoint file (default %s)" !cpfile);
     ("-sumfile", Arg.Set_string sumfile,
      sprintf "file name of summary file (default %s)" !sumfile);
     ("-heavyfrac", Arg.Set_float heavyfrac,
      sprintf "f fraction of bodies that are heavy (default %g)" !heavyfrac);
     ("-heavyfac", Arg.Set_float heavyfac,
      sprintf "f factor of increase in mass for heavy bodies (default %g)" !heavyfac);
     ("-massspec", Arg.Set massspec,
      sprintf " flag for mass spectrum (default %b)" !massspec);
     ("-alpha", Arg.Set_float alpha,
      sprintf "alpha power-law exponent for mass spectrum (default %g)" !alpha);
     ("-mmin", Arg.Set_float mmin,
      sprintf "mmin minimum mass for spectrum (default %g)" !mmin);
     ("-mmax", Arg.Set_float mmax,
      sprintf "mmax maximum mass for spectrum (default %g)" !mmax);
     ("-sf", Arg.Set_float sf,
      sprintf "sf safety factor for integrator timestep (default %g)" !sf)]

(* Should happen before we parse args *)
let _ = 
  let devr = open_in_bin "/dev/urandom" in 
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
    close_out out;
    bs
    
let checkpoint bs = 
  let cp_tmpname = !cpfile ^ ".tmp" in
  let cp_tmp = open_out_bin cp_tmpname in 
    Marshal.to_channel cp_tmp bs [];
    close_out cp_tmp;
    Unix.rename cp_tmpname !cpfile;
    bs

let write_states bs = 
  let out = open_out_gen [Open_wronly; Open_append; Open_creat; Open_binary] 0o644 !states_file in
    Marshal.to_channel out bs [];
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

let draw_mass_cont () = 
  let max = max ((!mmin)**(!alpha)) ((!mmax)**(!alpha)) in 
  let rec loop () = 
    let x = !mmin +. (Random.float (!mmax -. !mmin)) and 
        y = Random.float max in 
      if x**(!alpha) < y then 
        x
      else
        loop () in 
    loop ()

let draw_mass_disc () = 
  if Random.float 1.0 < !heavyfrac then 
    !heavyfac
  else
    1.0

let _ = 
  Arg.parse options (fun _ -> ()) "el_mass_segregation.native OPTIONS";
  let bs = 
    Ic.rescale_to_standard_units 
      (Ic.adjust_frame 
         (Ic.add_mass_spectrum (fun _ -> if !massspec then draw_mass_cont () else draw_mass_disc ())
            (Ic.make_plummer !n))) in 
  let _ = filter bs in
  let rec advance_loop bs = 
    if A.time bs > !tevolve then 
      ()
    else
      let new_bs = Adv.advance bs !dtcheck !sf in 
        advance_loop (filter new_bs) in 
    advance_loop bs
