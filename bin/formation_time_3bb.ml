(*  Calculate formation time of first hard binary in CMC potential.
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
module An = Analysis.Make(A)
module Ic = Ics.Make(A)
module E = An.E

let outdir = ref ""
let in_snapshot = ref ""
let bin_out_threshold = ref 0.1
let threshold = ref 5.0
let de_max = ref 0.01
let eta = ref 1e-3
let dt = ref 0.1
let nnbr = ref 6
let dtstop = ref infinity

let options = 
  [("-outpath", Arg.Set_string outdir, "DIR directory for output");
   ("-in", Arg.Set_string in_snapshot, "FILE input CMC snapshot file");
   ("-threshold", Arg.Set_float threshold, "E threshold for 'hard' binary (in units of kT)");
   ("-de", Arg.Set_float de_max, "DE max allowed fractional energy error");
   ("-eta", Arg.Set_float eta, "ETA accuracy parameter");
   ("-dt", Arg.Set_float dt, "DT checkpoint timestep (N-body time units)");
   ("-dtstop", Arg.Set_float dtstop, "DT (units are same as snapshot time)");
   ("-binthresh", Arg.Set_float bin_out_threshold, "E threshold for binary output (in units of kT)");
   ("-nnbr", Arg.Set_int nnbr, "N number of neighbors used to estimate density")]

let dump_body out b = 
  let m = b.A.m in 
    Printf.fprintf out "%g %g %g %g %g %g %g %g\n"
      b.A.t m 
      b.A.q.(0) b.A.q.(1) b.A.q.(2) 
      (b.A.p.(0)/.m) (b.A.p.(1)/.m) (b.A.p.(2)/.m)

let dump_snapshot out bs = 
  Array.iter (dump_body out) bs

let dump_and_get_binaries kt out bs = 
  let binaries = An.binaries_threshold (!bin_out_threshold *. kt) bs in 
    Array.iter (fun (b1, b2, e) -> Printf.fprintf out "%g %d %d\n" (e /. kt) b1.A.id b2.A.id) binaries;
    binaries

let tightest_binary_energy binaries = 
  if Array.length binaries = 0 then 
    0.0
  else
    Array.fold_left (fun eg (_,_,e) -> if e < eg then e else eg) infinity binaries

exception Tight_binary of A.b array * float
exception Energy_error of A.b array
exception Time_elapsed of A.b array

let tscale ms rs = rs**1.5 /. (sqrt ms)

let filter e0 energy rc_file mscale rscale bs = 
  let t = bs.(0).A.t in 
  let kt = An.body_temperature bs in 
  let e = energy bs in 
    Printf.eprintf "T = %g; E = %g; dE/E = %g\n%!" t e (abs_float ((e-.e0)/.e0));
    if abs_float ((e -. e0) /. e0) > !de_max then 
      raise (Energy_error bs)
    else begin
      let snap = open_out (Printf.sprintf "%s/snapshot-%g.dat" !outdir t) and 
          bin = open_out (Printf.sprintf "%s/binaries-%g.dat" !outdir t) and
          rho2s = An.density_squared_estimators !nnbr bs in 
      let rc = An.density_radius rho2s bs and 
          r10 = An.lagrange_radius bs 0.1 and 
          r50 = An.lagrange_radius bs 0.5 and 
          r90 = An.lagrange_radius bs 0.9 in 
        Printf.fprintf rc_file "%g %g %g %g %g %g %g %g %g %g\n%!" 
          t (t*.(tscale mscale rscale)) 
          rc (rc*.rscale)
          r10 (r10*.rscale)
          r50 (r50*.rscale)
          r90 (r90*.rscale);
        dump_snapshot snap bs;
        let binaries = dump_and_get_binaries kt bin bs in 
          close_out snap;
          close_out bin;
          let kT_threshold = !threshold *. kt in 
          let e_tight = tightest_binary_energy binaries in 
            if abs_float e_tight > kT_threshold then 
              raise (Tight_binary(bs, e_tight))
            else if t*.(tscale mscale rscale) >= !dtstop then 
              raise (Time_elapsed bs)
            else
              bs
    end

let _ = 
  Base.set_eps 0.0;
  Random.self_init ();
  Arg.parse (Arg.align options) (fun _ -> ()) "3bb_formation_time.native [OPTION ...]";
  if !in_snapshot = "" then begin
    Printf.eprintf "No input snapshot specified!\n";
    exit 1
  end;
  if !outdir = "" then begin
    Printf.eprintf "No output directory specified!\n";
    exit 1
  end;
  let (bs, v, gradV, mscale, rscale) = Ic.make_from_cmc_snapshot !in_snapshot in 
  let energy bs = Array.fold_left (fun e b -> e +. (v b)) (E.energy bs) bs in 
  let e0 = energy bs in 
  let rc_file = open_out (Printf.sprintf "%s/rc.dat" !outdir) in
    Printf.fprintf rc_file "# T_NB T_CMC RC_NB RC_CMC L10_NB L10_CMC L50_NB L50_CMC L90_NB L90_CMC\n";
    try 
      let rec loop bs = 
        loop (A.advance ~extpot:gradV (filter e0 energy rc_file mscale rscale bs) !dt !eta) in 
        loop bs
    with 
      | Energy_error(bs) -> 
        Printf.eprintf "Energy error exceeded limit at time %g\n" bs.(0).A.t;
        exit 1
      | Tight_binary(bs, e) -> 
        let kt = An.body_temperature bs in 
        let out = open_out (Printf.sprintf "%s/tight-binary.dat" !outdir) in 
        let nb_t = bs.(0).A.t in 
        let cmc_t = nb_t *. (tscale mscale rscale) in
          Printf.fprintf out "# T_NB T_CMC (E/kT)\n";
          Printf.fprintf out "%g %g %g\n" nb_t cmc_t (e /. kt);
          close_out out;
          close_out rc_file;
          exit 0
      | Time_elapsed(bs) -> 
        close_out rc_file;
        exit 0
