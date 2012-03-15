open Ocamlbuild_plugin

let oUnit_dir = "/Users/farr/Documents/code/oUnit"
let gsl_dir = "+site-lib/gsl"

let _ = dispatch begin function 
  | After_rules -> 
    ocaml_lib ~extern:true ~dir:oUnit_dir "oUnit";
    ocaml_lib ~extern:true ~dir:gsl_dir "gsl";
    flag ["native"; "compile"]
      (S[A"-inline"; A"1000"]);
    flag ["native"; "link"]
      (S[A"-inline"; A"1000"]);
    ocaml_lib ~extern:false "nbody"
  | _ -> ()
end
