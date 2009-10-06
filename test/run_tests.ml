open OUnit

let tests = "all tests" >:::
  ["ic.ml tests" >: Ic_test.tests]

let _ = 
  let results = run_test_tt_main tests in 
  let nfail = 
    List.fold_left
      (fun nfail res -> 
        match res with 
        | RSuccess(_) -> nfail
        | _ -> nfail + 1)
      0
      results in 
  exit nfail
