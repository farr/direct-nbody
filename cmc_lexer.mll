{ 
  open Cmc_parser
}

let eol = "\r\n" | '\n' | '\r'
let digit = ['0' - '9']
let unsigned_num = digit+
let int_num = unsigned_num | '+' unsigned_num | '-' unsigned_num
let real_base = int_num '.' unsigned_num | int_num '.'
let real_expt = int_num
let real_num = "nan" | "inf" | "-inf" | real_base | (real_base | int_num) ('e' | 'E' | 'd' | 'D') real_expt
let whitespace = [' ' '\t']
let begin_comment = '#'

rule tokenize = parse 
  | begin_comment { parse_comment lexbuf }
  | int_num as num { FLOAT (float_of_string num) }
  | real_num as num { FLOAT (float_of_string num) }
  | eol { EOL }
  | eof { EOF }
  | whitespace { tokenize lexbuf }
and 
  parse_comment = parse 
    | eol { tokenize lexbuf }
    | _ { parse_comment lexbuf }
