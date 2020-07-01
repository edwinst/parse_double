# Double-precision floating-point parser

This is a parser converting a decimal string representation of a number to
IEEE 754 double-precision floating-point.

Some properties of the code:

* The parser proper is roughly 400 lines of which some could still be eliminated.

* It does not do perfect unlimited-precision parsing. Instead the parser uses the
  maximum amount of significant mantissa digits that fit in a 64-bit integer and
  then it ignores the remaining digits. This leads to a worst-case rounding error
  that is slightly larger than that of implementations considering all digits.
  (See "Rounding" below.)

* The code has no dependencies except for the availability of a
  (64-bit)*(64-bit)->(128-bit) unsigned integer multiplication
  (in particular no floating-point functions are required).
  Note: Actually only the high 64 bits of the result are used. Therefore a
  64x64->64 "mulhi" operation as it is available on some cores can also be used.

* Faster than Visual C++ 2019 `atof` by a factor of 2.6-5x according to my crude measurements
  (also the standard deviation of the cycle times is much less than for Visual C++ atof).
  (Note: This is not a "fair" comparison in the sense that `atof` has slightly lower
   worst-case rounding error than this implementation; see "Rounding")

* Round-trip-safe in the sense that binary->string->binary (using Visual C++
  sprintf for the first step) is the identity if the sprintf precision is large
  enough; even for subnormal numbers.

* Round-trip-safe in the sense that string->binary->string (using Visual C++
  sprintf to generate the string and for the final step) is the identity in all cases I tested so far.

* The parser is written to run off a datastream with only single-byte-lookahead which makes
  it a bit more awkward than it would need to be for in-memory parsing. That would be easy to change.

The file parse_double.cpp contains the parser, round-trip testing and simple performance measurements.

histogram.* is only used to visualize the performance measurements. I included them
in case anyone wants to reproduce my results.

I put the code in the public domain (see "Unlicense" text in the files).

The code is certainly not polished and there are points I'd gladly discuss if anyone is interested in using it.

## Rounding

Because this parser considers only up to 20 decimal digits in the mantissa part,
it does not achieve the theoretically optimal rounding error in all cases.
More specifically, there is a range of numbers starting (exclusively) at each midway
point between any two adjacent exactly representable double-precision numbers and extending
towards larger absolute values in which this parser rounds towards zero while the
nearest representable number actually lies towards -/+ infinity of the exact result.

The included test program demonstrates such a case when invoked with the command line argument `--rounding-demo`.
It produces the following output at the time of this writing:

    number_a = 1.0000000000000000000000000000000000000000000000000000
    number_b = 1.0000000000000002220446049250313080847263336181640625
    number_a          : result = 1 (0x3ff0000000000000) atof_result = 1 (0x3ff0000000000000)
    number_b          : result = 1 (0x3ff0000000000001) atof_result = 1 (0x3ff0000000000001)
    halfpoint         : result = 1 (0x3ff0000000000000) atof_result = 1 (0x3ff0000000000000)
    
    1.00000000000000011102230246251565404236316680908203125: result = 1 (0x3ff0000000000000) atof_result = 1 (0x3ff0000000000000)
     our rounding error: -0.00000000000000011102230246251565404236316680908203125
    atof rounding error: -0.00000000000000011102230246251565404236316680908203125
    1.00000000000000011102230246251565404236316680908203126: result = 1 (0x3ff0000000000000) atof_result = 1 (0x3ff0000000000001)
     our rounding error: -0.00000000000000011102230246251565404236316680908203126
    atof rounding error: +0.00000000000000011102230246251565404236316680908203124
    1.00000000000000011102230246251565404236316680908203127: result = 1 (0x3ff0000000000000) atof_result = 1 (0x3ff0000000000001)
     our rounding error: -0.00000000000000011102230246251565404236316680908203127
    atof rounding error: +0.00000000000000011102230246251565404236316680908203123
    ...
    1.00000000000000011109999999999999999999999999999999998: result = 1 (0x3ff0000000000000) atof_result = 1 (0x3ff0000000000001)
     our rounding error: -0.00000000000000011109999999999999999999999999999999998
    atof rounding error: +0.00000000000000011094460492503130808472633361816406252
    1.00000000000000011109999999999999999999999999999999999: result = 1 (0x3ff0000000000000) atof_result = 1 (0x3ff0000000000001)
     our rounding error: -0.00000000000000011109999999999999999999999999999999999
    atof rounding error: +0.00000000000000011094460492503130808472633361816406251
    1.00000000000000011110000000000000000000000000000000000: result = 1 (0x3ff0000000000001) atof_result = 1 (0x3ff0000000000001)
     our rounding error: +0.00000000000000011094460492503130808472633361816406250
    atof rounding error: +0.00000000000000011094460492503130808472633361816406250
    
    largest seen absolute rounding errors:
        atof result: 0.00000000000000011102230246251565404236316680908203125
         our result: 0.00000000000000011109999999999999999999999999999999999 (0.0700% larger)


## Measurements

EXAMPLE measurements with scientific notation and low precision (short strings to parse):

    our  cycles> mean =       188.2
    our  cycles>  std =       806.6 (428.7%); n = 16'777'216
    our  cycles>  min =        56
    our  cycles>  max = 1'172'156
    our  cycles> q{01 25; 50; 75 99} = {149.5 149.5; 149.5; 149.5 749.5}
    our  cycles> [    0;        99]   0%   0%    32756|
    our  cycles> [  100;       199]  89%  88% 14819939|**********************************************************************
    our  cycles> [  200;       299]  97%   9%  1443126|*******
    our  cycles> [  300;       399]  98%   1%   158079|*
    our  cycles> [  400;       499]  98%   0%    66358|
    our  cycles> [  500;       599]  99%   0%    23682|
    our  cycles> [  600;       699]  99%   0%    34859|
    our  cycles> [  700;       799]  99%   0%    37860|
    our  cycles> [  800;       899]  99%   0%    61705|
    our  cycles> [  900;       999] 100%   0%    52744|
    our  cycles> [1'000;     1'099] 100%   0%    18761|
    our  cycles> [1'100;     1'199] 100%   0%     9250|
    our  cycles> [1'200;     1'299] 100%   0%     5183|
    our  cycles> [1'300;     1'399] 100%   0%     4518|
    our  cycles> [1'400;     1'499] 100%   0%     2138|
    our  cycles> [1'500;     1'599] 100%   0%      648|
    our  cycles> [1'600;     1'699] 100%   0%      291|
    our  cycles> [1'700;     1'799] 100%   0%      141|
    our  cycles> [1'800;     1'899] 100%   0%      105|
    our  cycles> [1'900;     1'999] 100%   0%       56|
    our  cycles> [2'000;     2'099] 100%   0%       62|
    our  cycles> [2'100;     2'199] 100%   0%       39|
    our  cycles> [2'200;     2'299] 100%   0%       23|
    our  cycles> [2'300;     2'399] 100%   0%       20|
    our  cycles> [2'400;     2'499] 100%   0%       11|
    our  cycles> [2'500;     2'599] 100%   0%       10|
    our  cycles> [2'600;     2'699] 100%   0%        4|
    our  cycles> [2'700;     2'799] 100%   0%        6|
    our  cycles> [2'800;     2'899] 100%   0%        3|
    our  cycles> [2'900;     2'999] 100%   0%        5|
    our  cycles> [3'000; 1'172'156] 100%   0%     4834>


    atof cycles> mean =       801.2
    atof cycles>  std =     3'059.1 (381.8%); n = 16'777'216
    atof cycles>  min =        84
    atof cycles>  max = 7'747'891
    atof cycles> q{01 25; 50; 75 99} = {249.5 549.5; 649.5; 949.5 nan  }
    atof cycles> [    0;        99]   0%   0%    4077|
    atof cycles> [  100;       199]   0%   0%    7538|
    atof cycles> [  200;       299]   2%   2%  263122|******
    atof cycles> [  300;       399]   8%   6%  992262|***********************
    atof cycles> [  400;       499]  23%  15% 2594007|*************************************************************
    atof cycles> [  500;       599]  41%  18% 2974849|**********************************************************************
    atof cycles> [  600;       699]  56%  15% 2504519|***********************************************************
    atof cycles> [  700;       799]  63%   7% 1159288|***************************
    atof cycles> [  800;       899]  68%   6%  954478|**********************
    atof cycles> [  900;       999]  76%   8% 1281171|******************************
    atof cycles> [1'000;     1'099]  84%   8% 1393911|*********************************
    atof cycles> [1'100;     1'199]  91%   7% 1178507|****************************
    atof cycles> [1'200;     1'299]  96%   5%  827418|*******************
    atof cycles> [1'300;     1'399]  97%   1%  187799|****
    atof cycles> [1'400;     1'499]  98%   0%   36099|*
    atof cycles> [1'500;     1'599]  98%   0%   25898|*
    atof cycles> [1'600;     1'699]  98%   0%   22783|*
    atof cycles> [1'700;     1'799]  98%   0%   20393|
    atof cycles> [1'800;     1'899]  98%   0%   19317|
    atof cycles> [1'900;     1'999]  98%   0%   18060|
    atof cycles> [2'000;     2'099]  98%   0%   16378|
    atof cycles> [2'100;     2'199]  98%   0%   15477|
    atof cycles> [2'200;     2'299]  98%   0%   14986|
    atof cycles> [2'300;     2'399]  99%   0%   14636|
    atof cycles> [2'400;     2'499]  99%   0%   15256|
    atof cycles> [2'500;     2'599]  99%   0%   16391|
    atof cycles> [2'600;     2'699]  99%   0%   15783|
    atof cycles> [2'700;     2'799]  99%   0%   13269|
    atof cycles> [2'800;     2'899]  99%   0%   11026|
    atof cycles> [2'900;     2'999]  99%   0%   10380|
    atof cycles> [3'000; 7'747'891] 100%   1% 168138>****
    Completed 16777216 tests, 0 failed.
    OK

EXAMPLE MEASUREMENTS without scientific notation and with very large precision (long strings) and bit-wise round-trip:

    our  cycles> mean =      2'171.3
    our  cycles>  std =     15'798.3 (727.6%); n = 16'777'216
    our  cycles>  min =         53
    our  cycles>  max = 41'892'541
    our  cycles> q{01 25; 50; 75 99} = {1'349.5 1'749.5; 1'949.5; 2'249.5 nan  }
    our  cycles> [    0;         99]   0%   0%    7747|
    our  cycles> [  100;        199]   0%   0%     182|
    our  cycles> [  200;        299]   0%   0%       8|
    our  cycles> [  300;        399]   0%   0%     178|
    our  cycles> [  400;        499]   0%   0%      51|
    our  cycles> [  500;        599]   0%   0%      18|
    our  cycles> [  600;        699]   0%   0%       4|
    our  cycles> [  700;        799]   0%   0%       4|
    our  cycles> [  800;        899]   0%   0%       0|
    our  cycles> [  900;        999]   0%   0%       0|
    our  cycles> [1'000;      1'099]   0%   0%       0|
    our  cycles> [1'100;      1'199]   0%   0%       0|
    our  cycles> [1'200;      1'299]   0%   0%   23072|*
    our  cycles> [1'300;      1'399]   3%   3%  471484|*******************
    our  cycles> [1'400;      1'499]   8%   5%  866938|***********************************
    our  cycles> [1'500;      1'599]  15%   7% 1211070|*************************************************
    our  cycles> [1'600;      1'699]  24%   8% 1422804|**********************************************************
    our  cycles> [1'700;      1'799]  34%  10% 1633439|******************************************************************
    our  cycles> [1'800;      1'899]  44%  10% 1672347|********************************************************************
    our  cycles> [1'900;      1'999]  53%  10% 1646542|*******************************************************************
    our  cycles> [2'000;      2'099]  64%  10% 1715046|*********************************************************************
    our  cycles> [2'100;      2'199]  74%  10% 1731136|**********************************************************************
    our  cycles> [2'200;      2'299]  83%   9% 1593589|****************************************************************
    our  cycles> [2'300;      2'399]  90%   6% 1054436|*******************************************
    our  cycles> [2'400;      2'499]  94%   4%  651631|**************************
    our  cycles> [2'500;      2'599]  94%   1%   96209|****
    our  cycles> [2'600;      2'699]  94%   0%   51133|**
    our  cycles> [2'700;      2'799]  95%   0%   53063|**
    our  cycles> [2'800;      2'899]  95%   0%   65873|***
    our  cycles> [2'900;      2'999]  96%   0%   61981|***
    our  cycles> [3'000; 41'892'541] 100%   4%  747231>******************************


    atof cycles> mean =      5'594.8
    atof cycles>  std =     12'478.6 (223.0%); n = 16'777'216
    atof cycles>  min =         81
    atof cycles>  max = 18'882'139
    atof cycles> q{01 25; 50; 75 99} = {2'449.5 nan  ; nan ; nan   nan  }
    atof cycles> [    0;         99]   0%   0%     3445|
    atof cycles> [  100;        199]   0%   0%     4471|
    atof cycles> [  200;        299]   0%   0%        5|
    atof cycles> [  300;        399]   0%   0%       18|
    atof cycles> [  400;        499]   0%   0%      120|
    atof cycles> [  500;        599]   0%   0%       89|
    atof cycles> [  600;        699]   0%   0%       30|
    atof cycles> [  700;        799]   0%   0%       10|
    atof cycles> [  800;        899]   0%   0%        3|
    atof cycles> [  900;        999]   0%   0%        1|
    atof cycles> [1'000;      1'099]   0%   0%        0|
    atof cycles> [1'100;      1'199]   0%   0%        1|
    atof cycles> [1'200;      1'299]   0%   0%        0|
    atof cycles> [1'300;      1'399]   0%   0%        0|
    atof cycles> [1'400;      1'499]   0%   0%        0|
    atof cycles> [1'500;      1'599]   0%   0%        0|
    atof cycles> [1'600;      1'699]   0%   0%        0|
    atof cycles> [1'700;      1'799]   0%   0%        0|
    atof cycles> [1'800;      1'899]   0%   0%        0|
    atof cycles> [1'900;      1'999]   0%   0%        2|
    atof cycles> [2'000;      2'099]   0%   0%       54|
    atof cycles> [2'100;      2'199]   0%   0%     9088|
    atof cycles> [2'200;      2'299]   0%   0%    38937|
    atof cycles> [2'300;      2'399]   1%   1%    93124|
    atof cycles> [2'400;      2'499]   2%   1%   147312|*
    atof cycles> [2'500;      2'599]   3%   1%   192214|*
    atof cycles> [2'600;      2'699]   4%   1%   197156|*
    atof cycles> [2'700;      2'799]   5%   1%   203257|*
    atof cycles> [2'800;      2'899]   7%   1%   223368|*
    atof cycles> [2'900;      2'999]   8%   2%   261199|*
    atof cycles> [3'000; 18'882'139] 100%  92% 15403312>**********************************************************************
    Completed 16777216 tests, 0 failed.
    OK  
