// decimal string to double-precision floating-point number parsing

// This is free and unencumbered software released into the public domain.
//
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
//
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// For more information, please refer to <http://unlicense.org/>


#include <cstdint>
#include <cassert>
#include <intrin.h>
#pragma intrinsic(_umul128)
#pragma intrinsic(__rdtsc)
#ifdef USE_INTRINSIC_BITSCANREVERSE64
#pragma intrinsic(_BitScanReverse64)
#endif

#include <cstdarg>   // for error reporting code

#include <cstdio>    // for test code only
#include <cstring>   // for test code only
#include <cmath>     // for test code only
#include <cinttypes> // for test code only

#include "histogram.h" // for profiling results

////////////////// FRAMEWORK AND ERROR REPORTING ////////////////////////////////////////////////

struct Status
{
    bool failed;
    char error_message[1024];
};

struct DataStream
{
    uint8_t *ptr;
    uint8_t *end;
};

void datastream_fetch(Status *st, DataStream *s)
{
    assert(s->ptr == s->end); // precondition (this function is only called if there is no more input buffered)
    // for a real stream implementation, this function would do its best to provide at least one more input byte
}

void fail_here(Status *st, DataStream *s, char *fmt, ...)
{
    st->failed = true;
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(st->error_message, sizeof(st->error_message), fmt, vl);
    va_end(vl);
}

////////////////////////////////// MAIN CODE ///////////////////////////////////////////////

/**
 * Software floating-point multiplication of two normalized unsigned numbers with
 * full 64-bit mantissae and separate binary exponents.
 */
inline void mulf64(uint64_t *mx, int *ex, uint64_t my, int ey)
{
    // We are multiplying two normalized numbers of the forms:
    //     X = 2^ex + mx * 2^(ex - 64)
    //     Y = 2^ey + my * 2^(ey - 64)
    // The mathematical result is
    //     XY = 2^(ex+ey) + m' * 2^(ex+ey-64)
    //     where m' = mx + my + mx*my*2^(-64)
    // With maximal mantissae (mx, my = 2^64 - 1), we get
    //     m'_max = 2^65 + 2^64 - 2^2 + 2^(-64) >= 2^65 (!)
    //
    // CAUTION: This means that floor(m') does not fit in a uint64_t!
    //
    // We will represent m' as follows:
    //     m' = overflow * 2^64 + m + lo * 2^(-64)
    // where
    //     overflow = floor(m' / (2^64)) is in {0, 1, 2}
    //     m = floor(m') mod 2^64        is in 0..2^64-1
    //     lo = (m' - floor(m')) * 2^64  is in 0..2^64-1
    // substituting for m' we have:
    //     m = (mx + my + hi) mod 2^64
    //     lo = mx*my - hi*2^64 = (mx*my mod 2^64)
    // where
    //     hi = floor(mx*my*2^(-64))
    //
    // The result of the multiplication written in these variables is
    //     XY = 2^(ex+ey) + (overflow * 2^64 + m + lo * 2^(-64)) * 2^(ex+ey-64)
    //        = 2^(ex+ey) * (overflow + 1) + m * 2^(ex+ey-64) + lo * 2^(ex+ey-128)
    //
    // Therefore, if overflow is non-zero we must renormalize the intermediate
    // result from exponent (ex+ey) to exponent (ex+ey+1) to get
    //
    //     XY = 2^(ex+ey+1) * {{floor((overflow + 1)/2)}}
    //        + ( ((overflow + 1) mod 2)*2^63 + floor(m/2) ) * 2^(ex+ey+1-64)
    //        + ( (m mod 2)*2^63 + floor(lo/2) ) * 2^(ex+ey+1-128)
    //        + (lo mod 2)*2^(ex+ey+1-129)
    //
    // Note that the part marked by {{...}} is always equal to 1 if overflow is non-zero.
    //
    int overflow = 0;
    uint64_t hi;
    _umul128(*mx, my, &hi);
    *ex += ey;
    if ((*mx += my) < my) overflow++;
    if ((*mx += hi) < hi) overflow++;
    if (overflow) {
        // renormalize result after overflow
        *ex += 1;
        *mx = (*mx >> 1) | ((uint64_t)(overflow + 1) << 63);
    }
    // XXX @Clarify consider rounding in this function
}

// The following tables were calculated using Mathematica and the command
// Table[ScientificForm[BaseForm[N[10^(s*(2^k)), 20], 2], 65], {s, {1, -1}}, {k, 8, 0, -1}] // Flatten // TableForm

// 10^(2^k) in normalized floating-point format with uint64_t mantissa and integer exponent
// the leading 1. bit is not stored as part of the mantissa
constexpr uint64_t mant1e256  = 0b0101010011111101110101111111011100111011111100111011110100011100; constexpr int exp1e256  = 850;
constexpr uint64_t mant1e128  = 0b0010011101110100100011111001001100000001110100110001100111000000; constexpr int exp1e128  = 425;
constexpr uint64_t mant1e64   = 0b1000010011110000001111101001001111111111100111110100110110101010; constexpr int exp1e64   = 212;
constexpr uint64_t mant1e32   = 0b0011101110001011010110110101000001010110111000010110101100111100; constexpr int exp1e32   = 106;
constexpr uint64_t mant1e16   = 0b0001110000110111100100110111111000001000000000000000000000000000; constexpr int exp1e16   =  53;
constexpr uint64_t mant1e8    = 0b0111110101111000010000000000000000000000000000000000000000000000; constexpr int exp1e8    =  26;
constexpr uint64_t mant1e4    = 0b0011100010000000000000000000000000000000000000000000000000000000; constexpr int exp1e4    =  13;
constexpr uint64_t mant1e2    = 0b1001000000000000000000000000000000000000000000000000000000000000; constexpr int exp1e2    =   6;
constexpr uint64_t mant1e1    = 0b0100000000000000000000000000000000000000000000000000000000000000; constexpr int exp1e1    =   3;

// 10^(-2^k) in normalized floating-point format with uint64_t mantissa and integer exponent
// the leading 1. bit is not stored as part of the mantissa
constexpr uint64_t mant1en256 = 0b1000000001100010100001100100101011000110111101000011001001110100; constexpr int exp1en256 = -851;
constexpr uint64_t mant1en128 = 0b1011101110100000100011001111100011001001011110011100100101000001; constexpr int exp1en128 = -426;
constexpr uint64_t mant1en64  = 0b0101000011111111110101000100111101001010011100111101001101001010; constexpr int exp1en64  = -213;
constexpr uint64_t mant1en32  = 0b1001111101100010001111010101101010001010011100110010100101110101; constexpr int exp1en32  = -107;
constexpr uint64_t mant1en16  = 0b1100110100101011001010010111110110001000100110111100001010110111; constexpr int exp1en16  =  -54;
constexpr uint64_t mant1en8   = 0b0101011110011000111011100010001100001000110000111001110111111010; constexpr int exp1en8   =  -27;
constexpr uint64_t mant1en4   = 0b1010001101101110001011101011000111000100001100101100101001011000; constexpr int exp1en4   =  -14;
constexpr uint64_t mant1en2   = 0b0100011110101110000101000111101011100001010001111010111000010100; constexpr int exp1en2   =   -7;
constexpr uint64_t mant1en1   = 0b1001100110011001100110011001100110011001100110011001100110011010; constexpr int exp1en1   =   -4;

// stupid helper for parsing "inf" and "nan" with our single-byte-lookahead data stream
void datastream_consume_word(Status *st, DataStream *s, char *word)
{
    while (*word) {
        if (s->ptr == s->end) {
            datastream_fetch(st, s);
            if (st->failed)
                return;
            if (s->ptr == s->end) {
                fail_here(st, s, "Unexpected end of stream; expected character '%c'.\n", *word);
                return;
            }
        }
        uint8_t ch = *s->ptr;
        if (ch != *word) {
            char format_buf[] = "'x' ";
            format_buf[1] = ch;
            fail_here(st, s, "Unexpected character (%s0x%02x); expected '%c'.\n",
                      (ch >= 0x20 && ch < 0x7f) ? format_buf : "", ch, *word);
        }
        s->ptr++;
        word++;
    }
}

double parse_double(Status *st, DataStream *s)
{
    if (s->ptr == s->end) {
        datastream_fetch(st, s);
        if (st->failed)
            goto failed;
        if (s->ptr == s->end) {
            fail_here(st, s, "Unexpected end-of-stream; expected decimal floating-point number.\n");
            goto failed;
        }
    }

    uint8_t ch = *s->ptr++;

    bool negative = false;
    uint64_t mantissa = 0;
    int digits = 0;
    int significant_digits = 0;
    // Note: fractional_digits and decimal_exponent should probably be combined into one signed decimal_exponent. See :ExponentAwkwardness
    int64_t fractional_digits = -1; // negative means no decimal point seen, yet
    uint64_t decimal_exponent = 0; // counts only towards positive decimal exponents

    switch (ch) {
        case '-': negative = true; break;
        case '+': /* ignore */   ; break;
        case '0': digits++;      ; break;
        case '1': case '2': case '3': case '4': case '5':
        case '6': case '7': case '8': case '9':
            mantissa = ch - '0';
            digits++;
            significant_digits++;
            break;
        case '.': fractional_digits = 0; break;
        case 'i':
                  s->ptr--; // put back character
                  goto parse_infinity;
        case 'n':
                  s->ptr--; // put back character
                  goto parse_nan;
        default:
                  {
                      char format_buf[] = "'x' ";
                      format_buf[1] = ch;
                      fail_here(st, s, "Unexpected character (%s0x%02x) in numeric literal.\n", (ch >= 0x20 && ch < 0x7f) ? format_buf : "", ch);
                      goto failed;
                  }
    }

    // XXX @Speed maybe split the following loop in a part for fractional_digits < 0 and one for fractional_digits >= 0 (after the decimal point)?
    do {
        while (s->ptr != s->end) {
            ch = *s->ptr;
            uint64_t digit_value = ch - '0';
            if (digit_value <= 9) {
                digits++;
                if (mantissa > (uint64_t)1844674407370955161 - (digit_value >= 6)) {
                    // :MantissaLimit
                    // NOTE: In this case, we cannot represent the 20th significant digit
                    // in the mantissa since we would overflow the uint64_t range.
                    // If we have not yet seen a decimal point (fractional_digits < 0),
                    // we simply let the mantissa be 19 digits and increase the decimal_exponent
                    // instead. If we have already seen a decimal point (fractional_digits >= 0),
                    // we also let the mantissa remain at 19 decimal digits and compensate
                    // for that by *not* incrementing fractional_digits.
                    // Note also that from this point on, the 'digits' value is mostly meaningless.
                    // Setting significant_digits to '20' is a lie in this case
                    // used to trigger the terminating condition below. The
                    // mantissa really remains at 19 decimal digits.
                    //
                    // A uint64_t can represent all 19-digit decimal unsigned integers.
                    // In addition it can represent 20-digit decimal unsigned integers
                    // up to and including 18446744073709551615.
                    assert(significant_digits == 19);
                    significant_digits = 20;
                    if (fractional_digits < 0)
                        decimal_exponent++;
                }
                else {
                    mantissa = 10 * mantissa + digit_value;
                    if (digit_value || significant_digits)
                        significant_digits++;
                    if (fractional_digits >= 0) fractional_digits++;
                }
            }
            else if (ch == '.') {
                if (fractional_digits >= 0)
                    goto end_of_number;
                fractional_digits = 0;
            }
            else if ((ch == 'e' || ch == 'E') && digits) {
parse_exponent:
                s->ptr++;
                if (s->ptr == s->end) {
                    datastream_fetch(st, s);
                    if (st->failed)
                        goto failed;
                    if (s->ptr == s->end) {
                        fail_here(st, s, "Incomplete exponent in decimal floating-point number.");
                        goto failed;
                    }
                }
                ch = *s->ptr++;
                bool negative_exponent = false;
                if (ch == '-')
                    negative_exponent = true;
                else if (ch == '+')
                    ; // ignore
                else
                    s->ptr--; // put back byte
                uint64_t absolute_exponent = 0;
                do {
                    while (s->ptr != s->end) {
                        ch = *s->ptr;
                        digit_value = ch - '0';
                        if (digit_value <= 9) {
                            if (absolute_exponent > (uint64_t)1844674407370955161 - (digit_value >= 6)) {
                                if (negative_exponent || significant_digits == 0)
                                    goto return_possibly_signed_zero;
                                else
                                    goto return_possibly_signed_infinity;
                            }
                            absolute_exponent = 10 * absolute_exponent + digit_value;
                            s->ptr++;
                        }
                        else
                            goto end_of_exponent;
                    }
                    datastream_fetch(st, s);
                    if (st->failed)
                        goto failed;
                } while (s->ptr != s->end);
end_of_exponent:
                // Note: OK, this is awkward. The reason for this mess is that I wrote this
                //       parser for PDF parsing where no explicit exponents occur.
                //       Now that I want to also handle the explicit exponents, the
                //       split handling of fractional_digits and decimal_exponent
                //       bites me. :(
                //       Maybe they should be combined again into a signed decimal_exponent.
                //       :ExponentAwkwardness
                if (negative_exponent && absolute_exponent) {
                    if (decimal_exponent >= absolute_exponent) {
                        decimal_exponent -= absolute_exponent;
                        absolute_exponent = 0;
                    }
                    else {
                        absolute_exponent -= decimal_exponent;
                        decimal_exponent = 0;
                    }
                    if (fractional_digits < 0)
                        fractional_digits = 0;
                    fractional_digits += absolute_exponent;
                    if (fractional_digits < 0 || (uint64_t)fractional_digits < absolute_exponent)
                        goto return_possibly_signed_zero; // underflow
                }
                else {
                    if (fractional_digits > 0) {
                        if ((uint64_t)fractional_digits >= absolute_exponent) {
                            fractional_digits -= (int64_t)absolute_exponent;
                            absolute_exponent = 0;
                        }
                        else {
                            absolute_exponent -= (uint64_t)fractional_digits;
                            fractional_digits = 0;
                        }
                    }
                    decimal_exponent += absolute_exponent;
                    if (decimal_exponent < absolute_exponent)
                        goto return_possibly_signed_infinity; // overflow
                }
                goto end_of_number;
            }
            else if (ch == 'i' && !digits) {
parse_infinity:
                datastream_consume_word(st, s, "inf");
                if (st->failed)
                    goto failed;
return_possibly_signed_infinity:
                uint64_t bits;
                if (negative)
                    bits = 0xFFF0000000000000; // -inf
                else
                    bits = 0x7FF0000000000000; // +inf
                return *reinterpret_cast<double*>(&bits);
            }
            else if (ch == 'n' && !digits) {
parse_nan:
                datastream_consume_word(st, s, "nan");
                if (st->failed)
                    goto failed;
                uint64_t bits;
                if (negative)
                    bits = 0xFFFFFFFFFFFFFFFF; // -nan
                else
                    bits = 0x7FFFFFFFFFFFFFFF; // +nan
                return *reinterpret_cast<double*>(&bits);
            }
            else {
                // we do not understand what follows; consider the number terminated here (s->ptr will point to the first unconsumed byte)
                goto end_of_number;
            }
            s->ptr++;
            if (significant_digits >= 20) // XXX @Cleanup can we come up with a better condition here? (See :MantissaLimit above)
            {
                // skip the remaining digits, counting the remaining
                // decimal places before the decimal point, if any, towards the decimal_exponent.
                do {
                    while (s->ptr != s->end) {
                        ch = *s->ptr;
                        if (ch == '.' && fractional_digits < 0)
                            fractional_digits = 0;
                        else if (ch == 'e' || ch == 'E')
                            goto parse_exponent;
                        else if (ch < '0' || ch > '9')
                            goto end_of_number;
                        if (fractional_digits < 0)
                            decimal_exponent++;
                        s->ptr++;
                    }
                    datastream_fetch(st, s);
                    if (st->failed)
                        goto failed;
                } while(s->ptr != s->end);
                goto end_of_number;
            }
        }
        datastream_fetch(st, s);
        if (st->failed)
            goto failed;
    } while(s->ptr != s->end);

end_of_number:
    // Note: 'digits' is not too meaningful from here on except for the (digits != 0) condition. (See also :MantissaLimit above.)
    assert((fractional_digits <= 0) || (decimal_exponent == 0));

    // XXX @Speed maybe put this check only in the less frequent code paths where it is needed?
    if (digits == 0) {
        fail_here(st, s, "Numeric literal without digits.\n");
        goto failed;
    }

    if (significant_digits == 0) {
return_possibly_signed_zero:
        if (negative) {
            uint64_t bits = 0x8000000000000000;
            return *reinterpret_cast<double*>(&bits);
        }
        return 0.0;
    }

    // printf("mantissa: %" PRIu64 " (0x%016" PRIx64 "), digits: %d significant_digits: %d decimal_exponent: %d fractional_digits: %d\n",
    //        mantissa, mantissa, digits, significant_digits, decimal_exponent, fractional_digits); // XXX DEBUG

    // normalize the mantissa such that: old_mantissa == 2^binexp + new_mantissa * 2^(binexp-64).
    // that is binexp = 63 - (number of leading zeroes in old_mantissa)
#ifdef USE_INTRINSIC_BITSCANREVERSE64
    unsigned long index;
    unsigned char bsr_result = _BitScanReverse64(&index, mantissa);
    assert(bsr_result);
    int binexp = index;
    // XXX @Speed Note: We could probably always just do mantissa <<= (64 - index). Strictly speaking
    //                  this is wrong for index == 0 (i.e. old mantissa == 1) because the shift is a
    //                  no-op then instead of clearing the mantissa. This error of +1 in mantissa
    //                  does not seem to affect the final result, though. I have not seen enough of
    //                  a benefit to do this yet.
    if (index)
        mantissa <<= (64 - index);
    else
        mantissa = 0; // instead of <<= 64.
#else
    int binexp = 63;
    if (!(mantissa & 0xffffffff00000000)) { mantissa <<= 32; binexp -= 32; }
    if (!(mantissa & 0xffff000000000000)) { mantissa <<= 16; binexp -= 16; }
    if (!(mantissa & 0xff00000000000000)) { mantissa <<=  8; binexp -=  8; }
    if (!(mantissa & 0xf000000000000000)) { mantissa <<=  4; binexp -=  4; }
    if (!(mantissa & 0xc000000000000000)) { mantissa <<=  2; binexp -=  2; }
    if (!(mantissa & 0x8000000000000000)) { mantissa <<=  1; binexp -=  1; }
    mantissa <<= 1; // Note: this shifts out the leading '1'-bit
#endif

    // printf("normalized: 0x%016" PRIx64 ", binexp: %d\n", mantissa, binexp); // XXX DEBUG

    // shift the decimal point to the correct position by multiplying with powers of ten
    if (fractional_digits > 0) {
        if (fractional_digits >= 512)
            goto return_possibly_signed_zero; // underflow
        if (fractional_digits & 256) mulf64(&mantissa, &binexp, mant1en256, exp1en256);
        if (fractional_digits & 128) mulf64(&mantissa, &binexp, mant1en128, exp1en128);
        if (fractional_digits &  64) mulf64(&mantissa, &binexp, mant1en64 , exp1en64 );
        if (fractional_digits &  32) mulf64(&mantissa, &binexp, mant1en32 , exp1en32 );
        if (fractional_digits &  16) mulf64(&mantissa, &binexp, mant1en16 , exp1en16 );
        if (fractional_digits &   8) mulf64(&mantissa, &binexp, mant1en8  , exp1en8  );
        if (fractional_digits &   4) mulf64(&mantissa, &binexp, mant1en4  , exp1en4  );
        if (fractional_digits &   2) mulf64(&mantissa, &binexp, mant1en2  , exp1en2  );
        if (fractional_digits &   1) mulf64(&mantissa, &binexp, mant1en1  , exp1en1  );
    }
    else if (decimal_exponent) {
        if (decimal_exponent >= 512)
            goto return_possibly_signed_infinity; // overflow
        if (decimal_exponent & 256) mulf64(&mantissa, &binexp, mant1e256, exp1e256);
        if (decimal_exponent & 128) mulf64(&mantissa, &binexp, mant1e128, exp1e128);
        if (decimal_exponent &  64) mulf64(&mantissa, &binexp, mant1e64 , exp1e64 );
        if (decimal_exponent &  32) mulf64(&mantissa, &binexp, mant1e32 , exp1e32 );
        if (decimal_exponent &  16) mulf64(&mantissa, &binexp, mant1e16 , exp1e16 );
        if (decimal_exponent &   8) mulf64(&mantissa, &binexp, mant1e8  , exp1e8  );
        if (decimal_exponent &   4) mulf64(&mantissa, &binexp, mant1e4  , exp1e4  );
        if (decimal_exponent &   2) mulf64(&mantissa, &binexp, mant1e2  , exp1e2  );
        if (decimal_exponent &   1) mulf64(&mantissa, &binexp, mant1e1  , exp1e1  );
    }

    if (binexp <= -1023) {
        // subnormal number
        if (binexp <= -1023 - 64) {
            // underflow XXX @Clarify @Diagnostics do we want an error in this case?
            goto return_possibly_signed_zero;
        }
        // re-insert the implicit leading '1.' (which actually needs to be explicit for a subnormal number)
        mantissa = (mantissa >> 1) | 0x8000000000000000;
        // shift the fraction to bring the binary exponent into the representable range
        mantissa >>= (-binexp - 1023);
        binexp = -1023;
    }

    // printf("before rounding: 0x%016" PRIx64 ", binexp: %d\n", mantissa, binexp); // XXX DEBUG

    // round the mantissa before truncating it to 52 bits
    if ((mantissa & 0xFFF) > 0x800) {
        mantissa += (1 << 11);
        if (mantissa < ((uint64_t)1 << 11)) {
            binexp++;
            mantissa >>= 1;
        }
    }

    if (binexp > 1023) {
        // overflow XXX @Clarify @Diagnostics do we want an error in this case?
        goto return_possibly_signed_infinity;
    }

    // printf("before truncation: 0x%016" PRIx64 ", binexp: %d\n", mantissa, binexp); // XXX DEBUG

    assert(binexp <= 1023);
    assert(binexp >= -1023);

    // truncate the mantissa to 52 bits and assemble the double-precision IEEE 754 format number
    mantissa >>= 12;
    mantissa |= (uint64_t)(binexp + 1023) << 52;
    if (negative)
        mantissa |= 0x8000000000000000;

    // fprintf(stderr, "==> result: 0x%016" PRIx64 "\n", mantissa); // XXX DEBUG

    return *reinterpret_cast<double*>(&mantissa);
failed:
    assert(st->failed);
    return 0.0;
}

///////////////////////////// TEST CODE BELOW //////////////////////////////////////////

// from https://en.wikipedia.org/wiki/Xorshift ; slightly modified
uint64_t xorshift64(uint64_t *state)
{
    uint64_t x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return (*state = x);
}

HISTOGRAM_DEFINE_DEFAULTS_AND_INIT(cycles, 0);
HISTOGRAM_DEFINE_DEFAULTS_AND_INIT(atof_cycles, 0);

int main(int argc, char **argv)
{
    // override automatic binwidth detection (it does not work well because we have only fast cases first)
    hist_cycles.binwidth = 100;
    hist_atof_cycles.binwidth = 100;

    if (argc >= 2) {
        Status status;
        status.failed = false;
        status.error_message[0] = 0;

        DataStream stream;
        stream.ptr = (uint8_t*)argv[1];
        stream.end = stream.ptr + strlen(argv[1]);

        double result = parse_double(&status, &stream);
        if (status.failed)
            printf("ERROR: parsing failed: %s\n", status.error_message);
        else
            printf("OK; result = %g (0x%016" PRIx64 ")\n", result, *reinterpret_cast<uint64_t*>(&result));
    }
    else {
        // round-trip test
        // we iterate over all possible sign- and exponent bits
        // we also iterate over all combinations of the first few and the last few bits in the fraction part
        // the rest of the bits in the middle of the fraction part run over all zeroes, all ones, and a few pseudo-random values
        constexpr uint32_t n_fraction_leading_bits_enumerated = 4;
        constexpr uint32_t n_fraction_trailing_bits_enumerated = 4;

        constexpr uint32_t test_precision = 340;
        constexpr uint32_t threshold_precision_for_bit_exactness = 340; // ceil((1023 + 52) / log2(10)) + a few that I dont fully understand (probably due to subnormal numbers)

        constexpr bool use_scientific_notation = false;

        uint64_t pseudo_random_state = 1; // must be non-zero

        char buffer[4096];
        char result_buffer[4096];
        uint32_t number_of_tests = 0;
        uint32_t number_of_fails = 0;

        for (uint64_t signexp = 0; signexp < 4096; ++signexp) {
            for (uint64_t fraction_leading_bits = 0; fraction_leading_bits < (1 << n_fraction_leading_bits_enumerated); ++fraction_leading_bits) {
                for (uint64_t fraction_trailing_bits = 0; fraction_trailing_bits < (1 << n_fraction_trailing_bits_enumerated); ++fraction_trailing_bits) {
                    for (uint32_t middle = 0; middle < 16; ++middle) {

                        uint64_t fraction_middle_bits;
                        if (middle == 0)
                            fraction_middle_bits = 0;
                        else if (middle == 1)
                            fraction_middle_bits = UINT64_MAX;
                        else
                            fraction_middle_bits = xorshift64(&pseudo_random_state);

                        // reduce to the number of pseudo-random bits we need
                        fraction_middle_bits >>= ((64 - 52) + n_fraction_leading_bits_enumerated + n_fraction_trailing_bits_enumerated);

                        // another variant:
                        // fraction_middle_bits &= 0xF << (4 * middle);

                        uint64_t fraction = (fraction_leading_bits << (52 - n_fraction_leading_bits_enumerated))
                                          | (fraction_middle_bits << n_fraction_trailing_bits_enumerated)
                                          | fraction_trailing_bits;

                        uint64_t original_bits = (signexp << 52) | fraction;
                        double original = *reinterpret_cast<double*>(&original_bits);

                        char *printf_format = use_scientific_notation ? "%.*g" : "%.*f";
                        int n_chars = snprintf(buffer, sizeof(buffer), printf_format , test_precision, original);
                        assert(n_chars >= 0 && n_chars < sizeof(buffer));

                        // set up framework for input and error reporting
                        Status status;
                        status.failed = false;
                        status.error_message[0] = 0;
                        DataStream stream;
                        stream.ptr = (uint8_t*)buffer;
                        stream.end = (uint8_t*)buffer + n_chars;

                        // DUT, go for it!
                        uint64_t start_time = __rdtsc();
                        double result = parse_double(&status, &stream);
                        uint64_t end_time = __rdtsc();

                        uint64_t atof_start_time = __rdtsc();
                        double atof_result = atof(buffer);
                        uint64_t atof_end_time = __rdtsc();

                        // printf("result = %f, atof_result = %f\n", result, atof_result); // XXX DEBUG

                        uint64_t result_bits = *reinterpret_cast<uint64_t*>(&result);
                        uint64_t atof_result_bits = *reinterpret_cast<uint64_t*>(&atof_result);

                        if (number_of_tests % 101001 < 16) {
                            printf("sample(%8u): original = %14g (0x%016" PRIx64 "), result = %14g (0x%016" PRIx64 "),"
                                    " atof_result = %14g (0x%016" PRIx64 ")\n",
                                    number_of_tests, original, original_bits, result, result_bits, atof_result, atof_result_bits);
                        }

                        bool pass = true;
                        if (status.failed) {
                            printf("ERROR(%u): parsing failed: %s\n(buffer was '%s')\n", number_of_tests, status.error_message, buffer);
                            pass = false;
                        }
                        else {
                            // printf("PARSED; result = %g (0x%016" PRIx64 ")\n", result, *reinterpret_cast<uint64_t*>(&result)); // XXX DEBUG

                            bool roundtrip_ok;
                            if (!use_scientific_notation && test_precision >= threshold_precision_for_bit_exactness)
                                roundtrip_ok = (result_bits == original_bits);
                            else {
                                int n_chars = snprintf(result_buffer, sizeof(result_buffer), printf_format, test_precision, result);
                                assert(n_chars >= 0 && n_chars < sizeof(result_buffer));
                                roundtrip_ok = (strcmp(result_buffer, buffer) == 0);

                                // special case: scientific notation prints numbers that are out of the double-precision range due to rounding
                                if (use_scientific_notation && std::isinf(result) && result == atof_result)
                                    roundtrip_ok = true;
                            }
                            if (!roundtrip_ok) {
                                // we only support one flavor of (+/-) "nan"; therefore our round-trip does not reproduce all bits in these cases
                                if ((signexp == 0x7FF || signexp == 0xFFF) && (fraction && fraction != 0xFFFFFFFFFFFFF)) {
                                    if (!std::isnan(result)) {
                                        printf("FAILED(%u): expected nan result; buffer = '%s', original = %f (0x%016" PRIx64 "), result = %f (0x%016" PRIx64 ")\n",
                                           number_of_tests, buffer, original, original_bits, result, result_bits);
                                        pass = false;
                                    }
                                }
                                else {
                                    uint32_t output_precision = test_precision >= threshold_precision_for_bit_exactness ? threshold_precision_for_bit_exactness : test_precision;
                                    printf("FAILED(%u): round-trip error; buffer = '%s', original = %.*g (0x%016" PRIx64 "), result = %.*g (0x%016" PRIx64 ")\n",
                                           number_of_tests, buffer, output_precision, original, original_bits, output_precision, result, result_bits);
                                    pass = false;
                                }
                            }

                            if (test_precision >= threshold_precision_for_bit_exactness && result_bits != atof_result_bits) {
                                if ((signexp == 0x7FF || signexp == 0xFFF) && (fraction && fraction != 0xFFFFFFFFFFFFF)) {
                                    if (!std::isnan(atof_result)) {
                                        printf("FAILED(%u): disagreement with atof which did not yield the expected nan result;buffer = '%s', original = %f (0x%016" PRIx64 "), result = %f (0x%016" PRIx64 "),"
                                                " atof_result = %f (0x%016" PRIx64 ")\n",
                                                number_of_tests, buffer, original, original_bits, result, result_bits, atof_result, atof_result_bits);
                                        pass = false;
                                    }
                                }
                                else {
                                    printf("FAILED(%u): disagreement with atof; buffer = '%s', original = %f (0x%016" PRIx64 "), result = %f (0x%016" PRIx64 "),"
                                            " atof_result = %f (0x%016" PRIx64 ")\n",
                                           number_of_tests, buffer, original, original_bits, result, result_bits, atof_result, atof_result_bits);
                                    pass = false;
                                }
                            }
                        }
                        number_of_tests++;
                        if (!pass)
                            number_of_fails++;

                        uint64_t cycles = end_time - start_time;
                        uint64_t atof_cycles = atof_end_time - atof_start_time;
                        histogram_add(&hist_cycles, cycles);
                        histogram_add(&hist_atof_cycles, atof_cycles);
                    }
                }
            }
        }

        histogram_show(&hist_cycles     , "our  cycles> ", 70, stdout);
        histogram_show(&hist_atof_cycles, "atof cycles> ", 70, stdout);

        printf("Completed %u tests, %u failed.\n", number_of_tests, number_of_fails);
        printf("%s\n", number_of_fails ? "FAILED" : "OK");
    }
    return 0;
}
