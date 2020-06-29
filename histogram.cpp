// nice ASCII histograms

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


#include <cstring>
#include <cassert>
#include <cinttypes>
#include <cmath>
// XXX @Todo @NiceToHave Add optional buffers for N smallest/largest samples (interesting for tail analysis)
#include <cstdlib>
#include <algorithm> // for std::min

#include "histogram.h"

namespace stats {
    namespace histogram {
        // XXX @Todo @NiceToHave Add optional buffers for N smallest samples (interesting for tail analysis)
        void histogram_init(Histogram *hist, uint64_t *bins, uint32_t nbins, int64_t low, uint64_t binwidth,
                            HistogramSampleN *sample_buf, uint32_t n_sample_buf_elements)
        {
            // XXX assert against overflow of low + binwidth * nbins
            assert(hist);
            assert(bins);
            assert(nbins >= 3); // at least one low, one regular, and one high bin
            memset(hist, 0, sizeof(*hist));
            memset(bins, 0, nbins * sizeof(bins[0]));
            assert(binwidth > 0 || (sample_buf && n_sample_buf_elements > 0));
            assert(n_sample_buf_elements < ((uint64_t)1 << HistogramSampleN::n_bits_sequence_index));
            hist->low = low;
            hist->bins = bins;
            hist->binwidth = binwidth;
            hist->min_sample = INT64_MAX;
            hist->max_sample = INT64_MIN;
            hist->nbins = nbins;
            hist->sample_buf = sample_buf;
            hist->n_sample_buf_allocated = n_sample_buf_elements;
            hist->n_sample_buf_used = 0;
            hist->n_sample_buf_overflow = false;
            hist->largest_samples_buf = nullptr;
            hist->n_largest_samples_buf_allocated = 0;
            hist->threshold_for_largest = INT64_MAX; // only go into slow add_to_largest code path for samples that are INT64_MAX
        }

        void histogram_record_largest_samples(Histogram *hist, HistogramSampleN *largest_samples_buf, uint32_t n_largest_samples_buf_allocated)
        {
            hist->largest_samples_buf = largest_samples_buf;
            hist->n_largest_samples_buf_allocated = n_largest_samples_buf_allocated;
            hist->threshold_for_largest = INT64_MIN;
            for (uint32_t i = 0; i < hist->n_largest_samples_buf_allocated; ++i) {
                HistogramSampleN *entry = hist->largest_samples_buf + i;
                entry->sample = INT64_MIN;
                entry->n = 0;
                entry->sequence_index = 0;
            }
        }

        namespace {
            /**
             * \pre hist->binwidth > 0, i.e. this function only works once the binwidth is known.
             */
            inline uint32_t histogram_bin_index(Histogram *hist, int64_t x)
            {
                assert(hist->binwidth > 0);
                uint32_t index;
                if (x < hist->low) {
                    index = 0;
                }
                else {
                    // calculate delta with well-defined modulo 2^64 semantics
                    uint64_t delta = (uint64_t)x - (uint64_t)hist->low;
                    delta /= hist->binwidth;
                    if (delta >= hist->nbins - 1)
                        index = hist->nbins - 1;
                    else
                        index =  1 + (uint32_t)delta;
                }
                return index;
            }

            /**
             * \pre hist->binwidth > 0, i.e. this function only works once the binwidth is known.
             */
            void histogram_add_to_bin(Histogram *hist, int64_t x, uint64_t n)
            {
                assert(hist->binwidth > 0);
                uint32_t index = histogram_bin_index(hist, x);

                assert(index < hist->nbins);
                if ((hist->bins[index] += n) < n) {
                    hist->bins[index] = UINT64_MAX;
                    hist->count_overflowed = true;
                }
            }

            uint64_t round_to_multiple_of(uint64_t factor, uint64_t x)
            {
                if (x > UINT64_MAX - (factor / 2))
                    return factor * (UINT64_MAX / factor);
                else
                    return factor * ((x + factor / 2) / factor);
            }

            uint64_t histogram_round_binwidth(uint64_t binwidth)
            {
                uint64_t step = 10000000000000000000ui64;
                constexpr uint64_t divisor = 10;
                while (step >= divisor) {
                    if (binwidth >= step) return round_to_multiple_of(step / divisor, binwidth);
                    step /= 10;
                }
                return binwidth;
            }

            /**
             * \pre hist->binwidth == 0, i.e. the binwidth is not yet known.
             */
            // XXX implement optional automatic determination of hist->low
            void histogram_estimate_and_commit_binwidth(Histogram *hist)
            {
                assert(!hist->binwidth);
                assert(!hist->n_sample_buf_overflow); // sample_buf must not have been overrun already
                if (!hist->nsamples) {
                    hist->binwidth = 1;
                    return;
                }
                double quantile = histogram_quantile(hist, 0.90);
                assert(!isnan(quantile));
                double span = (quantile - hist->low) / 0.90 * 1.20;
                double binwidth_fractional = span / (hist->nbins - 2); // we divide by the # of regular bins
                hist->binwidth = histogram_round_binwidth((uint64_t)ceil(binwidth_fractional));
                if (hist->binwidth < 1)
                    hist->binwidth = 1;
                for (uint32_t i = 0; i < hist->n_sample_buf_used; ++i) {
                    HistogramSampleN *sn = hist->sample_buf + i;
                    histogram_add_to_bin(hist, sn->sample, sn->n);
                }
                assert(hist->binwidth);

                // printf("ESTIMATED BINWIDTH %" PRIu64 " based on %u of %u sample/n pair(s) (quantile = %g, span = %g, binwidth_fractional = %g)\n",
                //        hist->binwidth, hist->n_sample_buf_used, hist->n_sample_buf_allocated,
                //        quantile, span, binwidth_fractional);
            }

            void histogram_add_to_largest(Histogram *hist, int64_t x, uint64_t n)
            {
                assert(x >= hist->threshold_for_largest);
                assert(!hist->n_largest_samples_buf_allocated || x >= hist->largest_samples_buf[hist->n_largest_samples_buf_allocated - 1].sample);
                for (uint32_t i = 0; i < hist->n_largest_samples_buf_allocated; ++i) {
                    HistogramSampleN *entry = hist->largest_samples_buf + i;
                    if (x > entry->sample) {
                        // shift the rest of the buffer down and create a new entry here
                        memmove(entry + 1, entry, sizeof(*entry) * (hist->n_largest_samples_buf_allocated - i - 1));
                        entry->sample = x;
                        entry->n = n; // XXX @Incomplete check whether n fits
                        hist->threshold_for_largest = hist->largest_samples_buf[hist->n_largest_samples_buf_allocated - 1].sample;
                        return;
                    }
                    else if (x == entry->sample) {
                        if ((entry->n += n) < n) {
                            entry->n = UINT64_MAX;
                            hist->count_overflowed = true;
                        }
                        return;
                    }
                }
            }
        } // anonymous namespace

        // XXX @NiceToHave add trend analysis (linear regression of #sequence number ~ sample) with correlation estimate
        void histogram_add(Histogram *hist, int64_t x, uint64_t n)
        {
            if (!n)
                return;

            // first, update statistics which do not depend on knowing the binwidth
            if (hist->min_sample > x)
                hist->min_sample = x;
            if (hist->max_sample < x)
                hist->max_sample = x;

            hist->sum_samples += (double)x * (double)n;
            hist->sum_square_samples += (double)x * (double)x * (double)n;
            if ((hist->nsamples += n) < n) {
                hist->nsamples = UINT64_MAX;
                hist->count_overflowed = true;
            }

            if (x >= hist->threshold_for_largest)
                histogram_add_to_largest(hist, x, n);

            // invariant:
            assert(hist->binwidth > 0 || (hist->n_sample_buf_used < hist->n_sample_buf_allocated));

            assert(n < ((uint64_t)1 << HistogramSampleN::n_bits_n));
            if (!hist->n_sample_buf_overflow && hist->n_sample_buf_used
                && hist->sample_buf[hist->n_sample_buf_used - 1].sample == x
                && hist->sample_buf[hist->n_sample_buf_used - 1].n < (((uint64_t)1 << HistogramSampleN::n_bits_n) - n)
                && hist->sample_buf[hist->n_sample_buf_used - 1].sequence_index == hist->n_sample_buf_used - 1) {
                // add to the current sample_buf entry
                hist->sample_buf[hist->n_sample_buf_used - 1].n += n;
            }
            else if (hist->n_sample_buf_used == hist->n_sample_buf_allocated)
                hist->n_sample_buf_overflow = true;
            else {
                assert(hist->n_sample_buf_used < hist->n_sample_buf_allocated);
                hist->sample_buf[hist->n_sample_buf_used].sample = x;
                hist->sample_buf[hist->n_sample_buf_used].n = n;
                hist->sample_buf[hist->n_sample_buf_used].sequence_index = hist->n_sample_buf_used;
                hist->n_sample_buf_used++;
            }

            assert(hist->n_sample_buf_used <= hist->n_sample_buf_allocated);
            if (!hist->binwidth && hist->n_sample_buf_used == hist->n_sample_buf_allocated) {
                // the binwidth needs to be automatically adapted
                histogram_estimate_and_commit_binwidth(hist);
            }

            if (hist->binwidth)
                histogram_add_to_bin(hist, x, n);

            // invariant:
            assert(hist->binwidth > 0 || (hist->n_sample_buf_used < hist->n_sample_buf_allocated));
        }

        namespace {
            uint32_t integer_display_width(int64_t num, uint32_t n_digits_per_group = 3)
            {
                assert(n_digits_per_group <= 19);
                uint32_t width = 1;
                uint64_t magnitude = (uint64_t)num;
                if (num < 0) {
                    width++;
                    #pragma warning(suppress : 4146) // don't sweat about unsigned negation
                    magnitude = -magnitude;
                }
                if (n_digits_per_group) {
                    uint64_t group_range = 10;
                    for (uint32_t i = 1; i < n_digits_per_group; ++i)
                        group_range *= 10;
                    while (magnitude >= group_range) {
                        width += n_digits_per_group + 1;
                        magnitude /= group_range;
                    }
                }
                while (magnitude >= 10) {
                    width++;
                    magnitude /= 10;
                }
                return width;
            }

            uint32_t double_display_width(Histogram *hist, uint32_t n_fract_digits, double num, uint32_t *integer_width)
            {
                if (isnan(num)) {
                    return 3;
                }
                if (num < INT64_MIN || num > INT64_MAX) {
                    return 0;
                }
                uint64_t integer_part = (uint64_t)floor(abs(num));
                double fractional_part = abs(num) - (double)integer_part;
                uint32_t intwidth = integer_display_width(integer_part);
                if (num < 0)
                    intwidth++;
                if (integer_width && *integer_width < intwidth)
                    *integer_width = intwidth;
                uint32_t width = intwidth;
                width += 1 + n_fract_digits;
                return width;
            }

            bool print_digit_group(FILE *file, uint32_t width, bool leading, uint64_t lsd, uint32_t lsd_exp, uint64_t multiple,
                                   uint32_t n_digits_per_group = 3, bool negative = false)
            {
                if (n_digits_per_group == 0) {
                    assert(lsd == 1);
                    assert(leading);
                    fprintf(file, "%*" PRId64, width, negative ? -(int64_t)multiple : multiple);
                    return false;
                }
                int32_t n_digits = (int32_t)width - ((int32_t)lsd_exp / (int32_t)n_digits_per_group) * (int32_t)(1 + n_digits_per_group);
                if (n_digits < 0)
                    n_digits = 0;
                else if ((uint32_t)n_digits > n_digits_per_group)
                    n_digits = n_digits_per_group;
                multiple /= lsd;
                if (!leading) {
                    uint64_t group_range = 10;
                    for (uint32_t i = 1; i < n_digits_per_group; ++i)
                        group_range *= 10;
                    multiple = multiple % group_range;
                    n_digits = n_digits_per_group;
                }
                if (!leading || multiple > 0 || lsd == 1) {
                    fprintf(file, leading ? "%*" PRId64 : "%0*" PRId64, n_digits, (leading && negative) ? -(int64_t)multiple : multiple);
                    leading = false;
                }
                else
                    fprintf(file, "%*s", n_digits, "");
                return leading;
            }

            void pretty_print_double(Histogram *hist, uint32_t n_fract_digits, double num, FILE *file, uint32_t integer_width = 0)
            {
                // XXX @Hack implement this in a more complete and correct way
                // XXX verify this function
                if (isnan(num)) {
                    fprintf(file, "%*s%*s", integer_width, "nan", n_fract_digits + 1, "");
                    return;
                }
                if (num < INT64_MIN || num > INT64_MAX) {
                    fprintf(file, "%*.*g", integer_width + 1 + n_fract_digits, n_fract_digits, num);
                    return;
                }
                if (num < 0) {
                    fputc('-', file);
                    if (integer_width) integer_width--;
                }
                uint64_t integer_part = (uint64_t)floor(abs(num));
                double fractional_part = abs(num) - (double)integer_part;
                uint32_t width = integer_display_width(integer_part);
                if (width < integer_width)
                    width = integer_width;
                pretty_print_integer(width, integer_part, 1, file);
                uint32_t factor = 1;
                for (uint32_t i = 0; i < n_fract_digits; ++i)
                    factor *= 10;
                fprintf(file, ".%0*" PRIu64, n_fract_digits, (uint64_t)round(fractional_part * factor));
            }

            void determine_display_widths(Histogram *hist,
                    uint64_t *max_count, uint32_t *width_low, uint32_t *width_high, uint32_t *width_count)
            {
                assert(hist->binwidth > 0);

                *max_count = 0;

                for (uint32_t i = 0; i < hist->nbins; ++i) {
                    if (*max_count < hist->bins[i])
                        *max_count = hist->bins[i];
                }

                *width_low  = *max_count ? integer_display_width(hist->min_sample) : 0;
                *width_high = *max_count ? integer_display_width(hist->max_sample) : 0;
                uint32_t w[4];
                w[0] = integer_display_width(hist->low - 1);
                w[1] = integer_display_width(hist->low);
                w[2] = integer_display_width(hist->low + hist->nbins * hist->binwidth - 1);
                w[3] = integer_display_width(hist->low + hist->nbins * hist->binwidth);
                for (uint32_t i = 0; i < sizeof(w)/sizeof(w[0]); ++i) {
                    if (*width_low < w[i])
                        *width_low = w[i];
                    if (*width_high < w[i])
                        *width_high = w[i];
                }
                *width_count = integer_display_width(*max_count, 0);
            }

            int compare_sample_and_n_and_index(const void *a, const void *b)
            {
                HistogramSampleN *sn_a = (HistogramSampleN*)a;
                HistogramSampleN *sn_b = (HistogramSampleN*)b;
                return (sn_a->sample  < sn_b->sample)
                       ? -1
                       : (sn_a->sample == sn_b->sample) ? (
                           (sn_a->n > sn_b->n)
                           ? -1
                           : (sn_a->n == sn_b->n) ? (
                                 (sn_a->sequence_index < sn_b->sequence_index)
                                 ? -1
                                 : (sn_a->sequence_index == sn_b->sequence_index) ? 0
                                 : +1 )
                           : +1 )
                       : +1;
            }

            int compare_index(const void *a, const void *b)
            {
                HistogramSampleN *sn_a = (HistogramSampleN*)a;
                HistogramSampleN *sn_b = (HistogramSampleN*)b;
                return  (sn_a->sequence_index < sn_b->sequence_index)
                        ? -1
                        : (sn_a->sequence_index == sn_b->sequence_index) ? 0
                        : +1;
            }
        } // anonymous namespace

        void pretty_print_integer(uint32_t width, int64_t num, int64_t factor, FILE *file, uint32_t n_digits_per_group)
        {
            assert(num % factor == 0); // we only allow precise printing
            int64_t multiple = num / factor;

            bool negative = (num < 0);
            if (negative)
                multiple = -multiple;

            if (n_digits_per_group) {
                // Note: INT64_MIN,INT64_MAX have 19 decimal digits each
                uint32_t max_n_groups = (19 + n_digits_per_group - 1) / n_digits_per_group;
                uint64_t group_range = 10;
                for (uint32_t i = 1; i < n_digits_per_group; ++i)
                    group_range *= 10;
                uint32_t lsd_exp = 0;
                uint64_t lsd = 1;
                for (uint32_t i = 1; i < max_n_groups; ++i) {
                    lsd_exp += n_digits_per_group;
                    lsd *= group_range;
                }

                bool leading = true;
                while (true) {
                    leading = print_digit_group(file, width, leading, lsd, lsd_exp, multiple, n_digits_per_group, negative);
                    if (!lsd_exp)
                        break;
                    if (leading) {
                        // handle negative '-' signs that fall into a group-separating column
                        if (negative && lsd >= group_range && multiple / (lsd / group_range) >= group_range / 10) {
                            fputc('-', file);
                            leading = false;
                            negative = false;
                        }
                        else if (width >= (lsd_exp / n_digits_per_group) * (1 + n_digits_per_group))
                            fputc(' ', file);
                    }
                    else
                        fputc('\'', file);
                    lsd /= group_range;
                    lsd_exp -= n_digits_per_group;
                }
            }
            else {
                print_digit_group(file, width, true, 1, 0, multiple, n_digits_per_group, negative);
            }
            fprintf(file, "%s",
                    factor == 1 ? "" :
                    factor == 1'000 ? "k" :
                    factor == 1'000'000 ? "M" :
                    factor == 1'000'000'000 ? "G" :
                    factor == 1'000'000'000'000 ? "T" :
                    factor == 1'000'000'000'000'000 ? "P" :
                    factor == 1'000'000'000'000'000'000 ? "E" :
                    "<INVALID FACTOR>"
                    );
        }

        /**
         * Add recorded samples to bins if that has not already been done.
         * \note This is needed when a histogram with automatic binwidth adaptation
         *     is going to be read out before the sample_buf has been completely filled.
         * \note This function is idempotent.
         */
        void histogram_commit(Histogram *hist)
        {
            // if we already know the binwidth, this function is a no-op
            if (hist->binwidth)
                return;

            histogram_estimate_and_commit_binwidth(hist);
        }

        /**
         * \warning This function will commit the histogram (see histogram_commit) if
         *     that has not already happened.
         */
        // XXX @NiceToHave option to print vertically reversed?
        void histogram_show(Histogram *hist, char *indent, uint32_t width, FILE *file)
        {
            histogram_commit(hist);

            assert(width >= 1);

            uint64_t max_count;
            uint32_t width_low;
            uint32_t width_high;
            uint32_t width_count;
            determine_display_widths(hist, &max_count, &width_low, &width_high, &width_count);
            // XXX DEBUG printf("display_widths: %u %u %u (max_count=%" PRIu64 ")\n", width_low, width_high, width_count, max_count);

            if (hist->count_overflowed)
                fprintf(file, "%s!!! COUNTS OVERFLOWED IN THE FOLLOWING HISTOGRAM !!!\n", indent);

            double sample_mean = histogram_mean(hist);
            double sample_std = histogram_std(hist);
            uint32_t mean_std_intwidth = 0;
            uint32_t mean_std_width = std::max(
                    double_display_width(hist, 1, sample_mean, &mean_std_intwidth),
                    double_display_width(hist, 1, sample_std , &mean_std_intwidth));
            if (hist->nsamples) {
                mean_std_intwidth = std::max(mean_std_intwidth,
                    std::max(integer_display_width(hist->min_sample),
                             integer_display_width(hist->max_sample)));
                mean_std_width = std::max(mean_std_width, mean_std_intwidth + 2);
            }
            fprintf(file, "%smean = ", indent);
            pretty_print_double(hist, 1, sample_mean, file, mean_std_intwidth);
            fprintf(file, "\n");
            fprintf(file, "%s std = ", indent);
            pretty_print_double(hist, 1, sample_std, file, mean_std_intwidth);
            fprintf(file, " (%.1f%%); n = ",
                    sample_mean ? sample_std / sample_mean * 100.0 : std::nan(""));
            pretty_print_integer(0, hist->nsamples, 1, file);
            fprintf(file, "\n");
            if (hist->nsamples) {
                fprintf(file, "%s min = ", indent);
                pretty_print_integer(mean_std_intwidth, hist->min_sample, 1, file);
                fprintf(file, "\n");
                fprintf(file, "%s max = ", indent);
                pretty_print_integer(mean_std_intwidth, hist->max_sample, 1, file);
                fprintf(file, "\n");
                fprintf(file, "%sq{01 25; 50; 75 99} = {", indent);
                pretty_print_double(hist, 1, histogram_quantile(hist, 0.01), file);
                fprintf(file, " ");
                pretty_print_double(hist, 1, histogram_quantile(hist, 0.25), file);
                fprintf(file, "; ");
                pretty_print_double(hist, 1, histogram_median(hist), file);
                fprintf(file, "; ");
                pretty_print_double(hist, 1, histogram_quantile(hist, 0.75), file);
                fprintf(file, " ");
                pretty_print_double(hist, 1, histogram_quantile(hist, 0.99), file);
                fprintf(file, "}\n");
            }

            int64_t low = hist->min_sample;
            int64_t high = hist->low - 1;
            uint64_t cumul_count = 0;
            uint64_t denom = hist->nsamples ? hist->nsamples : 1;
            uint64_t binwidth_factor = 1;
            // XXX @Incomplete this does not really work, yet, as the min/max sample
            //     are typically not multiples of the binwidth_factor.
#if 0
            while (hist->binwidth % (1000 * binwidth_factor) == 0) {
                binwidth_factor *= 1000;
                width_high -= 1 + 3;
                width_low  -= 1 + 3;
            }
#endif
            // XXX @NiceToHave (optionally) skip unused bins at the beginning of the histogram (maybe with indication "# bins unused" if # > 1)
            // XXX @NiceToHave (optionally) skip unused bins at the end of the histogram (maybe with indication "# bins unused" if # > 1)
            // XXX @Incomplete at least print empty bins in a visually more distinct way (think about the logarithmic histogram made up of the bin counts!)
            for (uint32_t i = 0; i < hist->nbins; ++i) {
                bool is_high_overflow = (i == hist->nbins - 1);
                bool is_regular = (i > 0) && !is_high_overflow;
                uint64_t count = hist->bins[i];
                cumul_count += count;
                if (count || is_regular) {
                    if (is_high_overflow) {
                        assert(hist->max_sample >= low);
                        high = hist->max_sample;
                    }
                    fprintf(file, "%s[", indent);
                    pretty_print_integer(width_low, low, binwidth_factor, file);
                    fprintf(file, "; ");
                    pretty_print_integer(width_high, high, binwidth_factor, file);
                    fprintf(file, "] %3.0f%% %3.0f%% ",
                            (double)cumul_count / (double)denom * 100.0,
                            (double)count / (double)denom * 100.0);
                    pretty_print_integer(width_count, count, 1, file, 0);
                    fprintf(file, "%c",
                            is_regular ? '|' : is_high_overflow ? '>' : '<');
#if 0
                    fprintf(file, "%s[%*" PRId64 "; %*" PRId64 "] %3.0f%% %3.0f%% %*" PRIu64 "%c",
                            indent, width_low, low, width_high, high,
                            (double)cumul_count / (double)denom * 100.0,
                            (double)count / (double)denom * 100.0,
                            width_count, count,
                            is_regular ? '|' : is_high_overflow ? '>' : '<');
#endif
                    if (max_count) {
                        uint32_t bar = (uint32_t)((double)count / (double)max_count * (double)width + 0.5);
                        if (bar > width)
                            bar = width;
                        for (uint32_t j = 0; j < bar; ++j)
                            fputc('*', file);
                    }
                    fputc('\n', file);
                }
                low = high + 1;
                high += hist->binwidth;
            }
        }

        void histogram_show_largest_samples(Histogram *hist, char *indent, uint32_t width, FILE *file)
        {
            uint64_t max_count = 0;
            uint32_t width_sample = 0;
            for (uint32_t i = 0; i < hist->n_largest_samples_buf_allocated; ++i) {
                HistogramSampleN *entry = hist->largest_samples_buf + i;
                if (!entry->n)
                    continue;
                if (max_count < entry->n)
                    max_count = entry->n;
                uint32_t w = integer_display_width(entry->sample);
                if (width_sample < w)
                    width_sample = w;
            }
            uint32_t width_count = integer_display_width(max_count);

            for (uint32_t i = 0; i < hist->n_largest_samples_buf_allocated; ++i) {
                HistogramSampleN *entry = hist->largest_samples_buf + i;
                if (entry->n) {
                    fprintf(file, "%s", indent);
                    pretty_print_integer(width_sample, entry->sample, 1, file);
                    fprintf(file, "  (");
                    pretty_print_integer(width_count, entry->n, 1, file);
                    fprintf(file, "x)\n");
                }
            }
        }

        void histogram_print_sample_buf(Histogram *hist, char *indent, FILE *file)
        {
            uint64_t n_samples_buffered = 0;
            for (uint32_t i = 0; i < hist->n_sample_buf_used; ++i)
                n_samples_buffered += hist->sample_buf[i].n;

            if (!hist->sample_buf || !hist->n_sample_buf_allocated) {
                fprintf(file, "%s<no sample buffer>\n", indent);
            }
            fprintf(file, "%s%" PRIu64 " sample%s of %" PRIu64 " buffered (in %u value-n-pair%s%s)",
                    indent, n_samples_buffered, n_samples_buffered == 1 ? "" : "s", hist->nsamples,
                    hist->n_sample_buf_used, hist->n_sample_buf_used == 1 ? "" : "s",
                    hist->n_sample_buf_overflow ? "; sample_buf does not contain all samples" : "; sample_buf contains all samples");
            if (hist->n_sample_buf_used) {
                fprintf(file, "; ordered by sample value:\n");
                qsort(hist->sample_buf, hist->n_sample_buf_used, sizeof(*hist->sample_buf),
                      compare_sample_and_n_and_index);
                for (uint32_t i = 0; i < hist->n_sample_buf_used; ++i)
                    fprintf(file, "%s    @%-4" PRIu64 ": %8" PRId64 " x %" PRIu64 "\n", indent,
                            hist->sample_buf[i].sequence_index, hist->sample_buf[i].sample, hist->sample_buf[i].n);

                fprintf(file, "%sordered by sequence of sampling:\n", indent);
                qsort(hist->sample_buf, hist->n_sample_buf_used, sizeof(*hist->sample_buf),
                      compare_index);
                for (uint32_t i = 0; i < hist->n_sample_buf_used; ++i)
                    fprintf(file, "%s    %8" PRId64 " x %" PRIu64 "\n", indent,
                            hist->sample_buf[i].sample, hist->sample_buf[i].n);
            }
            else
                fprintf(file, "\n");
        }

        /**
         * \warning This function will commit the histogram (see histogram_commit) if
         *     that has not already happened.
         */
        // XXX @NiceToHave option to print vertically reversed?
        // XXX @NiceToHave option to print integrated value and related percentages
        void histogram_show_sample_buf(Histogram *hist, char *indent, FILE *file)
        {
            histogram_commit(hist);

            uint64_t max_count;
            uint32_t width_low;
            uint32_t width_high;
            uint32_t width_count;
            determine_display_widths(hist, &max_count, &width_low, &width_high, &width_count);

            qsort(hist->sample_buf, hist->n_sample_buf_used, sizeof(*hist->sample_buf),
                  compare_index);

            int64_t low = hist->min_sample;
            int64_t high = hist->low - 1;
            uint64_t cumul_count = 0;
            uint64_t denom = hist->nsamples ? hist->nsamples : 1;
            for (uint32_t i = 0; i < hist->nbins; ++i) {
                bool is_high_overflow = (i == hist->nbins - 1);
                bool is_regular = (i > 0) && !is_high_overflow;
                uint64_t count = hist->bins[i];
                cumul_count += count;
                if (count || is_regular) {
                    if (is_high_overflow) {
                        assert(hist->max_sample >= low);
                        high = hist->max_sample;
                    }
                    // XXX refactor into common function
                    fprintf(file, "%s[%*" PRId64 "; %*" PRId64 "] %3.0f%% %3.0f%% %*" PRIu64 "%c",
                            indent, width_low, low, width_high, high,
                            (double)cumul_count / (double)denom * 100.0,
                            (double)count / (double)denom * 100.0,
                            width_count, count,
                            is_regular ? '|' : is_high_overflow ? '>' : '<');
                    uint64_t prev_sequence_index = 0;
                    uint64_t buffered_count = 0;
                    for (uint32_t j = 0; j < hist->n_sample_buf_used; ++j) {
                        HistogramSampleN *sn = hist->sample_buf + j;
                        assert(sn->sequence_index >= prev_sequence_index);
                        if (histogram_bin_index(hist, sn->sample) != i)
                            continue;
                        while (prev_sequence_index < sn->sequence_index) {
                            fputc('.', file);
                            prev_sequence_index++;
                        }
                        assert(prev_sequence_index == sn->sequence_index);
                        if (prev_sequence_index != sn->sequence_index) {
                            fflush(stdout);
                            fprintf(stderr, "INDEX!!!!!!"); fflush(stderr);
                            abort(); // XXX DEBUG
                        }
                        if (sn->n < 10)
                            fprintf(file, "%" PRIu64, sn->n);
                        else
                            fputc('#', file);
                        buffered_count += sn->n;
                        prev_sequence_index++;
                    }
                    if (buffered_count > count) {
                        fflush(stdout);
                        fprintf(stderr, "BUFFERED_COUNT"); fflush(stderr);
                        abort(); // XXX DEBUG
                    }
                    if (buffered_count < count)
                        fputc('>', file);
                    fputc('\n', file);
                }
                low = high + 1;
                high += hist->binwidth;
            }
        }

        /**
         * Find bins ordered by decreasing count, i.e. the first call finds the
         * bin with the highest count. Bins with the same count are returned in
         * increasing order of their bounds.
         * \param index: As an input, *index shall be the zero-based index of the previous
         *     bin found by this function. If *n == 0 (on the first call), *index is
         *     ignored as an input.
         *     As an output, *index is set to the zero-based index of the next-frequent
         *     bin if such a bin could be found. *index is undefined on output if
         *     this function returns false.
         * \param n: As an input, *n shall be the count of the previous bin found
         *     by this function, or zero if this is the first call.
         *     As an output, *n is set to the count of the next-frequent
         *     bin if such a bin could be found. *n is undefined on output if this
         *     function return false;
         * \returns true if a next-frequent bin could be found, false otherwise.
         */
        bool histogram_next_frequent(Histogram *hist, uint32_t *index, uint64_t *n)
        {
            assert(hist->binwidth > 0);
            assert(hist);
            assert(index);
            assert(n);
            uint64_t limit;
            if (!*n) {
                // this is the first call; find the largest count of all
                limit = UINT64_MAX;
            }
            else {
                // first, try to find another bin with the same count
                for (uint32_t i = *index + 1; i < hist->nbins; ++i) {
                    if (hist->bins[i] == *n) {
                        *index = i;
                        return true;
                    }
                }
                // There is no bin with the same count. Find the next-highest count.
                limit = *n - 1;
            }
            *n = 0;
            for (uint32_t i = 0; i < hist->nbins; ++i) {
                if (hist->bins[i] > *n && hist->bins[i] <= limit) {
                    *index = i;
                    *n = hist->bins[i];
                }
            }
            return (*n > 0);
        }

        double histogram_mean(Histogram *hist)
        {
            if (!hist->nsamples)
                return nan("");
            return hist->sum_samples / (double)hist->nsamples;
        }

        double histogram_std(Histogram *hist)
        {
            uint64_t n = hist->nsamples;
            if (n < 2)
                return nan("");
            double mean = histogram_mean(hist);
            double mean_square_corrected = hist->sum_square_samples / (n - 1);
            return sqrt(mean_square_corrected - n / (n-1) * mean * mean);
        }

        double histogram_bin_center(Histogram *hist, uint32_t index)
        {
            assert(hist->binwidth > 0);
            if (index >= hist->nbins)
                return nan("");

            int64_t low;
            int64_t high;

            if (index == 0) {
                if (!hist->nsamples || hist->low == INT64_MIN)
                    return nan("");

                low = hist->min_sample;
                high = hist->low - 1;
                if (high > hist->max_sample)
                    high = hist->max_sample;
            }
            else if (index == hist->nbins - 1) {
                if (!hist->nsamples)
                    return nan("");
                low = hist->low + (index - 1) * hist->binwidth;
                high = hist->max_sample;
                if (low < hist->min_sample)
                    low = hist->min_sample;
            }
            else {
                low = hist->low + (index - 1) * hist->binwidth;
                high = low + hist->binwidth - 1;
            }

            if ((uint64_t)high - (uint64_t)low >= hist->binwidth)
                return nan("");

            return ((double)low + (double)high) / 2.0;
        }

        /**
         * \returns an estimate of the given order statistic, i.e. an estimate
         *     of the order-th sample value (1-based) when sample values are
         *     sorted in increasing numerical order. If there are no samples or
         *     if the order statistic cannot be estimated with an uncertainty
         *     of at most hist->binwidth / 2, the return value is NaN.
         * \note If the number of (sample, repeat count) pairs added so far
         *     fits in the allocated sample_buf, then the returned order
         *     statistic will be exact. Otherwise, it will be an estimate based
         *     on histogram bin centers.
         */
        double histogram_order_statistic(Histogram *hist, uint64_t order)
        {
            if (order < 1 || order > hist->nsamples)
                return nan("");

            assert(hist->nsamples > 0);
            if (order == 1)
                return (double)hist->min_sample;
            if (order == hist->nsamples)
                return (double)hist->max_sample;

            if (!hist->n_sample_buf_overflow) {
                // find the exact order statistic
                assert(hist->sample_buf && hist->n_sample_buf_allocated > 0);
                qsort(hist->sample_buf, hist->n_sample_buf_used, sizeof(*hist->sample_buf),
                      compare_sample_and_n_and_index);
                HistogramSampleN *sn = hist->sample_buf;
                while (order > sn->n) {
                    order -= sn->n;
                    sn++;
                    assert(sn < hist->sample_buf + hist->n_sample_buf_used);
                }
                return (double)sn->sample;
            }

            // estimate the order statistic from the histogram
            assert(hist->binwidth); // if we do not have a known binwidth at this time: REEEEEEEEEEEEEEEEEEEEEEEEE
            uint32_t index = 0;
            assert(index < hist->nbins);
            while (order > hist->bins[index]) {
                order -= hist->bins[index];
                index++;
                assert(index < hist->nbins);
            }
            return histogram_bin_center(hist, index);
        }

        double histogram_quantile(Histogram *hist, double p)
        {
            if (hist->nsamples < 1)
                return nan("");

            // We implement the estimator Q^_8(p) defined in:
            // R. J. Hyndman, Y. Fan (1996) Sample Quantiles in Statistical Packages,
            // "Americian Statistician" 50, pp. 361--365.
            // doi:10.2307/2684934
            double continuous_j = hist->nsamples * p + (p + 1)/3;
            double j = floor(continuous_j);
            double gamma = continuous_j - j;
            if (j < 1)
                return (double)hist->min_sample;
            if (j >= hist->nsamples)
                return (double)hist->max_sample;

            assert(j <= UINT64_MAX - 1);
            uint64_t order = (uint64_t)j;
            assert(hist->nsamples >= 2);
            assert(order >= 1);
            assert(order <= hist->nsamples - 1);

            double X_j   = histogram_order_statistic(hist, order);
            double X_jp1 = histogram_order_statistic(hist, order + 1);

            return (1 - gamma) * X_j + gamma * X_jp1;
        }

        double histogram_median(Histogram *hist)
        {
            return histogram_quantile(hist, .5);
        }

        namespace {
            int compare_category_count(const void *a, const void *b)
            {
                CategoryAndCount *cat_a = (CategoryAndCount*)a;
                CategoryAndCount *cat_b = (CategoryAndCount*)b;
                if (cat_a->count > cat_b->count)
                    return -1;
                if (cat_a->count < cat_b->count)
                    return +1;
                return 0;
            }
        }

        void show_categorical_histogram(CategoryAndCount *data, uint32_t n, char *prefix, uint32_t width, FILE *file,
                                        double limit_cumulative)
        {
            assert(width >= 1);

            qsort(data, n, sizeof(CategoryAndCount), compare_category_count);

            uint64_t max_count = 0;
            uint64_t total_count = 0;
            uint64_t rest_count = 0;
            uint64_t cumul_count = 0;
            uint32_t width_name = 5; // strlen("total")
            constexpr uint32_t max_name_width = 80;

            // loop over *all* items (even non-displayed ones)
            for (uint32_t i = 0; i < n; ++i) {
                CategoryAndCount *cat = data + i;
                total_count += cat->count;
            }

            uint64_t denom = total_count ? total_count : 1;
            uint32_t n_displayed = n;
            uint32_t n_rest = 0;

            // loop again but consider which items will actually be displayed
            for (uint32_t i = 0; i < n; ++i) {
                CategoryAndCount *cat = data + i;
                double cumul_fraction_before = (double)cumul_count / (double)denom;
                // Note: We do not create an "others" line if only one item remains.
                if (limit_cumulative != 0.0 && cumul_fraction_before > limit_cumulative && i < n - 1) {
                    n_displayed = i;
                    n_rest = n - n_displayed;
                    assert(n_rest >= 2);
                    rest_count = total_count - cumul_count;
                    if (max_count < rest_count)
                        max_count = rest_count;
                    uint32_t rest_width = integer_display_width(n_rest, 0) + 9; // strlen("( others)")
                    if (width_name < rest_width)
                        width_name = rest_width;
                    break;
                }
                cumul_count += cat->count;
                if (max_count < cat->count)
                    max_count = cat->count;
                size_t len = strlen(cat->category_name);
                if (len > max_name_width)
                    len = max_name_width;
                if (width_name < len)
                    width_name = (uint32_t)len;
            }

            uint32_t width_rank = integer_display_width(n_displayed + (n_rest > 0), 0);
            if (width_rank < 4)
                width_rank = 4;
            uint32_t width_count = integer_display_width(total_count, 0);
            cumul_count = 0;

            fprintf(file, "%s%*s %*s ", prefix, width_rank + 1, "rank", width_name, "total");
            fprintf(file, "%3s%% %3s%% ", "cum", "ind");
            pretty_print_integer(width_count, total_count, 1, file, 0);
            fprintf(file, "\n%s", prefix);
            for (uint32_t j = 0; j < width_rank + 2 + width_name + 1 + 10 + width_count; ++j)
                fputc('-', file);
            fputc('\n', file);

            for (uint32_t i = 0; i < n; ++i) {
                uint64_t count;
                char name_buf[max_name_width + 1];
                if (i >= n_displayed) {
                    count = rest_count;
                    assert(sizeof(name_buf) >= 9 + integer_display_width(n_rest) + 1);
                    snprintf(name_buf, sizeof(name_buf), "(%u others)", n_rest);
                }
                else {
                    count = data[i].count;
                    #pragma warning(suppress : 4996)
                    strncpy(name_buf, data[i].category_name, max_name_width + 1);
                    if (name_buf[max_name_width]) {
                        // XXX @Cleanup removed for the parse_double test (do not want to require C++17): static_assert(max_name_width >= 6);
                        memcpy(name_buf + max_name_width + 1 - 6, "[...]", 6);
                    }
                }
                cumul_count += count;
                fprintf(file, "%s", prefix);
                pretty_print_integer(width_rank, 1 + i, 1, file, 0);
                fprintf(file, ". %*s ", width_name, name_buf);
                fprintf(file, "%3.0f%% %3.0f%% ",
                        (double)cumul_count / (double)denom * 100.0,
                        (double)count / (double)denom * 100.0);
                pretty_print_integer(width_count, count, 1, file, 0);
                fputc(i >= n_displayed ? '>' : '|', file);
                if (max_count) {
                    uint32_t bar = (uint32_t)((double)count / (double)max_count * (double)width + 0.5);
                    if (bar > width)
                        bar = width;
                    for (uint32_t j = 0; j < bar; ++j)
                        fputc('*', file);
                }
                fputc('\n', file);
                if (i >= n_displayed)
                    break;
            }
        }
    } // namespace histogram
} // namespace stats
