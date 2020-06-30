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


#pragma once

#include <cstdio>
#include <cstdint>

#define HISTOGRAM_DEFINE(name, nbins, nsamplepairs) \
    ::stats::histogram::Histogram hist_ ## name; \
    uint64_t bins_ ## name[(nbins)]; \
    ::stats::histogram::HistogramSampleN samplebuf_ ## name[(nsamplepairs)];

#define HISTOGRAM_DEFINE_DEFAULTS(name) \
    HISTOGRAM_DEFINE(name, ::stats::histogram::default_n_bins, ::stats::histogram::default_n_samplebuf_elements)

#define HISTOGRAM_INIT_AUTO(name, low) \
    ::stats::histogram::histogram_init(&hist_ ## name, bins_ ## name, sizeof(bins_ ## name)/sizeof(bins_ ## name[0]), \
            (low), 0, samplebuf_ ## name, sizeof(samplebuf_ ## name)/sizeof(samplebuf_ ## name[0]));

#define HISTOGRAM_DEFINE_AND_INIT(name, nbins, nsamplepairs, low) \
    HISTOGRAM_DEFINE(name, nbins, nsamplepairs); \
    namespace { \
        struct HISTINIT__ ## name { HISTINIT__ ## name() { HISTOGRAM_INIT_AUTO(name, low); } }; \
        HISTINIT__ ## name histinit__ ## name; \
    }

#define HISTOGRAM_DEFINE_AND_INIT_CODE(name, nbins, nsamplepairs, low) \
    HISTOGRAM_DEFINE(name, nbins, nsamplepairs); \
    HISTOGRAM_INIT_AUTO(name, low);

#define HISTOGRAM_DEFINE_DEFAULTS_AND_INIT(name, low) \
    HISTOGRAM_DEFINE_AND_INIT(name, ::stats::histogram::default_n_bins, ::stats::histogram::default_n_samplebuf_elements, low)

#define HISTOGRAM_DEFINE_DEFAULTS_AND_INIT_CODE(name, low) \
    HISTOGRAM_DEFINE_AND_INIT_CODE(name, ::stats::histogram::default_n_bins, ::stats::histogram::default_n_samplebuf_elements, low)

namespace stats {
    namespace histogram {
        constexpr uint32_t default_n_bins = 64;
        constexpr uint32_t default_n_samplebuf_elements = 32;

        struct HistogramSampleN {
            static constexpr uint8_t n_bits_n = 48;
            static constexpr uint8_t n_bits_sequence_index = 64 - n_bits_n;
            int64_t sample;
            uint64_t n              : n_bits_n;
            uint64_t sequence_index : n_bits_sequence_index;
        };

        // XXX @NiceToHave optional descriptive stats per bin (at least mean and std) --> possible mini-graphical box plot per bin [  --^-] or ([   | ] for very small std)
        struct Histogram {
            int64_t low;
            uint64_t *bins;
            uint64_t binwidth;
            uint64_t nsamples;
            int64_t min_sample;
            int64_t max_sample;
            double sum_samples;
            double sum_square_samples;
            HistogramSampleN *sample_buf;
            HistogramSampleN *largest_samples_buf;
            int64_t threshold_for_largest;
            uint32_t nbins;
            uint32_t n_sample_buf_allocated;
            uint32_t n_sample_buf_used;
            uint32_t n_largest_samples_buf_allocated;
            bool count_overflowed;
            bool n_sample_buf_overflow;
        };

        void histogram_init(Histogram *hist, uint64_t *bins, uint32_t nbins, int64_t low, uint64_t binwidth,
                            HistogramSampleN *sample_buf = nullptr, uint32_t n_sample_buf_elements = 0);
        void histogram_record_largest_samples(Histogram *hist, HistogramSampleN *largest_samples_buf, uint32_t n_largest_samples_buf_allocated);
        void histogram_add(Histogram *hist, int64_t x, uint64_t n = 1);
        void histogram_commit(Histogram *hist);
        void histogram_show(Histogram *hist, char *indent, uint32_t width, FILE *file);
        void histogram_show_largest_samples(Histogram *hist, char *indent, uint32_t width, FILE *file);
        void histogram_print_sample_buf(Histogram *hist, char *indent, FILE *file);
        void histogram_show_sample_buf(Histogram *hist, char *indent, FILE *file);
        bool histogram_next_frequent(Histogram *hist, uint32_t *index, uint64_t *n);
        double histogram_mean(Histogram *hist);
        double histogram_std(Histogram *hist);
        double histogram_bin_center(Histogram *hist, uint32_t index);
        double histogram_order_statistic(Histogram *hist, uint64_t order);
        double histogram_quantile(Histogram *hist, double p);
        double histogram_median(Histogram *hist);

        struct CategoryAndCount {
            char *category_name;
            uint64_t count;
        };

        void show_categorical_histogram(CategoryAndCount *data, uint32_t n, char *prefix, uint32_t width, FILE *file,
                                        double limit_cumulative = 0.0);

        void pretty_print_integer(uint32_t width, int64_t num, int64_t factor, FILE *file, uint32_t n_digits_per_group = 3);
    } // namespace histogram
} // namespace stats

