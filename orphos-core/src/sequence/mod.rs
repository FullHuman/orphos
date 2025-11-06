//! Sequence encoding and manipulation utilities.
//!
//! This module provides functions for encoding DNA sequences into compact bitmap
//! representations and performing sequence analysis operations.
//!
//! ## Overview
//!
//! DNA sequences are encoded using a 2-bit representation where:
//! - A (adenine): 00
//! - C (cytosine): 01
//! - G (guanine): 10
//! - T/U (thymine/uracil): 11
//!
//! This encoding reduces memory usage by 75% compared to ASCII representation
//! and enables fast bitwise operations for sequence analysis.
//!
//! ## Modules
//!
//! - [`encoded`]: Encoded sequence structures with forward and reverse-complement
//! - [`io`]: FASTA file reading and parsing
//! - [`processing`]: Sequence analysis functions (GC content, codon detection)
//!
//! ## Examples
//!
//! ### Encode a sequence
//!
//! ```rust
//! use orphos_core::sequence::encoded::EncodedSequence;
//!
//! let sequence = b"ATGAAACGCATTAGCACCACCATT";
//! let encoded = EncodedSequence::without_masking(sequence);
//!
//! println!("Length: {} bp", encoded.sequence_length);
//! println!("GC content: {:.2}%", encoded.gc_content * 100.0);
//! ```
//!
//! ### Test for specific nucleotides
//!
//! ```rust
//! use orphos_core::sequence::{is_a, is_gc};
//! use orphos_core::sequence::encoded::EncodedSequence;
//!
//! let sequence = b"ATGC";
//! let encoded = EncodedSequence::without_masking(sequence);
//!
//! assert!(is_a(&encoded.forward_sequence, 0)); // Position 0 is 'A'
//! assert!(is_gc(&encoded.forward_sequence, 3)); // Position 3 is 'C'
//! ```

use crate::bitmap;
use crate::constants::MASK_SIZE;
use crate::types::*;

pub mod encoded;
pub mod io;
pub mod processing;

use crate::bitmap::{set_bit, test_bit, toggle_bit};
use rayon::prelude::*;

pub use io::*;
pub use processing::*;
use wide::CmpEq;
use wide::u8x32;

/// Converts nucleotide character to 2-bit encoding for bitmap storage.
///
/// Maps nucleotides to compact 2-bit representation for efficient storage
/// and fast sequence operations.
///
/// # Encoding
///
/// - A: 00 (0)
/// - C: 01 (1)
/// - G: 10 (2)
/// - T/U: 11 (3)
/// - Other: 4 (invalid marker)
///
/// # Arguments
///
/// * `c` - ASCII nucleotide character (case-insensitive)
///
/// # Returns
///
/// A value 0-3 for valid nucleotides, 4 for invalid characters.
///
/// # Examples
///
/// ```rust
/// use orphos_core::sequence::char_to_nuc;
///
/// assert_eq!(char_to_nuc(b'A'), 0);
/// assert_eq!(char_to_nuc(b'a'), 0);
/// assert_eq!(char_to_nuc(b'C'), 1);
/// assert_eq!(char_to_nuc(b'G'), 2);
/// assert_eq!(char_to_nuc(b'T'), 3);
/// assert_eq!(char_to_nuc(b'N'), 4); // Invalid
/// ```
#[must_use]
pub const fn char_to_nuc(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

/// Test if nucleotide at given position is adenine (A)
#[must_use]
pub fn is_a(encoded_sequence: &[u8], n: usize) -> bool {
    let bit_index = n * 2;
    !(test_bit(encoded_sequence, bit_index) || test_bit(encoded_sequence, bit_index + 1))
}

/// Test if nucleotide at given position is cytosine (C)
#[must_use]
pub fn is_c(encoded_sequence: &[u8], n: usize) -> bool {
    let bit_index = n * 2;
    !test_bit(encoded_sequence, bit_index) && test_bit(encoded_sequence, bit_index + 1)
}

/// Test if nucleotide at given position is guanine (G)
#[must_use]
pub fn is_g(encoded_sequence: &[u8], n: usize) -> bool {
    let bit_index = n * 2;
    test_bit(encoded_sequence, bit_index) && !test_bit(encoded_sequence, bit_index + 1)
}

/// Test if nucleotide at given position is thymine (T)
#[must_use]
pub fn is_t(encoded_sequence: &[u8], n: usize) -> bool {
    let bit_index = n * 2;
    test_bit(encoded_sequence, bit_index) && test_bit(encoded_sequence, bit_index + 1)
}

/// Test if position contains an unknown nucleotide (N)
pub fn is_n(unknown_sequence: &[u8], n: usize) -> bool {
    if n >= unknown_sequence.len() * 8 {
        return false;
    }
    test_bit(unknown_sequence, n)
}

/// Test if nucleotide at given position is G or C (high GC content indicator)
pub fn is_gc(encoded_sequence: &[u8], n: usize) -> bool {
    let bit_index = n * 2;
    test_bit(encoded_sequence, bit_index) != test_bit(encoded_sequence, bit_index + 1)
}

/// Check if genetic code table uses only ATG as start codon
const fn uses_only_atg(trans_table: i32) -> bool {
    matches!(trans_table, 6 | 10 | 14 | 15 | 16 | 22)
}

/// Check if GTG is not used as start codon in given translation table
const fn gtg_not_start(trans_table: i32) -> bool {
    matches!(trans_table, 1 | 3 | 12 | 22)
}

/// Check if TTG is not used as start codon in given translation table
fn ttg_not_start(trans_table: i32) -> bool {
    trans_table < 4 || trans_table == 9 || (21..25).contains(&trans_table)
}

/// Test if codon at given position is a valid start codon
///
/// Checks ATG, GTG, and TTG based on the genetic code table rules.
/// ATG is universally accepted as a start codon across all tables.
pub fn is_start(encoded_sequence: &[u8], pos: usize, training: &Training) -> bool {
    // ATG is always a start
    if is_atg(encoded_sequence, pos) {
        return true;
    }

    // Tables that only use ATG
    if uses_only_atg(training.translation_table) {
        return false;
    }

    // GTG
    if is_gtg(encoded_sequence, pos) && !gtg_not_start(training.translation_table) {
        return true;
    }

    // TTG
    if is_ttg(encoded_sequence, pos) && !ttg_not_start(training.translation_table) {
        return true;
    }

    false
}

/// Test if codon at position is ATG (methionine start codon)
pub fn is_atg(encoded_sequence: &[u8], pos: usize) -> bool {
    is_a(encoded_sequence, pos)
        && is_t(encoded_sequence, pos + 1)
        && is_g(encoded_sequence, pos + 2)
}

/// Test if codon at position is GTG (valine start codon)
pub fn is_gtg(encoded_sequence: &[u8], pos: usize) -> bool {
    is_g(encoded_sequence, pos)
        && is_t(encoded_sequence, pos + 1)
        && is_g(encoded_sequence, pos + 2)
}

/// Test if codon at position is TTG (leucine start codon)
pub fn is_ttg(encoded_sequence: &[u8], pos: usize) -> bool {
    is_t(encoded_sequence, pos)
        && is_t(encoded_sequence, pos + 1)
        && is_g(encoded_sequence, pos + 2)
}

/// Check if TAG is recognized as stop codon in the given translation table
const fn is_tag_stop(trans_table: i32) -> bool {
    !matches!(trans_table, 6 | 15 | 16 | 22)
}

/// Check if TGA is recognized as stop codon in the given translation table
const fn is_tga_stop(trans_table: i32) -> bool {
    !matches!(trans_table, 2..=5 | 9 | 10 | 13 | 14 | 21 | 25)
}

/// Check if TAA is recognized as stop codon in the given translation table
const fn is_taa_stop(trans_table: i32) -> bool {
    !matches!(trans_table, 6 | 14)
}

/// Test if codon at given position is a stop codon
///
/// Checks for standard stop codons (TAA, TAG, TGA) and special cases
/// based on the genetic code translation table being used.
#[inline]
pub fn is_stop(encoded_sequence: &[u8], pos: usize, training: &Training) -> bool {
    if is_t(encoded_sequence, pos) {
        if is_a(encoded_sequence, pos + 1) {
            if is_g(encoded_sequence, pos + 2) {
                return is_tag_stop(training.translation_table);
            }
            if is_a(encoded_sequence, pos + 2) {
                return is_taa_stop(training.translation_table);
            }
        } else if is_g(encoded_sequence, pos + 1) && is_a(encoded_sequence, pos + 2) {
            return is_tga_stop(training.translation_table);
        }
    }

    // Special cases for different translation tables
    match training.translation_table {
        2 => {
            // AGA or AGG are stop codons in translation table 2
            is_a(encoded_sequence, pos)
                && is_g(encoded_sequence, pos + 1)
                && (is_a(encoded_sequence, pos + 2) || is_g(encoded_sequence, pos + 2))
        }
        22 => {
            // TCA is a stop codon in translation table 22
            is_t(encoded_sequence, pos)
                && is_c(encoded_sequence, pos + 1)
                && is_a(encoded_sequence, pos + 2)
        }
        23 => {
            // TTA is a stop codon in translation table 23
            is_t(encoded_sequence, pos)
                && is_t(encoded_sequence, pos + 1)
                && is_a(encoded_sequence, pos + 2)
        }
        _ => false,
    }
}

/// Calculate the GC content of a sequence region
///
/// Returns the fraction of nucleotides that are G or C within
/// the specified range (inclusive).
pub fn gc_content(encoded_sequence: &[u8], start: usize, end: usize) -> f64 {
    if start > end {
        return 0.0;
    }

    let (gc_count, total) = (start..=end)
        .map(|i| {
            if is_g(encoded_sequence, i) || is_c(encoded_sequence, i) {
                (1, 1)
            } else {
                (0, 1)
            }
        })
        .fold((0, 0), |(gc, tot), (g, t)| (gc + g, tot + t));

    if total == 0 {
        0.0
    } else {
        f64::from(gc_count) / f64::from(total)
    }
}

/// Convert a forward strand reading frame to its corresponding reverse strand frame
///
/// Maps reading frames between forward and reverse strands accounting for
/// sequence length and frame relationships.
pub const fn reverse_strand_reading_frame(forward_frame: usize, sequence_length: usize) -> usize {
    let frame_modulus = if sequence_length.is_multiple_of(3) {
        3
    } else {
        sequence_length % 3
    };
    (frame_modulus - 1 - forward_frame) % 3
}

/// Determine which of three reading frames has the highest score
///
/// Returns the index (0, 1, or 2) of the frame with maximum value.
pub const fn find_max_reading_frame(
    frame_0_value: i32,
    frame_1_value: i32,
    frame_2_value: i32,
) -> usize {
    if frame_0_value > frame_1_value {
        if frame_0_value > frame_2_value { 0 } else { 2 }
    } else if frame_1_value > frame_2_value {
        1
    } else {
        2
    }
}

/// Generate the reverse complement of an encoded DNA sequence
///
/// Creates a reverse complement sequence using 2-bit encoding,
/// handling both known and unknown nucleotides properly.
pub fn create_reverse_complement_sequence(
    forward_sequence: &[u8],
    unknown_sequence: &[u8],
    nucleotide_length: usize,
) -> Vec<u8> {
    let mut reverse_complement_encoded_sequence = vec![0; forward_sequence.len()];
    let sequence_length = nucleotide_length * 2;

    for i in 0..sequence_length {
        if !test_bit(forward_sequence, i) {
            let target_pos = if i % 2 == 0 {
                sequence_length - i - 2
            } else {
                sequence_length - i
            };
            if target_pos < sequence_length {
                set_bit(&mut reverse_complement_encoded_sequence, target_pos);
            }
        }
    }

    for i in 0..nucleotide_length {
        if test_bit(unknown_sequence, i) && sequence_length >= 2 + i * 2 {
            toggle_bit(
                &mut reverse_complement_encoded_sequence,
                sequence_length - 1 - i * 2,
            );
            toggle_bit(
                &mut reverse_complement_encoded_sequence,
                sequence_length - 2 - i * 2,
            );
        }
    }
    reverse_complement_encoded_sequence
}

/// Return the minimum of two integers (utility function)
#[inline]
pub fn min_of_two_integers(first_value: i32, second_value: i32) -> i32 {
    first_value.min(second_value)
}

/// Calculate k-mer index from sequence position for frequency analysis
///
/// Converts a sequence position to a numeric index representing the
/// k-mer pattern, used for codon usage and frequency calculations.
#[must_use]
pub fn calculate_kmer_index(kmer_length: usize, encoded_sequence: &[u8], position: usize) -> usize {
    let mut kmer_index = 0;
    for i in 0..(2 * kmer_length) {
        let bit_pos = position * 2 + i;
        kmer_index |= usize::from(test_bit(encoded_sequence, bit_pos)) << i;
    }
    kmer_index
}

/// Calculate background k-mer frequencies for both strands
///
/// Computes frequency distributions of k-mers across the entire sequence,
/// used for statistical modeling of codon usage patterns.
pub fn calculate_background_mer_frequencies(
    length: usize,
    encoded_sequence: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    bg: &mut [f64],
) {
    let mut size = 1usize;

    for _i in 1..=length {
        size *= 4;
    }

    // Use parallel processing to count k-mers
    let chunk_size = std::cmp::max(
        1000,
        (sequence_length - length + 1) / rayon::current_num_threads(),
    );
    let total_counts: Vec<i32> = (0..(sequence_length - length + 1))
        .into_par_iter()
        .chunks(chunk_size)
        .map(|chunk| {
            let mut local_counts = vec![0i32; size];

            for i in chunk {
                let seq_idx = calculate_kmer_index(length, encoded_sequence, i);
                if seq_idx < size {
                    local_counts[seq_idx] += 1;
                }

                let rseq_idx = calculate_kmer_index(length, reverse_complement_encoded_sequence, i);
                if rseq_idx < size {
                    local_counts[rseq_idx] += 1;
                }
            }

            local_counts
        })
        .reduce(
            || vec![0i32; size],
            |mut acc, local_counts| {
                for (i, &count) in local_counts.iter().enumerate() {
                    acc[i] += count;
                }
                acc
            },
        );

    let glob = (sequence_length - length + 1) * 2;

    bg.par_iter_mut()
        .enumerate()
        .take(size)
        .for_each(|(i, bg_val)| {
            *bg_val = f64::from(total_counts[i]) / (glob as f64);
        });
}

/// Convert k-mer index back to nucleotide sequence representation
///
/// Decodes a numeric k-mer index back to its original DNA sequence
/// for display and debugging purposes.
pub fn mer_text(len: usize, bit_index: usize) -> String {
    use crate::constants::NUCLEOTIDE_LETTERS;

    if len == 0 {
        return "None".to_string();
    }

    let mut result = String::with_capacity(len);
    let index = bit_index;

    for i in 0..len {
        // Extract 2 bits for position i
        let val = (index & (1 << (2 * i))) | (index & (1 << (2 * i + 1)));
        let val = val >> (i * 2);
        let base_idx = val & 0b11; // Ensure we only get 2 bits

        if base_idx < 4 {
            result.push(NUCLEOTIDE_LETTERS[base_idx]);
        } else {
            result.push('N'); // Fallback for invalid values
        }
    }

    result
}

/// Encode a DNA sequence into compact 2-bit representation
///
/// Converts raw DNA sequence to bitmap format for efficient storage and processing.
/// Also handles masking of low-complexity regions and tracks unknown nucleotides.
/// Returns the GC content of the encoded sequence.
pub fn encode_sequence(
    sequence: &[u8],
    encoded_sequence: &mut [u8],
    unknown_sequence: &mut [u8],
    masks: &mut Vec<Mask>,
    do_mask: bool,
) -> Result<f64, OrphosError> {
    let mut gc_count = 0;
    let mut total_count = 0;
    let mut mask_start: Option<usize> = None;

    for (i, &byte) in sequence.iter().enumerate() {
        if i * 2 + 1 >= encoded_sequence.len() * 8 {
            break;
        }

        // Handle masking for runs of N's
        if do_mask {
            if let Some(start) = mask_start {
                if byte != b'N' && byte != b'n' {
                    if i - start >= MASK_SIZE {
                        masks.push(Mask {
                            begin: start,
                            end: i - 1,
                        });
                    }
                    mask_start = None;
                }
            } else if byte == b'N' || byte == b'n' {
                mask_start = Some(i);
            }
        }

        // Encode nucleotide in bitmap format
        let bctr = i * 2;
        match byte.to_ascii_uppercase() {
            b'A' => {
                total_count += 1;
            }
            b'C' => {
                bitmap::set_bit(encoded_sequence, bctr + 1);
                gc_count += 1;
                total_count += 1;
            }
            b'G' => {
                bitmap::set_bit(encoded_sequence, bctr);
                gc_count += 1;
                total_count += 1;
            }
            b'T' | b'U' => {
                bitmap::set_bit(encoded_sequence, bctr);
                bitmap::set_bit(encoded_sequence, bctr + 1);
                total_count += 1;
            }
            _ => {
                bitmap::set_bit(encoded_sequence, bctr + 1);
                bitmap::set_bit(unknown_sequence, i);
                total_count += 1;
            }
        }
    }

    // Handle final mask if sequence ends with N's
    if do_mask
        && let Some(start) = mask_start
        && sequence.len() - start >= MASK_SIZE
    {
        masks.push(Mask {
            begin: start,
            end: sequence.len() - 1,
        });
    }

    let gc_content = if total_count > 0 {
        f64::from(gc_count) / f64::from(total_count)
    } else {
        0.0
    };

    Ok(gc_content)
}

/// SIMD-accelerated encoding using the `wide` crate with u8x32
///
/// Uses portable SIMD operations to process 32 nucleotides at once.
/// Returns the GC content of the encoded sequence.
pub fn encode_sequence_simd_wide(
    sequence: &[u8],
    encoded_sequence: &mut [u8],
    unknown_sequence: &mut [u8],
) -> Result<f64, OrphosError> {
    let mut gc_count = 0u32;
    let mut total_count = 0u32;

    // Process 32 bytes at a time with u8x32
    use crate::constants::CHUNK_SIZE;
    let chunks = sequence.len() / CHUNK_SIZE;

    // SIMD constants for nucleotide detection
    let a_upper = u8x32::splat(b'A');
    let c_upper = u8x32::splat(b'C');
    let g_upper = u8x32::splat(b'G');
    let t_upper = u8x32::splat(b'T');
    let u_upper = u8x32::splat(b'U');
    // let n_upper = u8x32::splat(b'N');

    let a_lower = u8x32::splat(b'a');
    let c_lower = u8x32::splat(b'c');
    let g_lower = u8x32::splat(b'g');
    let t_lower = u8x32::splat(b't');
    let u_lower = u8x32::splat(b'u');
    // let n_lower = u8x32::splat(b'n');

    for chunk_idx in 0..chunks {
        let chunk_start = chunk_idx * CHUNK_SIZE;

        // Safely load 32 bytes into SIMD vector
        let input_slice = &sequence[chunk_start..chunk_start + CHUNK_SIZE];

        // Convert slice to array for u8x32::from()
        let mut input_array = [0u8; 32];
        input_array.copy_from_slice(input_slice);
        let input = u8x32::from(input_array);

        // Create masks for each nucleotide type (case-insensitive)
        let is_a = input.cmp_eq(a_upper) | input.cmp_eq(a_lower);
        let is_c = input.cmp_eq(c_upper) | input.cmp_eq(c_lower);
        let is_g = input.cmp_eq(g_upper) | input.cmp_eq(g_lower);
        let is_t = input.cmp_eq(t_upper)
            | input.cmp_eq(t_lower)
            | input.cmp_eq(u_upper)
            | input.cmp_eq(u_lower);
        // let is_n = input.cmp_eq(n_upper) | input.cmp_eq(n_lower);

        let gc_mask = is_g | is_c;
        let valid_mask = is_a | is_c | is_g | is_t;

        // Convert to bitmask and count set bits
        gc_count += gc_mask.move_mask().count_ones();
        total_count += valid_mask.move_mask().count_ones();

        // Process each nucleotide for 2-bit encoding
        // let input_array: [u8; 32] = input.into();
        let is_a_array: [u8; 32] = is_a.into();
        let is_c_array: [u8; 32] = is_c.into();
        let is_g_array: [u8; 32] = is_g.into();
        let is_t_array: [u8; 32] = is_t.into();
        // let is_n_array: [u8; 32] = is_n.into();

        for i in 0..CHUNK_SIZE {
            let pos = chunk_start + i;
            if pos >= sequence.len() || pos * 2 + 1 >= encoded_sequence.len() * 8 {
                break;
            }

            let bit_pos = pos * 2;

            // Use SIMD results to determine nucleotide type
            if is_a_array[i] != 0 {
                // A = 00 (default, no bits to set)
            } else if is_c_array[i] != 0 {
                // C = 01
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
            } else if is_g_array[i] != 0 {
                // G = 10
                crate::bitmap::set_bit(encoded_sequence, bit_pos);
            } else if is_t_array[i] != 0 {
                // T/U = 11
                crate::bitmap::set_bit(encoded_sequence, bit_pos);
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
            } else {
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                crate::bitmap::set_bit(unknown_sequence, pos);
            }
        }
    }

    // Handle remaining bytes (scalar fallback)
    for (pos, byte) in sequence.iter().enumerate().skip(chunks * CHUNK_SIZE) {
        if pos * 2 + 1 >= encoded_sequence.len() * 8 {
            break;
        }

        // let byte = sequence[pos];
        let bit_pos = pos * 2;

        match byte.to_ascii_uppercase() {
            b'A' => {
                total_count += 1;
            }
            b'C' => {
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                gc_count += 1;
                total_count += 1;
            }
            b'G' => {
                crate::bitmap::set_bit(encoded_sequence, bit_pos);
                gc_count += 1;
                total_count += 1;
            }
            b'T' | b'U' => {
                crate::bitmap::set_bit(encoded_sequence, bit_pos);
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                total_count += 1;
            }
            _ => {
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                crate::bitmap::set_bit(unknown_sequence, pos);
                total_count += 1;
            }
        }
    }

    let gc_content = if total_count > 0 {
        gc_count as f64 / total_count as f64
    } else {
        0.0
    };

    Ok(gc_content)
}

/// Optimized packed encoding version with u8x32 and batch bit operations
pub fn encode_sequence_simd_wide_packed(
    sequence: &[u8],
    encoded_sequence: &mut [u8],
    unknown_sequence: &mut [u8],
) -> Result<f64, OrphosError> {
    let mut gc_count = 0u32;
    let mut total_count = 0u32;

    // Process 32 bytes at a time
    use crate::constants::CHUNK_SIZE;
    let chunks = sequence.len() / CHUNK_SIZE;

    for chunk_idx in 0..chunks {
        let chunk_start = chunk_idx * CHUNK_SIZE;

        // Load 32 bytes
        let input_slice = &sequence[chunk_start..chunk_start + CHUNK_SIZE];
        let mut input_array = [0u8; 32];
        input_array.copy_from_slice(input_slice);
        let input = u8x32::from(input_array);

        // SIMD nucleotide detection
        let a_upper = u8x32::splat(b'A');
        let c_upper = u8x32::splat(b'C');
        let g_upper = u8x32::splat(b'G');
        let t_upper = u8x32::splat(b'T');
        let u_upper = u8x32::splat(b'U');

        let a_lower = u8x32::splat(b'a');
        let c_lower = u8x32::splat(b'c');
        let g_lower = u8x32::splat(b'g');
        let t_lower = u8x32::splat(b't');
        let u_lower = u8x32::splat(b'u');

        let is_a = input.cmp_eq(a_upper) | input.cmp_eq(a_lower);
        let is_c = input.cmp_eq(c_upper) | input.cmp_eq(c_lower);
        let is_g = input.cmp_eq(g_upper) | input.cmp_eq(g_lower);
        let is_t = input.cmp_eq(t_upper)
            | input.cmp_eq(t_lower)
            | input.cmp_eq(u_upper)
            | input.cmp_eq(u_lower);

        let gc_mask = is_g | is_c;
        let valid_mask = is_a | is_c | is_g | is_t;

        gc_count += gc_mask.move_mask().count_ones();
        total_count += valid_mask.move_mask().count_ones();

        // Extract SIMD results for bit setting - move_mask() returns i32
        // let is_a_mask: i32 = is_a.move_mask();
        let is_c_mask: i32 = is_c.move_mask();
        let is_g_mask: i32 = is_g.move_mask();
        let is_t_mask: i32 = is_t.move_mask();
        let unknown_mask: i32 = !valid_mask.move_mask();

        // Batch process nucleotides using bit operations
        for i in 0..CHUNK_SIZE {
            let pos = chunk_start + i;
            if pos >= sequence.len() || pos * 2 + 1 >= encoded_sequence.len() * 8 {
                break;
            }

            let bit_pos = pos * 2;
            let bit_flag = 1i32 << i;

            // Use bit operations to check SIMD results
            if (is_c_mask & bit_flag) != 0 {
                // C = 01
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
            } else if (is_g_mask & bit_flag) != 0 {
                // G = 10
                crate::bitmap::set_bit(encoded_sequence, bit_pos);
            } else if (is_t_mask & bit_flag) != 0 {
                // T/U = 11
                crate::bitmap::set_bit(encoded_sequence, bit_pos);
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
            } else if (unknown_mask & bit_flag) != 0 {
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                crate::bitmap::set_bit(unknown_sequence, pos);
            }
            // A = 00 (default, no bits to set)
        }
    }

    // Handle remaining bytes with scalar fallback
    // for pos in (chunks * CHUNK_SIZE)..sequence.len() {
    for (pos, byte) in sequence.iter().enumerate().skip(chunks * CHUNK_SIZE) {
        if pos * 2 + 1 >= encoded_sequence.len() * 8 {
            break;
        }

        // let byte = sequence[pos];
        let bit_pos = pos * 2;

        match byte.to_ascii_uppercase() {
            b'A' => {
                total_count += 1;
            }
            b'C' => {
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                gc_count += 1;
                total_count += 1;
            }
            b'G' => {
                crate::bitmap::set_bit(encoded_sequence, bit_pos);
                gc_count += 1;
                total_count += 1;
            }
            b'T' | b'U' => {
                crate::bitmap::set_bit(encoded_sequence, bit_pos);
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                total_count += 1;
            }
            _ => {
                crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                crate::bitmap::set_bit(unknown_sequence, pos);
                total_count += 1;
            }
        }
    }

    let gc_content = if total_count > 0 {
        gc_count as f64 / total_count as f64
    } else {
        0.0
    };

    Ok(gc_content)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Training;

    #[test]
    fn test_char_to_nuc_valid_bases() {
        assert_eq!(char_to_nuc(b'A'), 0);
        assert_eq!(char_to_nuc(b'a'), 0);
        assert_eq!(char_to_nuc(b'C'), 1);
        assert_eq!(char_to_nuc(b'c'), 1);
        assert_eq!(char_to_nuc(b'G'), 2);
        assert_eq!(char_to_nuc(b'g'), 2);
        assert_eq!(char_to_nuc(b'T'), 3);
        assert_eq!(char_to_nuc(b't'), 3);
        assert_eq!(char_to_nuc(b'U'), 3);
        assert_eq!(char_to_nuc(b'u'), 3);
    }

    #[test]
    fn test_char_to_nuc_invalid_bases() {
        assert_eq!(char_to_nuc(b'N'), 4);
        assert_eq!(char_to_nuc(b'n'), 4);
        assert_eq!(char_to_nuc(b'X'), 4);
        assert_eq!(char_to_nuc(b'-'), 4);
        assert_eq!(char_to_nuc(b' '), 4);
    }

    #[test]
    fn test_nucleotide_check_functions() {
        let mut seq = vec![0u8; 10];

        // Encode ATCG at positions 0,1,2,3
        // A = 00 (default)
        crate::bitmap::set_bit(&mut seq, 2); // T = 11
        crate::bitmap::set_bit(&mut seq, 3);
        crate::bitmap::set_bit(&mut seq, 5); // C = 01
        crate::bitmap::set_bit(&mut seq, 6); // G = 10

        assert!(is_a(&seq, 0));
        assert!(!is_a(&seq, 1));
        assert!(!is_a(&seq, 2));
        assert!(!is_a(&seq, 3));

        assert!(!is_t(&seq, 0));
        assert!(is_t(&seq, 1));
        assert!(!is_t(&seq, 2));
        assert!(!is_t(&seq, 3));

        assert!(!is_c(&seq, 0));
        assert!(!is_c(&seq, 1));
        assert!(is_c(&seq, 2));
        assert!(!is_c(&seq, 3));

        assert!(!is_g(&seq, 0));
        assert!(!is_g(&seq, 1));
        assert!(!is_g(&seq, 2));
        assert!(is_g(&seq, 3));
    }

    #[test]
    fn test_start_codon_functions() {
        let mut seq = vec![0u8; 20];

        // Encode ATG at position 0: A(00) T(11) G(10)
        crate::bitmap::set_bit(&mut seq, 2); // T
        crate::bitmap::set_bit(&mut seq, 3);
        crate::bitmap::set_bit(&mut seq, 4); // G

        assert!(is_atg(&seq, 0));
        assert!(!is_gtg(&seq, 0));
        assert!(!is_ttg(&seq, 0));

        // Encode GTG at position 3: G(10) T(11) G(10)
        crate::bitmap::set_bit(&mut seq, 6); // G
        crate::bitmap::set_bit(&mut seq, 8); // T
        crate::bitmap::set_bit(&mut seq, 9);
        crate::bitmap::set_bit(&mut seq, 10); // G

        assert!(!is_atg(&seq, 3));
        assert!(is_gtg(&seq, 3));
        assert!(!is_ttg(&seq, 3));
    }

    #[test]
    fn test_stop_codon_functions() {
        let mut seq = vec![0u8; 20];
        let training = Training::default();

        // Encode TAA at position 0: T(11) A(00) A(00)
        crate::bitmap::set_bit(&mut seq, 0); // T
        crate::bitmap::set_bit(&mut seq, 1);
        // A and A are default (00)

        assert!(is_stop(&seq, 0, &training));

        // Encode TAG at position 3: T(11) A(00) G(10)
        crate::bitmap::set_bit(&mut seq, 6); // T
        crate::bitmap::set_bit(&mut seq, 7);
        crate::bitmap::set_bit(&mut seq, 10); // G

        assert!(is_stop(&seq, 3, &training));
    }

    #[test]
    fn test_gc_content_calculation() {
        let mut seq = vec![0u8; 20];

        // Encode ATCG: A(00) T(11) C(01) G(10)
        crate::bitmap::set_bit(&mut seq, 2); // T
        crate::bitmap::set_bit(&mut seq, 3);
        crate::bitmap::set_bit(&mut seq, 5); // C
        crate::bitmap::set_bit(&mut seq, 6); // G

        let gc = gc_content(&seq, 0, 3);
        assert!((gc - 0.5).abs() < 0.001); // 2 GC out of 4 = 50%
    }

    #[test]
    fn test_reverse_strand_reading_frame() {
        // Test with sequence length 9 (frame_modulus = 0, since 9 % 3 == 0, so frame_modulus = 3)
        assert_eq!(reverse_strand_reading_frame(0, 9), 2); // (3 - 1 - 0) % 3 = 2
        assert_eq!(reverse_strand_reading_frame(1, 9), 1); // (3 - 1 - 1) % 3 = 1
        assert_eq!(reverse_strand_reading_frame(2, 9), 0); // (3 - 1 - 2) % 3 = 0

        // Test with sequence length 10 (frame_modulus = 1, since 10 % 3 == 1)
        assert_eq!(reverse_strand_reading_frame(0, 10), 0); // (1 - 1 - 0) % 3 = 0

        // Test with sequence length 8 (frame_modulus = 2, since 8 % 3 == 2)
        assert_eq!(reverse_strand_reading_frame(0, 8), 1); // (2 - 1 - 0) % 3 = 1
        assert_eq!(reverse_strand_reading_frame(1, 8), 0); // (2 - 1 - 1) % 3 = 0
    }

    #[test]
    fn test_find_max_reading_frame() {
        assert_eq!(find_max_reading_frame(10, 5, 3), 0);
        assert_eq!(find_max_reading_frame(5, 10, 3), 1);
        assert_eq!(find_max_reading_frame(5, 3, 10), 2);
        assert_eq!(find_max_reading_frame(5, 5, 3), 1); // tie goes to second
    }

    #[test]
    fn test_min_of_two_integers() {
        assert_eq!(min_of_two_integers(5, 3), 3);
        assert_eq!(min_of_two_integers(3, 5), 3);
        assert_eq!(min_of_two_integers(-1, 5), -1);
        assert_eq!(min_of_two_integers(5, 5), 5);
    }

    #[test]
    fn test_calculate_kmer_index() {
        let mut seq = vec![0u8; 20];

        // Encode AC at position 0: A(00) C(01)
        crate::bitmap::set_bit(&mut seq, 1); // C

        let idx = calculate_kmer_index(2, &seq, 0);
        assert_eq!(idx, 2); // AC should give index 2
    }

    #[test]
    fn test_mer_text() {
        assert_eq!(mer_text(0, 0), "None");
        assert_eq!(mer_text(2, 0), "AA");
        assert_eq!(mer_text(2, 1), "GA");
        assert_eq!(mer_text(2, 2), "CA");
        assert_eq!(mer_text(2, 3), "TA");
    }

    #[test]
    fn test_encode_sequence_basic() {
        let sequence = b"ATCG";
        let mut encoded = vec![0u8; 10];
        let mut unknown_sequence = vec![0u8; 10];
        let mut masks = Vec::new();

        let gc = encode_sequence(
            sequence,
            &mut encoded,
            &mut unknown_sequence,
            &mut masks,
            false,
        )
        .unwrap();
        assert!((gc - 0.5).abs() < 0.001); // 2 GC out of 4 = 50%
    }

    #[test]
    fn test_encode_sequence_with_n() {
        let sequence = b"ATNG";
        let mut encoded = vec![0u8; 10];
        let mut unknown_sequence = vec![0u8; 10];
        let mut masks = Vec::new();

        let gc = encode_sequence(
            sequence,
            &mut encoded,
            &mut unknown_sequence,
            &mut masks,
            false,
        )
        .unwrap();
        assert!((gc - 0.25).abs() < 0.001); // 1 GC out of 4 = 25%
        assert!(crate::bitmap::test_bit(&unknown_sequence, 2)); // N should be marked in unknown_sequence
    }

    #[test]
    fn test_encode_sequence_masking() {
        // Create a sequence with 50+ N's to trigger masking (MASK_SIZE = 50)
        let mut sequence = b"ATC".to_vec();
        sequence.extend(vec![b'N'; 52]); // 52 N's should create a mask
        sequence.extend(b"GCG");

        let mut encoded = vec![0u8; 60];
        let mut unknown_sequence = vec![0u8; 60];
        let mut masks = Vec::new();

        let _gc = encode_sequence(
            &sequence,
            &mut encoded,
            &mut unknown_sequence,
            &mut masks,
            true,
        )
        .unwrap();
        assert!(!masks.is_empty()); // Should create at least one mask since we have 52 N's (> MASK_SIZE)
        assert_eq!(masks.len(), 1);
        assert_eq!(masks[0].begin, 3); // Start after "ATC"
        assert_eq!(masks[0].end, 54); // End at last N (3 + 52 - 1)
    }

    #[test]
    fn test_is_gc() {
        let mut seq = vec![0u8; 10];

        assert!(!is_gc(&seq, 0));

        crate::bitmap::set_bit(&mut seq, 1);
        assert!(is_gc(&seq, 0));

        let mut seq2 = vec![0u8; 10];
        crate::bitmap::set_bit(&mut seq2, 0);
        assert!(is_gc(&seq2, 0));

        let mut seq3 = vec![0u8; 10];
        crate::bitmap::set_bit(&mut seq3, 0);
        crate::bitmap::set_bit(&mut seq3, 1);
        assert!(!is_gc(&seq3, 0));
    }

    #[test]
    fn test_is_n() {
        let mut unknown_sequence = vec![0u8; 10];

        assert!(!is_n(&unknown_sequence, 0));
        assert!(!is_n(&unknown_sequence, 100));

        crate::bitmap::set_bit(&mut unknown_sequence, 5);
        assert!(is_n(&unknown_sequence, 5));
    }

    #[test]
    fn test_calculate_background_mer_frequencies() {
        let seq = vec![0u8; 20]; // All A's
        let rseq = vec![0u8; 20]; // All A's
        let mut bg = vec![0.0; 16]; // 4^2 = 16 possible 2-mers

        calculate_background_mer_frequencies(2, &seq, &rseq, 10, &mut bg);

        // Should have high frequency for AA (index 0) and low for others
        assert!(bg[0] > 0.5); // AA should be common
    }

    #[test]
    fn test_rcom_seq() {
        let seq = vec![0u8; 10];
        let unknown_sequence = vec![0u8; 10];

        // Encode A at position 0
        // A = 00, complement = T = 11

        let rseq = create_reverse_complement_sequence(&seq, &unknown_sequence, 2);

        assert!(is_t(&rseq, 1));
    }

    #[test]
    fn test_translation_table_functions() {
        assert!(uses_only_atg(6));
        assert!(uses_only_atg(10));
        assert!(!uses_only_atg(11));

        assert!(gtg_not_start(1));
        assert!(gtg_not_start(22));
        assert!(!gtg_not_start(11));

        assert!(ttg_not_start(1));
        assert!(ttg_not_start(9));
        assert!(!ttg_not_start(11));
    }

    #[test]
    fn test_start_codon_with_training() {
        let mut training = Training {
            translation_table: 11,
            ..Training::default()
        };

        let mut seq = vec![0u8; 20];

        // Encode ATG at position 0
        crate::bitmap::set_bit(&mut seq, 2); // T
        crate::bitmap::set_bit(&mut seq, 3);
        crate::bitmap::set_bit(&mut seq, 4); // G

        assert!(is_start(&seq, 0, &training));

        // Test with table that only uses ATG
        training.translation_table = 6;
        assert!(is_start(&seq, 0, &training)); // ATG still works

        // Encode GTG and test
        let mut seq2 = vec![0u8; 20];
        crate::bitmap::set_bit(&mut seq2, 0); // G
        crate::bitmap::set_bit(&mut seq2, 2); // T
        crate::bitmap::set_bit(&mut seq2, 3);
        crate::bitmap::set_bit(&mut seq2, 4); // G

        assert!(!is_start(&seq2, 0, &training)); // GTG not allowed in table 6
    }

    #[test]
    fn test_stop_codon_special_tables() {
        let mut training = Training::default();
        let mut seq = vec![0u8; 20];

        // Test AGA stop in table 2
        training.translation_table = 2;
        // Encode AGA: A(00) G(10) A(00)
        crate::bitmap::set_bit(&mut seq, 2); // G

        assert!(is_stop(&seq, 0, &training));

        // Test TCA stop in table 22
        training.translation_table = 22;
        let mut seq2 = vec![0u8; 20];
        // Encode TCA: T(11) C(01) A(00)
        crate::bitmap::set_bit(&mut seq2, 0); // T
        crate::bitmap::set_bit(&mut seq2, 1);
        crate::bitmap::set_bit(&mut seq2, 3); // C

        assert!(is_stop(&seq2, 0, &training));
    }
}
