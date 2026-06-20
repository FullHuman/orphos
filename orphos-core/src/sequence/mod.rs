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

use crate::bitmap::test_bit;
use rayon::prelude::*;

pub use io::*;
pub use processing::*;
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

const ENCODED_A: u8 = 0b00;
const ENCODED_G: u8 = 0b01;
const ENCODED_C: u8 = 0b10;
const ENCODED_T: u8 = 0b11;

/// Test if nucleotide at given position is adenine (A)
#[must_use]
#[inline]
pub fn is_a(encoded_sequence: &[u8], n: usize) -> bool {
    encoded_base_code(encoded_sequence, n) == ENCODED_A
}

/// Test if nucleotide at given position is cytosine (C)
#[must_use]
#[inline]
pub fn is_c(encoded_sequence: &[u8], n: usize) -> bool {
    encoded_base_code(encoded_sequence, n) == ENCODED_C
}

/// Test if nucleotide at given position is guanine (G)
#[must_use]
#[inline]
pub fn is_g(encoded_sequence: &[u8], n: usize) -> bool {
    encoded_base_code(encoded_sequence, n) == ENCODED_G
}

/// Test if nucleotide at given position is thymine (T)
#[must_use]
#[inline]
pub fn is_t(encoded_sequence: &[u8], n: usize) -> bool {
    encoded_base_code(encoded_sequence, n) == ENCODED_T
}

/// Test if position contains an unknown nucleotide (N)
pub fn is_n(unknown_sequence: &[u8], n: usize) -> bool {
    if n >= unknown_sequence.len() * 8 {
        return false;
    }
    test_bit(unknown_sequence, n)
}

/// Test if nucleotide at given position is G or C (high GC content indicator)
#[inline]
pub fn is_gc(encoded_sequence: &[u8], n: usize) -> bool {
    matches!(
        encoded_base_code(encoded_sequence, n),
        ENCODED_G | ENCODED_C
    )
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
#[inline]
pub fn is_atg(encoded_sequence: &[u8], pos: usize) -> bool {
    encoded_base_code(encoded_sequence, pos) == ENCODED_A
        && encoded_base_code(encoded_sequence, pos + 1) == ENCODED_T
        && encoded_base_code(encoded_sequence, pos + 2) == ENCODED_G
}

/// Test if codon at position is GTG (valine start codon)
#[inline]
pub fn is_gtg(encoded_sequence: &[u8], pos: usize) -> bool {
    encoded_base_code(encoded_sequence, pos) == ENCODED_G
        && encoded_base_code(encoded_sequence, pos + 1) == ENCODED_T
        && encoded_base_code(encoded_sequence, pos + 2) == ENCODED_G
}

/// Test if codon at position is TTG (leucine start codon)
#[inline]
pub fn is_ttg(encoded_sequence: &[u8], pos: usize) -> bool {
    encoded_base_code(encoded_sequence, pos) == ENCODED_T
        && encoded_base_code(encoded_sequence, pos + 1) == ENCODED_T
        && encoded_base_code(encoded_sequence, pos + 2) == ENCODED_G
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
    let first = encoded_base_code(encoded_sequence, pos);

    if first == ENCODED_T {
        let second = encoded_base_code(encoded_sequence, pos + 1);

        if second == ENCODED_A {
            let third = encoded_base_code(encoded_sequence, pos + 2);
            if third == ENCODED_G {
                return is_tag_stop(training.translation_table);
            }
            if third == ENCODED_A {
                return is_taa_stop(training.translation_table);
            }
        } else if second == ENCODED_G {
            return encoded_base_code(encoded_sequence, pos + 2) == ENCODED_A
                && is_tga_stop(training.translation_table);
        } else if (training.translation_table == 22 && second == ENCODED_C)
            || (training.translation_table == 23 && second == ENCODED_T)
        {
            return encoded_base_code(encoded_sequence, pos + 2) == ENCODED_A;
        }

        return false;
    }

    if training.translation_table == 2 && first == ENCODED_A {
        let second = encoded_base_code(encoded_sequence, pos + 1);
        if second != ENCODED_G {
            return false;
        }

        return matches!(
            encoded_base_code(encoded_sequence, pos + 2),
            ENCODED_A | ENCODED_G
        );
    }

    false
}

/// Calculate the GC content of a sequence region
///
/// Returns the fraction of nucleotides that are G or C within
/// the specified range (inclusive).
pub fn gc_content(encoded_sequence: &[u8], start: usize, end: usize) -> f64 {
    if start > end {
        return 0.0;
    }

    let total = end - start + 1;
    let gc_count = (start..=end)
        .filter(|&position| is_gc(encoded_sequence, position))
        .count();

    gc_count as f64 / total as f64
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

#[inline]
const fn encoded_base_code(encoded_sequence: &[u8], position: usize) -> u8 {
    let bit_index = position * 2;
    (encoded_sequence[bit_index >> 3] >> (bit_index & 0x07)) & 0b11
}

#[inline]
fn write_encoded_base_code(encoded_sequence: &mut [u8], position: usize, code: u8) {
    let bit_index = position * 2;
    let byte_index = bit_index >> 3;
    let shift = bit_index & 0x07;
    encoded_sequence[byte_index] =
        (encoded_sequence[byte_index] & !(0b11 << shift)) | ((code & 0b11) << shift);
}

#[inline]
const fn reverse_complement_byte(byte: u8) -> u8 {
    let base_0 = (!byte >> 6) & 0b11;
    let base_1 = (!byte >> 4) & 0b11;
    let base_2 = (!byte >> 2) & 0b11;
    let base_3 = !byte & 0b11;

    base_0 | (base_1 << 2) | (base_2 << 4) | (base_3 << 6)
}

fn copy_unknown_encoded_bases(
    reverse_complement_encoded_sequence: &mut [u8],
    forward_sequence: &[u8],
    unknown_sequence: &[u8],
    nucleotide_length: usize,
) {
    for (byte_index, &unknown_byte) in unknown_sequence.iter().enumerate() {
        let mut remaining_unknowns = unknown_byte;

        while remaining_unknowns != 0 {
            let bit_offset = remaining_unknowns.trailing_zeros() as usize;
            let source_position = byte_index * 8 + bit_offset;

            if source_position >= nucleotide_length {
                break;
            }

            let target_position = nucleotide_length - source_position - 1;
            let source_code = encoded_base_code(forward_sequence, source_position);
            write_encoded_base_code(
                reverse_complement_encoded_sequence,
                target_position,
                source_code,
            );

            remaining_unknowns &= remaining_unknowns - 1;
        }
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

    if nucleotide_length.is_multiple_of(4) {
        let encoded_bytes = nucleotide_length / 4;

        for (source_byte_index, &source_byte) in
            forward_sequence.iter().take(encoded_bytes).enumerate()
        {
            let target_byte_index = encoded_bytes - source_byte_index - 1;
            reverse_complement_encoded_sequence[target_byte_index] =
                reverse_complement_byte(source_byte);
        }

        copy_unknown_encoded_bases(
            &mut reverse_complement_encoded_sequence,
            forward_sequence,
            unknown_sequence,
            nucleotide_length,
        );
    } else {
        for source_position in 0..nucleotide_length {
            let target_position = nucleotide_length - source_position - 1;
            let source_code = encoded_base_code(forward_sequence, source_position);
            let target_code = if test_bit(unknown_sequence, source_position) {
                source_code
            } else {
                source_code ^ 0b11
            };

            write_encoded_base_code(
                &mut reverse_complement_encoded_sequence,
                target_position,
                target_code,
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
#[inline]
pub fn calculate_kmer_index(kmer_length: usize, encoded_sequence: &[u8], position: usize) -> usize {
    let mut kmer_index = 0;
    for offset in 0..kmer_length {
        kmer_index |=
            usize::from(encoded_base_code(encoded_sequence, position + offset)) << (offset * 2);
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

#[inline]
const fn expand_4_lanes_to_even_bits(bits: u32) -> u8 {
    ((bits & 0x1) | ((bits & 0x2) << 1) | ((bits & 0x4) << 2) | ((bits & 0x8) << 3)) as u8
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

    // SIMD splat constants hoisted outside the loop
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

    for chunk_idx in 0..chunks {
        let chunk_start = chunk_idx * CHUNK_SIZE;

        // Load 32 bytes
        let input_slice = &sequence[chunk_start..chunk_start + CHUNK_SIZE];
        let mut input_array = [0u8; 32];
        input_array.copy_from_slice(input_slice);
        let input = u8x32::from(input_array);

        // SIMD nucleotide detection
        let is_a = input.simd_eq(a_upper) | input.simd_eq(a_lower);
        let is_c = input.simd_eq(c_upper) | input.simd_eq(c_lower);
        let is_g = input.simd_eq(g_upper) | input.simd_eq(g_lower);
        let is_t = input.simd_eq(t_upper)
            | input.simd_eq(t_lower)
            | input.simd_eq(u_upper)
            | input.simd_eq(u_lower);

        let gc_mask = is_g | is_c;
        let valid_mask = is_a | is_c | is_g | is_t;

        gc_count += gc_mask.to_bitmask().count_ones();
        total_count += CHUNK_SIZE as u32;

        // Extract SIMD results for bit setting - to_bitmask() returns u32
        // let is_a_mask: u32 = is_a.to_bitmask();
        let is_c_mask: u32 = is_c.to_bitmask();
        let is_g_mask: u32 = is_g.to_bitmask();
        let is_t_mask: u32 = is_t.to_bitmask();
        let unknown_mask: u32 = !valid_mask.to_bitmask();

        let encoded_byte_start = chunk_start / 4;
        let unknown_byte_start = chunk_start / 8;

        if encoded_byte_start + 8 <= encoded_sequence.len()
            && unknown_byte_start + 4 <= unknown_sequence.len()
        {
            for byte_offset in 0..8 {
                let shift = byte_offset * 4;
                let even_bits =
                    expand_4_lanes_to_even_bits((is_g_mask | is_t_mask) >> shift & 0x0f);
                let odd_bits = expand_4_lanes_to_even_bits(
                    (is_c_mask | is_t_mask | unknown_mask) >> shift & 0x0f,
                ) << 1;
                encoded_sequence[encoded_byte_start + byte_offset] |= even_bits | odd_bits;
            }

            for byte_offset in 0..4 {
                unknown_sequence[unknown_byte_start + byte_offset] |=
                    (unknown_mask >> (byte_offset * 8)) as u8;
            }
        } else {
            for i in 0..CHUNK_SIZE {
                let pos = chunk_start + i;
                if pos >= sequence.len() || pos * 2 + 1 >= encoded_sequence.len() * 8 {
                    break;
                }

                let bit_pos = pos * 2;
                let bit_flag = 1u32 << i;

                if (is_c_mask & bit_flag) != 0 {
                    crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                } else if (is_g_mask & bit_flag) != 0 {
                    crate::bitmap::set_bit(encoded_sequence, bit_pos);
                } else if (is_t_mask & bit_flag) != 0 {
                    crate::bitmap::set_bit(encoded_sequence, bit_pos);
                    crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                } else if (unknown_mask & bit_flag) != 0 {
                    crate::bitmap::set_bit(encoded_sequence, bit_pos + 1);
                    crate::bitmap::set_bit(unknown_sequence, pos);
                }
            }
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

        // Encode C at position 0 using the physical bitmap representation.
        crate::bitmap::set_bit(&mut seq, 1); // C

        let idx = calculate_kmer_index(2, &seq, 0);
        assert_eq!(idx, 2);
    }

    fn calculate_kmer_index_reference(
        kmer_length: usize,
        encoded_sequence: &[u8],
        position: usize,
    ) -> usize {
        let mut kmer_index = 0;
        for i in 0..(2 * kmer_length) {
            let bit_pos = position * 2 + i;
            kmer_index |= usize::from(test_bit(encoded_sequence, bit_pos)) << i;
        }
        kmer_index
    }

    #[test]
    fn test_calculate_kmer_index_matches_bit_reference() {
        let sequence = b"ATCGNNGCATCGATGCGTACGATCGATCG";
        let nucleotide_length = sequence.len();
        let encoded_len = (nucleotide_length * 2).div_ceil(8);
        let unknown_len = nucleotide_length.div_ceil(8);
        let mut encoded = vec![0u8; encoded_len];
        let mut unknown_sequence = vec![0u8; unknown_len];
        let mut masks = Vec::new();

        encode_sequence(
            sequence,
            &mut encoded,
            &mut unknown_sequence,
            &mut masks,
            false,
        )
        .unwrap();

        for kmer_length in [1, 2, 3, 6] {
            for position in 0..=nucleotide_length - kmer_length {
                assert_eq!(
                    calculate_kmer_index(kmer_length, &encoded, position),
                    calculate_kmer_index_reference(kmer_length, &encoded, position),
                    "kmer_length={kmer_length}, position={position}",
                );
            }
        }
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
    fn test_simd_packed_encoding_matches_scalar_with_unknowns() {
        let sequence = b"ATCGNNNNGCATGCACTGACTNNATCGATCGXYZATCGATCGNNNNATCGATCGATCGATCG";
        let nucleotide_length = sequence.len();
        let encoded_len = (nucleotide_length * 2).div_ceil(8);
        let unknown_len = nucleotide_length.div_ceil(8);

        let mut scalar_encoded = vec![0u8; encoded_len];
        let mut scalar_unknown = vec![0u8; unknown_len];
        let mut scalar_masks = Vec::new();
        let scalar_gc = encode_sequence(
            sequence,
            &mut scalar_encoded,
            &mut scalar_unknown,
            &mut scalar_masks,
            false,
        )
        .unwrap();

        let mut simd_encoded = vec![0u8; encoded_len];
        let mut simd_unknown = vec![0u8; unknown_len];
        let simd_gc =
            encode_sequence_simd_wide_packed(sequence, &mut simd_encoded, &mut simd_unknown)
                .unwrap();

        assert_eq!(simd_encoded, scalar_encoded);
        assert_eq!(simd_unknown, scalar_unknown);
        assert_eq!(scalar_masks.len(), 0);
        assert!((simd_gc - scalar_gc).abs() < f64::EPSILON);
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

    fn create_reverse_complement_sequence_reference(
        forward_sequence: &[u8],
        unknown_sequence: &[u8],
        nucleotide_length: usize,
    ) -> Vec<u8> {
        let mut reverse_complement_encoded_sequence = vec![0; forward_sequence.len()];
        let sequence_length = nucleotide_length * 2;

        for i in 0..sequence_length {
            if !crate::bitmap::test_bit(forward_sequence, i) {
                let target_pos = if i % 2 == 0 {
                    sequence_length - i - 2
                } else {
                    sequence_length - i
                };
                if target_pos < sequence_length {
                    crate::bitmap::set_bit(&mut reverse_complement_encoded_sequence, target_pos);
                }
            }
        }

        for i in 0..nucleotide_length {
            if crate::bitmap::test_bit(unknown_sequence, i) && sequence_length >= 2 + i * 2 {
                crate::bitmap::toggle_bit(
                    &mut reverse_complement_encoded_sequence,
                    sequence_length - 1 - i * 2,
                );
                crate::bitmap::toggle_bit(
                    &mut reverse_complement_encoded_sequence,
                    sequence_length - 2 - i * 2,
                );
            }
        }

        reverse_complement_encoded_sequence
    }

    #[test]
    fn test_reverse_complement_matches_reference_with_unknowns() {
        for sequence in [
            b"ATCGNNGCATCG".as_slice(),
            b"ATCGNNGCATCGA".as_slice(),
            b"NNNNATCGXYZATCG".as_slice(),
        ] {
            let nucleotide_length = sequence.len();
            let encoded_len = (nucleotide_length * 2).div_ceil(8);
            let unknown_len = nucleotide_length.div_ceil(8);

            let mut encoded = vec![0u8; encoded_len];
            let mut unknown_sequence = vec![0u8; unknown_len];
            let mut masks = Vec::new();
            encode_sequence(
                sequence,
                &mut encoded,
                &mut unknown_sequence,
                &mut masks,
                false,
            )
            .unwrap();

            let expected = create_reverse_complement_sequence_reference(
                &encoded,
                &unknown_sequence,
                nucleotide_length,
            );
            let actual =
                create_reverse_complement_sequence(&encoded, &unknown_sequence, nucleotide_length);

            assert_eq!(actual, expected);
        }
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
